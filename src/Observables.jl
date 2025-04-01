
function ee(subrm::Matrix{ET}) where {ET}
    #  subrm=qi.ptrace(state*state',[2 for i in 1:N],[i for i in l+1:N])
    @assert ishermitian(subrm) "The reduced density matrix is not hermitian."
    spectrum=eigvals(subrm)
    EE=0
    for i in eachindex(spectrum)
        v=abs(spectrum[i])
            if v>1e-8
                EE+=-v*log(v)
            end
    end

    return EE
end

function ee_PXP_idx(N::Int64, splitlis::Vector{Int64}, idx::Int64) 
#only calculate half the EE list

    energy, states= eigen(PXP_Ham(BitStr{N, Int}))
    idx_state=states[:,idx]
    EE_lis=zeros(div(length(splitlis)+1,2))
    for m in 1:div(length(splitlis)+1,2)
        subidx=rdm_PXP(N, collect(1:splitlis[m]), idx_state)
        EE_lis[m]=ee(subidx)
    end
    EE_lis=[EE_lis; sort(EE_lis[1:div(length(splitlis)-1,2)],rev=true)]
    return EE_lis
end

function ee_PXP_state(N::Int64,splitlis::Vector{Int64},state::Vector{ET}) where {ET}
    EE_lis=zeros(length(splitlis))
    for m in eachindex(EE_lis)
        subrho=rdm_PXP(N, collect(1:splitlis[m]), state)
        EE_lis[m]=ee(subrho)
    end
    return EE_lis
end

function ee_PXP_scaling_fig(N::Int64, state::Vector{ET},fit::String) where {ET}
    splitlis=Vector(1:N-1)
    EElis=ee_PXP_state(N, splitlis, state)

    if fit=="CC" 
        cent, fig=fitCCEntEntScal(EElis; mincut=1, pbc=true)
    end

    if fit=="Page"
        cent, fig=fitpage_curve(EElis; mincut=1)
    end
    return cent, fig
end

function mutual_information(N::Int64, subsystems::Tuple{Vector{Int64}, Vector{Int64}}, state::Vector{ET}) where {ET}
    A, B = subsystems
    # MI formula defined as: I(A:B) = S_A + S_B - S_AB
    # Calculate the reduced density matrices
    ρ_A = rdm_PXP(N, A, state)
    ρ_B = rdm_PXP(N, B, state)
    ρ_AB = rdm_PXP(N, vcat(A, B), state)
    # Calculate the Von Neumann entropies
    S_A = ee(ρ_A)
    S_B = ee(ρ_B)
    S_AB = ee(ρ_AB)
    # Calculate the mutual information
    I_AB = S_A + S_B - S_AB
    return I_AB
    
end



function tri_mutual_information(N::Int64, subsystems::Tuple{Vector{Int64}, Vector{Int64}, Vector{Int64}}, state::Vector{ET}) where {ET}
    A, B, C = subsystems
    # TMI formula defined as: I(A:B:C) = S_A + S_B + S_C - S_AB - S_BC - S_AC + S_ABC
    
    ρ_A = rdm_PXP(N, A, state)
    ρ_B = rdm_PXP(N, B, state)
    ρ_C = rdm_PXP(N, C, state)

    ρ_AB = rdm_PXP(N, vcat(A,B), state)
    ρ_BC = rdm_PXP(N, vcat(B,C), state)
    ρ_AC = rdm_PXP(N, vcat(A,C), state)
    
    ρ_ABC = rdm_PXP(N, vcat(A,B,C), state)
    
    # Calculate the Von Neumann entropies
    
    S_A = ee(ρ_A)
    S_B = ee(ρ_B)
    S_C = ee(ρ_C)
    S_AB = ee(ρ_AB)
    S_BC = ee(ρ_BC)
    S_AC = ee(ρ_AC)
    S_ABC = ee(ρ_ABC)

    # Calculate the mutual information
    I_ABC = S_A + S_B + S_C - S_AB - S_BC - S_AC + S_ABC
    
    return I_ABC
end

function qfi(Ob::Vector{Float64}, state::Vector{T}) where {T <: Real}
    # Calculate the quantum fisher information, espeically for diagonal operators.
    DeltaOb=state'*(Ob.^2 .*state)-(state'*(Ob.*state))^2
    # Calculate the Quantum Fisher Information
    # For spin 1/2, w/o 4
    F_Q = 4*DeltaOb

    return F_Q
end

function qfi(Ob::Matrix{Float64}, state::Vector{T}) where {T <: Real}
    rho=state*state'
    DeltaOb=tr(rho*Ob^2)-tr(rho*Ob)^2
    # Calculate the Quantum Fisher Information
    # For spin 1/2, w/o 4
    F_Q = 4*DeltaOb

    return F_Q
end

function anti_ferro_order(::Type{T}, pbc::Bool=true) where {N, T <: BitStr{N}}
#param N: Number of sites
#return:  antiferromagnetic order diagonal elements
#The eigenvectors of this operator are going from -N to N, increasing by 2, totally N+1 eigenvectors. Number of each eigenvalues is N choose k, 
#where k is the number of domain walls when we consider total Hilbert space. Defined as sum_i Z_i =1/2 (-1)^(i+1) * Z_i, we aim for spin systems.(S_Z= 1/2 Pauli Z)
    basis = PXP_basis(T, pbc)
    l=length(basis)
    anti_ferro = zeros(l)

    mask = bmask(T, collect(2:2:N)...)
    for (idx, str) in enumerate(basis)
        masked_str = flip(str, mask)
        Zi=sum([masked_str...].-1/2)
        anti_ferro[idx] = Zi
    end

    return anti_ferro
end
anti_ferro_order(N::Int64, pbc::Bool=true) = anti_ferro_order(BitStr{N, Int}, pbc)

function domain_wall_density(::Type{T}, pbc::Bool=true) where {N, T <: BitStr{N}}
    # return domain_wall_density， defined as 1/N sum_i (1-Z_i*Z_{i+1})/2
    basis = PXP_basis(T, pbc)
    l=length(basis)
    anti_ferro = zeros(l)

    for (idx, str) in enumerate(basis)
        sum_walls = 0
        for i in 1:N-1
            # Convert bits to Z_i values: 0->1, 1->-1
            z_i = str[i] == 0 ? 1 : -1
            z_ip1 = str[i+1] == 0 ? 1 : -1
            # Add (1-Z_i*Z_{i+1})/2 which is 1 for a domain wall, 0 otherwise
            sum_walls += (1 - z_i * z_ip1) / 2
        end
        
        # Handle periodic boundary condition if needed
        if pbc
            z_1 = str[1] == 0 ? 1 : -1
            z_N = str[N] == 0 ? 1 : -1
            sum_walls += (1 - z_N * z_1) / 2
        end
        
        # Normalize by N
        anti_ferro[idx] = sum_walls / N
    end

    return anti_ferro
end
domain_wall_density(N::Int64, pbc::Bool=true) = domain_wall_density(BitStr{N, Int}, pbc)

function particlenumber(::Type{T},pbc::Bool=true) where {N, T <: BitStr{N}}
#param N: Number of sites,return: Particle number operator

    basis = PXP_basis(T, pbc)
    l=length(basis)
    P = zeros((l, l))
    for (idx, str) in enumerate(basis)
        P[idx, idx] = count_ones(str)
    end

    return P
end
particlenumber(N::Int64, pbc::Bool=true) = particlenumber(BitStr{N, Int}, pbc)

function on_siten(::Type{T}, i::Int64,pbc::Bool=true)  where {N, T <:BitStr{N}}
#param N: Number of sites,return: Particle number operator
    basis  = PXP_basis(T,pbc)
    l=length(basis)
    P = zeros((l, l))
    for (idx, str) in enumerate(basis)
        P[idx, idx] += str[N+1-i]
    end

    return P
    
end
on_siten(N::Int64, i::Int64, pbc::Bool=true) = on_siten(BitStr{N, Int}, i, pbc)

function ergotropy_PXP_idx(N::Int64, l::Int64, idx::Int64, pbc::Bool=true)
    HA=PXP_Ham(BitStr{l, Int}, false)
    energy, states= eigen(PXP_Ham(BitStr{N, Int}, pbc))
    subenergy, substates = eigen(HA)
    
    state = states[:,idx]
    subrho = rdm_PXP(N, collect(1:l), state, pbc) 
    GS_energy=tr(subrho*HA)

    spectrum=eigvals(subrho)
    sorted_spectrum=sort(spectrum, rev=true)
    passive_energy=dot(sorted_spectrum, subenergy)

    return GS_energy, subenergy[1], passive_energy
end

function ergotropy_PXP_state(N::Int64, l::Int64,  state::Vector{ET}, pbc::Bool=true) where {ET}
    HA=PXP_Ham(BitStr{l, Int}, false)
    subenergy, substates= eigen(HA)
    subrho = rdm_PXP(N, collect(1:l), state, pbc) 

    GS_energy=tr(subrho*HA)
    spectrum=eigvals(subrho)
    sorted_spectrum=sort(spectrum, rev=true)
    passive_energy=dot(sorted_spectrum, subenergy)

    return GS_energy, subenergy[1], passive_energy
end

function ergotropy_PXP_MSS_state(L::Int, l::Int, state::Vector{T}, k::Int=0) where T
    HA = PXP_Ham(l, false)
    subenergy, substates = eigen(HA)
    
    subrho = rdm_PXP_MSS(L, collect(1:l), state, k)
    
    GS_energy = tr(subrho * HA)
    spectrum = eigvals(subrho)
    sorted_spectrum = sort(spectrum, rev=true)
    passive_energy = dot(sorted_spectrum, subenergy)
    
    return GS_energy, subenergy[1], passive_energy
end
