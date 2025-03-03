
function ee(subrm::Matrix{Float64})
    #  subrm=qi.ptrace(state*state',[2 for i in 1:N],[i for i in l+1:N])
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

function ee_PXP_idx(::Type{T}, splitlis::Vector{Int64}, idx::Int64)  where {N, T <: BitStr{N}}
    """
    only calculate half the EE list
    """
    energy, states= eigen(PXP_Ham(T))
    idx_state=states[:,idx]
    EE_lis=zeros(div(length(splitlis)+1,2))
    for m in 1:div(length(splitlis)+1,2)
        subidx=rdm_PXP(BitStr{N, Int}, collect(1:splitlis[m]), idx_state)
        EE_lis[m]=ee(subidx)
    end
    EE_lis=[EE_lis; sort(EE_lis[1:div(length(splitlis)-1,2)],rev=true)]
    return EE_lis
end

function ee_PXP_state(::Type{T},splitlis::Vector{Int64},state::Vector{ET}) where {N, T <: BitStr{N}, ET <: Real}
    EE_lis=zeros(length(splitlis))
    for m in eachindex(EE_lis)
        subscar=rdm_PXP(T, collect(1:splitlis[m]), state)
        EE_lis[m]=ee(subscar)
    end
    return EE_lis
end

function ee_scaling_fig(::Type{T}, state::Vector{ET},fit::String) where {N, T <: BitStr{N}, ET <: Real}
    splitlis=Vector(1:N-1)
    EElis=ee_PXP_state(T,splitlis,state)

    if fit=="CC" 
        cent, fig=fitCCEntEntScal(EElis; mincut=1, pbc=true)
    end

    if fit=="Page"
        cent, fig=fitpage_curve(EElis; mincut=1)
    end
    return cent, fig
end

function mutual_information(::Type{T}, state::Vector{Float64}, subsystems::Tuple{Vector{Int64}, Vector{Int64}}) where {N, T <: BitStr{N}}
    A, B = subsystems
    # MI formula defined as: I(A:B) = S_A + S_B - S_AB
    # Calculate the reduced density matrices
    ρ_A = rdm_PXP(T, A, state)
    ρ_B = rdm_PXP(T, B, state)
    ρ_AB = rdm_PXP(T, vcat(A,B), state)
    # Calculate the Von Neumann entropies
    S_A = ee(ρ_A)
    S_B = ee(ρ_B)
    S_AB = ee(ρ_AB)
    # Calculate the mutual information
    I_AB = S_A + S_B - S_AB
    return I_AB
    
end

function tri_mutual_information(::Type{T}, state::Vector{Float64}, subsystems::Tuple{Vector{Int64}, Vector{Int64}, Vector{Int64}}) where {N, T <: BitStr{N}}
    A, B, C = subsystems
    # TMI formula defined as: I(A:B:C) = S_A + S_B + S_C - S_AB - S_BC - S_AC + S_ABC
    
    @time begin 
    ρ_A = rdm_PXP(T, A, state)
    ρ_B = rdm_PXP(T, B, state)
    ρ_C = rdm_PXP(T, C, state)

    ρ_AB = rdm_PXP(T, vcat(A,B), state)
    ρ_BC = rdm_PXP(T, vcat(B,C), state)
    ρ_AC = rdm_PXP(T, vcat(A,C), state)
    
    ρ_ABC = rdm_PXP(T, vcat(A,B,C), state)
    end
    myprint(stdout,"rho complete")
    
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

function domain_wall(::Type{T}, pbc::Bool=true) where {N, T <: BitStr{N}}
    """
    :param N: Number of sites
    :return: domain_wall_density diagonal elements
    The eigenvectors of this operator are going from -N to N, increasing by 2, totally N+1 eigenvectors. Number of each eigenvalues is N choose k, where k is the number of domain walls when we consider total Hilbert space. Defined as sum_i Z_i =1/2 (-1)^(i+1) * Z_i, we aim for spin systems.(S_Z= 1/2 Pauli Z)
    """
    
    basis = PXP_basis(T, pbc)
    l=length(basis)
    dw = zeros(l)

    mask = bmask(BitStr{N, Int}, collect(2:2:N)...)
    for (idx, str) in enumerate(basis)
        masked_str = flip(str, mask)
        Zi=sum([masked_str...].-1/2)
        dw[idx] = Zi
    end

    return dw
end


function particlenumber(::Type{T},pbc::Bool=true) where {N, T <: BitStr{N}}
    """
    :param N: Number of sites
    :return: Particle number operator
    """
    basis = PXP_basis(T, pbc)
    l=length(basis)
    P = zeros((l, l))
    for (idx, str) in enumerate(basis)
        P[idx, idx] = count_ones(str)
    end

    return P
    
end

function on_siten(::Type{T}, i::Int64,pbc::Bool=true)  where {N, T <:BitStr{N}}
    """
    :param N: Number of sites
    :return: Particle number operator
    """
    basis  = PXP_basis(N,pbc)
    l=length(basis)
    P = zeros((l, l))
    for (idx, str) in enumerate(basis)
        P[idx, idx] += (parse(Int, str[i], base=2))
    end

    return P
    
end

function swap_bits(bitstring::String)
    # 将比特字符串转换为字符数组
    bits = collect(bitstring)
    # 获取长度
    n = length(bits)
    
    # 创建一个新的数组来存储调换后的结果
    swapped_bits = copy(bits)
    
    # 遍历每两个元素进行调换
    for i in 1:2:n-1
        swapped_bits[i], swapped_bits[i+1] = swapped_bits[i+1], swapped_bits[i]
    end
    
    # 将字符数组转换回字符串并返回
    return String(swapped_bits)
end

function chiral(::Type{T}, pbc::Bool=true) where {N, T <:BitStr{N}}
    """
    :param N: Number of sites
    :return: Chiral operator
    N=6 found no index.
    """
    basis  = PXP_basis(N,pbc)
    l=length(basis)
    C = zeros((l, l))
    for (idx, str) in enumerate(basis)
        m=parse(Int, swap_bits(str), base=2)
        index=findfirst(x->x==m,basis)
        C[idx, index] +=1
    end

    return C
end

function ergotropy_PXP_idx(::Type{T}, l::Int64, idx::Int64) where {N, T <: BitStr{N}}
    HA=PXP_Ham(BitStr{l, Int},false)
    sub_basis=PXP_basis(BitStr{l, Int},false)
    energy, states= eigen(PXP_Ham(T))
    subenergy, substates = eigen(HA)
    
    scar_state = states[:,idx]
    subscarrho = rdm_PXP(T,collect(1:l),scar_state) 
    GS_energy=tr(subscarrho*HA)

    spectrum=eigvals(subscarrho)
    sorted_spectrum=sort(spectrum, rev=true)
    passive_energy=dot(sorted_spectrum, subenergy)

    return GS_energy, subenergy[1], passive_energy
end

function ergotropy_PXP_state(::Type{T}, l::Int64, state::Vector{ET}) where {N, T <: BitStr{N}, ET}
    HA=PXP_Ham(BitStr{l, Int},false)
    sub_basis=PXP_basis(BitStr{l, Int},false)
    subenergy, substates= eigen(HA)
    subscarrho = rdm_PXP(T,collect(1:l),scar_state) 

    GS_energy=tr(subscarrho*HA)
    spectrum=eigvals(subscarrho)
    sorted_spectrum=sort(spectrum, rev=true)
    passive_energy=dot(sorted_spectrum, subenergy)

    return GS_energy, subenergy[1], passive_energy
end

function ergotropy_PXP_idx_OBC(::Type{T}, l::Int64, idx::Int64) where {N, T <: BitStr{N}}
    HA=PXP_Ham(BitStr{l, Int},false)
    sub_basis=PXP_basis(BitStr{l, Int},false)
    energy, states= eigen(PXP_Ham(T,false))
    subenergy, substates= eigen(HA)

    state=states[:,idx]
    subrho = rdm_PXP(T, collect(1:l), state) 
     
    GS_energy=tr(subrho*HA)
    spectrum=eigvals(subrho)
    sorted_spectrum=sort(spectrum, rev=true)
    passive_energy=dot(sorted_spectrum, subenergy)

    return GS_energy, subenergy[1], passive_energy
end
