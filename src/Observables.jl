"""
    Observables.jl

Functions for computing physical observables in the PXP model.
This module provides tools for calculating entanglement measures, quantum information quantities,
and other physical observables from quantum states and reduced density matrices.
"""

"""
    ee(subrm::Matrix{ET}) where {ET}

Calculate the Von Neumann entanglement entropy of a reduced density matrix.

Computes S = -Tr(ρ log ρ) where ρ is the reduced density matrix.
Handles numerical precision by filtering out very small eigenvalues.

# Arguments
- `subrm::Matrix{ET}`: Reduced density matrix (must be Hermitian)

# Returns
- `Float64`: Von Neumann entropy

# Example
```julia
N=10
psi = zeros(length(PXP_basis(N)))  # Example state
rdm = rdm_PXP(10, collect(1:div(N,2)), psi)
entropy = ee(rdm)
```
"""
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

"""
    ee_PXP_idx(N::Int64, splitlis::Vector{Int64}, idx::Int64)

Calculate entanglement entropy profile for a specific eigenstate.

Computes the entanglement entropy between different bipartitions of the system
for the idx-th eigenstate of the PXP Hamiltonian.

# Arguments
- `N::Int64`: System size
- `splitlis::Vector{Int64}`: List of bipartition sizes to compute
- `idx::Int64`: Index of the eigenstate to analyze

# Returns
- `Vector{Float64}`: Entanglement entropies for each bipartition

# Example
```julia
# Compute EE profile for the ground state (idx=1)
ee_profile = ee_PXP_idx(12, [1,2,3,4,5,6], 1)
```
"""
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

"""
    ee_PXP_state(N::Int64, splitlis::Vector{Int64}, state::Vector{ET}, MSS::Bool=false) where {ET}

Calculate entanglement entropy profile for a given quantum state.

Computes the entanglement entropy for different subsystem sizes for an
arbitrary quantum state in the PXP model.

# Arguments
- `N::Int64`: System size
- `splitlis::Vector{Int64}`: List of subsystem sizes to compute
- `state::Vector{ET}`: Quantum state vector
- `MSS::Bool=false`: Whether the state is in maximum symmetry subspace

# Returns
- `Vector{Float64}`: Entanglement entropies for each subsystem size

# Example
```julia
# Compute EE profile for a custom state
N=10
psi = zeros(length(PXP_basis(N)))  # Example state
ee_profile = ee_PXP_state(10, [1,2,3,4,5,6], psi)
```
"""
function ee_PXP_state(N::Int64,splitlis::Vector{Int64},state::Vector{ET}, MSS::Bool=false) where {ET}
    EE_lis=zeros(length(splitlis))
    for m in eachindex(EE_lis)
        if MSS
            subrho = rdm_PXP_MSS(N, collect(1:splitlis[m]), state, 0)
        else
            subrho = rdm_PXP(N, collect(1:splitlis[m]), state)
        end
        EE_lis[m]=ee(subrho)
    end
    return EE_lis
end

"""
    mutual_information(N::Int64, subsystems::Tuple{Vector{Int64}, Vector{Int64}}, state::Vector{ET}) where {ET}

Calculate the mutual information between two subsystems.

Computes I(A:B) = S_A + S_B - S_AB where S represents the Von Neumann entropy.
Mutual information quantifies the total correlation between subsystems A and B.

# Arguments
- `N::Int64`: Total system size
- `subsystems::Tuple{Vector{Int64}, Vector{Int64}}`: Tuple of (A_sites, B_sites)
- `state::Vector{ET}`: Quantum state vector

# Returns
- `Float64`: Mutual information I(A:B)

# Example
```julia
A_sites = [1, 2, 3]
B_sites = [7, 8, 9]
mi = mutual_information(10, (A_sites, B_sites), psi)
```
"""
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



"""
    tri_mutual_information(N::Int64, subsystems::Tuple{Vector{Int64}, Vector{Int64}, Vector{Int64}}, state::Vector{ET}) where {ET}

Calculate the tripartite mutual information between three subsystems.

Computes I(A:B:C) = S_A + S_B + S_C - S_AB - S_BC - S_AC + S_ABC.
This measures genuine three-partie quantum correlations.

# Arguments
- `N::Int64`: Total system size
- `subsystems::Tuple{Vector{Int64}, Vector{Int64}, Vector{Int64}}`: Tuple of (A_sites, B_sites, C_sites)
- `state::Vector{ET}`: Quantum state vector

# Returns
- `Float64`: Tripartite mutual information I(A:B:C)

# Example
```julia
A_sites = [1, 2]
B_sites = [5, 6]
C_sites = [9, 10]
tmi = tri_mutual_information(10, (A_sites, B_sites, C_sites), psi)
```
"""
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

"""
    qfi(Ob::Vector{Float64}, state::Vector{T}) where T

Calculate the Quantum Fisher Information for a diagonal observable.

Computes F_Q = 4 * Var(O) where Var(O) is the variance of the observable O
in the given quantum state. The QFI quantifies the sensitivity of the state
to changes in a parameter encoded in the observable.

# Arguments
- `Ob::Vector{Float64}`: Diagonal observable (eigenvalues)
- `state::Vector{T}`: Quantum state vector

# Returns
- `Float64`: Quantum Fisher Information

# Example
```julia
# For a spin-1/2 observable
magnetization = magnetization = vcat(foldr(vcat, (fill([0.0, 1.0], 61))), 0.0)  # diagonal elements
qfi_val = qfi(magnetization, psi)
```
"""
function qfi(Ob::Vector{Float64}, state::Vector{T}) where T    
    # Calculate the quantum fisher information, espeically for diagonal operators.
    DeltaOb=state'*(Ob.^2 .*state)-(state'*(Ob.*state))^2
    # Calculate the Quantum Fisher Information
    # For spin 1/2, w/o 4
    F_Q = 4*DeltaOb

    return F_Q
end

function qfi(Ob::Matrix{Float64}, state::Vector{T}) where T
    rho=state*state'
    DeltaOb=tr(rho*Ob^2)-tr(rho*Ob)^2
    # Calculate the Quantum Fisher Information
    # For spin 1/2, w/o 4
    F_Q = 4*DeltaOb

    return F_Q
end

"""
    anti_ferro_order(::Type{T}, pbc::Bool=true) where {N, T <: BitStr{N}}
    anti_ferro_order(N::Int64, pbc::Bool=true)

Compute the antiferromagnetic order parameter for each basis state.

Calculates the staggered magnetization ∑ᵢ (-1)^(i+1) Zᵢ for each basis state,
where Zᵢ = 2nᵢ - 1 and nᵢ is the occupation number.

# Arguments
- `T::Type{BitStr{N}}` or `N::Int64`: System size specification
- `pbc::Bool=true`: Whether to use periodic boundary conditions

# Returns
- `Vector{Float64}`: Antiferromagnetic order for each basis state

# Example
```julia
anti_ferro = anti_ferro_order(8, true)
```
"""
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

"""
    domain_wall_density(::Type{T}, pbc::Bool=true) where {N, T <: BitStr{N}}

Compute the domain wall density for each basis state.

Calculates (1/N) ∑ᵢ (1-ZᵢZ_{i+1})/2 where Zᵢ are Pauli-Z eigenvalues.
Measures the fraction of nearest-neighbor pairs with opposite spins.

# Arguments
- `T::Type{BitStr{N}}`: System size specification
- `pbc::Bool=true`: Whether to use periodic boundary conditions

# Returns
- `Vector{Float64}`: Domain wall density for each basis state

# Example
```julia
dwd = domain_wall_density(BitStr{8, Int}, true)
```
"""
function domain_wall_density(::Type{T}, pbc::Bool=true) where {N, T <: BitStr{N}}
    # return domain_wall_density， defined as 1/N sum_i (1-Z_i*Z_{i+1})/2
    basis = PXP_basis(T, pbc)
    l=length(basis)
    dwd = zeros(l)

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
        dwd[idx] = sum_walls / N
    end

    return dwd
end
domain_wall_density(N::Int64, pbc::Bool=true) = domain_wall_density(BitStr{N, Int}, pbc)

"""
    particlenumber(::Type{T}, pbc::Bool=true) where {N, T <: BitStr{N}}

Compute the particle number for each basis state.

Simply counts the number of excited sites (1s) in each basis state.

# Arguments
- `T::Type{BitStr{N}}`: System size specification
- `pbc::Bool=true`: Whether to use periodic boundary conditions

# Returns
- `Vector{Float64}`: Particle number for each basis state

# Example
```julia
n_particles = particlenumber(BitStr{8, Int}, true)
```
"""
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

"""
    ergotropy_PXP_state(N::Int64, l::Int64, state::Vector{ET}, pbc::Bool=true) where {ET}

Calculate the ergotropy of a PXP state. Ergotropy is the maximum extractable work from a quantum state.

# Arguments
- `N::Int64`: Total system size
- `l::Int64`: Subsystem size (number of sites in the reduced density matrix)
- `state::Vector{ET}`: Quantum state vector
- `pbc::Bool=true`: Whether to use periodic boundary conditions

# Returns
- `Float64`: Ground state energy, first eigenvalue, and passive energy

# Example
```julia
# For a PXP state with 10 sites and reduced density matrix of size 5
ergotropy = ergotropy_PXP_state(10, 5, psi)
```
"""
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

"""
    ergotropy_PXP_MSS_state(L::Int, l::Int, state::Vector{T}, k::Int=0, inv::Int64=1) where T

Calculate the ergotropy of a PXP state in the maximum symmetry subspace (MSS). Here the reference energy is still the total Hilbert space OBC PXP_Ham.

# Arguments
- `L::Int`: Total system size
- `l::Int`: Subsystem size (number of sites in the reduced density matrix)
- `state::Vector{T}`: Quantum state vector in MSS basis
- `k::Int=0`: Momentum quantum number (default 0)
- `inv::Int64=1`: Inversion symmetry (default 1)

# Returns
- `Float64`: Ground state energy, first eigenvalue, and passive energy

# Example
```julia
# For a PXP state in MSS with 12 sites and reduced density matrix of size 6
psi_mss = zeros(length(PXP_MSS_basis(BitStr{12, Int}, 0)[1]))  # Example state in MSS basis
ergotropy = ergotropy_PXP_MSS_state(12, 6, psi_mss)
```
"""
function ergotropy_PXP_MSS_state(L::Int, l::Int, state::Vector{T}, k::Int=0, inv::Int64=1) where T
    HA = PXP_Ham(l, false)
    subenergy, substates = eigen(HA)
    
    subrho = rdm_PXP_MSS(L, collect(1:l), state, k, inv)
    
    GS_energy = tr(subrho * HA)
    spectrum = eigvals(subrho)
    sorted_spectrum = sort(spectrum, rev=true)
    passive_energy = dot(sorted_spectrum, subenergy)
    
    return GS_energy, subenergy[1], passive_energy
end

function inversion_matrix(::Type{T}) where {N, T <: BitStr{N}}
    basis=PXP_basis(T)
    l=length(basis)
    Imatrix=zeros((l,l))
    # reversed_basis = map(breflect, basis) # The optimization try of using map function and broadcast
    reversed_basis=similar(basis)
    for i in eachindex(basis)
        reversed_basis[i]=breflect(basis[i])
    end
    # Imatrix[CartesianIndex.(collect(1:length(basis)),searchsortedfirst.(Ref(basis), reversed_basis))].+=1.0
    for i in eachindex(basis)
        output=reversed_basis[i]
        j=searchsortedfirst(basis,output)
        Imatrix[i,j]+=1.0
    end
   
    return Imatrix
end
inversion_matrix(N::Int) = inversion_matrix(BitStr{N, Int})

function translation_matrix(::Type{T}) where {N, T <: BitStr{N}}
    basis=PXP_basis(T)  
    Mat=zeros(Float64,(length(basis),length(basis)))
    for (i,n) in enumerate(basis)
        m=cyclebits(n)
        j=searchsortedfirst(basis, m)
        Mat[i,j]=1.0
    end
    
    return Mat
end
translation_matrix(N::Int) = translation_matrix(BitStr{N, Int})
