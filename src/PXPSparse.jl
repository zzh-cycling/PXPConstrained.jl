"""
    PXPSparse.jl

Sparse matrix implementations for the PXP model Hamiltonians.
This module provides memory-efficient sparse representations of PXP Hamiltonians
in various symmetry sectors, enabling simulations of larger system sizes.
"""

"""
    PXP_Ham_sparse(::Type{T}, pbc::Bool=true) where {N, T <: BitStr{N}}
    PXP_Ham_sparse(N::Int64, pbc::Bool=true)

Construct sparse matrix representation of the PXP Hamiltonian.

Creates a memory-efficient sparse matrix version of the PXP Hamiltonian,
suitable for large system sizes where dense matrices become impractical.

# Arguments
- `T::Type{BitStr{N}}` or `N::Int64`: System size specification
- `pbc::Bool=true`: Whether to use periodic boundary conditions

# Returns
- `SparseMatrixCSC{Float64}`: Sparse PXP Hamiltonian matrix

# Example
```julia
H_sparse = PXP_Ham_sparse(16, true)  # 16-site sparse PXP Hamiltonian
eigenvals = eigvals(H_sparse)  # Compute eigenvalues efficiently
```
"""
function PXP_Ham_sparse(::Type{T}, pbc::Bool=true) where {N, T <: BitStr{N}}
    # Generate Hamiltonian for PXP model, automotically contain pbc or obc
    basis=PXP_basis(T,pbc)

    l=length(basis)
    # H=spzeros(Float64,(l,l))
    I, J, V = Int[], Int[], Float64[]
    for i in 1:l
        output=actingH_PXP(T, basis[i], pbc) 
        for m in output 
            j=searchsortedfirst(basis,m)
            # H[i, j] += 1
            push!(I, i); push!(J, j); push!(V, 1.0)
        end
    end

    H = sparse(I, J, V, l, l)
    
    return H
end
PXP_Ham_sparse(N::Int64, pbc::Bool=true) = PXP_Ham_sparse(BitStr{N, Int}, pbc)

"""
    PXP_K_Ham_sparse(::Type{T}, k::Int, Omega::Float64=1.0) where {N, T <: BitStr{N}}
    PXP_K_Ham_sparse(N::Int64, k::Int)

Construct sparse Hamiltonian in momentum K subspace.

Creates a sparse matrix representation of the PXP Hamiltonian projected
onto the translational eigenspace with momentum quantum number k.

# Arguments
- `T::Type{BitStr{N}}` or `N::Int64`: System size specification
- `k::Int`: Momentum quantum number (0 ≤ k ≤ N-1)
- `Omega::Float64=1.0`: Overall energy scale

# Returns
- `SparseMatrixCSC`: Sparse Hamiltonian in K subspace (real for k=0,π)

# Example
```julia
H_k0 = PXP_K_Ham_sparse(12, 0)  # Zero-momentum subspace
eigenvals = eigvals(H_k0)
```
"""
function PXP_K_Ham_sparse(::Type{T}, k::Int, Omega::Float64=1.0) where {N, T <: BitStr{N}}
#params: a int of lattice number, momentum of system and interaction strength of system which default to be 1
#return: the Hamiltonian matrix in given K space
    @assert k in 0:N-1 "k can only be between from 0 to $(N-1)"
    basisK, basis_dic = PXP_K_basis(T, k)
    l = length(basisK)
    omegak = exp(2im * π * k / N)
    # H = spzeros(ComplexF64, (l, l))
    I, J, V = Int[], Int[], ComplexF64[]

    for i in 1:l
        n=basisK[i]
        output = actingH_PXP(T, n, true)
        for m in output
            mbar, d = get_representative(m)
            if mbar ∈ basisK
                j=searchsortedfirst(basisK, mbar)
                Yn= sqrt(length(basis_dic[n])) / N
                Ym= sqrt(length(basis_dic[mbar])) / N
                # H[i, j] += Yn/Ym * omegak^d
                push!(I, i); push!(J, j); push!(V, Yn/Ym * omegak^d)
            end
        end
    end

    H = sparse(I, J, V, l, l)

    H=(H+H')/2
    if k==0 || k==div(N,2)
        H=real(H)
    end
    return H
end
PXP_K_Ham_sparse(N::Int64, k::Int) = PXP_K_Ham_sparse(BitStr{N, Int}, k)

"""
    PXP_MSS_Ham_sparse(::Type{T}, k::Int, inv::Int64=1) where {N, T <: BitStr{N}}

Construct sparse Hamiltonian in maximum symmetry subspace (MSS).

Creates a sparse matrix representation in the subspace that respects both
translational and inversion symmetries. Only available for k=0 or k=π.

# Arguments
- `T::Type{BitStr{N}}`: System size specification  
- `k::Int`: Momentum quantum number (must be 0 or N/2)
- `inv::Int64=1`: Inversion eigenvalue (±1)

# Returns
- `SparseMatrixCSC{Float64}`: Sparse Hamiltonian in MSS

# Example
```julia
H_mss = PXP_MSS_Ham_sparse(BitStr{12, Int}, 0, 1)  # k=0, even inversion
eigenvals, eigenvecs = eigen(H_mss)
```
"""
function PXP_MSS_Ham_sparse(::Type{T}, k::Int, inv::Int64=1) where {N, T <: BitStr{N}}
    #params: a int of lattice number, momentum of system and interaction strength of system which default to be 1, k is the momentum of system, only can take 0 or pi, inv is the inversion of the Hamiltonian, only 1 or -1.
    #return: the Hamiltonian matrix in given maximum symmetry space
    @assert k == 0 || k==div(N,2) "k is expected to be 0 or $(div(N,2)), but got $k"
    @assert inv ==1 || inv==-1 "inv is expected to be 1 or -1, but got $(inv)"

    basisK, basis_dic = PXP_K_basis(T, k)
    nK = length(basisK)
    η = k == 0 ? 1.0 : -1.0

    iso = iso_K2MSS_sparse(T, k, inv)
    nMSS = size(iso, 2)

    K2MSS_col = fill(0, nK)
    K2MSS_coef = zeros(Float64, nK)
    rows, cols, vals = findnz(iso)
    for idx in eachindex(rows)
        K2MSS_col[rows[idx]] = cols[idx]
        K2MSS_coef[rows[idx]] = vals[idx]
    end

    k_index = Dict{T, Int}()
    for i in eachindex(basisK)
        k_index[basisK[i]] = i
    end

    I = Int[]
    J = Int[]
    V = Float64[]
    for ket_idx in 1:nK
        ket_col = K2MSS_col[ket_idx]
        ket_col == 0 && continue

        ket_coeff = K2MSS_coef[ket_idx]
        ket_rep = basisK[ket_idx]
        Yn = sqrt(length(basis_dic[ket_rep])) / N

        output = actingH_PXP(T, ket_rep, true)
        for out_state in output
            bra_rep, d = get_representative(out_state)
            bra_idx = get(k_index, bra_rep, 0)
            bra_idx == 0 && continue

            bra_col = K2MSS_col[bra_idx]
            bra_col == 0 && continue

            Ym = sqrt(length(basis_dic[bra_rep])) / N
            hij = Yn / Ym * (η^d)
            push!(I, bra_col)
            push!(J, ket_col)
            push!(V, K2MSS_coef[bra_idx] * hij * ket_coeff)
        end
    end

    H = sparse(I, J, V, nMSS, nMSS)
    H = (H + transpose(H)) / 2
    return H
end
PXP_MSS_Ham_sparse(N::Int64, k::Int, inv::Int64=1) = PXP_MSS_Ham_sparse(BitStr{N, Int}, k, inv)

function iso_total2K_sparse(::Type{T}, k::Int64) where {N, T <: BitStr{N}}
#Function to map the total basis to the K space basis, actually is the isometry, defined as W'*W=I, W*W'=P, P^2=P
    @assert k in 0:N-1 "k can only be between from 0 to $(N-1)"
    basis = PXP_basis(T)
    k_dic = Dict{Int, Vector{Int64}}()
    basisK = Vector{T}(undef, 0)
    # Categorize basis states by their representative
    for i in eachindex(basis)
        state = basis[i]
        category = get_representative(state)[1]
        if haskey(k_dic, category)
            push!(k_dic[category], i)
        else
            k_dic[category] = [i]
        end
    end

    for j in eachindex(basis)
        n=basis[j]
        RS = get_representative(n)[1]
        if RS == n && (k * length(k_dic[RS])) % N == 0
            push!(basisK, n)
        end
    end

    # Initialize sparse matrix
    num_states = length(basis)
    num_categories = length(keys(basisK))
    rows = Vector{Int64}[]
    cols = Vector{Int64}[]
    vals = Vector{Float64}[]

    # Fill the sparse matrix with isometry values
    for (i, state) in enumerate(basisK)
        state_indices = k_dic[state]  
        l = length(state_indices)
        push!(rows, state_indices)
        push!(cols, fill(i, l))
        push!(vals, fill(1.0 / sqrt(l), l))
    end
    
    rows = vcat(rows...)
    cols = vcat(cols...)
    vals = vcat(vals...)
    # Create sparse matrix
    iso_sparse = sparse(rows, cols, vals, num_states, num_categories)

    return iso_sparse
end
iso_total2K_sparse(N::Int, k::Int64) = iso_total2K_sparse(BitStr{N, Int}, k)



function iso_K2MSS_sparse(::Type{T}, k::Int64, inv::Int64=1) where {N, T <: BitStr{N}}
    # Function to map the MSS basis to the K space basis.
    @assert k == 0 || k == div(N, 2) "k is expected to be 0 or $(div(N,2)), but got $k"
    @assert inv == 1 || inv == -1 "inv is expected to be 1 or -1, but got $inv"

    basisK, _ = PXP_K_basis(T, k)
    nK = length(basisK)
    index_of = Dict{T, Int}()
    for (i, n) in enumerate(basisK)
        index_of[n] = i
    end

    η = k == 0 ? 1 : -1
    visited = falses(nK)

    I = Int[]
    J = Int[]
    V = Float64[]
    col = 0

    for i in 1:nK
        visited[i] && continue
        n = basisK[i]
        nR, d = get_representative(breflect(n))
        j = index_of[nR]
        phase = η^d

        if i == j
            visited[i] = true
            if phase == inv
                col += 1
                push!(I, i); push!(J, col); push!(V, 1.0)
            end
        else
            n2 = basisK[j]
            n2R, d2 = get_representative(breflect(n2))
            i2 = index_of[n2R]
            phase2 = η^d2
            @assert i2 == i "inversion pairing inconsistency between K-basis representatives"
            @assert phase2 == phase "inversion phases in paired K states are inconsistent"

            visited[i] = true
            visited[j] = true
            col += 1
            push!(I, i); push!(J, col); push!(V, 1 / sqrt(2))
            push!(I, j); push!(J, col); push!(V, inv * phase / sqrt(2))
        end
    end

    return sparse(I, J, V, nK, col)
end
iso_K2MSS_sparse(N::Int, k::Int64, inv::Int64=1) = iso_K2MSS_sparse(BitStr{N, Int}, k, inv)

function iso_total2MSS_sparse(::Type{T}, k::Int64, inv::Int64=1) where {N, T <: BitStr{N}}
    # Function to map the total basis to the MSS space basis, k can only equal to 0 or N/2(pi)
    iso = iso_total2K_sparse(T, k) * iso_K2MSS_sparse(T, k, inv)

    return iso
end
iso_total2MSS_sparse(N::Int, k::Int64, inv::Int64=1) = iso_total2MSS_sparse(BitStr{N, Int}, k, inv)
