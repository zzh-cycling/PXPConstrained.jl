"""
    PXPSymmetry.jl

Functions for implementing translational and inversion symmetries in the PXP model.
This module provides tools for constructing momentum (K) and maximum symmetry subspace (MSS)
representations, along with the corresponding isometries and basis transformations.
"""

"""
    iso_total2K(::Type{T}, k::Int64) where {N, T <: BitStr{N}}
    iso_total2K(N::Int, k::Int64)

Construct the isometry mapping from total basis to momentum K subspace.

Creates the transformation matrix W such that W'*W=I and W*W'=P where P is the
projector onto the K-momentum subspace. This implements translational symmetry.

# Arguments
- `T::Type{BitStr{N}}` or `N::Int`: System size specification
- `k::Int64`: Momentum quantum number (0 ≤ k ≤ N-1)

# Returns
- `Matrix{Float64}`: Isometry matrix mapping total space to K subspace

# Example
```julia
total_ham = PXP_Ham(BitStr{10, Int})  # Total Hamiltonian for N=10
W = iso_total2K(10, 0)  # Map to zero-momentum subspace
K_ham = W' * total_ham * W  # Project Hamiltonian to K subspace
```
"""
function iso_total2K(::Type{T}, k::Int64) where {N, T <: BitStr{N}}
#Function to map the total basis to the K space basis, actually is the isometry, defined as W'*W=I, W*W'=P, P^2=P
    @assert 0<=k<=N-1 "k is expected to be in [0, $(N-1)], but got $k"

    basis = PXP_basis(T)

    k_dic = Dict{Int, Vector{Int64}}()
    basisK = Vector{T}(undef, 0)
    for i in eachindex(basis)
        state=basis[i]
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

    iso = zeros((length(basis), length(keys(basisK))))
    
    for (i, state) in enumerate(basisK)
        state_indices = k_dic[state]  
        l = length(state_indices)
        iso[state_indices, i] .= 1/sqrt(l)
    end

    return iso
end
iso_total2K(N::Int, k::Int64) = iso_total2K(BitStr{N, Int}, k)

"""
    mapstate_K2total(::Type{T}, state::Vector{ET}, k::Int64) where {N, T <: BitStr{N}, ET}
    mapstate_K2total(N::Int, state::Vector{ET}, k::Int64) where {ET}

Map a state from momentum K subspace to the total space.

Transforms a quantum state from the momentum-selected subspace back to the
full Hilbert space using translational symmetry.

# Arguments
- `T::Type{BitStr{N}}` or `N::Int`: System size specification
- `state::Vector{ET}`: State vector in K subspace
- `k::Int64`: Momentum quantum number (0 ≤ k ≤ N-1)

# Returns
- `Vector{ET}`: State vector in total space

# Example
```julia
k_state = zeros(length(PXP_K_basis(10,0)[1]))  # Some state in K=0 subspace
k_state[1] = 1.0  # Set first component to 1
total_state = mapstate_K2total(10, k_state, 0)
```
"""
function mapstate_K2total(::Type{T}, state::Vector{ET}, k::Int64) where {N, T <: BitStr{N}, ET}
    # Map the K space state to total space state
    @assert 0<=k<=N-1 "k is expected to be in [0, $(N-1)], but got $k"

    basis = PXP_basis(T)
    k_dic = Dict{Int, Vector{Int64}}()
    basisK = Vector{T}(undef, 0)
    for i in eachindex(basis)
        base=basis[i]
        category = get_representative(base)[1]
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
    
    total_state = zeros(ET, length(basis))
    for (i, basis) in enumerate(basisK)
        state_indices = k_dic[basis]  
        l = length(state_indices)
        total_state[state_indices] .+= 1/sqrt(l) * state[i]
    end

    return total_state
end
mapstate_K2total(N::Int, state::Vector{ET}, k::Int64) where {ET} = mapstate_K2total(BitStr{N, Int}, state, k)

"""
    rdm_PXP_K(::Type{T}, subsystems::Vector{Int64}, kstate::Vector{ET}, k::Int64) where {N,T <: BitStr{N}, ET}
    rdm_PXP_K(N::Int, subsystems::Vector{Int64}, state::Vector{ET}, k::Int64) where {ET}

Compute reduced density matrix for a state in momentum K subspace.

Transforms the state from K subspace to total space and then computes
the reduced density matrix for the specified subsystem.

# Arguments
- `T::Type{BitStr{N}}` or `N::Int`: System size specification
- `subsystems::Vector{Int64}`: Subsystem sites indices
- `kstate::Vector{ET}`: State vector in K subspace
- `k::Int64`: Momentum quantum number

# Returns
- `Matrix{ET}`: Reduced density matrix of the subsystem

# Example
```julia
rdm = rdm_PXP_K(10, [1,2,3], k_state, 0)
```
"""
function rdm_PXP_K(::Type{T}, subsystems::Vector{Int64},kstate::Vector{ET}, k::Int64) where {N,T <: BitStr{N}, ET}
    @assert length(kstate) == length(PXP_K_basis(T,k)[1]) "state length is expected to be $(length(PXP_K_basis(T, k)[1])), but got $(length(kstate))"
    state = mapstate_K2total(T, kstate, k)
    reduced_dm = rdm_PXP(T, subsystems, state)
    return reduced_dm
end
rdm_PXP_K(N::Int, subsystems::Vector{Int64},state::Vector{ET}, k::Int64) where {ET} = rdm_PXP_K(BitStr{N, Int}, subsystems, state, k)


"""
    iso_K2MSS(::Type{T}, k::Int64, inv::Int64=1) where {N, T <: BitStr{N}}
    iso_K2MSS(N::Int, k::Int64, inv::Int64=1)

Construct isometry from momentum K subspace to maximum symmetry subspace (MSS).

Maps from the momentum subspace to the subspace that respects both translational
and inversion symmetries. Only works for k=0 or k=π (N/2).

# Arguments
- `T::Type{BitStr{N}}` or `N::Int`: System size specification
- `k::Int64`: Momentum quantum number (must be 0 or N/2)
- `inv::Int64=1`: Inversion eigenvalue (±1)

# Returns
- `Matrix{Float64}`: Isometry matrix from K space to MSS

# Example
```julia
W_MSS = iso_K2MSS(10, 0, 1)  # Zero momentum, even inversion
MSS_Ham = W_MSS' * PXP_K_Ham(10, 0) * W_MSS  # Project Hamiltonian to MSS subspace
```
"""
function iso_K2MSS(::Type{T}, k::Int64, inv::Int64=1) where {N, T <: BitStr{N}}
#Function to map the MSS basis to the K space basis
    @assert k == 0 || k==div(N,2) "k is expected to be 0 or $(div(N,2)), but got $k"
    @assert inv ==1 || inv==-1 "inv is expected to be 1 or -1, but got $(inv)"
    basisK, k_dic = PXP_K_basis(T, k)

    MSS_dic = Dict{Int, Vector{Int64}}()
    # MSS_dic is a dictionary, the key is the representative state of the inversion of n, and the value is the index of the state in the basisK. NOTE that MSS_dic is not sorted, so we need to sort it later.
    qlist = Vector{Int}(undef, 0)
    # Below procedure is to collapse the extra basis in K space that can be converted mutually to MSS space.
    if inv==1 && k==0 || inv==-1 && k==div(N,2)
        for i in eachindex(basisK)
            n = basisK[i]
            # here we calculate the representative state of the inversion of n
            nR = get_representative(breflect(n))[1]
            # For example, n = 41, nR=37, then we only need to keep n=37, and n=41 will be removed.
            if n <= min(nR, n)
                push!(qlist, length(Set([n, nR])))
            end
            n = min(nR, n)
                if haskey(MSS_dic, n)
                    push!(MSS_dic[n], i)
                else
                    MSS_dic[n] = [i]
                end
        end

    else
        for i in eachindex(basisK)
            n = basisK[i]
            nR = get_representative(breflect(n))[1]
            if n != nR
                n = min(nR, n)
                if haskey(MSS_dic, n)
                    push!(MSS_dic[n], i)
                else
                    MSS_dic[n] = [i]
                end
                push!(qlist, 2)
            end     
        end    
    end

    iso = zeros((length(basisK), length(MSS_dic)))
    MSS_dic=sort(MSS_dic)
    for (i, state_index) in enumerate(values(MSS_dic))
        iso[state_index, i] .= 1/sqrt(qlist[i])
    end

    return iso
end
function iso_K2MSS(::Type{T}, k::Int64, inv::Int64=1) where {N, T <: BitStr{N}}
     # Only k=0 and k=pi(=N/2) are allowed for MSS construction.
     @assert k == 0 || k == div(N, 2) "k is expected to be 0 or $(div(N,2)),but got $k"
     # Inversion eigenvalue must be +1 (even) or -1 (odd).
     @assert inv == 1 || inv == -1 "inv is expected to be 1 or -1, but got $inv"

     basisK, basis_dic = PXP_K_basis(T, k)
     nK = length(basisK)

     # k=0: reflection carries no momentum phase and the representative-pair
     # construction gives an explicit deterministic basis.
     if k == 0
         rep_dict = Dict{T, Vector{Int64}}()
         for i in 1:nK
             n = basisK[i]
             nR = get_representative(breflect(n))[1]
             if inv == 1 || n != nR
                 rep = min(n, nR)
                 if haskey(rep_dict, rep)
                     push!(rep_dict[rep], i)
                 else
                     rep_dict[rep] = [i]
                 end
             end
         end

         reps = sort(collect(keys(rep_dict)))
         nMSS = length(reps)
         iso = zeros(nK, nMSS)

         for (col, rep) in enumerate(reps)
             indices = rep_dict[rep]
             if inv == 1
                 if length(indices) == 1
                     iso[indices[1], col] = 1.0
                 elseif length(indices) == 2
                     iso[indices[1], col] = 1 / sqrt(2)
                     iso[indices[2], col] = 1 / sqrt(2)
                 end
             else
                 @assert length(indices) == 2 "I=-1 sector: expected 2 indices, got $(length(indices))"
                 idx1, idx2 = indices[1], indices[2]
                 if basisK[idx1] == rep
                     iso[idx1, col] = 1 / sqrt(2)
                     iso[idx2, col] = -1 / sqrt(2)
                 else
                     iso[idx1, col] = -1 / sqrt(2)
                     iso[idx2, col] = 1 / sqrt(2)
                 end
             end
         end
         return iso
     end

     # k=pi: reflection picks up a state-dependent phase omega_k^dR.
     # Build the exact reflection operator in K basis and project to its ±1 eigenspaces.
     omegak = exp(2im * π * k / N)
     Ik = zeros(ComplexF64, nK, nK)
     k_index = Dict{T, Int}()
     for (i, st) in enumerate(basisK)
         k_index[st] = i
     end

     for i in 1:nK
         n = basisK[i]
         nR, dR = get_representative(breflect(n))
         j = get(k_index, nR, 0)
         j == 0 && continue
         Yn = sqrt(length(basis_dic[n])) / N
         YR = sqrt(length(basis_dic[nR])) / N
         Ik[j, i] += Yn / YR * omegak^dR
     end

     vals, vecs = eigen(Hermitian((Ik + Ik') / 2))
     idx = findall(x -> isapprox(x, inv, atol=1e-8), vals)
     return real(vecs[:, idx])
  end
iso_K2MSS(N::Int, k::Int64, inv::Int64=1) = iso_K2MSS(BitStr{N, Int}, k, inv)

"""
    mapstate_MSS2K(::Type{T}, state::Vector{ET}, k::Int64, inv::Int64=1) where {N, T <: BitStr{N}, ET}
    mapstate_MSS2K(N::Int, state::Vector{ET}, k::Int64, inv::Int64=1) where {ET}

Map a state from maximum symmetry subspace (MSS) to momentum K subspace.

Transforms a quantum state from the MSS back to the momentum-selected subspace.

# Arguments
- `T::Type{BitStr{N}}` or `N::Int`: System size specification
- `state::Vector{ET}`: State vector in MSS
- `k::Int64`: Momentum quantum number (0 ≤ k ≤ N-1)
- `inv::Int64`: Inversion flag (-1 or 1)

# Returns
- `Vector{ET}`: State vector in K subspace

# Example
```julia
mss_state = zeros(length(PXP_MSS_basis(BitStr{10, Int}, 0)[1]))  # Some state in MSS
k_state = mapstate_MSS2K(10, mss_state, 0, 1)
```
"""
function mapstate_MSS2K(::Type{T}, state::Vector{ET}, k::Int64, inv::Int64=1) where {N, T <: BitStr{N}, ET}
    @assert k == 0 || k==div(N,2) "k is expected to be 0 or $(div(N,2)), but got $k"
    @assert inv ==1 || inv==-1 "inv is expected to be 1 or -1, but got $(inv)"

    basisK, k_dic = PXP_K_basis(T, k)

    MSS_dic = Dict{Int, Vector{Int64}}()
    qlist = Vector{Int}(undef, 0)
   
    if inv==1 && k==0 || inv==-1 && k==div(N,2)
        for i in eachindex(basisK)
            n = basisK[i]
            nR = get_representative(breflect(n))[1]
            if n <= min(nR, n)
                push!(qlist, length(Set([n, nR])))
            end
            n = min(nR, n)
            if haskey(MSS_dic, n)
                push!(MSS_dic[n], i)
            else
                MSS_dic[n] = [i]
            end
        end

    elseif inv==1 && k==div(N,2) || inv==-1 && k==0
        for i in eachindex(basisK)
            n = basisK[i]
            nR = get_representative(breflect(n))[1]
            if n != nR
                n = min(nR, n)
                if haskey(MSS_dic, n)
                    push!(MSS_dic[n], i)
                else
                    MSS_dic[n] = [i]
                end
                push!(qlist, 2)
            end
        end    
    end


    total_state = zeros(ET, length(basisK))
    MSS_dic=sort(MSS_dic)
    for (i, state_index) in enumerate(values(MSS_dic))
        total_state[state_index] .= 1/sqrt(qlist[i])*state[i]
    end

    return total_state
end
mapstate_MSS2K(N::Int, state::Vector{ET}, k::Int64, inv::Int64=1) where {ET} = mapstate_MSS2K(BitStr{N, Int}, state, k, inv)

mapstate_MSS2total(N::Int64, state::Vector{ET}, k::Int64, inv::Int64=1) where {ET} = mapstate_K2total(N, mapstate_MSS2K(N, state, k, inv), k)

"""
    iso_total2MSS(::Type{T}, k::Int64, inv::Int64=1) where {N, T <: BitStr{N}}

Construct the isometry mapping from total basis to maximum symmetry subspace (MSS).

This is the combination of the total-to-K isometry and the K-to-MSS isometry.

# Arguments
- `T::Type{BitStr{N}}` or `N::Int`: System size specification
- `k::Int64`: Momentum quantum number (0 ≤ k ≤ N-1)
- `inv::Int64`: Inversion flag (-1 or 1)

# Returns
- `Matrix{Float64}`: Isometry matrix mapping total space to MSS

# Example
```julia
W_total2mss = iso_total2MSS(10, 0, 1)  # Map total space to MSS with inversion
MSS_Ham = W_total2mss' * PXP_Ham(BitStr{10, Int}) * W_total2mss  # Project Hamiltonian to MSS subspace
```
"""
function iso_total2MSS(::Type{T}, k::Int64, inv::Int64=1) where {N, T <: BitStr{N}}
    # Function to map the total basis to the MSS space basis, k can only equal to 0 or N/2(pi)
    iso = iso_total2K(T, k) * iso_K2MSS(T, k, inv)

    return iso
end
iso_total2MSS(N::Int, k::Int64, inv::Int64=1) = iso_total2MSS(BitStr{N, Int}, k, inv)

"""
    rdm_PXP_MSS(::Type{T}, subsystems::Vector{Int64}, mssstate::Vector{ET}, k::Int64, inv::Int64=1) where {N,T <: BitStr{N}, ET}
    rdm_PXP_MSS(N::Int64, subsystems::Vector{Int64}, state::Vector{ET}, k::Int64, inv::Int64=1) where {ET}

Compute reduced density matrix for a state in maximum symmetry subspace (MSS).

Transforms the state from MSS to total space and then computes
the reduced density matrix for the specified subsystem.

# Arguments
- `T::Type{BitStr{N}}` or `N::Int`: System size specification
- `subsystems::Vector{Int64}`: Subsystem sites indices
- `mssstate::Vector{ET}`: State vector in MSS
- `k::Int64`: Momentum quantum number
- `inv::Int64`: Inversion flag (-1 or 1)

# Returns
- `Matrix{ET}`: Reduced density matrix of the subsystem

# Example
```julia
rdm = rdm_PXP_MSS(10, [1,2,3], mss_state, 0, 1)
```
"""
function rdm_PXP_MSS(::Type{T}, subsystems::Vector{Int64}, mssstate::Vector{ET}, k::Int64, inv::Int64=1) where {N,T <: BitStr{N}, ET}
    @assert length(PXP_MSS_basis(T, k, inv)[1]) == length(mssstate) "state length is expected to be $(length(PXP_MSS_basis(T, k, inv)[1])), but got $(length(mssstate))"
    state=mapstate_MSS2total(N, mssstate, k, inv)
    reduced_dm = rdm_PXP(T, subsystems, state)
    return reduced_dm
end
rdm_PXP_MSS(N::Int64, subsystems::Vector{Int64},state::Vector{ET}, k::Int64, inv::Int64=1) where {ET} = rdm_PXP_MSS(BitStr{N, Int}, subsystems, state, k, inv)


function cyclebits(state::T) where {N, T <: BitStr{N}}
    #params: t is an integer, N is the length of the binary string
    #We also use this order: system size, state, circular shift bitstring 1 bit.
    # In case need to shift more than 1 bit, we can use a loop or recursion. or we leave a interface here  n_translations::Int
    mask = 1 << N - 1
    return ((state << 1) | (state >> (N - 1))) & mask
end

function get_representative(state::T) where {N, T <: BitStr{N}}
#Finds representative and representative translation for a state.
#State should be a decimal integer.

    representative = state
    translation = 0
    # cycle(bits) = (bits << 1) % (2^N - 1)  # Left shift the state by one position
    current = state
    for n_translation_sites in 1:N-1
        current = cyclebits(current)  # Cycle the bits
        if current < representative
            representative = current
            translation = n_translation_sites
        end
    end

    return representative, translation
end

function PXP_K_basis(::Type{T}, k::Int64) where {N, T <: BitStr{N}}
#params: a int of lattice number and momentum of system
#return: computational basis in given momentum kinetically constrained subspace with decimal int form in PXP model
    @assert 0<=k<=N-1 "k is expected to be in [0, $(N-1)], but got $k"

    basisK = Vector{T}(undef, 0)
    basis = PXP_basis(T)


    basis_dic = Dict{T, Vector{T}}()
    for i in basis
        category = get_representative(i)[1]
        if haskey(basis_dic, category)
            push!(basis_dic[category], i)
        else
            basis_dic[category] = [i]
        end
    end

    for j in eachindex(basis)
        n=basis[j]
        RS = get_representative(n)[1]
        if RS == n && (k * length(basis_dic[RS])) % N == 0
            push!(basisK, n)
        end
    end

    return basisK, basis_dic
end
PXP_K_basis(N::Int, k::Int64) = PXP_K_basis(BitStr{N, Int}, k)

function PXP_MSS_basis(::Type{T}, k::Int64, inv::Int64=1) where {N, T <: BitStr{N}}
#params: a int of lattice number and momentum of system, we have considered the inversion symmetry
#return: computational basis in given momentum inversion symmetry subspace with decimal int form
    @assert k == 0 || k==div(N,2) "k is expected to be 0 or $(div(N,2)), but got $k"
    @assert inv ==1 || inv==-1 "inv is expected to be 1 or -1, but got $(inv)"
    # MSS is the list of states in the maximum symmetry sector
    MSS = Vector{T}(undef, 0)
    basisK, basis_dic = PXP_K_basis(T, k)
    MSS_dic = Dict{T, Vector{T}}()

    # q is the number of states that are equivalent under inversion
    qlist = Vector{Int}(undef, 0)
    for i in eachindex(basisK)
        n = basisK[i]
        # here we calculate the representative state of the inversion of n
        nR, d = get_representative(breflect(n))
        q = length(Set([n, nR]))
        if n <= min(nR, n)
            if k == 0
                # At k=0, q=1 states are always even under inversion
                if inv == 1
                    push!(MSS, n)
                    MSS_dic[n] = basis_dic[n]
                    push!(qlist, q)
                elseif inv == -1 && q == 2
                    push!(MSS, n)
                    MSS_dic[n] = basis_dic[n]
                    push!(qlist, q)
                end
    # For k=0, the inversion invariant states (n=nR) belong to the +1 sector, while the states that are not invariant (n!=nR) form pairs that contribute to both +1 and -1 sectors. 
    # For k=π, the situation is reversed: the invariant states belong to the -1 sector, and the non-invariant pairs contribute to both sectors. Therefore, we need to filter MSS based on the inversion eigenvalue inv and update MSS_dic accordingly.
            else  # k == div(N, 2), i.e. k=π
                # At k=π, q=1 states have inversion eigenvalue (-1)^d
                if q == 1
                    parity = d % 2 == 0 ? 1 : -1
                    if parity == inv
                        push!(MSS, n)
                        MSS_dic[n] = basis_dic[n]
                        push!(qlist, q)
                    end
                else
                    push!(MSS, n)
                    MSS_dic[n] = basis_dic[n]
                    push!(qlist, q)
                end
            end
        end
    end

    return MSS, MSS_dic, qlist
end
PXP_MSS_basis(N::Int, k::Int64, inv::Int64=1) = PXP_MSS_basis(BitStr{N, Int}, k, inv)

function PXP_K_Ham(::Type{T}, k::Int, Omega::Float64=1.0) where {N, T <: BitStr{N}}
#params: a int of lattice number, momentum of system and interaction strength of system which default to be 1
#return: the Hamiltonian matrix in given K space

    @assert 0<=k<=N-1 "k is expected to be in [0, $(N-1)], but got $k"

    basisK, basis_dic = PXP_K_basis(T, k)
    l = length(basisK)
    omegak = exp(2im * π * k / N)
    H = zeros(ComplexF64, (l, l))

    for i in 1:l
        n=basisK[i]
        output = actingH_PXP(T, n, true)
        for m in output
            mbar, d = get_representative(m)
            if mbar ∈ basisK
                j=searchsortedfirst(basisK, mbar)
                Yn= sqrt(length(basis_dic[n])) / N
                Ym= sqrt(length(basis_dic[mbar])) / N
                H[i, j] += Yn/Ym * omegak^d
            end
        end
    end
    if k==0 || k==div(N,2)
        H=real(H)
    end
    H=(H+H')/2
    return H
end
PXP_K_Ham(N::Int, k::Int, Omega::Float64=1.0) = PXP_K_Ham(BitStr{N, Int}, k, Omega)

function PXP_MSS_Ham(::Type{T}, k::Int, inv::Int64=1) where {N, T <: BitStr{N}}
#params: a int of lattice number, momentum of system and interaction strength of system which default to be 1
#return: the Hamiltonian matrix in given maximum symmetry space 
    @assert k == 0 || k==div(N,2) "k is expected to be 0 or $(div(N,2)), but got $k"
    @assert inv ==1 || inv==-1 "inv is expected to be 1 or -1, but got $(inv)"

    omegak = exp(2im * π * k / N)
    
    MSS, MSS_dic, qlist = PXP_MSS_basis(T, k, inv)
    l = length(MSS)
    H = zeros(ComplexF64, (l, l))

    if inv==1 && k==0 || inv==-1 && k==div(N,2)
        for i in 1:l
            n = MSS[i]
            Zn = sqrt(qlist[i]) * sqrt(length(MSS_dic[n]))
            # Zn is the normalization factor for the state n in MSS, which is the square root of the number of states in K space that are equivalent to n under inversion, multiplied by the square root of the number of states in MSS that are equivalent to n under inversion. This is because when we map from K space to MSS space, we need to sum over all the states in K space that are equivalent to n under inversion, and each state in K space has a normalization factor of 1/sqrt(length(basis_dic[n])), and there are length(MSS_dic[n]) states in MSS that are equivalent to n under inversion, so we need to multiply by sqrt(length(MSS_dic[n])) to get the correct normalization factor for the state n in MSS.
            output = actingH_PXP(T, n, true)
            for m in output
                mbar, d = get_representative(m)
                inv_mbar = get_representative(breflect(mbar))[1]
                mtilde = min(mbar, inv_mbar)
                if mtilde ∈ MSS
                    j = searchsortedfirst(MSS, mtilde)
                    Zm = sqrt(qlist[j]) * sqrt(length(MSS_dic[mtilde])) 
                    H[i, j] +=  Zn / Zm*omegak^d
                end
            end
        end

    elseif inv==-1 && k==0 || inv==1 && k==div(N,2)
        for i in 1:l
            n = MSS[i]
            Zn = sqrt(length(MSS_dic[n]))
            output = actingH_PXP(T, n, true)
            for m in output
                mbar, d = get_representative(m)
                inv_mbar = get_representative(breflect(mbar))[1]
                mtilde = min(mbar, inv_mbar)
                @show m.buf, mbar.buf, inv_mbar.buf, mtilde.buf
                if mtilde ∈ MSS
                    j=searchsortedfirst(MSS, mtilde)
                    Zm = sqrt(length(MSS_dic[mtilde]))
                    H[i, j] +=  Zn / Zm*omegak^d
                end
            end
        end
    end
    
    H=real(H)
    H = (H + H') / 2  
    
    return H
end
PXP_MSS_Ham(N::Int, k::Int, inv::Int64=1) = PXP_MSS_Ham(BitStr{N, Int}, k, inv)
