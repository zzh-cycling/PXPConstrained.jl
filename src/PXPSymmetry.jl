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
    # Function to map K basis to MSS basis.
    # k can only be 0 or N/2 (pi), inv can only be +1 or -1.
    @assert k == 0 || k == div(N, 2) "k is expected to be 0 or $(div(N,2)),but got $k"
    @assert inv == 1 || inv == -1 "inv is expected to be 1 or -1, but got $inv"

    basisK, _ = PXP_K_basis(T, k)
    nK = length(basisK)
    basisK_dic = Dict{T, Int}()
    for i in eachindex(basisK)
        basisK_dic[basisK[i]] = i
    end

    omegak = k == 0 ? 1 : -1
    rep_dic = Dict{T, Vector{Int64}}()
    phase_dic = Dict{T, Int64}()

    for i in eachindex(basisK)
        n = basisK[i]
        nR, d = get_representative(breflect(n))
        rep = min(n, nR)
        haskey(rep_dic, rep) && continue

        if n == nR
            phase = omegak^d
            if phase == inv
                rep_dic[rep] = [i]
                phase_dic[rep] = phase
            end
        else
            repR, drep = get_representative(breflect(rep))
            idx_rep = basisK_dic[rep]
            idx_repR = basisK_dic[repR]
            rep_dic[rep] = [idx_rep, idx_repR]
            phase_dic[rep] = omegak^drep
        end
    end

    reps = sort(collect(keys(rep_dic)))
    nMSS = length(reps)
    iso = zeros(nK, nMSS)

    for (col, rep) in enumerate(reps)
        indices = rep_dic[rep]
        phase = phase_dic[rep]

        if length(indices) == 1
            iso[indices[1], col] = 1.0
        else
            @assert length(indices) == 2 "Unexpected number of K-indices for rep $rep: $(length(indices))"
            idx1 = indices[1]
            idx2 = indices[2]
            iso[idx1, col] = 1 / sqrt(2)
            iso[idx2, col] = inv * phase / sqrt(2)
        end
    end

    return iso
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

    basisK, basis_dic = PXP_K_basis(T, k)
    nK = length(basisK)
    basisK_dic = Dict{T, Int}()
    for i in eachindex(basisK)
        basisK_dic[basisK[i]] = i
    end

    omegak = k == 0 ? 1 : -1
    rep_dic = Dict{T, Vector{Int64}}()
    phase_dic = Dict{T, Int64}()

    for i in eachindex(basisK)
        n = basisK[i]
        nR, d = get_representative(breflect(n))
        rep = min(n, nR)
        haskey(rep_dic, rep) && continue

        if n == nR
            phase = omegak^d
            if phase == inv
                rep_dic[rep] = [i]
                phase_dic[rep] = phase
            end
        else
            repR, drep = get_representative(breflect(rep))
            idx_rep = basisK_dic[rep]
            idx_repR = basisK_dic[repR]
            rep_dic[rep] = [idx_rep, idx_repR]
            phase_dic[rep] = omegak^drep
        end
    end

    reps = sort(collect(keys(rep_dic)))
    nMSS = length(reps)
    @assert length(state) == nMSS "state length is expected to be $nMSS, but got $(length(state))"

    total_state = zeros(ET, nK)
    for (col, rep) in enumerate(reps)
        indices = rep_dic[rep]
        phase = phase_dic[rep]
        if length(indices) == 1
            total_state[indices[1]] = state[col]
        else
            idx1 = indices[1]
            idx2 = indices[2]
            total_state[idx1] = state[col] / sqrt(2)
            total_state[idx2] = inv * phase * state[col] / sqrt(2)
        end
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
        q = length(Set([n, nR])) # whether is self-conjugate under inversion, q=1 means self-conjugate, q=2 means not self-conjugate
        if n <= min(nR, n) # only consider one representative from each inversion pair, and for self-conjugate states, we only consider it once
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
    # Here d is the translation distance that relates n and its inversion nR. For k=π, the inversion eigenvalue of a self-conjugate state is determined by whether this translation distance is even or odd. If d is even, the state is even under inversion (parity=1); if d is odd, the state is odd under inversion (parity=-1). Therefore, we compare parity with inv to decide whether to include this state in MSS.
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

    basisK, basis_dic = PXP_K_basis(T, k)
    nK = length(basisK)
    η = k == 0 ? 1.0 : -1.0

    iso = iso_K2MSS(T, k, inv)
    nMSS = size(iso, 2)

    K2MSS_col = fill(0, nK)
    K2MSS_coef = zeros(Float64, nK)
    for col in 1:nMSS
        for row in 1:nK
            coeff = iso[row, col]
            if coeff != 0.0
                K2MSS_col[row] = col
                K2MSS_coef[row] = coeff
            end
        end
    end

    k_index = Dict{T, Int}()
    for (i, n) in enumerate(basisK)
        k_index[n] = i
    end

    H = zeros(Float64, nMSS, nMSS)
    for ket_idx in 1:nK
        ket_col = K2MSS_col[ket_idx]
        ket_col == 0 && continue

        ket_coeff = K2MSS_coef[ket_idx]
        ket_rep = basisK[ket_idx]
        Yn = sqrt(length(basis_dic[ket_rep])) / N

        for out_state in actingH_PXP(T, ket_rep, true)
            bra_rep, d = get_representative(out_state)
            bra_idx = get(k_index, bra_rep, 0)
            bra_idx == 0 && continue

            bra_col = K2MSS_col[bra_idx]
            bra_col == 0 && continue

            Ym = sqrt(length(basis_dic[bra_rep])) / N
            hij = Yn / Ym * (η^d)
            H[bra_col, ket_col] += K2MSS_coef[bra_idx] * hij * ket_coeff
        end
    end

    return (H + H') / 2
end
PXP_MSS_Ham(N::Int, k::Int, inv::Int64=1) = PXP_MSS_Ham(BitStr{N, Int}, k, inv)
