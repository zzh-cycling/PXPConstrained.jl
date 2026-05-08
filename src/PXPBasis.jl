"""
    PXPBasis.jl

Functions for generating basis states and Hamiltonians for the PXP model.
This module handles the construction of constrained Fibonacci chains that satisfy
the Rydberg blockade condition and provides utilities for reduced density matrices.
"""

"""
    Fibonacci_chain_OBC(::Type{T}) where {N, T <: BitStr{N}}

Generate the constrained basis for the PXP model with open boundary conditions.

Creates all valid bit configurations where no two adjacent sites can both be excited,
following the Fibonacci sequence growth pattern.

# Arguments
- `T::Type{BitStr{N}}`: Bit string type specifying system size N

# Returns
- `Vector{T}`: All valid basis states satisfying the constraint
"""
function Fibonacci_chain_OBC(::Type{T}) where {N, T <: BitStr{N}}
    # Generate Fibonacci chain for PXP model with open boundary condition
    fib_chain=[[T(0), T(1)],[T(0), T(1), T(2)]]
    for i in 3:N
        push!(fib_chain,vcat([s << 1 for s in fib_chain[i-1]],[(s << 2 | T(1)) for s in fib_chain[i-2]]))
    end
    # each push we add a bit"0" or bit"01" to the end of the bit_string, and finally return the last element of the fib_chain
    return fib_chain[N]
end

"""
    Fibonacci_chain_PBC(::Type{T}) where {N, T <: BitStr{N}}

Generate the constrained basis for the PXP model with periodic boundary conditions.

Creates all valid bit configurations for a ring geometry where the first and last
sites are also neighbors, adding additional constraint checks.

# Arguments
- `T::Type{BitStr{N}}`: Bit string type specifying system size N

# Returns
- `Vector{T}`: All valid basis states satisfying PBC constraints
"""
function Fibonacci_chain_PBC(::Type{T}) where {N, T <: BitStr{N}}
    # Generate Fibonacci chain  for PXP model with periodic boundary condition
    return filter(c -> iszero((c >> (N-1)) & (c & 1)), Fibonacci_chain_OBC(T))
end

"""
    actingH_PXP(::Type{T}, state::T, pbc::Bool=true) where {N, T <: BitStr{N}}

Apply the PXP Hamiltonian to a given basis state.

The PXP Hamiltonian flips spins at sites where both neighbors are in the ground state,
respecting the Rydberg blockade constraint. Returns all possible output states.

# Arguments
- `T::Type{BitStr{N}}`: Bit string type specifying system size N
- `state::T`: Input basis state to act upon
- `pbc::Bool=true`: Whether to use periodic boundary conditions

# Returns
- `Vector{T}`: All output states after applying the Hamiltonian

# Example
```julia
outputs = actingH_PXP(BitStr{6, Int}, BitStr{6}(5), true)
```
"""
function actingH_PXP(::Type{T}, state::T, pbc::Bool=true) where {N, T <: BitStr{N}}
    # The type of n is DitStr{D, N, Int}, which is a binary string with length N in D-ary form.
    # Acting Hamiltonian on a given state in bitstr and return the output states in bitstr
    # Here need to note that the order of the bitstr is from right to left, which is different from the normal order.
    mask=bmask(T, N, N-2)
    fl=bmask(T, N-1)
    # output = [flip(state, fl >> (i-1)) for i in 1:N-2 if state & (mask >> (i-1)) == 0] 
    output = map(i -> flip(state, fl >> (i-1)), filter(i -> state & (mask >> (i-1)) == 0, 1:N-2)) # faster

    if pbc
        if state[2]==0 && state[N]==0
            flip_str=flip(state,bmask(T, 1))
            push!(output,flip_str)
        end
        if state[1]==0 && state[N-1]==0
            flip_str=flip(state,bmask(T, N))
            push!(output,flip_str)
        end
    else
        if state[N-1]==0
            flip_str=flip(state,bmask(T, N))
            push!(output,flip_str)
        end
        if state[2]==0
            flip_str=flip(state,bmask(T, 1))
            push!(output,flip_str)
        end
    end
    return output
end

"""
    PXP_basis(::Type{T}, pbc::Bool=true) where {N, T <: BitStr{N}}
    PXP_basis(N::Int, pbc::Bool=true)

Generate the complete basis for the PXP model.

Creates all valid basis states that satisfy the Rydberg blockade constraint, sorted in ascending order for efficient searching.

# Arguments
- `T::Type{BitStr{N}}` or `N::Int`: System size specification
- `pbc::Bool=true`: Whether to use periodic boundary conditions

# Returns
- `Vector{BitStr{N}}`: Sorted list of all valid basis states

# Example
```julia
basis = PXP_basis(8, true)  # 8 sites with PBC
basis_obc = PXP_basis(8, false)  # 8 sites with OBC
```
"""
function PXP_basis(::Type{T},pbc::Bool=true) where {N, T <: BitStr{N}}
    # Generate basis for PXP model, return both decimal and binary form, where we both consider PBC and OBC
    if pbc
        basis=Fibonacci_chain_PBC(T)
    else
        basis=Fibonacci_chain_OBC(T)
    end
    sorted_basis=sort(basis)
    return sorted_basis
end
PXP_basis(N::Int, pbc::Bool=true) = PXP_basis(BitStr{N, Int}, pbc)

"""
    PXP_Ham(::Type{T}, pbc::Bool=true) where {N, T <: BitStr{N}}
    PXP_Ham(N::Int, pbc::Bool=true)

Construct the full Hamiltonian matrix for the PXP model.

Builds the Hamiltonian matrix in the constrained basis by applying the PXP
operator to each basis state and recording the matrix elements.

# Arguments
- `T::Type{BitStr{N}}` or `N::Int`: System size specification
- `pbc::Bool=true`: Whether to use periodic boundary conditions

# Returns
- `Matrix{Float64}`: The PXP Hamiltonian matrix

# Example
```julia
H = PXP_Ham(8, true)  # 8-site PXP Hamiltonian with PBC
eigenvals, eigenvecs = eigen(H)
```
"""
function PXP_Ham(::Type{T}, pbc::Bool=true) where {N, T <: BitStr{N}}
    # Generate Hamiltonian for PXP model, automotically contain pbc or obc
    basis=PXP_basis(T,pbc)

    l=length(basis)
    H=zeros(Float64,(l,l))
    for i in 1:l
        output=actingH_PXP(T, basis[i], pbc) 
        for m in output 
            j=searchsortedfirst(basis,m)
            H[i, j] += 1
        end
    end

    return H
end
PXP_Ham(N::Int, pbc::Bool=true) = PXP_Ham(BitStr{N, Int}, pbc)

"""
    process_join(a, b)

Join two lists by computing their Cartesian product.

Helper function for creating composite basis states from multiple subsystems.

# Arguments
- `a`, `b`: Lists to be joined

# Returns
- `Vector`: Vectorized Cartesian product of the two lists
"""
# join two lists of basis by make a product of two lists
function process_join(a, b)
    return vec([join(b, a) for a in a, b in b])
end

"""
    joint_pxp_basis(lengthlis::Vector{Int})

Create PXP basis for multiple disjoint sub-chains.

Generates the basis for a system composed of multiple disconnected chains,
each with open boundary conditions.

# Arguments
- `lengthlis::Vector{Int}`: Lengths of each sub-chain

# Returns
- `Vector`: Combined and sorted basis for the composite system

"""
# create pxp basis composed of multiple disjoint sub-chains
function joint_pxp_basis(lengthlis::Vector{Int})
    return sort(mapreduce(len -> PXP_basis(len, false), process_join, lengthlis))
end

"""
    connected_components(v::Vector{Int})

Find connected components in a sorted list of integers.

Groups consecutive integers into separate segments, useful for identifying
contiguous subsystems in quantum many-body calculations.

# Arguments
- `v::Vector{Int}`: Sorted vector of integers

# Returns
- `Vector{Vector{Int}}`: List of connected components (consecutive segments)

# Example
```julia
components = connected_components([1, 2, 3, 5, 6, 8])
# Returns: [[1, 2, 3], [5, 6], [8]]
```
"""
function connected_components(v::Vector{Int})
    if isempty(v)
        return []
    end

    sort!(v)

    result = []
    current_segment = [v[1]]

    for i in 2:length(v)
        if v[i] == v[i - 1] + 1
            push!(current_segment, v[i])
        else
            push!(result, current_segment)
            current_segment = [v[i]]
        end
    end

    push!(result, current_segment)

    return result
end

"""
    move_subsystem(::Type{BitStr{M, INT}}, basis::BitStr{N, INT}, subsystems::Vector{Int}) where {M, N, INT}

Move specified subsystem bits to the left of a larger bit string.

Rearranges bits to place the selected subsystem at the beginning of the bit string,
useful for constructing reduced density matrices.

# Arguments
- `BitStr{M, INT}`: Target bit string type with size M
- `basis::BitStr{N, INT}`: Source bit string with size N
- `subsystems::Vector{Int}`: Indices of subsystem sites

# Returns
- `BitStr{M}`: Rearranged bit string with subsystem bits at the left
"""
function move_subsystem(::Type{BitStr{M, INT}}, basis::BitStr{N, INT}, subsystems::Vector{Int}) where {M, N, INT}
    # Move the subsystem bits to the left of the basis, and return a new basis with the subsystem bits moved to the left, and move into the bigger basis.
    @assert length(subsystems) == N "subsystems length is expected to be $N, but got $(length(subsystems))"
    @assert M >= N "total length is expected to be greater than or equal to $N, but got $M"
    return sum(i -> BitStr{M}(readbit(basis.buf, i) << (M - subsystems[N-i+1])), 1:N)
end

"""
    takeenviron(x, mask::BitStr{l}) where {l}

Extract the environment part of a basis state.

Returns the bits that are not part of the specified subsystem.

# Arguments
- `x`: Full basis state
- `mask::BitStr{l}`: Mask specifying the subsystem bits

# Returns
- Environment part of the basis state
"""
# take environment part of a basis
takeenviron(x, mask::BitStr{l}) where {l} = x & (~mask)

"""
    takesystem(x, mask::BitStr{l}) where {l}

Extract the subsystem part of a basis state.

Returns only the bits that belong to the specified subsystem.

# Arguments
- `x`: Full basis state  
- `mask::BitStr{l}`: Mask specifying the subsystem bits

# Returns
- Subsystem part of the basis state
"""
# take system part of a basis
takesystem(x, mask::BitStr{l}) where {l} = (x & mask)

"""
    rdm_PXP(::Type{T}, subsystems::Vector{Int64}, state::Vector{ET}, pbc::Bool=true) where {N,T <: BitStr{N}, ET}
    rdm_PXP(N::Int, subsystems::Vector{Int64}, state::Vector{ET}, pbc::Bool=true) where {ET}

Compute the reduced density matrix for a subsystem of state in PXP basis.

Traces out the environment degrees of freedom to obtain the reduced density matrix
of the specified subsystem. Handles connected components automatically.

# Arguments
- `T::Type{BitStr{N}}` or `N::Int`: System size specification
- `subsystems::Vector{Int64}`: Indices of sites to include in the subsystem
- `state::Vector{ET}`: Input quantum state vector
- `pbc::Bool=true`: Whether to use periodic boundary conditions

# Returns
- `Matrix{ET}`: Reduced density matrix of the subsystem

# Example
```julia
# Get reduced density matrix for sites 1-4
psi = normalized_eigenstate  # some quantum state
rdm = rdm_PXP(12, collect(1:4), psi, true)
ee_value = ee(rdm)  # compute entanglement entropy
```
"""
function rdm_PXP(::Type{T}, subsystems::Vector{Int64}, state::Vector{ET}, pbc::Bool=true) where {N,T <: BitStr{N}, ET}
    # Usually subsystem indices count from the right of binary string. But this version we can count from the left, which is consistent with our intuition. This means that the systems we want to keep.
    # The function is to take common environment parts of the total basis, get the index of system parts in reduced basis, and then calculate the reduced density matrix.
    unsorted_basis = PXP_basis(T, pbc)
    @assert length(unsorted_basis) == length(state) "state length is expected to be $(length(unsorted_basis)), but got $(length(state))"
    
    subsystems=connected_components(subsystems)
    lengthlis=length.(subsystems)
    subsystems=vcat(subsystems...)
    mask = bmask(T, (N .-subsystems .+1)...)

    
    order = sortperm(unsorted_basis, by = x -> (takeenviron(x, mask), takesystem(x, mask))) #first sort by environment, then by system. The order of environment doesn't matter.
    basis, state = unsorted_basis[order], state[order]
    
    reduced_basis = move_subsystem.(T, joint_pxp_basis(lengthlis), Ref(subsystems))
    len = length(reduced_basis)
    # Initialize the reduced density matrix
    reduced_dm = zeros(ET, (len, len))

    # Keep track of indices where the key changes
    result_indices = Int[]
    current_key = -1
    for (idx, i) in enumerate(basis)
        key = takeenviron(i, mask)  # Get environment l bits
        if key != current_key
            @assert key > current_key "key is expected to be greater than $current_key, but got $key"
            push!(result_indices, idx)
            current_key = key
        end
    end
    # Add the final index to get complete ranges
    push!(result_indices, length(basis) + 1)

    for i in 1:length(result_indices)-1
        range = result_indices[i]:result_indices[i+1]-1         
        # Get indices in the reduced basis
        indices = searchsortedfirst.(Ref(reduced_basis), takesystem.(basis[range], mask))
        view(reduced_dm, indices, indices) .+= view(state, range) .* view(state, range)'
    end

    return reduced_dm
end
rdm_PXP(N::Int, subsystems::Vector{Int64}, state::Vector{ET}, pbc::Bool=true) where {ET} = rdm_PXP(BitStr{N, Int}, subsystems, state, pbc)


"""
    myprint(io::IO, xs...)

Enhanced print function with automatic flushing.

Prints arguments to the specified IO stream with an extra newline and
flushes the buffer to ensure immediate output, useful for real-time monitoring.

# Arguments
- `io::IO`: Output stream
- `xs...`: Arguments to print

# Example
```julia
myprint(stdout, "Results:", value1, value2)
```
"""
function myprint(io::IO, xs...)
    println(io, xs..., '\n')
    flush(io)
end

"""
    iso_full2cons(::Type{T}, pbc::Bool) where {N, T <: BitStr{N}}
    
    Construct the isometry matrix from the full Hilbert space to the constrained PXP basis.
    This function creates a mapping matrix that projects states from the full 2^N-dimensional Hilbert space onto the subspace defined by the PXP constraints, where only valid configurations are retained.
# Arguments
- `T::Type{BitStr{N}}`: Bit string type specifying system size N
- `pbc::Bool`: Whether to use periodic boundary conditions for generating the constrained basis
# Returns
- `Matrix{Int64}`: Isometry matrix mapping full Hilbert space to constrained PXP basis
# Example
```julia
iso_matrix = iso_full2cons(BitStr{6, Int}, true)
```
"""
function iso_full2cons(::Type{T}, pbc::Bool) where {N, T <: BitStr{N}}
    basis = PXP_basis(T, pbc)
    l = length(basis)
    map_matrix = zeros(Int64, (2^N, l))
    for i in 1:l
        base = basis[i]
        map_matrix[base.buf + 1, i] = 1
    end
    return map_matrix
end
iso_full2cons(N::Int, pbc::Bool) = iso_full2cons(BitStr{N, Int}, pbc)