#This PXP_basis package is used to generate the basis for PXP model, and also the Hamiltonian matrix for PXP model according to their basis.  
#Totally we have three kinds of functions, basis, Ham, reduced_dm. We consider both PBC and OBC in the basis and Hamiltonian matrix, and we also consider the translational symmetry and inversion symmetry, and FSA subspace.


function Fibonacci_chain_OBC(::Type{T}) where {N, T <: BitStr{N}}
    # Generate Fibonacci chain for PXP model with open boundary condition
    fib_chain=[[T(0), T(1)],[T(0), T(1), T(2)]]
    for i in 3:N
        push!(fib_chain,vcat([s << 1 for s in fib_chain[i-1]],[(s << 2 | T(1)) for s in fib_chain[i-2]]))
    end
    # each push we add a bit"0" or bit"01" to the end of the bit_string, and finally return the last element of the fib_chain
    return fib_chain[N]
end

function Fibonacci_chain_PBC(::Type{T}) where {N, T <: BitStr{N}}
    # Generate Fibonacci chain  for PXP model with periodic boundary condition
    return filter(c -> iszero((c >> (N-1)) & (c & 1)), Fibonacci_chain_OBC(T))
end

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

# join two lists of basis by make a product of two lists
function process_join(a, b)
    return vec([join(b, a) for a in a, b in b])
end

# create pxp basis composed of multiple disjoint sub-chains
function joint_pxp_basis(lengthlis::Vector{Int})
    return sort(mapreduce(len -> PXP_basis(len, false), process_join, lengthlis))
end

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

function move_subsystem(::Type{BitStr{M, INT}}, basis::BitStr{N, INT}, subsystems::Vector{Int}) where {M, N, INT}
    # Move the subsystem bits to the left of the basis, and return a new basis with the subsystem bits moved to the left, and move into the bigger basis.
    @assert length(subsystems) == N "subsystems length is expected to be $N, but got $(length(subsystems))"
    @assert M >= N "total length is expected to be greater than or equal to $N, but got $M"
    return sum(i -> BitStr{M}(readbit(basis.buf, i) << (M - subsystems[N-i+1])), 1:N)
end

# take environment part of a basis
takeenviron(x, mask::BitStr{l}) where {l} = x & (~mask)
# take system part of a basis
takesystem(x, mask::BitStr{l}) where {l} = (x & mask)

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


function myprint(io::IO, xs...)
    println(io, xs..., '\n')
    flush(io)
end

