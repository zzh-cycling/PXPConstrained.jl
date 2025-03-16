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
    output = [flip(state, fl >> (i-1)) for i in 1:N-2 if state & (mask >> (i-1)) == 0] # faster method

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

function find_indices(reduced_basis, basis, mask)
    indices = Int[]           # 用于存储 basis 的索引
    bit_vector = falses(length(basis)) 

    for r in eachindex(basis)
        value = takesystem(basis[r], mask)
        index = searchsortedfirst(reduced_basis, value)

        if index <= length(reduced_basis) && reduced_basis[index] == value
            push!(indices, index) 
            bit_vector[r] = true    
        end
    end

    return indices, bit_vector 
end

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
    @assert length(subsystems) == N "subsystems length is expected to be $N, but got $(length(subsystems))"
    @assert M >= N "total length is expected to be greater than or equal to $N, but got $M"
    return sum(i -> BitStr{M}(readbit(basis.buf, i) << (subsystems[i] - 1)), 1:N)
end

# take environment part of a basis
takeenviron(x, mask::BitStr{l}) where {l} = x & (~mask)
# take system part of a basis
takesystem(x, mask::BitStr{l}) where {l} = (x & mask)

function rdm_PXP(::Type{T}, subsystems::Vector{Int64}, state::Vector{ET}, pbc::Bool=true) where {N,T <: BitStr{N}, ET}
    # Usually subsystem indices count from the right of binary string.
    # The function is to take common environment parts of the total basis, get the index of system parts in reduced basis, and then calculate the reduced density matrix.
    subsystems=connected_components(subsystems)
    lengthlis=length.(subsystems)
    subsystems=vcat(subsystems...)
    mask = bmask(T, subsystems...)

    unsorted_basis = PXP_basis(T, pbc)
    order = sortperm(unsorted_basis, by = x -> (takeenviron(x, mask), takesystem(x, mask))) #first sort by environment, then by system. The order of environment doesn't matter.
    basis, state = unsorted_basis[order], state[order]
    @assert length(basis) == length(state) "state length is expected to be $(length(basis)), but got $(length(state))"
    
    reduced_basis = move_subsystem.(T, joint_pxp_basis(lengthlis), Ref(subsystems))
    len = length(reduced_basis)
    # Initialize the reduced density matrix
    reduced_dm = zeros(ET, (len, len))

    # Keep track of indices where the key changes
    result_indices = Int[]
    current_key = -1
    for (idx, i) in enumerate(basis)
        key = takeenviron(i, mask)  # Get system l bits
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
        indices, _ = find_indices(reduced_basis, basis[range], mask)
        # s = view(state, range[exists])
        view(reduced_dm, indices, indices) .+= view(state, range) .* view(state, range)'
    end

    return reduced_dm
end
rdm_PXP(N::Int, subsystems::Vector{Int64}, state::Vector{ET}, pbc::Bool=true) where {ET} = rdm_PXP(BitStr{N, Int}, subsystems, state, pbc)

function iso_total2K(::Type{T}, k::Int64) where {N, T <: BitStr{N}}
#Function to map the total basis to the K space basis, actually is the isometry, defined as W'*W=I, W*W'=P, P^2=P

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

function rdm_PXP_K(::Type{T}, subsystems::Vector{Int64},state::Vector{ET}, k::Int64) where {N,T <: BitStr{N}, ET}
    iso = iso_total2K(T, k)
    @assert size(iso)[2] == length(state) "state length is expected to be $(size(iso)[2]), but got $(length(state))"
    state=iso*state
    reduced_dm = rdm_PXP(T, subsystems, state)
    return reduced_dm
end
rdm_PXP_K(N::Int, subsystems::Vector{Int64},state::Vector{ET}, k::Int64) where {ET} = rdm_PXP_K(BitStr{N, Int}, subsystems, state, k)


function iso_K2MSS(::Type{T}, k::Int64, inv::Int64=1) where {N, T <: BitStr{N}}
#Function to map the MSS basis to the K space basis

    basis = PXP_basis(T)
    basisK, k_dic = PXP_K_basis(T, k)

    MSS_dic = Dict{Int, Vector{Int64}}()
    # MSS_dic is a dictionary, the key is the representative state of the inversion of n, and the value is the index of the state in the basisK. NOTE that MSS_dic is not sorted, so we need to sort it later.
    qlist = Vector{Int}(undef, 0)
    # Below procedure is to collapse the extra basis in K space that can be converted mutually to MSS space.
    if inv==1
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
iso_K2MSS(N::Int, k::Int64, inv::Int64=1) = iso_K2MSS(BitStr{N, Int}, k, inv)

function iso_total2MSS(::Type{T}, k::Int64, inv::Int64=1) where {N, T <: BitStr{N}}
    # Function to map the total basis to the MSS space basis, k can only equal to 0 or N/2(pi)
    iso = iso_total2K(T, k) * iso_K2MSS(T, k, inv)

    return iso
end
iso_total2MSS(N::Int, k::Int64, inv::Int64=1) = iso_total2MSS(BitStr{N, Int}, k, inv)

function rdm_PXP_MSS(::Type{T}, subsystems::Vector{Int64}, state::Vector{ET}, k::Int64, inv::Int64=1) where {N,T <: BitStr{N}, ET}
    iso = iso_total2MSS(T, k, inv)
    @assert size(iso)[2] == length(state) "state length is expected to be $(size(iso)[2]), but got $(length(state))"
    state=iso*state
    reduced_dm = rdm_PXP(T, subsystems, state)
    return reduced_dm
end
rdm_PXP_MSS(N::Int64, subsystems::Vector{Int64},state::Vector{ET}, k::Int64, inv::Int64=1) where {ET} = rdm_PXP_MSS(BitStr{N, Int}, subsystems, state, k, inv)


function inversion_matrix(::Type{T}) where {N, T <: BitStr{N}}
    basis=PXP_basis(T)
    l=length(basis)
    Imatrix=zeros((l,l))
    reversed_basis=similar(basis)
    for i in eachindex(basis)
        reversed_basis[i]=breflect(basis[i])
    end
    for i in eachindex(basis)
        output=reversed_basis[i]
        j=searchsortedfirst(basis,output)
        Imatrix[i,j]+=1.0
    end
   
    return Imatrix
end
inversion_matrix(N::Int) = inversion_matrix(BitStr{N, Int})

function cyclebits(state::T, n_translations::Int) where {N, T <: BitStr{N}}
#params: t is an integer, N is the length of the binary string
#n_translations: number of positions to shift
#return: the binary left shift
#We also use this order: system size, state, number of translations

    return (state << n_translations) % (2^N - 1)
end

function get_representative(state::T) where {N, T <: BitStr{N}}
#Finds representative and representative translation for a state.
#State should be a decimal integer.

    representative = state
    translation = 0
    for n_translation_sites in 0:N-1
        new_state = cyclebits(state, n_translation_sites)
        if new_state < representative
            representative = new_state
            translation = n_translation_sites
        end
    end
    return representative, translation
end

function PXP_K_basis(::Type{T}, k::Int64) where {N, T <: BitStr{N}}
#params: a int of lattice number and momentum of system
#return: computational basis in given momentum kinetically constrained subspace with decimal int form in PXP model

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


function PXP_MSS_basis(::Type{T}, k::Int64,inv::Int64=1) where {N, T <: BitStr{N}}
#params: a int of lattice number and momentum of system, we have considered the inversion symmetry
#return: computational basis in given momentum inversion symmetry subspace with decimal int form

    # MSS is the list of states in the maximum symmetry sector
    MSS = Vector{T}(undef, 0)
    basisK, basis_dic = PXP_K_basis(T, k)
    MSS_dic = Dict{T, Vector{T}}()

    # q is the number of states that are equivalent under inversion
    qlist = Vector{Int}(undef, 0)
    if inv==1
        for i in eachindex(basisK)
            n = basisK[i]
            # here we calculate the representative state of the inversion of n
            nR = get_representative(breflect(n))[1]
            if n <= min(nR, n)
                push!(MSS, n)
                MSS_dic[n] = basis_dic[n]
                push!(qlist, length(Set([n, nR])))
            end
        end

        return MSS, MSS_dic, qlist

    else
        for i in eachindex(basisK)
            n = basisK[i]
            nR = get_representative(breflect(n))[1]
            if n <= min(nR, n)
                push!(MSS, n)
                MSS_dic[n] = basis_dic[n]
                push!(qlist, length(Set([n, nR])))
            end
        end    
        index=findall(x -> x==2, qlist)
        new_MSS_dic = Dict(k => v for k in MSS[index] for v in [MSS_dic[k]])
        return MSS[index], new_MSS_dic
    end
    
end
PXP_MSS_basis(N::Int, k::Int64, inv::Int64=1) = PXP_MSS_basis(BitStr{N, Int}, k, inv)

function PXP_K_Ham(::Type{T}, k::Int, Omega::Float64=1.0) where {N, T <: BitStr{N}}
#params: a int of lattice number, momentum of system and interaction strength of system which default to be 1
#return: the Hamiltonian matrix in given K space

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
    H=real(H)
    H=(H+H')/2
    return H
end
PXP_K_Ham(N::Int, k::Int, Omega::Float64=1.0) = PXP_K_Ham(BitStr{N, Int}, k, Omega)

function PXP_MSS_Ham(::Type{T}, k::Int, inv::Int64=1) where {N, T <: BitStr{N}}
#params: a int of lattice number, momentum of system and interaction strength of system which default to be 1
#return: the Hamiltonian matrix in given maximum symmetry space
  
    omegak = exp(2im * π * k / N)
    
    if inv==1
        MSS, MSS_dic, qlist = PXP_MSS_basis(T, k)

        l = length(MSS)
        H = zeros(ComplexF64, (l, l))
        for i in 1:l
            n = MSS[i]
            Zn = sqrt(qlist[i]) / 4 * sqrt(length(MSS_dic[n])) / N
            output = actingH_PXP(T, n, true)
            for m in output
                mbar, d = get_representative(m)
                inv_mbar = get_representative(breflect(mbar))[1]
                mtilde = min(mbar, inv_mbar)
                if mtilde ∈ MSS
                    j = searchsortedfirst(MSS, mtilde)
                    Zm = sqrt(qlist[j]) / 4 * sqrt(length(MSS_dic[mtilde])) / N
                    H[i, j] +=  Zn / Zm*omegak^d
                end
            end
        end
        H=real(H)
        H = (H + H') / 2
    else
        MSS, MSS_dic = PXP_MSS_basis(T, k, -1)

        l = length(MSS)
        H = zeros(ComplexF64, (l, l))
        for i in 1:l
            n = MSS[i]
            Zn = 1 / 4 * sqrt(length(MSS_dic[n])) / N
            output = actingH_PXP(T, n, true)
            for m in output
                mbar, d = get_representative(m)
                if mbar ∈ MSS
                    j=searchsortedfirst(MSS, mbar)
                    Zm = 1 / 4 * sqrt(length(MSS_dic[mbar])) / N
                    H[i, j] +=  Zn / Zm*omegak^d
                end
            end
        end
        H=real(H)
        H = (H + H') / 2  
    end

    return H
end
PXP_MSS_Ham(N::Int, k::Int, inv::Int64=1) = PXP_MSS_Ham(BitStr{N, Int}, k, inv) 

function myprint(io::IO, xs...)
    println(io, xs..., '\n')
    flush(io)
end


function translation_matrix(::Type{T}) where {N, T <: BitStr{N}}
    basis=PXP_basis(T)  
    Mat=zeros(Float64,(length(basis),length(basis)))
    for (i,n) in enumerate(basis)
        m=cyclebits(n, 1)
        j=searchsortedfirst(basis, m)
        Mat[i,j]=1.0
    end
    
    return Mat
end
translation_matrix(N::Int) = translation_matrix(BitStr{N, Int})

function actingHplus_PXP(::Type{T}, rowstate::T) where {N, T <: BitStr{N}}
    mask = bmask(T, N, N-2)
    fl = bmask(T, N-1)


    output = [
    flip(rowstate, fl >> (j-1))
    for j in 1:N-2
    if isodd(j) && rowstate & (mask >> (j-1)) == 0 && rowstate[N-j] == 0 ||
       iseven(j) && rowstate & (mask >> (j-1)) == 0 && rowstate[N-j] == 1
]

    if rowstate[2] == 0 && rowstate[N] == 0
        if iseven(N) && rowstate[1] == 0
            flip_str = flip(rowstate, bmask(T, 1))
            push!(output, flip_str)
        elseif isodd(N) && rowstate[1] == 1
            flip_str = flip(rowstate, bmask(T, 1))
            push!(output, flip_str)
        end
    end
    if rowstate[1] == 0 && rowstate[N-1] == 0 && rowstate[N] == 1
        flip_str = flip(rowstate, bmask(T, N))
        push!(output, flip_str)
    end
    return output
end


function actingHminus_PXP(::Type{T}, rowstate::T) where {N, T <: BitStr{N}}
    output = T[]
    
    mask = bmask(T, N, N-2)
    fl = bmask(T, N-1)

    for j in 1:N-2
        if isodd(j)  # odd site, then apply σ⁻ operator
            if rowstate & (mask >> (j-1)) ==0 && rowstate[N-j] == 1
                flip_str = flip(rowstate, fl >> (j-1))
                push!(output, flip_str)
            end
        else  # even site, then apply σ⁺ operator
            if rowstate & (mask >> (j-1)) ==0 && rowstate[N-j] == 0
                flip_str = flip(rowstate, fl >> (j-1))
                push!(output, flip_str)
            end
        end
    end
    if rowstate[2] == 0 && rowstate[N] == 0
        if iseven(N) && rowstate[1] == 1
            flip_str = flip(rowstate, bmask(T, 1))
            push!(output, flip_str)
        elseif isodd(N) && rowstate[1] == 0
            flip_str = flip(rowstate, bmask(T, 1))
            push!(output, flip_str)
        end
    end
    if rowstate[1] == 0 && rowstate[N-1] == 0 && rowstate[N] == 0
        flip_str = flip(rowstate, bmask(T, N))
        push!(output, flip_str)
    end
    return output
end

function iso_total2FSA(::Type{T}) where {N, T <: BitStr{N}}
    # Once you have isometry, you can use it to map the total basis to the target basis. So you do not need to write the PXP_FSA_basis function.
    basis= PXP_basis(T)
    l = length(basis)
    actingH_PXP_matr = zeros(Int64, (l, l))
    for i in 1:l
        output = actingHplus_PXP(T, basis[i])
        for m in output
            j = searchsortedfirst(basis, m)
            actingH_PXP_matr[j, i] += 1
        end
    end

    initial_state = bmask(T, 2:2:N)
    final_state = initial_state >> 1

    statelis = Matrix{Float64}(undef, l, N+1)
    statelis[:, 1] .= 0; statelis[end, 1] = 1
    statelis[:, end] .= 0; statelis[searchsortedfirst(basis, final_state), end] = 1

    for i in 2:N
        state = actingH_PXP_matr * statelis[:, i-1]
        state /= norm(state)
        statelis[:, i] = state
    end
    return statelis
end
iso_total2FSA(N::Int) = iso_total2FSA(BitStr{N, Int})

function PXP_FSA_Ham(::Type{T}) where {N, T <: BitStr{N}}
#This function is based on Forward Scattering Approximation, utilizes the function PXP_FSA_basis to build the basis, 
#and project the PXP Hamiltonian to its FSA subspace. It has an input parameter N, the system size, and outputs the FSA matrix.
# Arguments
#N::Int: The system size.
#Returns
#`Matrix{Float64}`: The FSA Hamiltonian matrix.
# Examples

    Ham = PXP_Ham(T, true)
    file_path = "a/Users/cycling/Documents/projects/big_data/scar_thermal_FSA/iso_FSA/iso_total2FSA$(N).jld"

    if isfile(file_path)
        iso = load(file_path, "iso")
    else
        iso = iso_total2FSA(T)
    end
    
    H = iso' * Ham * iso
    return H
end
PXP_FSA_Ham(n::Int) = PXP_FSA_Ham(BitStr{n, Int})
# Such behavior is allowed, it is not a performance bottleneck(additional cost is deserved).