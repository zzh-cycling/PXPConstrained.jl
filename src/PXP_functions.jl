"""
This PXP_basis package is used to generate the basis for PXP model, and also the Hamiltonian matrix for PXP model according to their basis.  Totally we have three kinds of functions, basis, Ham, reduced_dm. We consider both PBC and OBC in the basis and Hamiltonian matrix, and we also consider the translational symmetry and inversion symmetry, and FSA subspace.
"""

function Fibonacci_seq(N::Int64)
    # Generate Fibonacci sequence
    if N==0
        return []
    elseif  N==1
        return [1]
    else
        fib=[0,1]
        for i in 3:N
            push!(fib,fib[i-1]+fib[i-2])
        end
        return fib
    end

end


function Fibonacci_chain_OBC(::Type{T}) where {N, T <: BitStr{N}}
    # Generate Fibonacci chain for PXP model with open boundary condition
    fib_chain=[[T(0), T(1)],[T(0), T(1), T(2)]]
    for i in 3:N
        push!(fib_chain,vcat([s << 1 for s in fib_chain[i-1]],[(s << 2 |  T(1)) for s in fib_chain[i-2]]))
    end
    # each push we add a bit"0" or bit"01" to the end of the bit_string, and finally return the last element of the fib_chain
    return fib_chain[N]
end

function Fibonacci_chain_PBC(::Type{T}) where {N, T <: BitStr{N}}
    # Generate Fibonacci chain  for PXP model with periodic boundary condition
    chain=Fibonacci_chain_OBC(T)
    filtered_fib_chain = T[]
    mask = bmask(T, 1, N)
    for s in chain
        if allone(s, mask)
            continue
        else
            push!(filtered_fib_chain, s)
        end
    end
    return filtered_fib_chain
end

function actingH_PXP(::Type{T}, n::T, pbc::Bool=true) where {N, T <: BitStr{N}}
    # The type of n is DitStr{D, N, Int}, which is a binary string with length N in D-ary form.
    # Acting Hamiltonian on a given state in bitstr and return the output states in bitstr
    # Here need to note that the order of the bitstr is from right to left, which is different from the normal order.
    mask=bmask(T, N, N-2)
    fl=bmask(T, N-1)
    output = [flip(n, fl >> (i-1)) for i in 1:N-2 if n & (mask >> (i-1)) == 0] # faster method

    if pbc
        if n[2]==0 && n[N]==0
            flip_str=flip(n,bmask(T, 1))
            push!(output,flip_str)
        end
        if n[1]==0 && n[N-1]==0
            flip_str=flip(n,bmask(T, N))
            push!(output,flip_str)
        end
    else
        if n[N-1]==0
            flip_str=flip(n,bmask(T, N))
            push!(output,flip_str)
        end
        if n[2]==0
            flip_str=flip(n,bmask(T, 1))
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

function find_common_elements(list1, list2)
    # Find common elements in two lists
    common_elements = Set(list1) ∩ Set(list2)
    indices1 = findall(x -> x in common_elements, list1)
    indices2 = findall(x -> x in common_elements, list2)
    return indices1, indices2
end

function reduced_dm_PXP(N::Int64, l::Int64, state::Vector{Float64})
    # At least faster than tn contraction, and idea is similar to the above function
    """
    Function to compute a reduced density matrix out of a given wave function. Its principle is to divide the system into left and right two parts, and use formula rho_reduced=sum_k I * <k|rho|k> * I, then compare how many states in the left part which are the same, while the right part is still the same.
    """
    basis=PXP_basis(N)
    left_sys=[readbit(i, l+1:N...) for i in basis]
    right_sys=[readbit(i, 1:l...) for i in basis]
    m=length(left_sys)
    
    # Here we use a dictionary to store the category of left part of the basis, which is the basis of subsystem.
    sub_basis = Dict{BitStr{N, Int}, Vector{BitStr{N, Int}}}()

    for i in 1:m
        string= left_sys[i]
        if haskey(sub_basis, string)
            push!(sub_basis[string], i)
        else
            sub_basis[string] = [i]
        end
    end
    """
    For example, if we have a basis of 4 qubits, then sub_basis is like
    Dict{DitStr{2, 4, Int64}, Vector{DitStr{2, 4, Int64}}} with 3 entries:
    0000 ₍₂₎ => [0001 ₍₂₎, 0010 ₍₂₎, 0011 ₍₂₎]
    0010 ₍₂₎ => [0110 ₍₂₎, 0111 ₍₂₎]
    0001 ₍₂₎ => [0100 ₍₂₎, 0101 ₍₂₎]

    ordered_sub_basis is
    3-element Vector{DitStr{2, 4, Int64}}:
    0000 ₍₂₎
    0001 ₍₂₎
    0010 ₍₂₎

    3-element Vector{Vector{Int64}}:
    [1, 2, 3]
    [4, 5]
    [6, 7]
    """

    # Sort the sub_basis in natural order, and obtain the corresponding index in PXP_basis of sub_basis
    ordered_sub_basis = sort(collect(keys(sub_basis)), by = x -> x)
    sub_basis_index=[Int.(sub_basis[k]) for k in ordered_sub_basis]
    size=length(ordered_sub_basis)

    # for given basis of subsystem, we first obtain which state in total system will give the target subsystem basis.
    reduced_dm=zeros((size,size))
    for i in 1:size
        for j in i:size
            row=right_sys[sub_basis_index[i]]
            column=right_sys[sub_basis_index[j]]
            idx1,idx2=find_common_elements(row,column)
            """
            Actually, in i,j=1,2 case, basis[index1],basis[index2]=(["0000", "0001"], ["0100", "0101"]), row is right part of |a> will give "00" for k <k| * I|a>, column is right part of |a> will give "01" for k <a|I * |k>, find_common_elements is try to keep k is the same. 
            """
            # @show sub_basis_index, row
            
            index1=sub_basis_index[i][idx1]
            index2=sub_basis_index[j][idx2]
            reduced_dm[i,j]= dot(state[index1],state[index2])
            reduced_dm[j,i]=reduced_dm[i,j]
        end
    end

    return reduced_dm, ordered_sub_basis
end

function myreadbit(bit::BitStr, index::Vector{Int64})
    # read the indexed bit of a bit string
    mapreduce(el -> readbit(bit, el[2]) << (el[1]-1), |, enumerate(index))
end

function rdm_PXP(::Type{T}, subsystems::Vector{Int64}, state::Vector{ET}) where {N, T <: BitStr{N}, ET}
    # At least faster than tn contraction
    # Function to compute a reduced density matrix out of a given wave function. Its principle is to divide the system into left and right two parts, and use formula rho_reduced=sum_k I * <k|rho|k> * I, then compare how many states in the left part which are the same, while the right part is still the same.

    basis = PXP_basis(T)
    system = [myreadbit(i, subsystems) for i in basis]
    environ = [myreadbit(i, setdiff(1:N, subsystems)) for i in basis]
    
    # Here we use a dictionary to store the category of left part of the basis, which is the basis of subsystem.
    sub_basis = Dict{T, Vector{Int}}()  # maps subsystem basis to its indices in `state`.

    for (i, string) in enumerate(system)
        if haskey(sub_basis, string)
            push!(sub_basis[string], i)
        else
            sub_basis[string] = [i]
        end
    end

    # Sort the sub_basis in natural order, and obtain the corresponding index in PXP_basis of sub_basis
    ordered_sub_basis = sort(collect(keys(sub_basis)))
    sub_basis_index = [sub_basis[k] for k in ordered_sub_basis]
    size=length(ordered_sub_basis)

    # for given basis of subsystem, we first obtain which state in total system will give the target subsystem basis.
    reduced_dm = zeros(ET, (size, size))
    for i in 1:size, j in i:size
        idrs = sub_basis_index[i]
        idcs = sub_basis_index[j]
        val = zero(ET)
        for idr in idrs
            r = environ[idr]
            idc = findfirst(==(r), idcs)
            if idc !== nothing
                val += state[idr]' * state[idc]
            end
        end
        reduced_dm[i, j] = val
        reduced_dm[j, i] = reduced_dm[i, j]
    end
    
    return reduced_dm
end

function iso_total2K(::Type{T}, k::Int64) where {N, T <: BitStr{N}}
    """
    Function to map the total basis to the K space basis
    """
    basis = PXP_basis(N)

    indices = Dict{Int, Vector{Int64}}()

    for i in eachindex(basis)
        state=basis[i]
        category = get_representative(state)[1]
        if haskey(indices, category)
            push!(indices[category], i)
        else
            indices[category] = [i]
        end
    end
    basisK=sort(collect(keys(indices)), by = x -> x)
    iso = zeros((length(basis), length(basisK)))
    
    # 预分配 temp_lis
    temp_lis = zeros(length(basis))

    for i in eachindex(basisK)
        state_index = indices[basisK[i]]
        l = length(state_index)

        # 重置 temp_lis
        fill!(temp_lis, 0)  # 先将 temp_lis 清零
        temp_lis[state_index] .= 1/sqrt(l)
        iso[:, i] = temp_lis
    end
    return iso
end

function rdm_PXP_K(N::Int64, subsystems::Vector{Int64}, state::Vector{Float64}, k::Int64)
    state=iso_total2K(N,k)*state
    reduced_dm = rdm_PXP(N, subsystems, state)
    return reduced_dm
end

function Inv_proj_matrix(N::Int64)
    basis=PXP_basis(N)
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

function cyclebits(state::T, n_translations::Int) where {N, T <: BitStr{N}}
    """
    :params: t is an integer, N is the length of the binary string
    :n_translations: number of positions to shift
    :return: the binary left shift
    We also use this order: system size, state, number of translations
    """
    return (state << n_translations) % (2^N - 1)

    # n_translations = mod(n_translations, N)
    # shifted_right = state >> n_translations
    # shifted_left = (state << (N - n_translations)) & ((1 << N) - 1)
    
    # rotated_n = shifted_right | shifted_left
end

function get_representative(state::T) where {N, T <: BitStr{N}}
    """
    Finds representative and representative translation for a state.
    State should be a decimal integer.
    """
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
    """
    :params: a int of lattice number and momentum of system
    :return: computational basis in given momentum kinetically constrained subspace with decimal int form in PXP model
    """
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

function PXP_MSS_basis(::Type{T}, k::Int64,inv::Int64=1) where {N, T <: BitStr{N}}
    """
    :params: a int of lattice number and momentum of system, we have considered the inversion symmetry
    :return: computational basis in given momentum inversion symmetry subspace with decimal int form
    """
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


function PXP_K_Ham(::Type{T}, k::Int, Omega::Float64=1.0) where {N, T <: BitStr{N}}
    """
    :params: a int of lattice number, momentum of system and interaction strength of system which default to be 1
    :return: the Hamiltonian matrix in given K space
    """
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

function PXP_MSS_Ham(::Type{T}, k::Int, inv::Int64=1) where {N, T <: BitStr{N}}
    """
    :params: a int of lattice number, momentum of system and interaction strength of system which default to be 1
    :return: the Hamiltonian matrix in given maximum symmetry space
    """    
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
        return H
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
        return H
    end
end

function wf_time_evolution(psi0::Union{Vector{Float64}, Vector{ComplexF64}}, times::Vector{Float64}, energy::Vector{Float64},states::Matrix{Float64})
    T=typeof(psi0)
    wflis=Vector{T}(undef,length(times))
    for (i,t) in enumerate(times)
        wf=similar(psi0)
        for (j, state) in enumerate(eachcol(states))
            wf+=exp(-1im*t*energy[j])*dot(psi0,state)*state
        end
        wflis[i]=wf
        # state[end]
    end
    return wflis
end

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


function sum_of_powers(N::Int)
    total_sum = 0
    current_power = N - 1
    while current_power >= 1
        total_sum += 2^current_power
        current_power -= 2
    end
    return total_sum
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

function PXP_FSA_Ham(::Type{T}) where {N, T <: BitStr{N}}
    """
    This function is based on Forward Scattering Approximation, utilizes the function PXP_FSA_basis to build the basis, and project the PXP Hamiltonian to its FSA subspace. It has an input parameter N, the system size, and outputs the FSA matrix.

    # Arguments
    - `N::Int`: The system size.

    # Returns
    - `Matrix{Float64}`: The FSA Hamiltonian matrix.

    # Examples
    """
    Ham = PXP_Ham(T, true)
    iso = load("/Users/cycling/Documents/projects/big_data/scar_thermal_FSA/iso_FSA/iso_total2FSA$(N).jld", "iso")
    # iso = iso_total2FSA(N)
    H = iso' * Ham * iso
    return H
end