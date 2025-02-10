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

function Fibonacci_chain_OBC(N::Int64)
    # Generate Fibonacci chain for PXP model with open boundary condition
    fib_chain=[["0", "1"],["00", "01","10"]]
    for i in 3:N
        push!(fib_chain,vcat([s * "0" for s in fib_chain[i-1]],[s * "01" for s in fib_chain[i-2]]))
    end
    return fib_chain[N]
end


function Fibonacci_chain_PBC(N::Int64)
    # Generate Fibonacci chain  for PXP model with periodic boundary condition
    if N == 0
        return []
    elseif N == 1
        return ["0", "1"]
    else
        fib_chain=[["0", "1"],["00", "01","10"]]
        for i in 3:N
            push!(fib_chain,vcat([s * "0" for s in fib_chain[i-1]],[s * "01" for s in fib_chain[i-2]]))
        end
        filtered_fib_chain = []
        for s in fib_chain[N]
            if s[1] == '1' && s[end] == '1'
                continue
            else
                push!(filtered_fib_chain, s)
            end
        end
        return filtered_fib_chain
    end
end


function actingH_PXP(N::Int64, n::String,pbc::Bool=true)
    # Acting Hamiltonian on a given state in str and return the output states in str
    output=[]
    for i in 1:N-2
        if n[i]=='0' && n[i+2]=='0'
            column=n[1:i]*string(xor(1,parse(Int,n[i+1])))*n[i+2:end]
            push!(output,column)
        end
    end
    if pbc
        if n[N-1]=='0' && n[1]=='0'
            column=n[1:N-1]*string(xor(1,parse(Int,n[N])))
            push!(output,column)
        end
        if n[N]=='0' && n[2]=='0'
            column=string(xor(1,parse(Int,n[1])))*n[2:N]
            push!(output,column)
        end
    else
        if n[2]=='0'
            column=string(xor(1,parse(Int,n[1])))*n[2:end]
            push!(output,column)
        end
        if n[end-1]=='0'
            column=n[1:end-1]*string(xor(1,parse(Int,n[end])))
            push!(output,column)
        end
    end
    return output
end

function PXP_basis(N::Int64,pbc::Bool=true)
    # Generate basis for PXP model, return both decimal and binary form, where we both consider PBC and OBC
    if pbc
        basis_string=Fibonacci_chain_PBC(N)
        basis_int=[parse(Int,i,base=2) for i in basis_string]
    else
        basis_string=Fibonacci_chain_OBC(N)
        basis_int=[parse(Int,i,base=2) for i in basis_string]
    end
    zipped_lis=zip(basis_int,basis_string)
    sorted_lis=sort(collect(zipped_lis),by=x -> x[1])
    unzip(a) = map(x->getfield.(a, x), fieldnames(eltype(a)))
    basis_int,basis_string=unzip(sorted_lis)
    
    return basis_int, basis_string
end

function PXP_Ham(N::Int64, pbc::Bool=true)
        # Generate Hamiltonian for PXP model, automotically contain pbc or obc
        basis_int, basis_string=PXP_basis(N,pbc)

        l=length(basis_int)
        H=zeros(Float64,(l,l))
        for i in 1:l
            output=actingH_PXP(N,basis_string[i],pbc) 
            for m in output 
                j=searchsortedfirst(basis_int,parse(Int, m, base=2))
                H[i,j]+=1
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

function reduced_dm(psi::Union{Vector{Float64}, Vector{ComplexF64}}, n_left::Int, n_total::Int)
    """
    计算约化密度矩阵
    参数:
    - psi: 纯态的波函数，表示为一个复数向量
    - n_left: 左边子系统的量子比特数
    - n_total: 总的量子比特数
    """
    # 计算右边子系统的量子比特数
    n_right = n_total - n_left

    # 计算子系统的维度
    dim_left = 2^n_left
    dim_right = 2^n_right

    # 初始化约化密度矩阵
    rho_reduced = zeros(ComplexF64, dim_left, dim_left)

    # 遍历所有可能的子系统状态
    for i in 1:dim_left
        for j in 1:dim_left
            for k in 1:dim_right
                # 计算全局态的索引
                index_i = (i - 1) * dim_right + k
                index_j = (j - 1) * dim_right + k

                # 累加到约化密度矩阵
                rho_reduced[i, j] += psi[index_i] * conj(psi[index_j])
            end
        end
    end

    return rho_reduced
end

function reduced_dm_PXP(N::Int64, l::Int64, state::Vector{Float64})
    # At least faster than tn contraction, and idea is similar to the above function
    """
    Function to compute a reduced density matrix out of a given wave function. Its principle is to divide the system into left and right two parts, and use formula rho_reduced=sum_k I * <k|rho|k> * I, then compare how many states in the left part which are the same, while the right part is still the same.
    """
    basis_int, basis_string=PXP_basis(N)
    left_sys_string=[i[1:l] for i in basis_string]
    right_sys_string=[i[l+1:N] for i in basis_string]
    m=length(left_sys_string)
    
    # Here we use a dictionary to store the category of left part of the basis, which is the basis of subsystem.
    sub_basis = Dict{String, Vector{Int64}}()

    for i in 1:m
        string= left_sys_string[i]
        if haskey(sub_basis, string)
            push!(sub_basis[string], i)
        else
            sub_basis[string] = [i]
        end
    end
    """
    For example, if we have a basis of 4 qubits, then sub_basis is like
    Dict{String, Vector{Int64}} with 3 entries:
     "00" => [1, 2, 3]
     "10" => [6, 7]
     "01" => [4, 5]

    ordered_sub_basis is 
    3-element Vector{String}:
     "00"
     "01"
     "10"

    sub_basis_number is 
    3-element Vector{Vector{Int64}}:
     [1, 2, 3]
     [4, 5]
     [6, 7]
    """

    # Sort the sub_basis in natural order, and obtain the corresponding index in PXP_basis of sub_basis
    ordered_sub_basis = sort(collect(keys(sub_basis)), by = x -> parse(Int, x, base = 2))
    sub_basis_number=[sub_basis[k] for k in ordered_sub_basis]
    size=length(ordered_sub_basis)

    # for given basis of subsystem, we first obtain which state in total system will give the target subsystem basis.
    reduced_dm=zeros((size,size))
    for i in 1:size
        for j in i:size
            row=right_sys_string[sub_basis_number[i]]
            column=right_sys_string[sub_basis_number[j]]
            idx1,idx2=find_common_elements(row,column)
            """
            Actually, in i,j=1,2 case, basis_string[index1],basis_string[index2]=(["0000", "0001"], ["0100", "0101"]), row is right part of |a> will give "00" for k <k| * I|a>, column is right part of |a> will give "01" for k <a|I * |k>, find_common_elements is try to keep k is the same. 
            """
            # @show sub_basis_number, row
            
            index1=sub_basis_number[i][idx1]
            index2=sub_basis_number[j][idx2]
            reduced_dm[i,j]= dot(state[index1],state[index2])
            reduced_dm[j,i]=reduced_dm[i,j]
        end
    end

    return reduced_dm, ordered_sub_basis
end

function rdm_PXP(N::Int64, subsystems::Vector{Int64}, state::Vector{Float64})
    # At least faster than tn contraction
    """
    Function to compute a reduced density matrix out of a given wave function. Its principle is to divide the system into left and right two parts, and use formula rho_reduced=sum_k I * <k|rho|k> * I, then compare how many states in the left part which are the same, while the right part is still the same.
    """
    basis_int, basis_string = PXP_basis(N)
    system_string = [i[subsystems] for i in basis_string]
    environ_string = [i[setdiff(1:N, subsystems)] for i in basis_string]
    m = length(system_string)
    l = length(subsystems)
    
    # Here we use a dictionary to store the category of left part of the basis, which is the basis of subsystem.
    sub_basis = Dict{String, Vector{Int64}}()

    for i in 1:m
        string = system_string[i]
        if haskey(sub_basis, string)
            push!(sub_basis[string], i)
        else
            sub_basis[string] = [i]
        end
    end

    """
    For example, if we have a basis of 4 qubits, then sub_basis might look like:
    Dict{String, Vector{Int64}} with 3 entries:
     "00" => [1, 2, 3]
     "10" => [6, 7]
     "01" => [4, 5]

    sorted_keys is 
    3-element Vector{String}:
     "00"
     "01"
     "10"

    valuelis is 
    3-element Vector{Vector{Int64}}:
     [1, 2, 3]
     [4, 5]
     [6, 7]
    """

    # Sort the sub_basis in natural order, and obtain the corresponding index in PXP_basis of sub_basis
    ordered_sub_basis = sort(collect(keys(sub_basis)), by = x -> parse(Int, x, base = 2))
    sub_basis_number=[sub_basis[k] for k in ordered_sub_basis]
    size=length(ordered_sub_basis)

    # for given basis of subsystem, we first obtain which state in total system will give the target subsystem basis.
    reduced_dm = zeros((size, size))
    for i in 1:size
        for j in i:size
            row = environ_string[sub_basis_number[i]]
            column = environ_string[sub_basis_number[j]]
            idx1, idx2 = find_common_elements(row, column)
            """
            Actually, in i,j=1,2 case, basis_string[index1],basis_string[index2]=(["0000", "0001"], ["0100", "0101"]), row is right part of |a> will give "00" for k <k| * I|a>, column is right part of |a> will give "01" for k <a|I * |k>, find_common_elements is try to keep k is the same. 
            """
            index1 = sub_basis_number[i][idx1]
            index2 = sub_basis_number[j][idx2]
            reduced_dm[i, j] = dot(state[index1], state[index2])
            reduced_dm[j, i] = reduced_dm[i, j]
        end
    end
    
    return reduced_dm
end


function iso_total2K(N::Int64, k::Int64)
    """
    Function to map the total basis to the K space basis
    """
    basis, basis_string = PXP_basis(N)

    indices = Dict{Int, Vector{Int64}}()

    for i in eachindex(basis)
        state=basis[i]
        category = get_representative(N, state)[1]
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


function delete_at(arr::AbstractMatrix{T}, index::Int, axis::Int) where T
    if axis == 1  # 删除行
        return vcat(arr[1:index-1, :], arr[index+1:end, :])
    elseif axis == 2  # 删除列
        return hcat(arr[:, 1:index-1], arr[:, index+1:end])
    else
        throw(ArgumentError("axis must be 1 (row) or 2 (column)"))
    end
end

function delete_rowcolumn(arr, indices)
    # Sort indices in descending order to avoid shifting issues
    sorted_indices = sort(indices, rev=true)
    result = copy(arr)
    for idx in sorted_indices
        result = delete_at(result, idx, 1)  # 删除行
    end
    for idx in sorted_indices
        result = delete_at(result, idx, 2)  # 删除列
    end
    return result
end


function inversion_strlis(strlis::Vector{String})
    """
    - `strlis`: input a vector of strings like ["0000", "0100", "1000", "0010", "1010", "0001", "0101", "1001"]
    Return the inversion symmetry states in the string list 
    """
    inversion(str)=@. reverse(str)
    inversed_lis = similar(strlis)  # Initialize an empty array for the inversed strings

    for (idx,i) in enumerate(strlis)
        inversed_lis[idx]=inversion(i)  # Append the inverted string to the list
    end

    return inversed_lis
end

function Inv_proj_matrix(N::Int64, inv::Int64=1)
    basis_int, basis_string=PXP_basis(N)
    l=length(basis_int)
    Imatrix=zeros((l,l))
    reversed_basis_string=inversion_strlis(basis_string)
    for i in eachindex(basis_int)
        output=reversed_basis_string[i]
        j=searchsortedfirst(basis_int,parse(Int, output, base=2))
        Imatrix[i,j]+=1.0
    end
    # positive_inversion_states=Inv_proj_matrix(16,1)*states[:,1075]
    # negative_inversion_states=Inv_proj_matrix(16,-1)*states[:,1075]
    # if inv==1
    #     return (I(l)+Imatrix)/2
    # end
    # if inv==-1
    #     return (I(l)-Imatrix)/2
    # end
    return Imatrix
end

function cyclebits(N::Int,state::Int64, n_translations::Int)
    """
    :params: t is an integer, N is the length of the binary string
    :n_translations: number of positions to shift
    :return: the binary left shift
    We also use this order: system size, state, number of translations
    """
    return (state << n_translations) % (2^N - 1)
end

function cyclestring(state::String, n_translations::Int=1)
    """
    :params: str is a string
    :n_translations: >=1
    :return: the binary left shift
    """
    return state[n_translations+1:end] * state[1:n_translations]
end

function bin(N::Int,n::Int)
    return lpad(string(n, base=2), N, '0')
end

function get_representative(N::Int, state::Int)
    """
    Finds representative and representative translation for a state.
    State should be a decimal integer.
    """
    representative = state
    translation = 0
    for n_translation_sites in 0:N-1
        new_state = cyclebits(N,state, n_translation_sites)
        if new_state < representative
            representative = new_state
            translation = n_translation_sites
        end
    end
    return representative, translation
end

function PXP_K_basis(N::Int64, k::Int64)
    """
    :params: a int of lattice number and momentum of system
    :return: computational basis in given momentum kinetically constrained subspace with decimal int form in PXP model
    """
    basisK_string = Vector{String}(undef, 0)
    basisK_int = Vector{Int}(undef, 0)
    basis_int, basis_string = PXP_basis(N)


    basis_dic = Dict{Int, Vector{Int}}()
    for i in basis_int
        category = get_representative(N, i)[1]
        if haskey(basis_dic, category)
            push!(basis_dic[category], i)
        else
            basis_dic[category] = [i]
        end
    end

    for j in eachindex(basis_int)
        n=basis_int[j]
        n_str=basis_string[j]
        RS = get_representative(N, n)[1]
        if RS == n && (k * length(basis_dic[RS])) % N == 0
            push!(basisK_int, n)
            push!(basisK_string, n_str)
        end
    end

    return basisK_int, basisK_string, basis_dic
end

function PXP_MSS_basis(N::Int64, k::Int64,inv::Int64=1)
    """
    :params: a int of lattice number and momentum of system, we have considered the inversion symmetry
    :return: computational basis in given momentum inversion symmetry subspace with decimal int form
    """
    # MSS_int is the list of states in the maximum symmetry sector
    MSS_int = Vector{Int}(undef, 0)
    basisK_int, basisK_string, basis_dic = PXP_K_basis(N, k)
    MSS_dic = Dict{Int, Vector{Int}}()

    # q is the number of states that are equivalent under inversion
    qlist = Vector{Int}(undef, 0)
    if inv==1
        for i in eachindex(basisK_int)
            n = basisK_int[i]
            n_str=basisK_string[i]
            # here we calculate the representative state of the inversion of n
            nR = get_representative(N, parse(Int, reverse(n_str), base=2))[1]
            if n <= min(nR, n)
                push!(MSS_int, n)
                MSS_dic[n] = basis_dic[n]
                push!(qlist, length(Set([n, nR])))
            end
        end

        return MSS_int, MSS_dic, qlist

    else
        for i in eachindex(basisK_int)
            n = basisK_int[i]
            n_str=basisK_string[i]
            # here we calculate the representative state of the inversion of n
            nR = get_representative(N, parse(Int, reverse(n_str), base=2))[1]
            if n <= min(nR, n)
                push!(MSS_int, n)
                MSS_dic[n] = basis_dic[n]
                push!(qlist, length(Set([n, nR])))
            end
        end    
        index=findall(x -> x==2, qlist)
        new_MSS_dic = Dict(k => v for k in MSS_int[index] for v in [MSS_dic[k]])
        return MSS_int[index], new_MSS_dic
    end
    
end


function PXP_K_Ham(N::Int, k::Int, Omega::Float64=1.0)
    """
    :params: a int of lattice number, momentum of system and interaction strength of system which default to be 1
    :return: the Hamiltonian matrix in given K space
    """
    basisK_int, basisK_string, basis_dic = PXP_K_basis(N, k)
    l = length(basisK_int)
    omegak = exp(2im * π * k / N)
    H = zeros(ComplexF64, (l, l))

    for i in 1:l
        n=basisK_int[i]
        n_str=basisK_string[i]
        output = actingH_PXP(N, n_str, true)
        for m in output
            mbar, d = get_representative(N, parse(Int, m, base=2))
            if mbar ∈ basisK_int
                j=searchsortedfirst(basisK_int,mbar)
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

function PXP_MSS_Ham(N::Int, k::Int, inv::Int64=1)
    """
    :params: a int of lattice number, momentum of system and interaction strength of system which default to be 1
    :return: the Hamiltonian matrix in given maximum symmetry space
    """    
    omegak = exp(2im * π * k / N)
    
    if inv==1
        MSS_int, MSS_dic, qlist = PXP_MSS_basis(N, k)
        MSS_string=[bin(N,i) for i in MSS_int]

        l = length(MSS_int)
        H = zeros(ComplexF64, (l, l))
        for i in 1:l
            n = MSS_int[i]
            n_str = bin(N, n)
            Zn = sqrt(qlist[i]) / 4 * sqrt(length(MSS_dic[n])) / N
            output = actingH_PXP(N, n_str, true)
            for m in output
                mbar, d = get_representative(N, parse(Int, m, base=2))
                inv_mbar = get_representative(N, parse(Int, reverse(bin(N, mbar)), base=2))[1]
                mtilde = min(mbar, inv_mbar)
                if mtilde ∈ MSS_int
                    j=searchsortedfirst(MSS_int,mtilde)
                    Zm = sqrt(qlist[j]) / 4 * sqrt(length(MSS_dic[mtilde])) / N
                    H[i, j] +=  Zn / Zm*omegak^d
                end
            end
        end
        H = (H + H') / 2
        return H
    else
        MSS_int, MSS_dic = PXP_MSS_basis(N, k, -1)
        MSS_string=[bin(N,i) for i in MSS_int]
        l = length(MSS_int)
        H = zeros(ComplexF64, (l, l))
        for i in 1:l
            n = MSS_int[i]
            n_str = bin(N, n)
            Zn = 1 / 4 * sqrt(length(MSS_dic[n])) / N
            output = actingH_PXP(N, n_str, true)
            for m in output
                mbar, d = get_representative(N, parse(Int, m, base=2))
                if mbar ∈ MSS_int
                    j=searchsortedfirst(MSS_int,mbar)
                    Zm = 1 / 4 * sqrt(length(MSS_dic[mbar])) / N
                    H[i, j] +=  Zn / Zm*omegak^d
                end
            end
        end
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


function translation_matrix(N::Int64)
    basis_int, basis_string=PXP_basis(N)
    Mat=zeros(Float64,(length(basis_int),length(basis_int)))
    for (i,n_str) in enumerate(basis_string)
        m_str=cyclestring(n_str,1)
        m=parse(Int,m_str,base=2)
        j=searchsortedfirst(basis_int,m)
        Mat[i,j]=1.0
    end
    
    return Mat
end


function actingHplus_PXP(N::Int, state::Int)
    output = Int64[]
    rowstate = lpad(string(state, base=2), N, '0')
    for j in 1:N-2
        if isodd(j)  # even site, then apply σ⁺ operator
            if rowstate[j] == '0' && rowstate[j+2] == '0' && rowstate[j+1] == '0'
                column = rowstate[1:j] * '1' * rowstate[j+2:end]
                push!(output, parse(Int, column, base=2))
            end
        else  # odd site, then apply σ⁻ operator
            if rowstate[j] == '0' && rowstate[j+2] == '0' && rowstate[j+1] == '1'
                column = rowstate[1:j] * '0' * rowstate[j+2:end]
                push!(output, parse(Int, column, base=2))
            end
        end
    end
    if rowstate[N-1] == '0' && rowstate[1] == '0'
        if iseven(N) && rowstate[N] == '0'
            column = rowstate[1:N-1] * '1'
            push!(output, parse(Int, column, base=2))
        elseif isodd(N) && rowstate[N] == '1'
            column = rowstate[1:N-1] * '0'
            push!(output, parse(Int, column, base=2))
        end
    end
    if rowstate[N] == '0' && rowstate[2] == '0' && rowstate[1] == '1'
        column = '0' * rowstate[2:end]
        push!(output, parse(Int, column, base=2))
    end
    return output
end

function actingHminus_PXP(N::Int, state::Int)
    output = Int64[]
    rowstate = lpad(string(state, base=2), N, '0')
    for j in 1:N-2
        if isodd(j)  # odd site, then apply σ⁻ operator
            if rowstate[j] == '0' && rowstate[j+2] == '0' && rowstate[j+1] == '1'
                column = rowstate[1:j] * '0' * rowstate[j+2:end]
                push!(output, parse(Int, column, base=2))
            end
        else  # even site, then apply σ⁺ operator
            if rowstate[j] == '0' && rowstate[j+2] == '0' && rowstate[j+1] == '0'
                column = rowstate[1:j] * '1' * rowstate[j+2:end]
                push!(output, parse(Int, column, base=2))
            end
        end
    end
    if rowstate[N-1] == '0' && rowstate[1] == '0'
        if iseven(N) && rowstate[N] == '1'
            column = rowstate[1:N-1] * '0'
            push!(output, parse(Int, column, base=2))
        elseif isodd(N) && rowstate[N] == '0'
            column = rowstate[1:N-1] * '1'
            push!(output, parse(Int, column, base=2))
        end
    end
    if rowstate[N] == '0' && rowstate[2] == '0' && rowstate[1] == '0'
        column = '1' * rowstate[2:end]
        push!(output, parse(Int, column, base=2))
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

function iso_total2FSA_Matr(N::Int64)
    basis_int, basis_string = PXP_basis(N)
    l=length(basis_int)
    actingH_PXP_matr=zeros(Int64,(l,l))
    for i in 1:l
        output=actingHplus_PXP(N,basis_int[i])
        for m in output
            j=searchsortedfirst(basis_int,m)
            actingH_PXP_matr[j,i]+=1
        end
    end

    initial_value = sum_of_powers(N)
    initial_state = zeros(l)
    final_state = zeros(l)
    initial_state[end]=1.0
    final_state[searchsortedfirst(basis_int, div(initial_value,2))]=1.0
    statelis = Vector{Vector{Float64}}(undef, N+1)
    statelis[1] = initial_state
    statelis[N+1] = final_state

    for i in 2:N
        state=actingH_PXP_matr*statelis[i-1]
        state/=norm(state)
        statelis[i]=state
    end

    Matrix=zeros(Float64,(l,N+1))
    for i in 1:N+1
        Matrix[:,i]=statelis[i]
    end
    return Matrix
end 

function PXP_FSA_basis(N::Int)
    """
        PXP_FSA_basis(N::Int) -> statelis_comp, statelen

    This function takes one integer arguments, `N`, and returns two integer Vector result. 
    It generates the FSA basis for N particles PXP model. Starting from very initial state, it generates the basis by acting the generating operator on the previous states. It returns the FSA basis in total basis as state list and the list containing repeated elements of the basis.

    # Arguments
    - `N::Int64`: The system size.

    # Returns
    - `statelis_comp`: The compressed statelis eliminating the repeated elements.
    - `statelen`: The list containing the number of times each element is repeated.
    # Examples
    - N=6
    - statelis_comp 7-element Vector{Vector{Int64}}:
        [42]
        [34, 40, 10]
        [32, 2, 8]
        [36, 0, 18, 9]
        [4, 16, 1]
        [20, 5, 17]
        [21]
    - statelen 7-element Vector{Vector{Int64}}:
        [1]
        [1, 1, 1]
        [2, 2, 2]
        [2, 6, 2, 2]
        [8, 8, 8]
        [16, 16, 16]
        [48]
    """
    initial_state = sum_of_powers(N)
    statelis = Vector{Vector{Int64}}(undef, N+1)
    statelis[1] = [initial_state]
    statelen = Vector{Vector{Int64}}(undef, N+1)
    statelen[1] = [1]
    statelis_comp=Vector{Vector{Int64}}(undef, N+1)
    statelis_comp[1] = [initial_state]
    for i in 2:N
        # temp = Set{Int64}() #这里不能直接用Set，因为并不能保证每一个元素重复的次数一样多。
        temp=Vector{Int64}(undef, 0)
        @time for j in eachindex(statelis_comp[i-1])
            temp_state=statelis_comp[i-1][j]
            mapped_states = actingHplus_PXP(N, temp_state)
            for state in mapped_states
                push!(temp, fill(state,statelen[i-1][j])...)
            end
        end
        statelis[i] = temp

        counts = Dict{Int, Int}()
        for element in temp
            counts[element] = get(counts, element, 0) + 1
        end

        statelis_comp[i] = collect(keys(counts))
        statelen[i]=[counts[key] for key in keys(counts)]
        
    end
    #     statelis_comp[i] = unique_temp #这里不能直接赋值statelis，因为每一次H^= action需要的元素是给定的，先set就减少了很多。
    #     # len_temp=[count(x-> x== element, temp) for element in unique_temp]
    #     # statelen[i]=div.(len_temp,minimum(len_temp)) #这里也不能这样除
    #     statelen[i]=[counts[key] for key in keys(counts)]
    #     statelis[i] = temp  

    # end

    statelis[N+1] = [div(initial_state, 2)]
    statelen[N+1] = [1]
    statelis_comp[N+1] = [div(initial_state, 2)]
    return statelis_comp, statelen
end


function iso_total2FSA(N::Int64)
    statelis, statelen = PXP_FSA_basis(N)
    l = length(statelis)

    basis_int, basis_string = PXP_basis(N, true)
    basis_dict = Dict(basis_int[i] => i for i in eachindex(basis_int))


    totalspace_state = zeros(Float64, length(basis_int), l)

    for i in eachindex(statelis)
        statelis[i] = [basis_dict[state] for state in statelis[i]]
    end

    for (i, state) in enumerate(statelis)
        total_state = zeros(Float64, length(basis_int))
        total_state[state] .= statelen[i]

        state = total_state / norm(total_state)
        totalspace_state[:, i] = state
    end

    return totalspace_state
end

function PXP_FSA_Ham(N::Int)
    """
    This function is based on Forward Scattering Approximation, utilizes the function PXP_FSA_basis to build the basis, and project the PXP Hamiltonian to its FSA subspace. It has an input parameter N, the system size, and outputs the FSA matrix.

    # Arguments
    - `N::Int`: The system size.

    # Returns
    - `Matrix{Float64}`: The FSA Hamiltonian matrix.

    # Examples
    """
   
    Ham = PXP_Ham(N, true)
    iso=load("/Users/cycling/Documents/projects/big_data/scar_thermal_FSA/iso_FSA/iso_total2FSA$(N).jld", "iso")
    # iso=iso_total2FSA(N)
    H=iso'*Ham*iso
    return H
end