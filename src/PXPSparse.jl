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

function PXP_MSS_Ham_sparse(::Type{T}, k::Int, inv::Int64=1) where {N, T <: BitStr{N}}
    #params: a int of lattice number, momentum of system and interaction strength of system which default to be 1, k is the momentum of system, only can take 0 or pi, inv is the inversion of the Hamiltonian, only 1 or -1.
    #return: the Hamiltonian matrix in given maximum symmetry space   
    # omegak = exp(2im * π * k / N)
    @assert k == 0 || k==div(N,2) "k is expected to be 0 or $(div(N,2)), but got $k"
    @assert inv ==1 || inv==-1 "inv is expected to be 1 or -1, but got $(inv)"

    omegak= k == 0 ? 1 : -1
    
    if inv==1 && k==0 || k==div(N,2) && inv==-1
        MSS, MSS_dic, qlist = PXP_MSS_basis(T, k, inv)

        l = length(MSS)
        # H = spzeros(Float64, (l, l))
        # matrix_elements = Dict{Tuple{Int,Int}, Float64}()
        I, J, V = Int[], Int[], Float64[]
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
                    # H[i, j] += Zn / Zm*omegak^d 
                    push!(I, i); push!(J, j); push!(V, Zn / Zm*omegak^d)
                end
            end
        end
    
        H = sparse(I, J, V, l, l)
        H = (H + H') / 2
    elseif inv==1 && k==div(N,2) || inv==-1 && k==0
        MSS, MSS_dic = PXP_MSS_basis(T, k, inv)

        l = length(MSS)
        I, J, V = Int[], Int[], Float64[]
        for i in 1:l
            n = MSS[i]
            Zn = 1 / 4 * sqrt(length(MSS_dic[n])) / N
            output = actingH_PXP(T, n, true)
            for m in output
                mbar, d = get_representative(m)
                if mbar ∈ MSS
                    j=searchsortedfirst(MSS, mbar)
                    Zm = 1 / 4 * sqrt(length(MSS_dic[mbar])) / N
                    # H[i, j] +=  Zn / Zm*omegak^d
                    push!(I, i); push!(J, j); push!(V, Zn / Zm*omegak^d)
                end
            end
        end
    
        H = sparse(I, J, V, l, l)
        H = (H + H') / 2
    end

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
#Function to map the MSS basis to the K space basis
    @assert k == 0 || k==div(N,2) "k is expected to be 0 or $(div(N,2)), but got $k"
    @assert inv ==1 || inv==-1 "inv is expected to be 1 or -1, but got $(inv)"

    basis = PXP_basis(T)
    basisK, k_dic = PXP_K_basis(T, k)

    MSS_dic = Dict{Int, Vector{Int64}}()
    qlist = Vector{Int}(undef, 0)
    # Below procedure is to collapse the extra basis in K space that can be converted mutually to MSS space.
    if inv==1 && k==0 || k==div(N,2) && inv==-1
        for i in eachindex(basisK)
            n = basisK[i]
            # here we calculate the representative state of the inversion of n
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

    num_states = length(basisK)
    num_categories = length(keys(MSS_dic))
    rows = Vector{Int64}[]
    cols = Vector{Int64}[]
    vals = Vector{Float64}[]

    MSS_dic=sort(MSS_dic)
    for (i, state_indices) in enumerate(values(MSS_dic))
        l = qlist[i]
        push!(rows, state_indices)
        push!(cols, fill(i, l))
        push!(vals, fill(1.0 / sqrt(l), l))
    end

    rows = vcat(rows...)
    cols = vcat(cols...)
    vals = vcat(vals...)

    iso_sparse = sparse(rows, cols, vals, num_states, num_categories)

    return iso_sparse
end
iso_K2MSS_sparse(N::Int, k::Int64, inv::Int64=1) = iso_K2MSS_sparse(BitStr{N, Int}, k, inv)

function iso_total2MSS_sparse(::Type{T}, k::Int64, inv::Int64=1) where {N, T <: BitStr{N}}
    # Function to map the total basis to the MSS space basis, k can only equal to 0 or N/2(pi)
    iso = iso_total2K_sparse(T, k) * iso_K2MSS_sparse(T, k, inv)

    return iso
end
iso_total2MSS_sparse(N::Int, k::Int64, inv::Int64=1) = iso_total2MSS_sparse(BitStr{N, Int}, k, inv)



