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

function rdm_PXP_K(::Type{T}, subsystems::Vector{Int64},kstate::Vector{ET}, k::Int64) where {N,T <: BitStr{N}, ET}
    @assert length(kstate) == length(PXP_K_basis(T,k)[1]) "state length is expected to be $(length(PXP_K_basis(T, k)[1])), but got $(length(kstate))"
    state = mapstate_K2total(T, kstate, k)
    reduced_dm = rdm_PXP(T, subsystems, state)
    return reduced_dm
end
rdm_PXP_K(N::Int, subsystems::Vector{Int64},state::Vector{ET}, k::Int64) where {ET} = rdm_PXP_K(BitStr{N, Int}, subsystems, state, k)


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
iso_K2MSS(N::Int, k::Int64, inv::Int64=1) = iso_K2MSS(BitStr{N, Int}, k, inv)

function mapstate_MSS2K(::Type{T}, state::Vector{ET}, k::Int64, inv::Int64=1) where {N, T <: BitStr{N}, ET}
    @assert k == 0 || k==div(N,2) "k is expected to be 0 or $(div(N,2)), but got $k"
    @assert inv ==1 || inv==-1 "inv is expected to be 1 or -1, but got $(inv)"

    basisK, k_dic = PXP_K_basis(T, k)

    MSS_dic = Dict{Int, Vector{Int64}}()
    qlist = Vector{Int}(undef, 0)
   
    if inv==1
        for i in eachindex(basisK)
            n = basisK[i]
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

mapstate_MSS2total(::Type{T}, state::Vector{ET}, k::Int64, inv::Int64=1) where {N, T <: BitStr{N}, ET} = mapstate_K2total(T, mapstate_MSS2K(T, state, k, inv), k)
mapstate_MSS2total(N::Int64, state::Vector{ET}, k::Int64, inv::Int64=1) where {ET} = mapstate_MSS2total(BitStr{N, Int}, state, k, inv)

function iso_total2MSS(::Type{T}, k::Int64, inv::Int64=1) where {N, T <: BitStr{N}}
    # Function to map the total basis to the MSS space basis, k can only equal to 0 or N/2(pi)
    iso = iso_total2K(T, k) * iso_K2MSS(T, k, inv)

    return iso
end
iso_total2MSS(N::Int, k::Int64, inv::Int64=1) = iso_total2MSS(BitStr{N, Int}, k, inv)

function rdm_PXP_MSS(::Type{T}, subsystems::Vector{Int64}, mssstate::Vector{ET}, k::Int64, inv::Int64=1) where {N,T <: BitStr{N}, ET}
    @assert length(PXP_MSS_basis(T, k, inv)[1]) == length(mssstate) "state length is expected to be $(length(PXP_MSS_basis(T, k, inv)[1])), but got $(length(mssstate))"
    state=mapstate_MSS2total(T, mssstate, k, inv)
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


function PXP_MSS_basis(::Type{T}, k::Int64,inv::Int64=1) where {N, T <: BitStr{N}}
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
    if inv==1 && k==0 || inv==-1 && k==div(N,2)
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
        
    elseif inv==1 && k==div(N,2) || inv==-1 && k==0
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
            return MSS[index], new_MSS_dic, qlist
    end
          
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
    
    if inv==1 && k==0 || inv==-1 && k==div(N,2)
        MSS, MSS_dic, qlist = PXP_MSS_basis(T, k, inv)

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
    elseif inv==1 && k==div(N,2) || inv==-1 && k==0
        MSS, MSS_dic, _ = PXP_MSS_basis(T, k, inv)

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
PXP_MSS_Ham(N::Int, k::Int, inv::Int64=1) = PXP_MSS_Ham(BitStr{N, Int}, k, inv) 
