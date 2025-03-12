function PXP_Ham_sparse(::Type{T}, pbc::Bool=true) where {N, T <: BitStr{N}}
    # Generate Hamiltonian for PXP model, automotically contain pbc or obc
    basis=PXP_basis(T,pbc)

    l=length(basis)
    H=spzeros(Float64,(l,l))
    for i in 1:l
        output=actingH_PXP(T, basis[i], pbc) 
        for m in output 
            j=searchsortedfirst(basis,m)
            H[i, j] += 1
        end
    end

    return H
end

function PXP_K_Ham_sparse(::Type{T}, k::Int, Omega::Float64=1.0) where {N, T <: BitStr{N}}
    """
    :params: a int of lattice number, momentum of system and interaction strength of system which default to be 1
    :return: the Hamiltonian matrix in given K space
    """
    basisK, basis_dic = PXP_K_basis(T, k)
    l = length(basisK)
    omegak = exp(2im * π * k / N)
    H = spzeros(ComplexF64, (l, l))

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

function PXP_MSS_Ham_sparse(::Type{T}, k::Int, inv::Int64=1) where {N, T <: BitStr{N}}
    """
    :params: a int of lattice number, momentum of system and interaction strength of system which default to be 1
    :return: the Hamiltonian matrix in given maximum symmetry space
    """    
    omegak = exp(2im * π * k / N) 
    
    if inv==1
        MSS, MSS_dic, qlist = PXP_MSS_basis(T, k)

        l = length(MSS)
        H = spzeros(ComplexF64, (l, l))
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
                    H[i, j] += Zn / Zm*omegak^d 
                end
            end
        end
        H=real(H)
        H = (H + H') / 2
    else
        MSS, MSS_dic = PXP_MSS_basis(T, k, -1)

        l = length(MSS)
        H = spzeros(ComplexF64, (l, l))
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
PXP_MSS_Ham_sparse(N::Int64, k::Int) = PXP_MSS_Ham_sparse(BitStr{N, Int}, k)

function wf_time_evolution_mss(L::Int, k::Int64, psi0::Vector{ET}, t_values::Vector{Float64}) where {ET}
    # 预先计算时间步长
    dt = t_values[2] - t_values[1]
    H = PXP_MSS_Ham_sparse(L, k)
    # 初始化存储波函数的数组
    @assert length(psi0) == size(H, 1) "state length is expected to be $(size(H, 1)), but got $(length(psi0))"
    wavefunctions = Vector{Vector{ComplexF64}}(undef, length(t_values))
    wavefunctions[1] = copy(psi0)
        
    initial_norm = norm(psi0)
    if !isapprox(initial_norm, 1.0, atol=1e-10)
        @warn "The initial state isn't normalized, with norm: $initial_norm"
        psi0 = psi0 / initial_norm
    end
        
    psi_t = copy(psi0)
    for i in 2:length(t_values)
        # 使用 expv 计算波函数的时间演化
        psi_t = expv(-im * dt, H, psi_t; ishermitian=true)   
        wavefunctions[i] = psi_t
    end
    
    return wavefunctions
end

function ergotropy_PXP_MSS_state(L::Int, l::Int, state::Vector{T}, k::Int=0) where T
    HA = PXP_Ham(l, false)
    subenergy, substates = eigen(HA)
    
    subrho = rdm_PXP_MSS(L, collect(1:l), state, k)
    
    GS_energy = tr(subrho * HA)
    spectrum = eigvals(subrho)
    sorted_spectrum = sort(spectrum, rev=true)
    passive_energy = dot(sorted_spectrum, subenergy)
    
    return GS_energy, subenergy[1], passive_energy
end