function wf_time_evolution(psi0::Vector{T}, times::Vector{Float64}, energy::Vector{Float64},states::Matrix{Float64}) where {T <: Real}
    wflis=Vector{Vector{ComplexF64}}(undef,length(times))
    c = states'*psi0
    exp_factors = [exp.(-1im * t * energy) for t in times]
    
    # Use multi-threading for parallel computation
    Threads.@threads for i in eachindex(times)
        wflis[i] = states * (c .* exp_factors[i])
    end
    return wflis
end

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

function rotated_psi_state(::Type{T}, θ::Float64) where {N, T<: BitStr{N}}
    """
    :params: the particlenumber of the space, and rotation angle θ for the Z2 state
    :return: the state rotated by on site rotation exp(i θ/2 Y)
    """

    basis=PXP_basis(T)
    γ = tan(θ/2)
    rotated_state = zeros(length(basis))
    for (i,base) in enumerate(basis)
        # base #to be filled
        # rotated_state[i]=base # to be filled
    end
    return rotated_state.*cos(θ/2)^N
end
rotated_psi_state(N::Int64, θ::Float64) = rotated_psi_state(BitStr{N, Int}, θ)

function rotated_psi_state_mss(::Type{T}, θ::Float64) where {N, T<: BitStr{N}}
    """
    :params: a state in maximum symmetry space, and the momentum of the state
    :return: the state in total space
    """
    MSS, MSS_dic = PXP_MSS_basis(T, k)
    γ = tan(θ/2)
    rotated_state = zeros(length(MSS))
    for (i,base) in enumerate(MSS)
        # base #to be filled
        # rotated_state[i]=base # to be filled
    end
    return rotated_state.*cos(θ/2)^N
    
end
rotated_psi_state_mss(N::Int64, k::Int64, θ::Float64) = rotated_psi_state_mss(BitStr{N, Int}, k, θ)