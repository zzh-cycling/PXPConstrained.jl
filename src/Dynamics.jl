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

function rotated_psi_state(::Type{T}, θ::Real) where {N, T<: BitStr{N}}
#params: the particlenumber of the space, and rotation angle θ for the Z2 state
#return: the state rotated by on site rotation exp(i θ/2 Y)
    basis = PXP_basis(T)
    γ = tan(θ/2)
    rotated_state = zeros(length(basis))
    
    for (i, base) in enumerate(basis)
        rotated_state[i] = even_zeros(base, γ)
    end
    
    return rotated_state .* cos(θ/2)^N
end
rotated_psi_state(N::Int64, θ::Real) = rotated_psi_state(BitStr{N, Int}, θ)

function even_zeros(base::BitStr{N}, γ::Float64) where {N}
    # 计算奇数位上1的数量
    odd_ones = 0
    # 计算偶数位上0的数量
    even_zeros = 0
    
    for j in 1:N
        if j % 2 == 1  # 奇数位
            if base[j] == 1
                odd_ones += 1
            end
        else  # 偶数位
            if base[j] == 0
                even_zeros += 1
            end
        end
    end
    
    # 计算振幅：γ^(奇数位1的数量) * (-γ)^(偶数位0的数量)
    return γ^(odd_ones + even_zeros) * (-1)^even_zeros
end

function even_ones(base::BitStr{N}, γ::Float64) where {N}
    # 计算奇数位上1的数量
    odd_zeros = 0
    # 计算偶数位上0的数量
    even_ones = 0
    
    for j in 1:N
        if j % 2 == 1  # 奇数位
            if base[j] == 0
                odd_zeros += 1
            end
        else  # 偶数位
            if base[j] == 1
                even_ones += 1
            end
        end
    end
    
    # 计算振幅：γ^(奇数位1的数量 + 偶数位0的数量) * (-1)^(偶数位0的数量)
    return γ^(odd_zeros + even_ones) * (-1)^even_ones
end


function rotated_psi_state_mss(::Type{T}, k::Int64, θ::Real) where {N, T<: BitStr{N}}
    # params: a state in maximum symmetry space, and the momentum of the state
    # return: the state in total space
    MSS, MSS_dic, qlist = PXP_MSS_basis(T, k)
    γ = tan(θ/2)
    basisK, k_dic = PXP_K_basis(T, k)
    
    rotated_state = zeros(Float64, length(MSS))

    for (i, base) in enumerate(MSS)
        Y= sqrt(length(k_dic[base]))/N
        Z= sqrt(qlist[i])*Y/2
        @show base,length(k_dic[base]), Y, Z, qlist[i]
        amp1 = even_zeros(base, γ)
        amp2 = even_ones(base, γ)
        
        rotated_state[i]=Z*N*(amp1+amp2)
    end
    
    return rotated_state.* cos(θ/2)^N
end
rotated_psi_state_mss(N::Int64, k::Int64, θ::Real) = rotated_psi_state_mss(BitStr{N, Int}, k, θ)
