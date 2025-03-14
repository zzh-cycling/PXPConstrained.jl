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

# function rotated_psi_state(::Type{T}, θ::Float64) where {N, T<: BitStr{N}}
# #params: the particlenumber of the space, and rotation angle θ for the Z2 state
# #return: the state rotated by on site rotation exp(i θ/2 Y)
#     basis=PXP_basis(T)
#     γ = tan(θ/2)
#     rotated_state = zeros(length(basis))
#     for (i,base) in enumerate(basis)
#         # base #to be filled
#         # rotated_state[i]=base # to be filled
#         # 计算旋转后的态在基态上的振幅
#         amplitude = 1.0
#         for j in 1:N
#             if base[j] == 0  # 如果第j位是0
#                 # 对应于 |0⟩ 变为 |0⟩ + γ|1⟩
#                 if j % 2 == 1  # 奇数位
#                     amplitude *= 1.0
#                 else  # 偶数位
#                     amplitude *= -γ
#                 end
#             else  # 如果第j位是1
#                 # 对应于 |1⟩ 变为 -γ|0⟩ + |1⟩
#                 if j % 2 == 1  # 奇数位
#                     amplitude *= γ
#                 else  # 偶数位
#                     amplitude *= 1.0
#                 end
#             end
#         end
#         rotated_state[i] = amplitude
#     end
#     return rotated_state.*cos(θ/2)^N
# end
# rotated_psi_state(N::Int64, θ::Float64) = rotated_psi_state(BitStr{N, Int}, θ)

# function rotated_psi_state_mss(::Type{T},θ::Float64) where {N, T<: BitStr{N}}
# #params: a state in maximum symmetry space, and the momentum of the state
# #return: the state in total space
#     MSS, MSS_dic = PXP_MSS_basis(T, k)
#     γ = tan(θ/2)
#     rotated_state = zeros(length(MSS))
#     for (i,base) in enumerate(MSS)
#         # 计算旋转后的态在基态上的振幅
#         amplitude = 1.0
#         for j in 1:N
#             if base[j] == 0  # 如果第j位是0
#                 # 对应于 |0⟩ 变为 |0⟩ + γ|1⟩
#                 if j % 2 == 1  # 奇数位
#                     amplitude *= 1.0
#                 else  # 偶数位
#                     amplitude *= -γ
#                 end
#             else  # 如果第j位是1
#                 # 对应于 |1⟩ 变为 -γ|0⟩ + |1⟩
#                 if j % 2 == 1  # 奇数位
#                     amplitude *= γ
#                 else  # 偶数位
#                     amplitude *= 1.0
#                 end
#             end
#         end
#         rotated_state[i] = amplitude
#     end
#     return rotated_state.*cos(θ/2)^N
# end


function rotated_psi_state(::Type{T}, θ::Real) where {N, T<: BitStr{N}}
#params: the particlenumber of the space, and rotation angle θ for the Z2 state
#return: the state rotated by on site rotation exp(i θ/2 Y)
    basis = PXP_basis(T)
    γ = tan(θ/2)
    rotated_state = zeros(length(basis))
    
    for (i, base) in enumerate(basis)
        rotated_state[i] = amplitude_rotated(base, γ, N)
    end
    
    return rotated_state .* cos(θ/2)^N
end

function amplitude_rotated(base::BitStr{N}, γ::Float64, L::Int) where {N}
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
    
    # 计算振幅：γ^(奇数位1的数量 + 偶数位0的数量) * (-1)^(偶数位0的数量)
    return γ^(odd_ones + even_zeros) * (-1)^even_zeros
end

rotated_psi_state(N::Int64, θ::Real) = rotated_psi_state(BitStr{N, Int}, θ)

# function rotated_psi_state_mss(::Type{T}, k::Int64, θ::Real) where {N, T<: BitStr{N}}
# #params: a state in maximum symmetry space, and the momentum of the state
# #return: the state in total space
#     MSS, MSS_dic = PXP_MSS_basis(T, k)
#     γ = tan(θ/2)
#     rotated_state = zeros(length(MSS))
    
#     for (i, base) in enumerate(MSS)
#         rotated_state[i] = amplitude_rotated(base, γ, N)
#     end
    
#     return rotated_state .* cos(θ/2)^N
# end




function rotated_psi_state_mss(::Type{T}, k::Int64, θ::Real) where {N, T<: BitStr{N}}
    # params: a state in maximum symmetry space, and the momentum of the state
    # return: the state in total space
    MSS, MSS_dic = PXP_MSS_basis(T, k)
    γ = tan(θ/2)
    rotated_state = zeros(ComplexF64, length(MSS))
    
    # 对于k=0的情况，需要特殊处理
        for (i, base) in enumerate(MSS)
            # 计算原始振幅
            amp1 = amplitude_rotated(base, γ, N)
            
            # 计算反演后的振幅
            inverted_base = flip(base, bmask(T, 1:N))
            amp2 = amplitude_rotated(inverted_base, γ, N)
            
            # 根据公式，对称约束空间的振幅是两者的和
            rotated_state[i] = (amp1 + amp2) * cos(θ/2)^N
        end
        
        rotated_state = rotated_state 
    # else
    #     # 对于其他动量值，使用原始方法
    #     for (i, base) in enumerate(MSS)
    #         rotated_state[i] = amplitude_rotated(base, γ, N) * cos(θ/2)^N
    #     end
        
    #     # 归一化
    #     rotated_state = rotated_state / norm(rotated_state)
    # end
    
    return rotated_state/= norm(rotated_state)
end
rotated_psi_state_mss(N::Int64, k::Int64, θ::Real) = rotated_psi_state_mss(BitStr{N, Int}, k, θ)
