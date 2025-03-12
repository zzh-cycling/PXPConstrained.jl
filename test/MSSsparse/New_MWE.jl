using Test
using Yao
using SparseArrays
using ManyBodyScar
using LinearAlgebra
using Plots
using LaTeXStrings

#writing test to get the correct sparse Hamiltonian
@testset "ishermitian for PXP_MSS_Ham_sparse" begin
    for L in 4:28
        @test ishermitian(PXP_MSS_Ham_sparse(L, 0))
        @test ishermitian(PXP_MSS_Ham_sparse(L, 12))
    end
end


@testset "fibonacci_chain_pbc" begin
    for L in 4:28
        @test Set(fibonacci_chain_pbc_test(L)) == Set(fbcp(L))
    end
end


# function test_equivalence(n::Int)
#     for n in 2:n  
#         result1 = fibonacci_chain_pbc(n)
#         result2 = fbcp(n)
        
#         result1_int = parse.(Int, result1, base=2)
        
#         if sort(result1_int) == sort(result2)
#             println("n = $n: equivalent")
#         else
#             println("n = $n: not equivalent")
#             println("fibonacci_chain_pbc: ", result1_int)
#             println("fbcp: ", result2)
#         end
#     end
# end

# test_equivalence(10)

# @testset "PXP_MSS_Ham_sparse" begin
#     @test Set(PXP_MSS_Ham(24, 0)) == Set(PXP_MSS_Ham_sparse(24, 0))
#     @test Set(PXP_MSS_Ham(24, 12)) == Set(PXP_MSS_Ham_sparse(24, 12))
#     @test Set(PXP_MSS_Ham(18, 12)) == Set(PXP_MSS_Ham_sparse(18, 12))
#     @test Set(PXP_MSS_Ham(18, 0)) == Set(PXP_MSS_Ham_sparse(18, 0))
# end


@testset "reverse bits" begin
    @test ManyBodyScar.reverse_bits(Int(0b10100), 5) == 5
    @test ManyBodyScar.reverse_bits(1<<39, 50) == 1<<10
end

@testset "get_representative_sparse" begin
    @test ManyBodyScar.cyclebits_sparse(34, Int(typemax(Int32)), 36) > 0
end

@testset "the normalization of neel_state_mss" begin
    for L in 4:32
        @test norm(neel_state_mss(L, 0)) == 1
        @test norm(neel_state_mss(L, L÷2)) == 1
    end
end

@testset "PXP_MSS_Ham_sparse size" begin
    for L in 4:28
        println("L = $L")
        # Test for correct dimensions
        @test size(PXP_MSS_Ham_sparse(L, 0)) == (length(PXP_MSS_basis_sparse(L, 0)[1]), length(PXP_MSS_basis_sparse(L, 0)[1]))
        @test size(PXP_MSS_Ham_sparse(L, L÷2)) == (length(PXP_MSS_basis_sparse(L, L÷2)[1]), length(PXP_MSS_basis_sparse(L, L÷2)[1]))
    end
end

# @testset "time_evolution_mss" begin
#     @test norm(time_evolution_mss(16, collect(0:0.1:20))[1]) ≈ 1
#     @test norm(time_evolution_mss(16, collect(0:0.1:20))[end]) ≈ 1
# end


# g=time_evolution_mss(18, collect(0:0.1:20));
# echo=zeros(201)
# for i in 1:201
#     echo[i]=norm(g[i][end])
# end

# Plots.plot(echo)

#目前需要解决sparse版本实现与dense版本实现的差异问题(已解决)
#sparse1 has some problems with respect to the eigenstate,sparse3 is the correct version
h = PXP_MSS_Ham_sparse(24,0)
#sparsity indicates the degree of connectivity/interaction between different quantum states in the system
sparsity = count(!iszero, h)/length(h)
energy, states=eigen(Matrix(h))
Plots.scatter(energy, log.(states[end,:].^2),color=:orange,ylim=(-17,1))
Plots.scatter(energy, log.(states[:,end].^2),ylim=(-17,1));


# visualize the connectivity of the Hamiltonian matrix
#the input need to be tuned manually
function visualize_connectivity(h)
    nonzero_vals = filter(!iszero, h[:])
    # 使用 spy 绘制非零元素的位置
    spy(h, 
        markersize = 2,
        title = L"Hamiltonian\ Matrix\ Connectivity\ Pattern\ (N=24)",
        xlabel = L"state\ index\ j",
        ylabel = L"state\ index\ i",
        color = :plasma,
        colorbar = true,  # 添加颜色条
        markerstrokewidth = 0,  # 移除标记边框
        clims = (minimum(abs.(nonzero_vals)), maximum(abs.(nonzero_vals)))  # 设置颜色范围
    )
end

visualize_connectivity(h)
savefig("ScarHconnectivity(N=26).pdf")

# h = PXP_MSS_Ham(24,0)
# energy, states=eigen(Matrix(h))
# using Plots
# Plots.scatter(energy, log.(states[end,:].^2),ylim=(-17,1))
# Plots.scatter(energy, log.(states[:,end].^2))
# histogram(energy,bins=100)
g = PXP_MSS_Ham(24,0)
energy2, states2=eigen(g)
Plots.scatter(energy2, log.(states2[end,:].^2),color=:red,ylim=(-17,1))
Plots.scatter(energy2, log.(states2[:,end].^2))


function neel_state_mss(L::Int, k::Int)
    # Generate the Neel state in the MSS basis
    # params L: System size
    # params k: Momentum number
    # return: Neel state vector in the MSS basis
    
    # get MSS basis
    MSS_int, MSS_dic, qlist = PXP_MSS_basis_sparse(L, k)
    
    # directly generate the integer representation of the neel state
    neel_int = 0
    for i in 0:(L÷2 - 1)
        neel_int |= (1 << (2*i))
    end
    
    # find the representative of the neel state
    neel_rep, _ = get_representative_sparse(L, neel_int)
    neel_inv_rep = get_representative_sparse(L, reverse_bits(neel_rep, L))[1]
    neel_mss_rep = min(neel_rep, neel_inv_rep)
    
    # find the index of the neel state in the MSS basis
    idx = findfirst(x -> x == neel_mss_rep, MSS_int)
    
    # create the state vector
    state = zeros(ComplexF64, length(MSS_int))
    if idx !== nothing
        # for k=0 case, we only need to set a suitable amplitude
        # in the MSS basis, the translation invariance is already considered
        # so we only need to set a normalized amplitude
        state[idx] = 1.0
        
        println("State construction details:")
        println("k = ", k)
        println("MSS basis size = ", length(MSS_int))
        # println("Representative state index = ", idx)
        # println("Final normalization = ", sqrt(sum(abs2.(state))))
    end
    
    return state
end


function time_evolution_mss(L::Int, t_values::Vector{Float64})
    # 预先计算时间步长
    dt = t_values[2] - t_values[1]
    psi0 = neel_state_mss(L, 0)
    H = PXP_MSS_Ham_sparse(L, 0)
    # 初始化存储波函数的数组
    wavefunctions = Vector{Vector{ComplexF64}}(undef, length(t_values))
    wavefunctions[1] = copy(psi0)
        
    initial_norm = norm(psi0)
    if !isapprox(initial_norm, 1.0, atol=1e-10)
        @warn "初始态未归一化! 初始范数: $initial_norm"
        psi0 = psi0 / initial_norm
    end
        
    psi_t = copy(psi0)
    for i in 2:length(t_values)
        # 使用 expv 计算波函数的时间演化
        psi_t = expv(-im * dt, H, psi_t; ishermitian=true)
        
        # 检查演化后的归一化
        current_norm = norm(psi_t)
        if !isapprox(current_norm, 1.0, atol=1e-10)
            @warn "时间步 $i 的波函数未归一化! 当前范数: $current_norm"
            psi_t = psi_t / current_norm 
        end
                   
        wavefunctions[i] = copy(psi_t)
    end
    
    return wavefunctions
end

function time_evolution_mss_rotation(L::Int, t_values::Vector{Float64}, θ::Real)
    dt = t_values[2] - t_values[1]
    psi0 = neel_state_mss(L, 0)
    
    # Apply rotation to the initial state
    rotation_matrix = exp(-im * θ / 2 * mat(Y))
    psi0 = rotation_matrix * psi0
    
    H = PXP_MSS_Ham_sparse(L, 0)
    wavefunctions = Vector{Vector{ComplexF64}}(undef, length(t_values))
    wavefunctions[1] = copy(psi0)
    
    # Initialize the wavefunction for time evolution
    psi_t = copy(psi0)
    for i in 2:length(t_values)
        # Use expv to compute the time evolution of the wavefunction
        psi_t = expv(-im * dt, H, psi_t; ishermitian=true)
        
        # Check normalization of the evolved wavefunction
        current_norm = norm(psi_t)
        if !isapprox(current_norm, 1.0, atol=1e-10)
            @warn "Wavefunction at time step $i is not normalized! Current norm: $current_norm"
            psi_t = psi_t / current_norm 
        end
        
        wavefunctions[i] = copy(psi_t)
    end
    
    return wavefunctions
end

time_evolution_mss_rotation(24, collect(0:0.1:20), π/2)
