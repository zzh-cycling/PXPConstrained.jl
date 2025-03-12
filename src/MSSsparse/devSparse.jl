#dev is cooking hard!
#construct the sparse Hamiltonian matrix for the PXP model and calculate the dynamics evolution
using Test
using Yao
using SparseArrays
using LinearAlgebra
using Plots
using LaTeXStrings
using ExponentialUtilities

function fbco(n::Int)
    n == 1 && return [0, 1]
    subconfigs = fbco(n - 1)
    return append!(subconfigs, [c + 1<<(n-1) for c in subconfigs if iszero(c >> (n-2))])
end

# Apply the periodic boundary condition to the Fibonacci chain
function fbcp(n::Int)
    return filter(c -> iszero((c >> (n-1)) & (c & 1)), fbco(n))
end

function apply_pxp_flip_sparse(N::Int64, state::Int, pbc::Bool=true)
    output = Int[]
    # Handle flips for middle positions (from position 1 to N-2)
    for i in 0:N-3
        # Check if surrounded by 0s
        if ((state >> i) & 1) == 0 && ((state >> (i+2)) & 1) == 0
            # Flip middle bit
            flipped = xor(state, 1 << (i+1))
            push!(output, flipped)
        end
    end

    if pbc  # Periodic boundary conditions
        # Handle flips for last position
        if ((state >> (N-2)) & 1) == 0 && (state & 1) == 0
            flipped = xor(state, 1 << (N-1))
            push!(output, flipped)
        end
        # Handle flips for first position
        if ((state >> (N-1)) & 1) == 0 && ((state >> 1) & 1) == 0
            flipped = xor(state, 1)
            push!(output, flipped)
        end
    else  # Open boundary
        # Handle flips for first position
        if ((state >> 1) & 1) == 0
            flipped = xor(state, 1)
            push!(output, flipped)
        end
        # Handle flips for last position
        if ((state >> (N-2)) & 1) == 0
            flipped = xor(state, 1 << (N-1))
            push!(output, flipped)
        end
    end
    return output
end


# Binary cyclic shift function for quantum state translation
function cyclebits_sparse(L::Int, state::Int64, n_translations::Int)
    n_translations = n_translations % L  
    mask = (1 << L) - 1  # Generate an all-1s mask to clear the highest bit and keep the last L bits
                         # Ensures the result only retains L bits, enabling cyclic shift
    return ((state << n_translations) & mask) | (state >> (L - n_translations)) # Cycle overflow bits to the right
end

function get_representative_sparse(L::Int, state::Int)
    # Initialize representative state as input state
    representative = state
    # Initialize translation amount to 0
    translation = 0
    
    # Try all possible translations (0 to L-1 positions)
    for n_translation_sites in 0:L-1
        # Use cyclebits to translate the state
        new_state = cyclebits_sparse(L,state, n_translation_sites)
        # If a smaller state is found, update representative and translation amount
        if new_state < representative
            representative = new_state
            translation = n_translation_sites
        end
    end
    # Return the translation-equivalent state and corresponding translation amount
    return representative, translation
end

function PXP_K_basis_sparse(L::Int, k::Int)
    basisK_int = Vector{Int}(undef, 0)
    
    # Generate basis states using periodic boundary conditions
    basis_int = fbcp(L)
    
    # Create a dictionary to store representative states and their translations
    # Key: representative state, Value: array of all translations of this state
    basis_dic = Dict{Int, Vector{Int}}()
    for i in basis_int
        category = get_representative_sparse(L, i)[1]
        if haskey(basis_dic, category)
            push!(basis_dic[category], i)
        else
            basis_dic[category] = [i]
        end
    end

    # Filter states that satisfy momentum conditions:
    # 1. State must be its own representative
    # 2. k * orbit_length must be divisible by L
    for i in basis_int
        RS = get_representative_sparse(L, i)[1]
        if RS == i && (k * length(basis_dic[RS])) % L == 0
            push!(basisK_int, i)
        end
    end

    return basisK_int, basis_dic
end

# Reverse the bits of a binary number
# Example: for L=4, 1010 becomes 0101
function reverse_bits(n::Int, L::Int)
    result = 0
    for i in 0:L-1
        if (n & (1 << i)) != 0
            result |= (1 << (L-1-i))
        end
    end
    return result
end

function PXP_MSS_basis_sparse(L::Int, k::Int)
    MSS_int = Vector{Int}(undef, 0)
    
    # Get momentum basis states and their translation dictionary
    basisK_int, basis_dic = PXP_K_basis_sparse(L, k)
    
    # Dictionary to store MSS states and their translations
    MSS_dic = Dict{Int, Vector{Int}}()
    # List to store the number of unique states in each symmetry sector
    qlist = Vector{Int}(undef, 0)
    
    # Iterate through momentum basis states to construct MSS basis
    for i in eachindex(basisK_int)
        n = basisK_int[i]
        
        # Get spatial inversion of state n
        reversed_n = reverse_bits(n, L)
        # Find representative of inverted state
        nR = get_representative_sparse(L, reversed_n)[1]
        
        # Include state if it's the smallest among itself and its spatial inversion
        if n <= min(nR, n)
            # @show n
            push!(MSS_int, n)
            MSS_dic[n] = basis_dic[n]
            # Store number of unique states (1 if state = its inversion, 2 otherwise)
            push!(qlist, length(Set([n, nR])))
        end
    end

    return MSS_int, MSS_dic, qlist
end
#如果两个态属于同一个key的value数组，说明它们是平移等价的
#value数组的长度告诉我们这个态有多少个平移等价态 
#for MSS_Dic,key：代表态（最小的那个等价态,value：所有通过平移操作可以得到的等价态（包括代表态本身）

#correct sparse Hamiltonian
function PXP_MSS_Ham_sparse2(L::Int, k::Int, Omega::Float64=1.0)
    MSS_int, MSS_dic, qlist = PXP_MSS_basis_sparse(L, k)
    l = length(MSS_int)
    matrix_elements = Dict{Tuple{Int,Int}, Float64}()
    
    for i in 1:l
        n = MSS_int[i]
        Zn = sqrt(qlist[i]) / 4 * sqrt(length(MSS_dic[n])) / L
        
        output = apply_pxp_flip_sparse(L, n, true)
        
        for m in output
            mbar, d = get_representative_sparse(L, m)
            inv_mbar = get_representative_sparse(L, reverse_bits(mbar, L))[1]
            mtilde = min(mbar, inv_mbar)
            j = findfirst(x -> x == mtilde, MSS_int)
            
            if j !== nothing
                Zm = sqrt(qlist[j]) / 4 * sqrt(length(MSS_dic[mtilde])) / L
                value = Omega * Zn / Zm / 2  # 在这里除以2来修正
                
                idx_pair = minmax(i, j)
                matrix_elements[idx_pair] = get(matrix_elements, idx_pair, 0.0) + value
            end
        end
    end
    
    I, J, V = Int[], Int[], Float64[]
    for ((i, j), value) in matrix_elements
        push!(I, i); push!(J, j); push!(V, value)
        if i != j  # 确保厄米性
            push!(I, j); push!(J, i); push!(V, value)
        end
    end
    
    return sparse(I, J, V, l, l)
end


h = PXP_MSS_Ham_sparse(20,0)
#sparsity indicates the degree of connectivity/interaction between different quantum states in the system
sparsity = count(!iszero, h)/length(h)
energy, states=eigen(Matrix(h))
Plots.scatter(energy, log.(states[end,:].^2),color=:orange)
Plots.scatter(energy, log.(states[:,end].^2),ylim=(-17,1));

# visualize the connectivity of the Hamiltonian matrix
#the input need to be tuned manually
function visualize_connectivity(h)
    nonzero_vals = filter(!iszero, h[:])
    # 使用 spy 绘制非零元素的位置
    spy(h, 
        markersize = 2,
        title = L"Hamiltonian\ Matrix\ Connectivity\ Pattern\ (N=20)",
        xlabel = L"state\ index\ j",
        ylabel = L"state\ index\ i",
        color = :plasma,
        colorbar = true,  
        markerstrokewidth = 0,  # 移除标记边框
        clims = (minimum(abs.(nonzero_vals)), maximum(abs.(nonzero_vals)))  # 设置颜色范围
    )
end

visualize_connectivity(h)
savefig("connectivity_pattern_N=20(sparse).pdf")

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

function rotation_neel(L::Int, k::Int, θ::Real)
    # 获取MSS基底维度
    MSS_int, _, _ = PXP_MSS_basis_sparse(L, k)
    dim = length(MSS_int)
    
    # 构建旋转矩阵中的元素
    reverse_matrix = [
        if i == dim - j + 1
            i > j ? 1.0 : -1.0  # 如果在反对角线上，根据位置决定符号
        else
            0.0  # 非反对角线位置为0
        end
        for i in 1:dim, j in 1:dim
    ]
    
    # 构建旋转矩阵：对角线上是cos(θ/2)，反对角线上是sin(θ/2)
    rotation_matrix = cos(θ/2) * Matrix{ComplexF64}(I(dim)) + sin(θ/2) * reverse_matrix
    
    return rotation_matrix * neel_state_mss(L, k)
end

b = rotation_neel(20,0,0)
norm(b)
a == b
a ≈ b

h = PXP_MSS_Ham_sparse(20,0)
f = Matrix(h)
energy2, states2=eigen(f)
Plots.scatter(energy2, log.(states2[end,:].^2),color=:red,ylim=(-17,1))

# 计算旋转Neel态与本征态的重叠
h = PXP_MSS_Ham_sparse(20, 0)
energy, states = eigen(Matrix(h))
rotated_neel = rotation_neel(20, 0, π/3)  # 可以调整旋转角度
rotated_neel2 = rotation_neel(20, 0, π/4)
rotated_neel3 = rotation_neel(20, 0, π/2)
rotated_neel4 = rotation_neel(20, 0, π)
rotated_neel5 = rotation_neel(20, 0, 2π)
norm(rotated_neel2)
rotated_neel6 = rotation_neel(20, 0, 0)
overlaps_0 = abs2.(adjoint(states) * rotated_neel6)
overlaps_2pi = abs2.(adjoint(states) * rotated_neel5)
overlaps_0 ≈ overlaps_2pi
#看rotation后overlap的差异


# 绘制重叠图
Plots.scatter(energy, log.(overlaps), 
    color=:pink,
    ylim=(-17,1),
    xlabel=L"Energy\ (E_n)",
    ylabel=L"\log(|\langle E_n|\psi_{\theta}\rangle|^2)",
    title=L"Rotated\ Neel\ State\ Overlaps\ (\theta=0)",
    legend=false
)

# #计算overlap
# function calculate_overlaps(L::Int, k::Int, thetas::Vector{Float64})
#     # 获取哈密顿量和本征态
#     h = PXP_MSS_Ham_sparse(L, k)
#     energy, eigenstates = eigen(Matrix(h))
    
#     # 初始化存储重叠的矩阵
#     # 行对应不同的theta值，列对应不同的本征态
#     overlaps = zeros(ComplexF64, length(thetas), length(energy))
    
#     for (i, theta) in enumerate(thetas)
#         # 计算旋转后的Neel态
#         psi_theta = rotation_neel(L, k, theta)
        
#         # 计算与每个本征态的重叠
#         for j in 1:length(energy)
#             overlaps[i,j] = dot(eigenstates[:,j], psi_theta)
#         end
#     end
    
#     return energy, overlaps
# end

# # 绘制重叠图
# function plot_overlaps(L::Int, k::Int, thetas::AbstractVector{<:Real})
#     energy, overlaps = calculate_overlaps(L, k, collect(thetas))
    
#     # 创建热图
#     heatmap(energy, thetas, abs2.(overlaps),
#         xlabel = L"Energy\ (E_n)",
#         ylabel = L"\theta",
#         title = L"|\langle E_n|\psi_0(\theta)\rangle|^2",
#         color = :viridis,
#         colorbar_title = L"Overlap",
#     )
# end

# thetas = range(0, 2π, length=100)
# plot_overlaps(20, 0, thetas)
# savefig("overlaps_heatmap.pdf")




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
    # Generate Fibonacci chain for PXP model with periodic boundary condition
    chain = Fibonacci_chain_OBC(T)
    filtered_fib_chain = T[]
    # 使用位操作直接创建掩码，而不是使用bmask函数
    mask = (1 << 1) | (1 << N)  # 创建第1位和第N位为1的掩码
    for s in chain
        # 检查s的第1位和第N位是否同时为1
        if (Int(s) & mask) == mask
            continue
        else
            push!(filtered_fib_chain, s)
        end
    end
    return filtered_fib_chain
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

function iso_total2K(::Type{T}, k::Int64) where {N, T <: BitStr{N}}
#Function to map the total basis to the K space basis, actually is the isometry, defined as W'*W=I, W*W'=P, P^2=P
    basis = PXP_basis(T)

    k_dic = Dict{Int, Vector{Int64}}()

    for i in eachindex(basis)
        state=basis[i]
        category = get_representative_sparse(L, state)[1]
        if haskey(k_dic, category)
            push!(k_dic[category], i)
        else
            k_dic[category] = [i]
        end
    end
    
    iso = zeros((length(basis), length(keys(k_dic))))
    
    for (i, state_index) in enumerate(values(k_dic))
        l = length(state_index)
        iso[state_index, i] .= 1/sqrt(l)
    end

    return iso
end

function cyclebits_sparse(L::Int, state::Int64, n_translations::Int)
    n_translations = n_translations % L  
    mask = (1 << L) - 1  # Generate an all-1s mask to clear the highest bit and keep the last L bits
                         # Ensures the result only retains L bits, enabling cyclic shift
    return ((state << n_translations) & mask) | (state >> (L - n_translations)) # Cycle overflow bits to the right
end

function get_representative_sparse(L::Int, state::Int)
    # Initialize representative state as input state
    representative = state
    # Initialize translation amount to 0
    translation = 0
    
    # Try all possible translations (0 to L-1 positions)
    for n_translation_sites in 0:L-1
        # Use cyclebits to translate the state
        new_state = cyclebits_sparse(L,state, n_translation_sites)
        # If a smaller state is found, update representative and translation amount
        if new_state < representative
            representative = new_state
            translation = n_translation_sites
        end
    end
    # Return the translation-equivalent state and corresponding translation amount
    return representative, translation
end

function rdm_PXP_K(::Type{T}, subsystems::Vector{Int64}, state::Vector{Float64}, k::Int64) where {N, T <: BitStr{N}}

    state=iso_total2K(T,k)*state
    reduced_dm = rdm_PXP(T, subsystems, state)
    return reduced_dm
end

function iso_K2MSS(::Type{T}, k::Int64, inv::Int64=1) where {N, T <: BitStr{N}}
#Function to map the MSS basis to the K space basis
    basis = PXP_basis(T)
    basisK, k_dic = PXP_K_basis(T, k)

    MSS_dic = Dict{Int, Vector{Int64}}()

    # Below procedure is to collapse the extra basis in K space that can be converted mutually to MSS space.
    if inv==1
        for i in eachindex(basisK)
            n = basisK[i]
            # here we calculate the representative state of the inversion of n
            nR = get_representative(breflect(n))[1]
            if n <= min(nR, n)
                if haskey(MSS_dic, nR)
                    push!(MSS_dic[nR], i)
                else
                    MSS_dic[nR] = [i]
                end
            end
            
        end
    else
        for i in eachindex(basisK)
            n = basisK[i]
            nR = get_representative(breflect(n))[1]
            if n != nR
                if n <= min(nR, n)
                    MSS_dic[n] = [i]
                end
            end
        end    
        
    
    end

    iso = zeros((length(basisK), length(MSS_dic)))

    for (i, state_index) in enumerate(values(MSS_dic))
        l = length(state_index)
        iso[state_index, i] .= 1/sqrt(l)
    end

    return iso
end

function iso_total2MSS(::Type{T}, k::Int64, inv::Int64=1) where {N, T <: BitStr{N}}
#Function to map the total basis to the MSS space basis, k can only equal to 0 or N/2(pi)
    iso = iso_total2K(T, k) * iso_K2MSS(T, k, inv)

    return iso
end

function rdm_PXP_MSS(::Type{T}, subsystems::Vector{Int64}, state::Vector{Float64}, k::Int64, inv::Int64=1) where {N, T <: BitStr{N}}

    state=iso_total2MSS(T, k, inv)*state
    reduced_dm = rdm_PXP(T, subsystems, state)
    return reduced_dm
end


function dynamics_ergotropy(N::Int64, l::Int64, wf_tlis::Vector{Vector{T}}) where {T <: Real}
    step=length(wf_tlis)
    ergotropy_lis=Vector{Float64}(undef,step)
    bound_energy_lis=Vector{Float64}(undef,step)
    DeltaE_lis=Vector{Float64}(undef,step)

    HA=PXP_Ham(l,false)
    sub_basis=PXP_basis(l,false)
    energy, states= eigen(PXP_Ham(N))
    subenergy, substates= eigen(HA)

    for (i,wf) in enumerate(wf_tlis)
        subscarrho, reduced_basis=rdm_PXP(N,l,wf) 
        println(tr(subscarrho))
        indices = findall(x -> !(x in sub_basis[2]), reduced_basis)
        subscarrho=delete_rowcolumn(subscarrho,indices) #Here we can also truncate the dim of reduce_rho(which with additional basis)
        GS_energy=tr(subscarrho*HA)
    
        spectrum=eigvals(subscarrho)
        sorted_spectrum=sort(spectrum, rev=true)
        passive_energy=dot(sorted_spectrum, subenergy)
        
        ergotropy=GS_energy-passive_energy
        bound_energy=passive_energy-subenergy[1]
        DeltaE=GS_energy-subenergy[1]

        ergotropy_lis[i]=ergotropy
        bound_energy_lis[i]=bound_energy
        DeltaE_lis[i]=DeltaE
    end

    return ergotropy_lis, bound_energy_lis, DeltaE_lis
end



plot(tlis, ergotropy_lis, label=L"W_A", marker=:circle, markersize=1.5, xlabel=L"t", ylabel=L"E")
plot!(tlis,bound_energy_lis, marker=:x, label=L"Q_A",markersize=1.5)
plot!(tlis,DeltaE_lis, marker=:+, label=L"\Delta E_A",markersize=1.5,color="green")


function ergotropy_PXP_idx(::Type{T}, l::Int64, idx::Int64) where {N, T <: BitStr{N}}
    HA=PXP_Ham(BitStr{l, Int},false)
    sub_basis=PXP_basis(BitStr{l, Int},false)
    energy, states= eigen(PXP_Ham(T))
    subenergy, substates = eigen(HA)
    
    state = states[:,idx]
    subrho = rdm_PXP(T,collect(1:l), state) 
    GS_energy=tr(subrho*HA)

    spectrum=eigvals(subrho)
    sorted_spectrum=sort(spectrum, rev=true)
    passive_energy=dot(sorted_spectrum, subenergy)

    return GS_energy, subenergy[1], passive_energy
end

function ergotropy_PXP_state(::Type{T}, l::Int64, state::Vector{ET}) where {N, T <: BitStr{N}, ET}
    HA=PXP_Ham(BitStr{l, Int},false)
    sub_basis=PXP_basis(BitStr{l, Int},false)
    subenergy, substates= eigen(HA)
    subrho = rdm_PXP(T,collect(1:l), state) 

    GS_energy=tr(subrho*HA)
    spectrum=eigvals(subrho)
    sorted_spectrum=sort(spectrum, rev=true)
    passive_energy=dot(sorted_spectrum, subenergy)

    return GS_energy, subenergy[1], passive_energy
end



#MSS下的ergotropy
function ergotropy_PXP_MSS_idx(L::Int, l::Int, idx::Int, k::Int=0)
    # 构建子系统哈密顿量
    HA = PXP_MSS_Ham_sparse(l, k)
    # 获取全系统哈密顿量的本征值和本征态
    h = PXP_MSS_Ham_sparse(L, k)
    energy, states = eigen(Matrix(h))
    # 获取子系统哈密顿量的本征值
    subenergy, substates = eigen(Matrix(HA))
    
    # 获取指定索引的本征态
    state = states[:, idx]
    # 将复数状态转换为实数状态（如果状态是实数的）
    state_real = real.(state)
    
    # 计算约化密度矩阵
    subrho = rdm_PXP_MSS(DitStr{2, L, Int64}, collect(1:l), state_real, k)
    
    # 计算基态能量
    GS_energy = tr(subrho * Matrix(HA))
    
    # 计算被动能量
    spectrum = eigvals(subrho)
    sorted_spectrum = sort(spectrum, rev=true)
    passive_energy = dot(sorted_spectrum, subenergy)
    
    return GS_energy, subenergy[1], passive_energy
end

function ergotropy_PXP_MSS_state(L::Int, l::Int, state::Vector{T}, k::Int=0) where T
    # 构建子系统哈密顿量
    HA = PXP_MSS_Ham_sparse(l, k)
    # 获取子系统哈密顿量的本征值
    subenergy, substates = eigen(Matrix(HA))
    
    # 将复数状态转换为实数状态（如果状态是实数的）
    if T <: Complex
        state_real = real.(state)
    else
        state_real = state
    end
    
    # 计算约化密度矩阵
    subrho = rdm_PXP_MSS(DitStr{2, L, Int64}, collect(1:l), state_real, k)
    
    # 计算基态能量
    GS_energy = tr(subrho * Matrix(HA))
    
    # 计算被动能量
    spectrum = eigvals(subrho)
    sorted_spectrum = sort(spectrum, rev=true)
    passive_energy = dot(sorted_spectrum, subenergy)
    
    return GS_energy, subenergy[1], passive_energy
end

# 设置参数
L = 18  # 总系统大小
l = 6   # 子系统大小
tmax = 10.0
dt = 0.1
t_values = collect(0:dt:tmax)

# 计算初始态（Neel态）
psi0 = neel_state_mss(N, 0)

# 计算时间演化
wavefunctions = time_evolution_mss(N, t_values)

# 初始化存储结果的数组
ergotropy_lis = Vector{Float64}(undef, length(t_values))
bound_energy_lis = Vector{Float64}(undef, length(t_values))
DeltaE_lis = Vector{Float64}(undef, length(t_values))

# 对每个时间点计算ergotropy
for (i, t) in enumerate(t_values)
    # 获取当前时间点的波函数
    psi_t = wavefunctions[i]
    
    # 直接调用ergotropy_PXP_MSS_state计算能量
    GS_energy, min_energy, passive_energy = ergotropy_PXP_MSS_state(N, l, psi_t, 0)
    
    # 计算各种能量
    ergotropy_lis[i] = GS_energy - passive_energy
    bound_energy_lis[i] = passive_energy - min_energy
    DeltaE_lis[i] = GS_energy - min_energy
    
    println("完成时间点 $t ($(i)/$(length(t_values)))")
end
