using PXPConstrained
using LinearAlgebra
using Plots
using LaTeXStrings

# Part 1: Calculate ergotropy and entanglement entropy for rotated states at different angles
# This section computes energy properties and entanglement entropy for a 16-qubit PXP model
# as a function of rotation angle θ
N=16
θlis=collect(0.0:0.1:2π)
ergolis=similar(θlis)
eelis=similar(θlis)
gslis=similar(θlis)
palsis=similar(θlis)
for i in eachindex(θlis)
    θ=θlis[i]
    rotated_state=rotated_psi_state(N, θ)
    GS_energy, subenergy, passive_energy=ergotropy_PXP_state(N ,div(N,2), rotated_state)
    # @show GS_energy, subenergy, passive_energy
    ergolis[i]= GS_energy-passive_energy  # Calculate ergotropy (maximum extractable work)
    subrho=rdm_PXP(N, collect(1:div(N,2)), rotated_state)  # Compute reduced density matrix
    eelis[i]=ee(subrho)  # Calculate entanglement entropy
    gslis[i]=GS_energy  # Store ground state energy
    palsis[i]=passive_energy  # Store passive energy
end

# Part 2: Plot the relationship between rotation angle, ergotropy and entanglement entropy
# This visualization shows how these quantum properties change with the rotation angle θ
fig = plot(θlis, ergolis, 
    title=L"Ergotropy\ and\ Entanglement\ Entropy\ (N=%$N)", 
    xlabel=L"\theta", 
    ylabel=L"Value", 
    label=L"Ergotropy", 
    c=:red,
    xlims=(0, π),
    xticks=(0:π/8:π, [L"0", L"\pi/4", L"\pi/2", L"3\pi/4", L"\pi", L"5\pi/4", L"3\pi/2", L"7\pi/4", L"2\pi"]))
plot!(θlis, eelis, 
    label=L"Entanglement\ Entropy", 
    c=:purple)

#plot entanglement entropy
# plot(subplot =2, θlis, eelis, title="Eigenenergy of rotated state", xlabel="θ", ylabel="Eigenenergy", legend=:outertopright)

# Part 3: Time evolution analysis of entanglement entropy
# Here we compute the time evolution of a specific rotated state (θ=π/2)
# and track how its entanglement entropy changes over time
energy, states = eigen(PXP_Ham(N))  # Diagonalize the PXP Hamiltonian
timelis=collect(0.0:0.2:200)  # Create time points from 0 to 200 with step 0.2
st=rotated_psi_state(N, π/2)  # Initialize a rotated state at θ=π/2
wflis= wf_time_evolution(st, timelis, energy, states)  # Perform time evolution
steelis=zeros(length(timelis))  # Array to store entanglement entropy values

for i in eachindex(timelis)
    subrho=rdm_PXP(N, collect(1:div(N,2)), wflis[i])  # Calculate reduced density matrix at each time point
    steelis[i]=ee(subrho)  # Compute entanglement entropy at each time point
end

# Part 4: Visualization of time evolution and Page curve for the final state
# Plot time evolution of entanglement entropy and analyze scaling properties of the final state
fig2=plot(timelis, steelis, title="Time evolution of z2 component", xlabel="Time", ylabel="Eigenenergy", label="z2 component", c=:red)

finalst=wflis[end]  # Get the final state from time evolution
figPage=ee_PXP_scaling_fig(N, finalst, "Page")[2]  # Generate Page scaling figure for final state




#new
# 计算N=16时，旋转pi/2后，演化(0.0:0.2:200)的ergotropy和ee，画在一张图上
N = 16
energy, states = eigen(PXP_Ham(N))  # 对角化PXP哈密顿量
timelis = collect(0.0:0.2:200)  # 创建时间点，从0到200，步长为0.2
st = rotated_psi_state(N, 0)  # 初始化旋转角度为π/2的量子态

# 执行时间演化
wflis = wf_time_evolution(st, timelis, energy, states)

# 创建数组存储纠缠熵和ergotropy
steelis = zeros(length(timelis))  # 存储纠缠熵
ergolis = zeros(length(timelis))  # 存储ergotropy

# 对每个时间点计算纠缠熵和能量功率性
for i in eachindex(timelis)
    # 计算纠缠熵
    subrho = rdm_PXP(N, collect(1:div(N,2)), wflis[i])
    steelis[i] = ee(subrho)
    
    # 计算ergotropy
    GS_energy, subenergy, passive_energy = ergotropy_PXP_state(N, div(N,2), wflis[i])
    ergolis[i] = real(GS_energy - passive_energy)
end

# 在一张图上绘制两个量
fig = plot(timelis, ergolis, 
    title=L"Evolution\ of\ Ergotropy\ & Entanglement\ Entropy\ (N=%$N,\ \theta=0)",
    xlabel=L"Time",
    ylabel=L"Value",
    label=L"Ergotropy",
    c=:red,
    legend=:topright,
    titlefont=font(12),    # 将标题字体大小设置为 8
    # guidefont=font(8),    # 坐标轴标签字体大小
    # tickfont=font(6)      # 坐标轴刻度字体大小
)

plot!(timelis, steelis,
    label=L"Entanglement\ Entropy",
    c=:purple)

# 可选：保存图像
# savefig(fig, "ergotropy_ee_evolution.pdf")

# 显示图像
display(fig)



using Statistics
#非mss情况计算dynamics
# 设置系统参数
N_values = [16]  
timelis = collect(0.0:1:1000)  # 时间列表

# 创建存储结果的矩阵
# 行：时间点，列：不同的N值

ergotropy_matrix = zeros(length(timelis), length(N_values))
entropy_matrix = zeros(length(timelis), length(N_values))
mean(ergotropy_matrix[800:1000,:])
mean(ergotropy_matrix[500:1000,:])

# 循环计算每个N值
for (idx, N) in enumerate(N_values)
    # 对角化完整PXP哈密顿量
    energy, states = eigen(PXP_Ham(N))
    
    # 生成初始态：旋转角度为0的态
    rotated_state = rotated_psi_state(N, π/2)
    
    wflis = wf_time_evolution(rotated_state, timelis, energy, states)
    
    # 对每个时间点计算ee和ergotropy
    for i in eachindex(timelis)
        # 获取时间演化后的态
        psi_t = wflis[i]
        
        # 计算纠缠熵
        subrho = rdm_PXP(N, collect(1:div(N,2)), psi_t)
        entropy_matrix[i, idx] = ee(subrho)
        
        GS_energy, subenergy, passive_energy = ergotropy_PXP_state(N, div(N,2), psi_t)
        ergotropy_matrix[i, idx] = real(GS_energy - passive_energy)
    end
    
    println("完成计算 N = $N")
end

# 设置主题
theme(:blues)

# 绘制组合图
fig_combined = plot(
    title = L"Evolution\ of\ Ergotropy\ & Entanglement\ Entropy\ (Non-MSS,\ \theta=\pi/2)",
    grid = false,
    layout = (2, 1),  # 两行一列的布局
    size = (800, 600)  # 设置图表大小
);

# 顶部子图：ergotropy
plot!(timelis, ergotropy_matrix, subplot = 1,
    label = permutedims(["N=$n" for n in sort(N_values, rev=true)]),
    marker = :none,  # 只使用线，不显示标记点
    line_z = permutedims(N_values),
    c = cgrad(:blues, length(N_values), rev=true),
    linewidth = 2.0,
    linealpha = 0.8,
    colorbar = false,
    ylim = (2, 5),  # 设置y轴范围
    title = L"Evolution\ of\ Ergotropy\ (Non-MSS,\ \theta=\pi/2)",
    legend = :topright,
    legendfontsize = 6,  # 图例字体大小
    legend_background_color = RGBA(1, 1, 1, 0.5),  # 半透明背景
    legendfontfamily = "serif"  # 使用衬线字体
);

# 底部子图：纠缠熵
plot!(timelis, entropy_matrix, subplot = 2,
    label = nothing,  # 不显示图例
    marker = :none,  # 只使用线，不显示标记点
    line_z = permutedims(N_values),
    c = cgrad(:blues, length(N_values), rev=true),
    linewidth = 2.0,
    linealpha = 0.8,
    colorbar = false,
    title = L"Evolution\ of\ Entanglement\ Entropy\ (Non-MSS,\ \theta=\pi/2)"
)

display(fig_combined)
# savefig(fig_combined, "ergotropy_entropy_evolution_nonMSS.pdf")




#非mss情况计算平均ergotropy
using Statistics
# 计算不同 N 和 θ 值下的平均 ergotropy (800-1000 时间段)
# 设置系统参数
N_values = [10, 12, 14, 16]  
theta_values = collect(0:0.05:π/2)  # θ 从 0 到 π，步长 0.1
timelis = collect(0.0:5.0:1000)  # 时间列表，从 0 到 1000，步长 1

# 确定用于计算平均值的时间范围（800-1000）
avg_time_range = 800:1000
avg_time_idx = findfirst(x -> x >= avg_time_range.start, timelis):findfirst(x -> x >= avg_time_range.stop, timelis)

# 创建存储结果的矩阵
# 行：θ值，列：不同的N值
avg_ergotropy_matrix = zeros(length(theta_values), length(N_values))

# 循环计算每个 N 值
for (n_idx, N) in enumerate(N_values)
    println("开始计算 N = $N")
    
    # 对角化完整 PXP 哈密顿量（对每个 N 只需计算一次）
    energy, states = eigen(PXP_Ham(N))
    
    # 循环计算每个 θ 值
    for (theta_idx, theta) in enumerate(theta_values)
        println("  计算 θ = $(round(theta, digits=2))")
        
        # 生成初始态：旋转角度为 θ 的态
        rotated_state = rotated_psi_state(N, theta)
        
        # 执行时间演化
        wflis = wf_time_evolution(rotated_state, timelis, energy, states)
        
        # 临时存储每个时间点的 ergotropy
        ergotropy_values = zeros(length(timelis))
        
        # 计算每个时间点的 ergotropy
        for i in eachindex(timelis)
            psi_t = wflis[i]
            GS_energy, subenergy, passive_energy = ergotropy_PXP_state(N, div(N,2), psi_t)
            ergotropy_values[i] = real(GS_energy - passive_energy)
        end
        
        # 计算时间段 800-1000 的平均 ergotropy
        avg_ergotropy = mean(ergotropy_values[avg_time_idx])
        avg_ergotropy_matrix[theta_idx, n_idx] = avg_ergotropy
    end
    
    println("完成计算 N = $N")
end

# 保存计算结果
# using DelimitedFiles
# writedlm("avg_ergotropy_data.csv", [theta_values avg_ergotropy_matrix], ',')
a=avg_ergotropy_matrix./N_values'
# 绘制平均 ergotropy 随 θ 变化的图表，使用不同蓝色渐变
fig_avg_ergotropy = plot(
    # title = L"Average\ Ergotropy\ vs\ \theta\ (Time\ Average:\ 800-1000)",
    xlabel = L"\theta/\pi",
    ylabel = L"\bar{W_A}/N",
    grid = false,
    legend = :topright,
    legendfontsize = 8,
    legend_background_color = RGBA(1, 1, 1, 0.7),
    legendfontfamily = "serif",
    size = (800, 600)
);

# blue_colors = cgrad(:blues, length(N_values) + 2)[1:end-1]; # 使用blues色阶中间的部分

# 为每个N值单独绘制一条线
plot!(theta_values./π, a, label = N_values', marker = :circle,markersize = 3, colorbar = false,
marker_z = N_values', line_z = N_values', linewidth = 2.0, color =cgrad(:viridis, length(N_values)))
# for (i, N) in enumerate(N_values)
#     plot!(
#         fig_avg_ergotropy,
#         theta_values./π, 
#         a[:, i],
#         label = "N=$N",
#         marker = :circle,  
#         markersize = 3,    
#         linewidth = 2.0,
#         color = :viridis
#     )
# end

display(fig_avg_ergotropy)
savefig(fig_avg_ergotropy, "avg_ergoden_thetalis.pdf")
