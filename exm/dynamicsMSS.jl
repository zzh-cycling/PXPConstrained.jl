using PXPConstrained
using LinearAlgebra
using Plots
using LaTeXStrings

# ================================
#  计算 MSS 下的演化、ergotropy 与纠缠熵
# ================================

# 设置主题为 plasma
# theme(:plasma)

# define the system size
N_values = [10, 12, 14, 16]  
k = 0  # 设置动量量子数
timelis = collect(0.0:0.2:150)  # timelist

# 创建存储结果的矩阵
# 行：时间点，列：不同的N值
ergotropy_matrix = zeros(length(timelis), length(N_values))
entropy_matrix = zeros(length(timelis), length(N_values))

# 循环计算每个N值
for (idx, N) in enumerate(N_values)
    # 生成初始态：在 MSS 中的旋转态 (θ=π/2)
    mss_rotated_state = rotated_psi_state_mss(N, k, 0)
    
    # 在 MSS 中执行时间演化
    mss_wflis = wf_time_evolution_mss(N, k, mss_rotated_state, timelis)
    
    # 对每个时间点计算ee和ergotropy
    for i in eachindex(timelis)
        # 获取时间演化后的态
        psi_t = mss_wflis[i]
        
        # 计算ee
        subrho = rdm_PXP_MSS(N, collect(1:div(N,2)), psi_t, k)
        entropy_matrix[i, idx] = ee(subrho)
        
        # 计算ergotropy
        GS_energy, subenergy, passive_energy = ergotropy_PXP_MSS_state(N, div(N,2), psi_t)
        ergotropy_matrix[i, idx] = real(GS_energy - passive_energy)
    end
    
    println("calculation done for N = $N")
end


#####画图部分
#1，单独画图
# 给出fig白板和基本设置
fig_ergotropy = plot(
    title = L"Evolution\ of\ Ergotropy\ in\ MSS\ (\theta=\pi/2,\ k=%$k)",
    xlabel = L"Time",
    ylabel = L"Ergotropy",
    legend = :topright,
    titlefont = font(10),
    ylim = (0, 12),
    grid = false
);

# 绘制ergotropy
plot!(timelis, ergotropy_matrix,
    label = permutedims(["N=$n" for n in sort(N_values, rev=true)]),
    marker = :circle,
    markersize = 2.5,
    markerstrokewidth = 0.4,
    markeralpha = 0.4,
    markerstrokealpha = 0.6,
    line_z = permutedims(N_values),
    marker_z = permutedims(N_values),
    c = cgrad(:viridis, length(N_values), rev=true),
    linewidth = 1.5,
    linealpha = 0.7,
    markerindex = 1:10:length(timelis),
    colorbar = false
)

display(fig_ergotropy)
#savefig(fig_ergotropy, "ergotropy_evolution_MSS.pdf")

# 绘制纠缠熵图
fig_entropy = plot(
    title = L"Evolution\ of\ Entanglement\ Entropy\ in\ MSS\ (\theta=\pi/2,\ k=%$k)",
    xlabel = L"Time",
    ylabel = L"Entanglement\ Entropy",
    legend = :topright,
    titlefont = font(10),
    ylim = (0, 8),
    grid = false
);

# 绘制纠缠熵
plot!(timelis, entropy_matrix,
    label = permutedims(["N=$n" for n in sort(N_values, rev=true)]),
    marker = :circle,
    markersize = 3.5,
    markerstrokewidth = 0.4,
    markeralpha = 0.4,
    markerstrokealpha = 0.6,
    line_z = permutedims(N_values),
    marker_z = permutedims(N_values),
    c = cgrad(:viridis, length(N_values), rev=true),
    linewidth = 1.5,
    linealpha = 0.7,
    markerindex = 1:10:length(timelis),
    colorbar = false
)

display(fig_entropy)
savefig(fig_entropy, "entropy_evolution_MSS.pdf")


######### 2.绘制组合图
fig_combined = plot(
    title = L"Evolution\ of\ Ergotropy\ & Entanglement\ Entropy\ in\ MSS\ (\theta=\pi,\ k=%$k)",
    xlabel = L"Time",
    ylabel = L"Value",
    titlefont = font(10),
    grid = false,
    layout = (2, 1),  # 两行一列的布局
    size = (800, 600)  # 设置图表大小
);

# 顶部子图：ergotropy
plot!(timelis, ergotropy_matrix, subplot = 1,
    label = permutedims(["N=$n" for n in sort(N_values, rev=true)]),
    marker = :circle,
    markersize = 2.5,
    markerstrokewidth = 0.4,
    markeralpha = 0.5,
    markerstrokealpha = 0.6,
    line_z = permutedims(N_values),
    marker_z = permutedims(N_values),
    c = cgrad(:blues, length(N_values), rev=true),
    linewidth = 1.5,
    linealpha = 0.7,
    markerindex = 1:10:length(timelis),
    colorbar = false,
    ylim = (4, 8),
    title = L"Evolution\ of\ Ergotropy\ (\theta=\pi,\ k=%$k)",
    legend = :topright,
    legendfontsize = 4,  # 缩小图例字体大小
    legend_background_color = RGBA(1, 1, 1, 0.5),  # 半透明背景
    legendfontfamily = "serif"  # 使用衬线字体
);

# 底部子图：纠缠熵
plot!(timelis, entropy_matrix, subplot = 2,
    label = nothing,  # 设置为nothing以不显示图例
    marker = :circle,
    markersize = 2.5,
    markerstrokewidth = 0.4,
    markeralpha = 0.9,
    markerstrokealpha = 0.7,
    line_z = permutedims(N_values),
    marker_z = permutedims(N_values),
    c = cgrad(:blues, length(N_values), rev=true),
    linewidth = 1,
    linealpha = 0.9,
    markerindex = 1:10:length(timelis),
    colorbar = false,
    ylim = (0.5, 2.2),
    title = L"Evolution\ of\ Entanglement\ Entropy\ (\theta=\pi,\ k=%$k)",
)

#保存图片
# display(fig_combined)
# savefig(fig_combined, "combined_evolution_MSS.pdf")


##### 3. 组合图，但是用直线的版本，不用circle
fig_combined = plot(
    title = L"Evolution\ of\ Ergotropy\ & Entanglement\ Entropy\ in\ MSS\ (\theta=\pi,\ k=%$k)",
    xlabel = L"Time",
    ylabel = L"Value",
    titlefont = font(10),
    grid = false,
    layout = (2, 1),  # 两行一列的布局
    size = (800, 600)  # 设置图表大小
);

plot!(timelis, ergotropy_matrix, subplot = 1,
    label = permutedims(["N=$n" for n in sort(N_values, rev=true)]),
    marker = :none,  # 不显示标记点，只显示线
    line_z = permutedims(N_values),
    c = cgrad(:blues, length(N_values), rev=true),
    linewidth = 2.0,  # 增加线宽以提高可见性
    linealpha = 0.8,  # 增加线的不透明度
    colorbar = false,
    ylim = (2, 5),
    title = L"Evolution\ of\ Ergotropy\ (Mss\ \theta=\pi,\ k=%$k)",
    legend = :topright,
    legendfontsize = 4,  # 缩小图例字体大小
    legend_background_color = RGBA(1, 1, 1, 0.5),  # 半透明背景
    legendfontfamily = "serif"  # 使用衬线字体
);

# 底部子图：纠缠熵
plot!(timelis, entropy_matrix, subplot = 2,
    label = nothing,  # 设置为nothing以不显示图例
    marker = :none,  # 不显示标记点，只显示线
    line_z = permutedims(N_values),
    c = cgrad(:blues, length(N_values), rev=true),
    linewidth = 2.0,  # 增加线宽以提高可见性
    linealpha = 0.8,  # 线的不透明度
    colorbar = false,
    ylim = (0.7, 2.2),
    title = L"Evolution\ of\ Entanglement\ Entropy\ (Mss\ \theta=0,\ k=%$k)",
)

###计算不同theta的ergotropy和entanglement entropy
N_values = 18 
k = 0  # 设置动量量子数
timelis = collect(0.0:0.2:150)  # timelist

# 创建存储结果的矩阵
# 行：时间点，列：不同的N值
ergotropy_matrix = zeros(length(timelis), length(N_values))
entropy_matrix = zeros(length(timelis), length(N_values))

# 循环计算每个N值
for (idx, N) in enumerate(N_values)
    # 生成初始态：在 MSS 中的旋转态 (θ=π/2)
    mss_rotated_state = rotated_psi_state_mss(N, k, π/2)
    
    # 在 MSS 中执行时间演化
    mss_wflis = wf_time_evolution_mss(N, k, mss_rotated_state, timelis)
    
    # 对每个时间点计算ee和ergotropy
    for i in eachindex(timelis)
        # 获取时间演化后的态
        psi_t = mss_wflis[i]
        
        # 计算ee
        subrho = rdm_PXP_MSS(N, collect(1:div(N,2)), psi_t, k)
        entropy_matrix[i, idx] = ee(subrho)
        
        # 计算ergotropy
        GS_energy, subenergy, passive_energy = ergotropy_PXP_MSS_state(N, div(N,2), psi_t)
        ergotropy_matrix[i, idx] = real(GS_energy - passive_energy)
    end
    
    println("calculation done for N = $N")
end





# 画不同theta的动力学演化
theme(:blues)

# 系统参数
N = 18  # 固定系统大小为18
k = 0   # 动量子空间
timelis = collect(0.0:0.2:150)  # 时间范围

# 定义不同的θ值
theta_values = collect(0:0.31:π)  # θ从0到π，分11个点

# 创建存储结果的矩阵
# 行：时间点，列：不同的θ值
ergotropy_matrix = zeros(length(timelis), length(theta_values))

# 循环计算每个θ值
for (idx, theta) in enumerate(theta_values)
    # 生成在MSS中的旋转态(对应角度θ)
    mss_rotated_state = rotated_psi_state_mss(N, k, theta)
    
    # 在MSS中执行时间演化
    mss_wflis = wf_time_evolution_mss(N, k, mss_rotated_state, timelis)
    
    # 对每个时间点计算ergotropy
    for i in eachindex(timelis)
        # 获取时间演化后的态
        psi_t = mss_wflis[i]
        
        # 计算ergotropy
        GS_energy, subenergy, passive_energy = ergotropy_PXP_MSS_state(N, div(N,2), psi_t)
        ergotropy_matrix[i, idx] = real(GS_energy - passive_energy)
    end
    
    println("计算完成 θ = $(round(theta, digits=2))")
end

# 绘制ergotropy随时间演化的图
fig_ergotropy = plot(
    title = L"Evolution\ of\ Ergotropy\ for\ Different\ Initial\ States\ (N=%$N,\ k=%$k)",
    xlabel = L"θ/\pi",
    ylabel = L"Ergotropy",
    legend = :topright,
    titlefont = font(10),
    grid = false,
    size = (800, 500)
)

# 为每个θ值绘制一条线
c_gradient = cgrad(:blues, length(theta_values), rev=true)

# 为每个θ绘制一条线
for (idx, theta) in enumerate(theta_values)
    theta_label = round(theta/π, digits=2)
    if theta_label == 0.0
        theta_str = "0"
    elseif theta_label == 1.0
        theta_str = "\\pi"
    else
        theta_str = "$(theta_label)\\pi"
    end
    
    plot!(
        fig_ergotropy,
        timelis, 
        ergotropy_matrix[:, idx],
        label = L"\theta = %$theta_str",
        linewidth = 2.0,
        linealpha = 0.8,
        c = c_gradient[idx]
    )
end

# 图例设置
plot!(
    fig_ergotropy,
    legendfontsize = 4,
    legend_background_color = RGBA(1, 1, 1, 0.7),
    legendfontfamily = "serif"
)

display(fig_ergotropy)
savefig(fig_ergotropy, "ergotropy_evolution_N18_different_theta.png")

# 创建热图（时间-θ-ergotropy）
fig_heatmap = heatmap(
    theta_values./π, 
    timelis, 
    ergotropy_matrix,
    title = L"Ergotropy\ vs\ Time\ and\ \theta\ (N=%$N,\ k=%$k)",
    xlabel = L"\theta/\pi",
    ylabel = L"Time",
    c = :thermal,
    colorbar_title = L"Ergotropy",
    size = (1000, 800),
    titlefont = font(10)
)

display(fig_heatmap)
savefig(fig_heatmap, "ergotropy_heatmap_N18_theta_time.pdf")




#### 创建三维表面图（时间-θ-ergotropy）
fig_3d = surface(
    timelis,                # x轴：时间
    theta_values./π,        # y轴：θ/π
    ergotropy_matrix',      # z轴：ergotrop
    title = L"Ergotropy\ vs\ Time\ and\ \theta\ (N=%$N,\ k=%$k)",
    xlabel = L"Time",
    ylabel = L"\theta/\pi",
    zlabel = L"Ergotropy",
    c = :thermal,         
    camera = (30, 30),     # 设置视角 (角度1, 角度2)
    size = (800, 600),    
    titlefont = font(10),
    grid = false,
    colorbar_title = L"Ergotropy",
)

display(fig_3d)
# savefig(fig_3d, "ergotropy_3d_N18_theta_time.pdf")

# 从不同角度创建多个视图
angles = [(30, 30), (60, 30), (90, 60), (20, 90)]
for (i, angle) in enumerate(angles)
    fig_view = surface(
        timelis,
        theta_values./π,
        ergotropy_matrix',
        title = L"Ergotropy\ 3D\ View\ %$i\ (N=%$N,\ k=%$k)",
        xlabel = L"Time",
        ylabel = L"\theta/\pi",
        zlabel = L"Ergotropy",
        c = :inferno,
        camera = angle,
        size = (800, 600),
        titlefont = font(10),
        grid = false
    )
    display(fig_view)
    # savefig(fig_view, "ergotropy_3d_view$(i)_N18.pdf")
end


