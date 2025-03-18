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
st = rotated_psi_state(N, π/2)  # 初始化旋转角度为π/2的量子态

# 执行时间演化
wflis = wf_time_evolution(st, timelis, energy, states)

# 创建数组存储纠缠熵和能量功率性
steelis = zeros(length(timelis))  # 存储纠缠熵
ergolis = zeros(length(timelis))  # 存储能量功率性

# 对每个时间点计算纠缠熵和能量功率性
for i in eachindex(timelis)
    # 计算纠缠熵
    subrho = rdm_PXP(N, collect(1:div(N,2)), wflis[i])
    steelis[i] = ee(subrho)
    
    # 计算能量功率性
    GS_energy, subenergy, passive_energy = ergotropy_PXP_state(N, div(N,2), wflis[i])
    ergolis[i] = real(GS_energy - passive_energy)
end

# 在一张图上绘制两个量
fig = plot(timelis, ergolis, 
    title=L"Evolution\ of\ Ergotropy\ & Entanglement\ Entropy\ (N=%$N,\ \theta=\pi/2)",
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
# savefig(fig, "ergotropy_ee_evolution.png")

# 显示图像
display(fig)

