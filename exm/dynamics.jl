using PXPConstrained
using LinearAlgebra
using BenchmarkTools
using Plots
using LaTeXStrings

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
    ergolis[i]= GS_energy-passive_energy
    subrho=rdm_PXP(N, collect(1:div(N,2)), rotated_state)
    eelis[i]=ee(subrho)
    gslis[i]=GS_energy
    palsis[i]=passive_energy
end

rotated_state=rotated_psi_state(N, π)

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

energy, states = eigen(PXP_Ham(N))
timelis=collect(0.0:0.2:500)
st=zeros(2207);st[end]=1
wflis= wf_time_evolution(st, timelis, energy, states)
z2eelis=zeros(length(timelis))

for i in eachindex(timelis)
    subrho=rdm_PXP(N, collect(1:div(N,2)), wflis[i])
    z2eelis[i]=ee(subrho)
end

fig2=plot(timelis, z2eelis, title="Time evolution of z2 component", xlabel="Time", ylabel="Eigenenergy", label="z2 component", c=:red)



