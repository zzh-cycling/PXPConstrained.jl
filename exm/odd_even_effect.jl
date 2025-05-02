using PXPConstrained
using LinearAlgebra
using Plots

N=20
energy, states=eigen(PXP_MSS_Ham(N,0))
st=rotated_psi_state_mss(N,0,0)
timelis=collect(0:1.0:1000);
wflis=wf_time_evolution(st, timelis, energy, states)
plot(timelis, [norm(i[end]) for i in wflis])
Slis=[ee(rdm_PXP_MSS(N, collect(1:div(N,2)),i,0)) for i in wflis]
plot(timelis, Slis)

energy, states=eigen(PXP_MSS_Ham(N,div(N,2),-1))
st=rotated_psi_state_mss(N,0,0,-1)
timelis=collect(0:1.0:1000);
wflis=wf_time_evolution(st, timelis, energy, states)
plot(timelis, [norm(i[end]) for i in wflis])
Slis=[ee(rdm_PXP_MSS(N, collect(1:div(N,2)),i,0)) for i in wflis]
plot(timelis, Slis)