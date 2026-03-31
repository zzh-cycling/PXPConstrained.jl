using PXPConstrained
using BitBasis
using LinearAlgebra
using Plots
using Optim
include("FitEntEntScal.jl")

L = 24
H = PXP_MSS_Ham(L, 0)
energy, states = eigen(H)
psi = rotated_psi_state_mss(L, 0, π/2)
tlis = collect(0.1:200:2000)
wflis = wf_time_evolution(psi, tlis, energy, states)

st = wflis[end]
"""
    find_β(ρ, H)

Finds the inverse temperature `β` for a given density matrix `ρ` and Hamiltonian `H`,
such that the energy of the corresponding thermal state matches the energy of `ρ`.

The energy of the state `ρ` is `E = tr(ρ * H)`.
The energy of a thermal state at inverse temperature `β` is `E(β) = tr(H * exp(-β*H)) / tr(exp(-β*H))`.
This function solves `E(β) = E` for `β` using a numerical root-finding algorithm.
"""
function find_β(ρ, H)
    # It's computationally better to work in the eigenbasis of ρ and H
    ϵ = reverse(eigvals(H))
    λ = eigvals(ρ)

    # Define the function whose root we want to find
    # E(β) - target_energy = 0
    function energy_diff(β)
        exp_vals = exp.(-β .* ϵ)
        partition_func = sum(exp_vals)

        return sum((exp_vals/partition_func .- λ).^2)
    end

    # Use automatic differentiation with optimization
    result = optimize(energy_diff, 1e-6, 100.0, Brent())
    
    return Optim.minimizer(result)
end


ρ = rdm_PXP_MSS(L, collect(1:4), st, 0)
eigvals(ρ)
HA = PXP_Ham(4, false)

β = find_β(ρ, HA)
subρ = exp(-β*HA)
subρ = subρ / tr(subρ)


GSlis = zeros(length(tlis))
subGSlis = zeros(length(tlis))
passive_energylis = zeros(length(tlis))

for (i, wf) in enumerate(wflis)
    @show i
    GSlis[i], subGSlis[i], passive_energylis[i] = real.(ergotropy_PXP_MSS_state(L, div(L,2), wf, 0))
end

Wlis = GSlis .- passive_energylis


plot(tlis, Wlis,xscale=:log10,  xlabel="Time", ylabel="Ergotropy", title="Ergotropy Dynamics in PXP Model", legend=false)

sum(Wlis[751:end])/length(Wlis[751:end])

L=16
H = PXP_Ham(L)
energy, states = eigen(H)
pseudo_scar = states[:,1075]
Inv = inversion_matrix(L)
T = translation_matrix(L)
pseudo_scarIp = (I(2207) .+ Inv)/2*pseudo_scar
pseudo_scarIm = (I(2207) .- Inv)/2*pseudo_scar

eelisIp = ee_PXP_state(L, collect(1:L-1), pseudo_scarIp)
plot(collect(1:L-1), eelisIp, xlabel="Subsystem Size", ylabel="Entanglement Entropy", title="EE Profile of Pseudo-Scar State", legend=false)
cent, fig = fitCCEntEntScal(eelisIp, mincut=2, pbc=true)
eelisIm = ee_PXP_state(L, collect(1:L-1), pseudo_scarIm)
cent, fig = fitCCEntEntScal(eelisIm, mincut=2, pbc=true)

scar1, scar2, thermal = sep_scar_FSA(L, energy, states)
thermal2 = thermal[:, end-1]
E_GS1, subE1, passive1 = ergotropy_PXP_state(L, div(L,2), scar1)
E_GS2, subE2, passive2 = ergotropy_PXP_state(L, div(L,2), scar2)
E_GS3, subE3, passive3 = ergotropy_PXP_state(L, div(L,2), thermal2)
W1 = E_GS1 - passive1
W2 = E_GS2 - passive2
W3 = E_GS3 - passive3
scar1'*Inv*scar1
scar2'*Inv*scar2
scar1'*T*scar1
scar2'*T*scar2

H_kpi = PXP_MSS_Ham(L, div(L,2), -1)
energy_kpi, states_kpi = eigen(H_kpi)
Slis = zeros(size(states_kpi, 2))

for i in 1:size(states_kpi, 2)
    st = states_kpi[:, i]
    ρ = rdm_PXP_MSS(L, collect(1:div(L,2)), st, div(L,2), -1)
    S = ee(ρ)
    Slis[i] = S
end

scatter(energy_kpi, Slis, xlabel="Energy", ylabel="Entanglement Entropy", title="EE vs Energy in PXP Model with k=π", legend=false)


# newst = scar1 + scar2 # or scar1 - scar2, both contribute less ergotropy than only scar1
# newst = newst / norm(newst)
# E_GS4, subE4, passive4 = ergotropy_PXP_state(L, div(L,2), newst)
# W4 = E_GS4 - passive4