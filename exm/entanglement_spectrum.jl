using PXPConstrained
using BitBasis
using LinearAlgebra
using Plots
using LaTeXStrings
include("theoretic_formula.jl")

# Numerical: entanglement entropy and spectrum vs θ
N = 20
θlis = range(0, π, length=49)
ψlis = [rotated_psi_state(N, θ) for θ in θlis]
ρlis = [rdm_PXP(N, collect(1:div(N,2)), ψ) for ψ in ψlis]
eelis = [ee(ρ) for ρ in ρlis]
spectrum = [sort(eigvals(ρ); rev=true) for ρ in ρlis]
nz_spectrum = [eigs[1:4] for eigs in spectrum]
nz_spectrum_matrix = hcat(nz_spectrum...)'

N16 = 16
ψlis16 = [rotated_psi_state(N16, θ) for θ in θlis]
ρlis16 = [rdm_PXP(N16, collect(1:div(N16,2)), ψ) for ψ in ψlis16]
eelis16 = [ee(ρ) for ρ in ρlis16]
spectrum16 = [sort(eigvals(ρ); rev=true) for ρ in ρlis16]
nz_spectrum16 = [eigs[1:4] for eigs in spectrum16]
nz_spectrum_matrix16 = hcat(nz_spectrum16...)'

N10 = 10
ψlis10 = [rotated_psi_state(N10, θ) for θ in θlis]
ρlis10 = [rdm_PXP(N10, collect(1:div(N10,2)), ψ) for ψ in ψlis10]
eelis10 = [ee(ρ) for ρ in ρlis10]
spectrum10 = [sort(eigvals(ρ); rev=true) for ρ in ρlis10]
nz_spectrum10 = [eigs[1:4] for eigs in spectrum10]
nz_spectrum_matrix10 = hcat(nz_spectrum10...)'

fig18 = plot(θlis, nz_spectrum_matrix,  lw=2, label=false,
    legend_background_color=nothing,
    legend_foreground_color=nothing,
    color = :blues,
    xlabel = L"θ",
    ylabel = L"ρ_i",
    title = "Entanglement Spectrum vs θ (N=20)",
)
fig12 = plot(θlis, nz_spectrum_matrix10, lw=2,
    color = :reds,
    label=false,
    legend_background_color=nothing,
    legend_foreground_color=nothing,
    xlabel = L"θ",
    ylabel = L"ρ_i",
    title = "Entanglement Spectrum vs θ (N=10)",
)


# Theoretical: entanglement entropy and spectrum vs θ
formula_eelis = zeros(length(θlis))
matrix_eelis = zeros(length(θlis))
for (i, θ) in enumerate(θlis)
    λ1_val = ρ3(θ)
    λ2_val = ρ3(θ)
    λ3_val = 1 - λ1_val - λ2_val
    formula_eelis[i] = -λ1_val*log(λ1_val) - λ2_val*log(λ2_val) - λ3_val*log(λ3_val)
end

spectrum = zeros(4, length(θlis))
⊗(a::AbstractArray, b::AbstractArray) = kron(a, b)
for (i, θ) in enumerate(θlis)
    # v2LABv2RAB = reshape(normvec2LAB(θ),2,2)*reshape(normvec2RAB(θ),2,2)
    # v2LBAv2RBA = reshape(normvec2LBA(θ),2,2)*reshape(normvec2RBA(θ),2,2)
    v2LABv2RBA = reshape(normvec2LAB(θ),2,2)*reshape(normvec2RBA(θ),2,2)
    v2LBAv2RAB = reshape(normvec2LBA(θ),2,2)*reshape(normvec2RAB(θ),2,2)
    M = v2LABv2RBA ⊗  v2LBAv2RAB
    spectrum[:, i] = diag(M)
    M /= tr(M)
    matrix_eelis[i] = -sum(spectrum[:, i].*log.(spectrum[:, i]))
end

plot(θlis[1:end-1], [normvec2LAB(θ)[1] for θ in θlis[1:end-1]], lw=2, label="LAB",
    legend_background_color=nothing,
    legend_foreground_color=nothing,
    xlabel = L"θ",
    ylabel = "Element1",
)
plot!(θlis[1:end-1], [normvec2RAB(θ)[1] for θ in θlis[1:end-1]], lw=2, label="RAB")
plot!(θlis[1:end-1], [normvec2LBA(θ)[1] for θ in θlis[1:end-1]], lw=2, label="LBA")
plot!(θlis[1:end-1], [normvec2RBA(θ)[1] for θ in θlis[1:end-1]], lw=2, label="RBA")

plot(θlis[1:end-1], [normvec2LAB(θ)[4] for θ in θlis[1:end-1]], lw=2, label="LAB",
    legend_background_color=nothing,
    legend_foreground_color=nothing,
    xlabel = L"θ",
    ylabel = "Element4",
)
plot!(θlis[1:end-1], [normvec2RAB(θ)[4] for θ in θlis[1:end-1]], lw=2, label="RAB")
plot!(θlis[1:end-1], [normvec2LBA(θ)[4] for θ in θlis[1:end-1]], lw=2, label="LBA")
plot!(θlis[1:end-1], [normvec2RBA(θ)[4] for θ in θlis[1:end-1]], lw=2, label="RBA")

plot(θlis[1:end-1], [normvec2LAB(θ)[1] for θ in θlis[1:end-1]], lw=2, label="LAB",
    legend_background_color=nothing,
    legend_foreground_color=nothing,
    xlabel = L"θ",
    ylabel = "Element1",
)
plot!(θlis[1:end-1], [normvec2RAB(θ)[1] for θ in θlis[1:end-1]], lw=2, label="RAB")
plot!(θlis[1:end-1], [normvec2LBA(θ)[1] for θ in θlis[1:end-1]], lw=2, label="LBA")
plot!(θlis[1:end-1], [normvec2RBA(θ)[1] for θ in θlis[1:end-1]], lw=2, label="RBA")


fig_eelis = plot(θlis, eelis, lw=2, label="Numerical",
    legend_background_color=nothing,
    legend_foreground_color=nothing,
    xlabel = L"θ",
    ylabel = L"S",
    title = "Entanglement Entropy vs θ (N=20)",
)
plot!(fig_eelis, θlis, formula_eelis, lw=2, ls=:dash, label="Theoretical")
plot!(fig_eelis, θlis, matrix_eelis, lw=2, ls=:dot, label="Theoretical (Matrix Method)")