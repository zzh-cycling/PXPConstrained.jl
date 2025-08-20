using PXPConstrained, BitBasis
using JLD
using Plots
using LinearAlgebra
using LaTeXStrings
ITensors.set_warn_order(60)

function plot_fit_CC_page()
    data=load("/Users/cycling/Documents/projects/big_data/scar_thermal_FSA/scar_thermal_FSA/scar_thermal_N16.jld")
    scar=data["scar"]
    the=data["thermal_ensemble"]

    N=16
    splitlis=collect(1:15)
    ee_scar=EE_PXP_state(BitStr{N,Int}, splitlis, scar)
    ee_the=[EE_PXP_state(BitStr{N,Int}, splitlis, the[:,i]) for i in 1:size(the)[2]]
    fig_scar=fitCCEntEntScal(ee_scar, pbc=true)
    savefig(fig_scar[2], "/Users/cycling/Documents/projects/quantumErgotropy/figs/PXP_N16_scar_scaling/fig0.pdf")

    sum=0
    for (i,ee) in enumerate(ee_the)
        cent, fig=fitpage_curve(ee)
        sum+=cent
        savefig(fig, "/Users/cycling/Documents/projects/quantumErgotropy/figs/PXP_N16_scar_scaling/fig$(i).pdf")
    end

    sum/size(the)[2]
   
end


function gene_2scar(N::Int)
    T=BitStr{N,Int}
    basis= PXP_basis(T)
    basis_int = [i.buf for i in basis]
    Trans = translation_matrix(T)
    scar=gene_scar(N)
    scar1=storage(scar)

    scar1=scar1[basis_int.+1]
    scar1/=norm(scar1)
    scar2 = Trans*scar1
    @test norm(scar2) ≈ 1

    scark0=scar1+scar2
    scark0/=norm(scark0)
    scarkpi=scar1-scar2
    scarkpi/=norm(scarkpi)

    return scark0, scarkpi
end

function ee_2scar(::Type{T}) where {N, T <: BitStr{N}}
    scark0, scarkpi = gene_2scar(N)
    Trans = translation_matrix(T)
    @test isapprox(scark0'*Trans*scark0, 1)
    rhok0, rhokpi = rdm_PXP(T, collect(1:div(N,2)), scark0), rdm_PXP(T, collect(1:div(N,2)), scarkpi)
    sk0, skpi = ee(rhok0), ee(rhokpi)
    return sk0, skpi
end

Nlis = collect(8:2:20)
deltaS = zeros(length(Nlis))

for i in eachindex(Nlis)
    sk0, skpi = ee_2scar(BitStr{Nlis[i],Int})
    deltaS[i] = skpi - sk0
end

Plots.plot(1 ./Nlis, deltaS, seriestype=:scatter, xlabel=L"1/N", ylabel=L"\Delta S= S(|\Psi, k=0 \rangle)- S(|\Psi, k=\pi \rangle)", legend=false)
savefig("/Users/cycling/Documents/projects/quantumErgotropy/figs/exact_scar/exact_wosymmetry/EE/DeltaSk0kpiinvN.pdf")

Plots.plot(Nlis, deltaS, seriestype=:scatter, xlabel=L"N", ylabel=L"\Delta S= S(|\Psi, k=0 \rangle)- S(|\Psi, k=\pi \rangle)", legend=false)
savefig("/Users/cycling/Documents/projects/quantumErgotropy/figs/exact_scar/exact_wosymmetry/EE/DeltaSk0kpiN.pdf")

N=16
thetalis = collect(range(0, π, 49))
stlis = rotated_psi_state.(N, thetalis)
ρlis = [rdm_PXP(N, collect(1:div(N,2)-1), st) for st in stlis]
eslis = [eigvals(ρ) for ρ in ρlis]
truncate_Slis = [sum(x-> -x*log(x), es[end-2:end]) for es in eslis]
numer_Slis = [ee(rho) for rho in ρlis]

plot(thetalis, numer_Slis, seriestype=:scatter, xlabel=L"\theta", ylabel=L"S(\theta)", label="", legend=false)
plot!(thetalis, truncate_Slis)
esmat = hcat(eslis...)

fig = plot(thetalis, esmat[end-4:end-1, :]', seriestype=:scatter, xlabel=L"\theta", ylabel=L"eigenvalues", label="", legend=false)

f(t) = sqrt(2)*sqrt(44*cos(2*t)-3*cos(4*t)+87)
h(t) = -16cos(t)- 2cos(2t) + 2
ρ1(t) = 1 + (4h(t)^2 - f(t)^2)/(2f(t)^2 + 4h(t)*f(t))
ρ2(t) = 1 + (8*(sin(2t) + 2sin(t))^2)/(f(t)^2 + 2h(t)*f(t))
plot!(thetalis, ρ1.(thetalis), label=L"\rho_1(\theta)", color=:red)
plot!(thetalis, ρ2.(thetalis), label=L"\rho_2(\theta)", color=:blue)

Δ(θ) = sqrt(174 + 88*cos(2θ) - 6*cos(4θ))

den(θ) = (-87 - cos(2θ)*(44 + Δ(θ)) + Δ(θ) + 8cos(θ)*Δ(θ) + 3cos(4θ))^2

# 第一项
expr1(θ) = 4096 * sin(θ/2)^8 * sin(θ)^4 / den(θ)

# 第二项
numer2(θ) = (-2 - 16cos(θ) + 2cos(2θ) + Δ(θ))^4
expr2(θ) = numer2(θ) / (16 * den(θ))

plot(thetalis, expr1.(thetalis), label=L"e_1(\theta)", color=:green)
plot!(thetalis, expr2.(thetalis), label=L"e_2(\theta)", color=:orange)

λ1lis = expr1.(thetalis)
λ2lis = expr2.(thetalis)
Slis = zeros(length(thetalis))

for i in eachindex(thetalis)
    Slis[i] = -λ1lis[i]*log(λ1lis[i]) - λ2lis[i]*log(λ2lis[i])
end

plot(thetalis, Slis, seriestype=:scatter, xlabel=L"\theta", ylabel=L"S(\theta)", label="", legend=false)