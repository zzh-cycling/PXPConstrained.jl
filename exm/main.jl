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