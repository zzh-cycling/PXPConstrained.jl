module PXPConstrained

using BitBasis
using LinearAlgebra, LaTeXStrings,Printf, LsqFit, Measurements, ITensors, JLD, Plots

export PXP_Ham, PXP_basis, PXP_K_basis, PXP_MSS_basis, PXP_K_Ham, PXP_MSS_Ham, iso_total2FSA, iso_total2K, rdm_PXP, rdm_PXP_K, iso_total2MSS, rdm_PXP_MSS, inversion, wf_time_evolution, myprint, PXP_FSA_Ham, translation, actingH_PXP

export EE, EE_PXP_idx, EE_scaling_fig, EE_PXP_state, Tri_mutual_information, Mutual_information, QFI, domain_wall, domain_wall_density, particlenumber

export proj_Z2, proj_invZ2, sep_scar_FSA, sep_scar_Ob, proj_Ob, proj_FSA, proj_FSA2total, sep_scar_FSA, sep_scar_exact, gene_scar

export fitCCEntEntScal, fitpage_curve, fit_both

include("PXP_functions.jl")
include("Observables.jl")
include("scar_seperate.jl")
include("fitCCEntEntScal.jl")

end
