module PXPConstrained

using BitBasis
using LinearAlgebra, LaTeXStrings,Printf, LsqFit, Measurements, ITensors, JLD, Plots

export PXP_Ham, PXP_basis, PXP_K_basis, PXP_MSS_basis, PXP_K_Ham, PXP_MSS_Ham, iso_total2FSA, iso_total2K, rdm_PXP, rdm_PXP_K, iso_total2MSS, iso_K2MSS, rdm_PXP_MSS, wf_time_evolution, myprint, PXP_FSA_Ham, translation_matrix, inversion_matrix, actingH_PXP

export ee, ee_PXP_idx, ee_PXP_scaling_fig, ee_PXP_state, tri_mutual_information, mutual_information, qfi, domain_wall, particlenumber, on_siten, ergotropy_PXP_state, ergotropy_PXP_idx, ergotropy_PXP_idx_OBC

export proj_Z2, proj_invZ2, sep_scar_FSA, sep_scar_Ob, proj_Ob, proj_FSA, proj_FSA2total, sep_scar_FSA, sep_scar_exact, gene_scar

export fitCCEntEntScal, fitpage_curve, fit_both

include("PXPFunctions.jl")
include("Observables.jl")
include("ScarSeparate.jl")
include("FitEntEntScal.jl")

end
