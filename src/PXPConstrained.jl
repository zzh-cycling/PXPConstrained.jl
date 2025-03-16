module PXPConstrained

using BitBasis
using SparseArrays, ExponentialUtilities
using LinearAlgebra, LaTeXStrings,Printf, LsqFit, Measurements, ITensors, JLD, Plots

export PXP_Ham, PXP_basis, PXP_K_basis, PXP_MSS_basis, PXP_K_Ham, PXP_MSS_Ham, iso_total2FSA, iso_total2K, rdm_PXP, rdm_PXP_K, iso_total2MSS, iso_K2MSS, rdm_PXP_MSS, myprint, PXP_FSA_Ham, translation_matrix, inversion_matrix, actingH_PXP

export ee, ee_PXP_idx, ee_PXP_scaling_fig, ee_PXP_state, tri_mutual_information, mutual_information, qfi, domain_wall, particlenumber, on_siten, ergotropy_PXP_state, ergotropy_PXP_idx, ergotropy_PXP_MSS_state

export proj_Z2, proj_invZ2, sep_scar_FSA, sep_scar_Ob, proj_Ob, proj_FSA, proj_FSA2total, sep_scar_FSA, sep_scar_exact, gene_scar

export fitCCEntEntScal, fitpage_curve, fit_both

export PXP_Ham_sparse, PXP_K_Ham_sparse, PXP_MSS_Ham_sparse, iso_total2K_sparse, iso_total2MSS_sparse, iso_K2MSS_sparse

export wf_time_evolution, wf_time_evolution_mss, rotated_psi_state, rotated_psi_state_mss

include("PXPFunctions.jl")
include("Observables.jl")
include("ScarSeparate.jl")
include("FitEntEntScal.jl")
include("PXPSparse.jl")
include("Dynamics.jl")
end
