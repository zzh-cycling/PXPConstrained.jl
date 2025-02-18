module PXPConstrained

using LinearAlgebra
using BitBasis

export PXP_Ham, PXP_basis, PXP_K_basis, PXP_MSS_basis, PXP_K_Ham, PXP_MSS_Ham, iso_total2FSA, iso_total2K, rdm_PXP, rdm_PXP_K, rdm_PXP_MSS, Inv_proj_matrix, wf_time_evolution, myprint, PXP_FSA_Ham, translation_matrix

export EE, EE_PXP_idx, EE_scaling_fig, EE_PXP_state, Tri_mutual_information, Mutual_information, Quantumfisherinfo, domain_wall_density, particlenumber

export proj_Z2, sep_scar_state, superposition_scar,sep_scar_FSA2, proj_invZ2, sep_scar_FSA, sep_scar_Ob

export fitCCEntEntScal, fitpage_curve, fit_both

include("PXP_functions.jl")
include("Observables.jl")
include("scar_seperate.jl")
inclue("fitCCEntEntScal.jl")

end
