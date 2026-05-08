"""
    PXPConstrained

A Julia package for simulating the 1D PXP chain with Rydberg blockade constraints.

This module provides comprehensive tools for:
- Generating constrained basis states (Fibonacci chains)
- Constructing PXP Hamiltonians with periodic/open boundary conditions
- Implementing translational and inversion symmetries
- Computing observables and entanglement measures
- Time evolution dynamics
- Sparse matrix representations for large systems

The PXP model describes Rydberg atom arrays with nearest-neighbor blockade constraints,
where atoms cannot be simultaneously excited if they are adjacent.
"""
module PXPConstrained

using BitBasis
using SparseArrays, ExponentialUtilities
using LinearAlgebra, ITensors

export actingH_PXP, PXP_Ham, PXP_basis, myprint, iso_full2cons

export PXP_K_basis, PXP_MSS_basis, PXP_K_Ham, PXP_MSS_Ham,  iso_total2K, rdm_PXP, rdm_PXP_K, iso_total2MSS, iso_K2MSS, rdm_PXP_MSS,   mapstate_K2total, mapstate_MSS2K, mapstate_MSS2total

export ee, ee_PXP_idx, ee_PXP_state, tri_mutual_information, mutual_information, qfi, domain_wall_density, particlenumber, on_siten, ergotropy_PXP_state, ergotropy_PXP_idx, ergotropy_PXP_MSS_state, anti_ferro_order, translation_matrix, inversion_matrix, OTOC, OTOC_map, apply_operator_map, Z_map, X_map

export iso_total2FSA, PXP_FSA_Ham, proj_Z2, proj_invZ2, sep_scar_FSA, sep_scar_Ob, proj_Ob, proj_FSA, proj_FSA2total, sep_scar_FSA, sep_scar_exact, gene_scar

export PXP_Ham_sparse, PXP_K_Ham_sparse, PXP_MSS_Ham_sparse, iso_total2K_sparse, iso_total2MSS_sparse, iso_K2MSS_sparse

export wf_time_evolution, wf_time_evolution_sparse, rotated_psi_state, rotated_psi_state_mss

# QEC (Quantum Error Correction) exports - constrained basis version
export coherent_information_constrained, apply_Z_dephasing_extended
export prepare_scar_encoding_constrained, prepare_thermal_encoding_constrained
export knill_laflamme_coefficient_Z, reference_code_state
export von_neumann_entropy, partial_trace_R
export neel_state_bitstr, build_extended_basis
export apply_Z_dephasing_constrained

include("PXPBasis.jl")
include("PXPSymmetry.jl")
include("Observables.jl")
include("ScarSeparate.jl")
include("PXPSparse.jl")
include("Dynamics.jl")
include("QEC.jl")
end
