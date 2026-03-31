#!/usr/bin/env julia
"""
    run_qec_scan.jl

Main driver for QEC coherent information scans.
Computes I_c(p) for scar encoding under Z-dephasing noise in constrained basis.

Usage:
    julia --project=. exm/scar_qec/run_qec_scan.jl
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using PXPConstrained
using LinearAlgebra
using JLD2
using Printf
using BitBasis

"""
    run_coherent_info_scan(L::Int; kwargs...)

Run coherent information scan for a given system size.

# Arguments
- `L::Int`: System size (number of sites, must be even)

# Keyword Arguments
- `p_values`: Range of dephasing strengths (default: 0.0:0.02:1.0)
- `save_data::Bool`: Whether to save results (default: true)

# Returns
- `results::Vector{Float64}`: I_c(p) values
- `b_coeff::Float64`: b(L) Knill-Laflamme coefficient
"""
function run_coherent_info_scan(L::Int; 
                                 p_values=0.0:0.02:1.0, 
                                 save_data=true)
    
    println("="^60)
    println("Scar QEC Coherent Information Scan (Constrained Basis)")
    println("System size L = $L")
    println("="^60)
    
    # Get constrained basis
    basis = PXP_basis(L, true)
    d_Q = length(basis)
    println("Constrained basis dimension: $d_Q")
    
    # Prepare scar encoding
    println("\nPreparing scar encoding...")
    ψ_scar, ext_basis = prepare_scar_encoding_constrained(L, basis)
    ρ_scar_0 = reference_code_state(ψ_scar)
    
    # Initial coherent information (should be log(2) ≈ 0.693)
    I_c_0 = coherent_information_constrained(ρ_scar_0, d_Q)
    @printf("Initial I_c = %.6f (expected: %.6f = log(2))\n", I_c_0, log(2))
    
    # Z-dephasing scan
    println("\nProcessing Z-dephasing...")
    I_c_values = Float64[]
    
    for (i, p) in enumerate(p_values)
        ρ_noisy = apply_Z_dephasing_extended(ρ_scar_0, p, L, ext_basis)
        I_c = coherent_information_constrained(ρ_noisy, d_Q)
        push!(I_c_values, I_c)
        
        # Print progress for key values
        if p ≈ 0.0 || p ≈ 0.25 || p ≈ 0.5 || p ≈ 0.75 || p ≈ 1.0
            @printf("  p = %.2f: I_c = %.6f\n", p, I_c)
        end
    end
    
    # Find approximate threshold (where I_c crosses 0)
    for i in 2:length(I_c_values)
        if I_c_values[i-1] > 0 && I_c_values[i] ≤ 0
            p_c = collect(p_values)[i-1]
            @printf("  Approximate threshold p_c ≈ %.2f\n", p_c)
            break
        end
    end
    
    # Compute b(L) coefficient
    println("\nKnill-Laflamme Analysis:")
    b_coeff = knill_laflamme_coefficient_Z(L, basis)
    @printf("  b(L=%d, Z) = %.6f\n", L, b_coeff)
    
    if save_data
        datadir = joinpath(@__DIR__, "data")
        mkpath(datadir)
        filename = joinpath(datadir, "coherent_info_scar_L$(L).jld2")
        p_vals_vec = collect(p_values)
        @save filename L p_vals_vec I_c_values b_coeff d_Q
        println("\nData saved to: $filename")
    end
    
    return I_c_values, b_coeff
end

"""
    main()

Run coherent information scans for multiple system sizes.
"""
function main()
    # System sizes to scan (must be even for PBC)
    L_values = [6, 8, 10, 12]
    
    println("QEC Scar Code Analysis (Constrained Basis)")
    println("System sizes: $L_values")
    println()
    
    all_results = Dict{Int, Vector{Float64}}()
    all_b_coeffs = Dict{Int, Float64}()
    
    for L in L_values
        println("\n" * "="^60)
        @time results, b_coeff = run_coherent_info_scan(L)
        all_results[L] = results
        all_b_coeffs[L] = b_coeff
    end
    
    # Summary
    println("\n" * "="^60)
    println("SUMMARY: b(L) Coefficients")
    println("="^60)
    println("L\td_Q\tb(L,Z)")
    for L in L_values
        d_Q = length(PXP_basis(L, true))
        @printf("%d\t%d\t%.6f\n", L, d_Q, all_b_coeffs[L])
    end
    
    return all_results, all_b_coeffs
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
