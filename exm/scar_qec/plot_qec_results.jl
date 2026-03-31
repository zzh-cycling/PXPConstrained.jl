#!/usr/bin/env julia
"""
    plot_qec_results.jl

Analysis and plotting for QEC coherent information results.
Loads scan data, computes scaling exponents, and generates summary.

Usage:
    julia --project=. exm/scar_qec/plot_qec_results.jl
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using JLD2
using Printf
using LinearAlgebra

"""
    load_and_analyze(L_values::Vector{Int})

Load QEC scan results and perform scaling analysis.

# Arguments
- `L_values::Vector{Int}`: System sizes to analyze

# Returns
- `all_results::Dict`: Coherent information data for each L
- `all_b_coeffs::Dict`: b(L) coefficients for each L
"""
function load_and_analyze(L_values::Vector{Int})
    datadir = joinpath(@__DIR__, "data")
    
    all_results = Dict{Int, Any}()
    all_b_coeffs = Dict{Int, Float64}()
    all_d_Q = Dict{Int, Int}()
    
    println("Loading data files...")
    for L in L_values
        filename = joinpath(datadir, "coherent_info_scar_L$(L).jld2")
        if isfile(filename)
            data = load(filename)
            all_results[L] = (p_values=data["p_vals_vec"], I_c=data["I_c_values"])
            all_b_coeffs[L] = data["b_coeff"]
            all_d_Q[L] = data["d_Q"]
            println("  ✓ Loaded L=$L (d_Q=$(data["d_Q"]))")
        else
            println("  ✗ Warning: $filename not found")
        end
    end
    
    if isempty(all_b_coeffs)
        println("\nNo data found. Run run_qec_scan.jl first.")
        return nothing, nothing
    end
    
    # Print b(L) summary
    println("\n" * "="^60)
    println("b(L) Coefficient Summary (Z-dephasing)")
    println("="^60)
    println("L\td_Q\tb(L)")
    println("-"^40)
    for L in sort(collect(keys(all_b_coeffs)))
        @printf("%d\t%d\t%.6f\n", L, all_d_Q[L], all_b_coeffs[L])
    end
    
    # Scaling analysis: b(L) = L for Z-dephasing on Néel states
    println("\n" * "="^60)
    println("Scaling Analysis: b(L) for Z-dephasing")
    println("="^60)
    println("Expected: b(L) = L (linear scaling)")
    
    L_vals = sort(collect(keys(all_b_coeffs)))
    b_vals = [all_b_coeffs[L] for L in L_vals]
    
    # Check linear relationship b(L) = L
    residuals = b_vals .- L_vals
    max_residual = maximum(abs.(residuals))
    @printf("Max deviation from b(L)=L: %.2e\n", max_residual)
    
    if max_residual < 1e-8
        println("✓ Confirmed: b(L) = L exactly")
    else
        # Linear regression
        n = length(L_vals)
        x = Float64.(L_vals)
        y = b_vals
        slope = (n * sum(x .* y) - sum(x) * sum(y)) / (n * sum(x.^2) - sum(x)^2)
        intercept = (sum(y) - slope * sum(x)) / n
        @printf("Linear fit: b(L) ≈ %.4f * L + %.4f\n", slope, intercept)
    end
    
    # Threshold analysis
    println("\n" * "="^60)
    println("Threshold Analysis: p_c where I_c → 0")
    println("="^60)
    
    for L in sort(collect(keys(all_results)))
        data = all_results[L]
        p_vals = data.p_values
        I_c_vals = data.I_c
        
        # Find threshold
        p_c = NaN
        for i in 2:length(I_c_vals)
            if I_c_vals[i-1] > 0 && I_c_vals[i] ≤ 0
                # Linear interpolation
                p_c = p_vals[i-1] + (p_vals[i] - p_vals[i-1]) * 
                      I_c_vals[i-1] / (I_c_vals[i-1] - I_c_vals[i])
                break
            end
        end
        
        if isnan(p_c)
            if all(I_c_vals .> 0)
                @printf("L=%d: I_c > 0 for all p (robust protection)\n", L)
            else
                @printf("L=%d: I_c ≤ 0 for all p (no protection)\n", L)
            end
        else
            @printf("L=%d: p_c ≈ %.3f\n", L, p_c)
        end
    end
    
    return all_results, all_b_coeffs
end

"""
    print_ascii_plot(p_vals, I_c_vals; width=50, height=12)

Print a simple ASCII plot of I_c(p) for visualization.
"""
function print_ascii_plot(p_vals, I_c_vals; width=50, height=12)
    I_max = maximum(I_c_vals)
    I_min = min(minimum(I_c_vals), 0.0)  # Include 0 line
    I_range = I_max - I_min
    
    if I_range < 1e-10
        println("  (Constant value, cannot plot)")
        return
    end
    
    for row in height:-1:0
        I_level = I_min + (row / height) * I_range
        line = @sprintf("%6.3f |", I_level)
        
        for col in 0:width-1
            idx = 1 + round(Int, col * (length(p_vals) - 1) / (width - 1))
            I_normalized = (I_c_vals[idx] - I_min) / I_range * height
            
            if abs(I_normalized - row) < 0.5
                line *= "*"
            elseif abs(I_level) < I_range / height / 2
                line *= "-"  # Zero line
            else
                line *= " "
            end
        end
        println(line)
    end
    
    println("       +" * "-"^width)
    @printf("        0%s1\n", " "^(width-2))
    println("        " * " "^(div(width,2)-1) * "p")
end

"""
    main()

Main analysis routine.
"""
function main()
    println("QEC Scar Code Analysis Results (Constrained Basis)")
    println("="^60)
    
    # Try to load all available sizes
    L_values = [6, 8, 10, 12, 14]
    
    all_results, all_b_coeffs = load_and_analyze(L_values)
    
    if all_results !== nothing && !isempty(all_results)
        # Print ASCII plots for the largest system size
        L_max = maximum(keys(all_results))
        data = all_results[L_max]
        
        println("\n" * "="^60)
        println("I_c(p) for L = $L_max (Z-dephasing)")
        println("="^60)
        print_ascii_plot(data.p_values, data.I_c)
    end
    
    println("\n" * "="^60)
    println("Analysis complete.")
    println("="^60)
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
