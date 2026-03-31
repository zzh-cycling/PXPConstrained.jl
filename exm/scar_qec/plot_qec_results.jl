#!/usr/bin/env julia
"""
    plot_qec_results.jl

Analysis and plotting for QEC coherent information results.
Loads scan data, computes scaling exponents, and generates PDF figures.

Usage:
    julia --project=. exm/scar_qec/plot_qec_results.jl
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using JLD2
using Printf
using LinearAlgebra
using Plots
using LaTeXStrings

# Set plot defaults
default(fontfamily="Computer Modern", framestyle=:box, grid=false, 
        legendfontsize=10, tickfontsize=10, guidefontsize=12)

"""
    load_and_analyze(L_values::Vector{Int})

Load QEC scan results and perform scaling analysis.
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
    
    return all_results, all_b_coeffs, all_d_Q
end

"""
    plot_coherent_info(all_results; savepath=nothing)

Plot coherent information I_c(p) for all system sizes.
"""
function plot_coherent_info(all_results; savepath=nothing)
    plt = plot(xlabel=L"p", ylabel=L"I_c(R \rangle Q)",
               title="Coherent Information vs Z-Dephasing",
               legend=:topright, size=(600, 450))
    
    colors = [:blue, :red, :green, :orange, :purple, :brown]
    markers = [:circle, :square, :diamond, :utriangle, :dtriangle, :pentagon]
    
    L_vals = sort(collect(keys(all_results)))
    for (i, L) in enumerate(L_vals)
        data = all_results[L]
        plot!(plt, data.p_values, data.I_c, 
              label="L=$L", color=colors[mod1(i, length(colors))],
              marker=markers[mod1(i, length(markers))], markersize=3,
              linewidth=1.5, markerstrokewidth=0.5)
    end
    
    # Add horizontal line at I_c = 0
    hline!(plt, [0.0], linestyle=:dash, color=:gray, label="", linewidth=1)
    
    # Add log(2) reference
    hline!(plt, [log(2)], linestyle=:dot, color=:black, label=L"\log(2)", linewidth=1)
    
    if savepath !== nothing
        savefig(plt, savepath)
        println("Saved: $savepath")
    end
    
    return plt
end

"""
    plot_threshold_scaling(all_results; savepath=nothing)

Plot threshold p_c vs system size L.
"""
function plot_threshold_scaling(all_results; savepath=nothing)
    L_vals = Int[]
    p_c_vals = Float64[]
    
    for L in sort(collect(keys(all_results)))
        data = all_results[L]
        p_vals = data.p_values
        I_c_vals = data.I_c
        
        # Find threshold via linear interpolation
        for i in 2:length(I_c_vals)
            if I_c_vals[i-1] > 0 && I_c_vals[i] ≤ 0
                p_c = p_vals[i-1] + (p_vals[i] - p_vals[i-1]) * 
                      I_c_vals[i-1] / (I_c_vals[i-1] - I_c_vals[i])
                push!(L_vals, L)
                push!(p_c_vals, p_c)
                break
            end
        end
    end
    
    plt = plot(L_vals, p_c_vals, 
               xlabel=L"L", ylabel=L"p_c",
               title="Threshold vs System Size",
               marker=:circle, markersize=6, linewidth=2,
               color=:blue, legend=false, size=(500, 400))
    
    if savepath !== nothing
        savefig(plt, savepath)
        println("Saved: $savepath")
    end
    
    return plt, L_vals, p_c_vals
end

"""
    plot_b_coefficient(all_b_coeffs; savepath=nothing)

Plot b(L) coefficient vs L showing linear scaling.
"""
function plot_b_coefficient(all_b_coeffs; savepath=nothing)
    L_vals = sort(collect(keys(all_b_coeffs)))
    b_vals = [all_b_coeffs[L] for L in L_vals]
    
    plt = plot(L_vals, b_vals,
               xlabel=L"L", ylabel=L"b(L)",
               title="Knill-Laflamme Coefficient",
               marker=:circle, markersize=6, linewidth=2,
               color=:red, label=L"b(L)", size=(500, 400))
    
    # Add b(L) = L reference line
    L_range = range(minimum(L_vals)-0.5, maximum(L_vals)+0.5, length=100)
    plot!(plt, L_range, L_range, linestyle=:dash, color=:black, 
          label=L"b(L) = L", linewidth=1.5)
    
    if savepath !== nothing
        savefig(plt, savepath)
        println("Saved: $savepath")
    end
    
    return plt
end

"""
    main()

Main analysis routine.
"""
function main()
    println("QEC Scar Code Analysis Results (Constrained Basis)")
    println("="^60)
    
    # Load data
    L_values = [6, 8, 10, 12, 14, 16]
    all_results, all_b_coeffs, all_d_Q = load_and_analyze(L_values)
    
    if isempty(all_results)
        println("\nNo data found. Run run_qec_scan.jl first.")
        return
    end
    
    # Print summary
    println("\n" * "="^60)
    println("b(L) Coefficient Summary")
    println("="^60)
    println("L\td_Q\tb(L)\tp_c")
    println("-"^40)
    
    for L in sort(collect(keys(all_b_coeffs)))
        data = all_results[L]
        p_vals = data.p_values
        I_c_vals = data.I_c
        
        # Find threshold
        p_c = NaN
        for i in 2:length(I_c_vals)
            if I_c_vals[i-1] > 0 && I_c_vals[i] ≤ 0
                p_c = p_vals[i-1] + (p_vals[i] - p_vals[i-1]) * 
                      I_c_vals[i-1] / (I_c_vals[i-1] - I_c_vals[i])
                break
            end
        end
        @printf("%d\t%d\t%.1f\t%.3f\n", L, all_d_Q[L], all_b_coeffs[L], p_c)
    end
    
    # Generate plots
    figdir = joinpath(@__DIR__, "figures")
    mkpath(figdir)
    
    println("\n" * "="^60)
    println("Generating plots...")
    println("="^60)
    
    plot_coherent_info(all_results; savepath=joinpath(figdir, "coherent_info_vs_p.pdf"))
    plot_threshold_scaling(all_results; savepath=joinpath(figdir, "threshold_vs_L.pdf"))
    plot_b_coefficient(all_b_coeffs; savepath=joinpath(figdir, "b_coefficient_vs_L.pdf"))
    
    println("\n" * "="^60)
    println("Analysis complete. Figures saved to: $figdir")
    println("="^60)
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
