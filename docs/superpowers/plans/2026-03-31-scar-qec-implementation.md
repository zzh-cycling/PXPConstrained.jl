# Scar QEC Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement numerical computation of coherent information for scar vs thermal encodings under local dephasing, demonstrating that quantum many-body scars form approximate quantum error correcting codes.

**Architecture:** Build on existing PXPConstrained package infrastructure. Add QEC module to `src/` for core functionality, calculation scripts in `exm/scar_qec/`.

**Tech Stack:** Julia, JLD2 for data storage, existing PXPConstrained package (PXP_basis, PXP_Ham, ee, rdm_PXP)

---

## File Structure

```
src/
├── PXPConstrained.jl              # Main module (add QEC exports)
├── QEC.jl                         # NEW: Core QEC functionality
│   ├── coherent_information       # I_c computation using existing ee()
│   ├── dephasing_channel          # Noise channel implementation
│   ├── scar_thermal_encoding      # Code space preparation
│   └── knill_laflamme             # KL condition & b(L) computation
exm/
├── scar_qec/
│   ├── run_qec_scan.jl            # Parameter sweep driver
│   ├── plot_qec_results.jl        # Visualization
│   └── data/                      # Output data directory
│       ├── coherent_info_scar_L{N}.jld2
│       └── coherent_info_thermal_L{N}.jld2
```

---

## Task 1: QEC Core Module in src/

**Files:**
- Create: `src/QEC.jl`
- Modify: `src/PXPConstrained.jl` (add include and exports)

- [ ] **Step 1: Create QEC.jl module**

```julia
# src/QEC.jl
"""
    QEC.jl

Quantum Error Correction analysis for PXP scar states.
Implements coherent information, dephasing channels, and Knill-Laflamme analysis.
"""

using LinearAlgebra
using SparseArrays
using BitBasis

export coherent_information, apply_dephasing_channel
export prepare_scar_encoding, prepare_thermal_encoding
export knill_laflamme_coefficient, reference_code_state

#=============================================================================
# Dephasing Channel
=============================================================================#

"""
    pauli_matrices()

Return Pauli matrices as sparse matrices: (σx, σy, σz, I2)
"""
function pauli_matrices()
    σx = sparse(ComplexF64[0 1; 1 0])
    σy = sparse(ComplexF64[0 -im; im 0])
    σz = sparse(ComplexF64[1 0; 0 -1])
    I2 = sparse(ComplexF64[1 0; 0 1])
    return σx, σy, σz, I2
end

"""
    single_site_pauli(α::Symbol, site::Int, L::Int)

Build single-site Pauli operator σ_α^[site] in full 2^L Hilbert space.
"""
function single_site_pauli(α::Symbol, site::Int, L::Int)
    σx, σy, σz, I2 = pauli_matrices()
    pauli = α == :X ? σx : α == :Y ? σy : α == :Z ? σz : error("Unknown Pauli: $α")
    ops = [i == site ? pauli : I2 for i in 1:L]
    return reduce(kron, ops)
end

"""
    apply_dephasing_channel(ρ::AbstractMatrix, p::Real, α::Symbol, L::Int)

Apply uniform dephasing channel to density matrix.

N_{p,α}(ρ) = ∏_j [(1-p/2)ρ + (p/2) σ_α^[j] ρ σ_α^[j]]

Args:
    ρ: density matrix in full 2^L space
    p: dephasing strength ∈ [0,1]
    α: Pauli type (:X, :Y, or :Z)
    L: number of sites
"""
function apply_dephasing_channel(ρ::AbstractMatrix, p::Real, α::Symbol, L::Int)
    ρ_out = copy(ρ)
    for j in 1:L
        σj = single_site_pauli(α, j, L)
        ρ_out = (1 - p/2) * ρ_out + (p/2) * σj * ρ_out * σj'
    end
    return ρ_out
end

#=============================================================================
# Coherent Information (using existing ee function pattern)
=============================================================================#

"""
    partial_trace_first(ρ_RQ::AbstractMatrix, d_R::Int, d_Q::Int)

Partial trace over the FIRST subsystem (R), keeping second (Q).
ρ_RQ is in R⊗Q ordering.
"""
function partial_trace_first(ρ_RQ::AbstractMatrix, d_R::Int, d_Q::Int)
    @assert size(ρ_RQ, 1) == d_R * d_Q "Dimension mismatch: expected $(d_R * d_Q), got $(size(ρ_RQ, 1))"
    
    ρ_Q = zeros(ComplexF64, d_Q, d_Q)
    for i in 1:d_Q, j in 1:d_Q
        for k in 1:d_R
            row_idx = (k-1)*d_Q + i
            col_idx = (k-1)*d_Q + j
            ρ_Q[i, j] += ρ_RQ[row_idx, col_idx]
        end
    end
    return ρ_Q
end

"""
    von_neumann_entropy(ρ::AbstractMatrix; tol::Real=1e-12)

Compute von Neumann entropy S(ρ) = -Tr(ρ log ρ).
Similar to existing ee() but works on any density matrix.
"""
function von_neumann_entropy(ρ::AbstractMatrix; tol::Real=1e-12)
    ρ_herm = (ρ + ρ') / 2
    λs = eigvals(Hermitian(Matrix(ρ_herm)))
    λs_pos = filter(λ -> real(λ) > tol, λs)
    isempty(λs_pos) && return 0.0
    return -sum(λ -> real(λ) * log(real(λ)), λs_pos)
end

"""
    coherent_information(ρ_RQ::AbstractMatrix, d_R::Int, d_Q::Int)

Compute coherent information I_c(R⟩Q) = S(ρ_Q) - S(ρ_RQ).

Args:
    ρ_RQ: joint density matrix (d_R * d_Q, d_R * d_Q)
    d_R: reference system dimension
    d_Q: code system dimension
"""
function coherent_information(ρ_RQ::AbstractMatrix, d_R::Int, d_Q::Int)
    ρ_Q = partial_trace_first(ρ_RQ, d_R, d_Q)
    S_Q = von_neumann_entropy(ρ_Q)
    S_RQ = von_neumann_entropy(ρ_RQ)
    return S_Q - S_RQ
end

#=============================================================================
# Scar and Thermal Encoding (using existing PXP_basis)
=============================================================================#

"""
    neel_state_indices(L::Int)

Get indices of Néel states |Z2⟩ and |Z2'⟩ in full 2^L computational basis.
|Z2⟩ = |101010...⟩, |Z2'⟩ = |010101...⟩
"""
function neel_state_indices(L::Int)
    z2 = sum(1 << i for i in 0:2:(L-1))      # |101010...⟩
    z2p = sum(1 << i for i in 1:2:(L-1))     # |010101...⟩
    return z2 + 1, z2p + 1  # +1 for Julia 1-based indexing
end

"""
    prepare_scar_encoding(L::Int)

Prepare reference-code entangled state for scar encoding:
|ψ_RQ⟩ = (1/√2)(|0⟩_R ⊗ |Z2⟩_Q + |1⟩_R ⊗ |Z2'⟩_Q)

Returns: (ψ_RQ, d_R, d_Q)
"""
function prepare_scar_encoding(L::Int)
    d_R = 2
    d_Q = 2^L  # Full Hilbert space dimension
    
    idx_z2, idx_z2p = neel_state_indices(L)
    
    ψ_RQ = zeros(ComplexF64, d_R * d_Q)
    ψ_RQ[idx_z2] = 1/sqrt(2)           # |0⟩_R ⊗ |Z2⟩_Q
    ψ_RQ[d_Q + idx_z2p] = 1/sqrt(2)    # |1⟩_R ⊗ |Z2'⟩_Q
    
    return ψ_RQ, d_R, d_Q
end

"""
    prepare_thermal_encoding(L::Int, th1::Vector, th2::Vector)

Prepare reference-code entangled state for thermal encoding:
|ψ_RQ^th⟩ = (1/√2)(|0⟩_R ⊗ |th1⟩_Q + |1⟩_R ⊗ |th2⟩_Q)

Args:
    L: system size
    th1, th2: orthogonal thermal eigenstates at E≈0

Returns: (ψ_RQ, d_R, d_Q)
"""
function prepare_thermal_encoding(L::Int, th1::Vector, th2::Vector)
    d_R = 2
    d_Q = length(th1)
    @assert d_Q == 2^L "Thermal states should be in full 2^L space"
    
    ψ_RQ = zeros(ComplexF64, d_R * d_Q)
    ψ_RQ[1:d_Q] = th1 / sqrt(2)
    ψ_RQ[d_Q+1:2*d_Q] = th2 / sqrt(2)
    
    return ψ_RQ, d_R, d_Q
end

"""
    reference_code_state(ψ_RQ::Vector)

Create density matrix from pure reference-code state.
"""
reference_code_state(ψ_RQ::Vector) = ψ_RQ * ψ_RQ'

#=============================================================================
# Knill-Laflamme Analysis
=============================================================================#

"""
    knill_laflamme_coefficient(L::Int, α::Symbol)

Compute the b(L) coefficient measuring approximate KL condition violation.

b(L) = (1/D²) Σ_j [D·Tr(σ_j P σ_j P) - Tr(σ_j P)²]

where P = |Z2⟩⟨Z2| + |Z2'⟩⟨Z2'| is the code projector.
"""
function knill_laflamme_coefficient(L::Int, α::Symbol)
    idx_z2, idx_z2p = neel_state_indices(L)
    D = 2  # Code dimension
    d = 2^L
    
    b_sum = 0.0
    for j in 1:L
        σj = Matrix(single_site_pauli(α, j, L))
        
        # P = |Z2⟩⟨Z2| + |Z2'⟩⟨Z2'|
        # σ_j P σ_j elements
        σP_z2 = σj[:, idx_z2]
        σP_z2p = σj[:, idx_z2p]
        
        # Tr(σ_j P σ_j P) = ⟨Z2|σ_j P σ_j|Z2⟩ + ⟨Z2'|σ_j P σ_j|Z2'⟩
        tr_σPσP = abs2(σP_z2[idx_z2]) + abs2(σP_z2[idx_z2p]) + 
                  abs2(σP_z2p[idx_z2]) + abs2(σP_z2p[idx_z2p])
        
        # Tr(σ_j P) = ⟨Z2|σ_j|Z2⟩ + ⟨Z2'|σ_j|Z2'⟩
        tr_σP = σj[idx_z2, idx_z2] + σj[idx_z2p, idx_z2p]
        
        b_sum += D * tr_σPσP - abs2(tr_σP)
    end
    
    return real(b_sum) / D^2
end
```

- [ ] **Step 2: Update PXPConstrained.jl to include QEC**

Add to `src/PXPConstrained.jl`:

```julia
# After other includes
include("QEC.jl")

# Add to exports
export coherent_information, apply_dephasing_channel
export prepare_scar_encoding, prepare_thermal_encoding
export knill_laflamme_coefficient, reference_code_state
export von_neumann_entropy, partial_trace_first
```

- [ ] **Step 3: Commit**

```bash
git add src/QEC.jl src/PXPConstrained.jl
git commit -m "feat(qec): add quantum error correction module

Implements coherent information, dephasing channels, and Knill-Laflamme
analysis for scar QEC codes. Built on existing PXPConstrained infrastructure."
```

---

## Task 2: Main Scan Driver

**Files:**
- Create: `exm/scar_qec/run_qec_scan.jl`

- [ ] **Step 1: Create parameter sweep script**

```julia
# exm/scar_qec/run_qec_scan.jl
"""
Main driver for QEC coherent information scans.
"""

using PXPConstrained
using LinearAlgebra
using JLD2
using Printf

function run_coherent_info_scan(L::Int; 
                                 p_values=0.0:0.02:1.0, 
                                 noise_types=[:X, :Y, :Z],
                                 save_data=true)
    
    println("="^60)
    println("Scar QEC Coherent Information Scan")
    println("System size L = $L")
    println("="^60)
    
    # Prepare scar encoding
    ψ_scar, d_R, d_Q = prepare_scar_encoding(L)
    ρ_scar_0 = reference_code_state(ψ_scar)
    
    # Results storage
    results = Dict{Symbol, Vector{Float64}}()
    
    for α in noise_types
        println("\nProcessing $α-dephasing...")
        I_c_values = Float64[]
        
        for p in p_values
            ρ_noisy = apply_dephasing_channel(ρ_scar_0, p, α, L)
            I_c = coherent_information(ρ_noisy, d_R, d_Q)
            push!(I_c_values, I_c)
            
            if p ≈ 0.0 || p ≈ 0.5 || p ≈ 1.0
                @printf("  p = %.2f: I_c = %.6f\n", p, I_c)
            end
        end
        
        results[α] = I_c_values
    end
    
    # Compute b(L) coefficient
    b_coeffs = Dict{Symbol, Float64}()
    for α in noise_types
        b_coeffs[α] = knill_laflamme_coefficient(L, α)
        @printf("\nb(L=%d, %s) = %.6f\n", L, α, b_coeffs[α])
    end
    
    if save_data
        datadir = joinpath(@__DIR__, "data")
        mkpath(datadir)
        filename = joinpath(datadir, "coherent_info_scar_L$(L).jld2")
        @save filename L p_values results b_coeffs
        println("\nData saved to: $filename")
    end
    
    return results, b_coeffs
end

# Run for multiple system sizes
function main()
    for L in [4, 6, 8, 10]  # Start small for testing
        @time run_coherent_info_scan(L)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
```

- [ ] **Step 2: Commit**

```bash
git add exm/scar_qec/run_qec_scan.jl
git commit -m "feat(qec): add coherent information scan driver

Parameter sweep over dephasing strength for scar encoding.
Computes I_c(p) and b(L) coefficients."
```

---

## Task 3: Plotting and Analysis

**Files:**
- Create: `exm/scar_qec/plot_qec_results.jl`

- [ ] **Step 1: Create plotting script**

```julia
# exm/scar_qec/plot_qec_results.jl
"""
Plotting and analysis for QEC results.
"""

using JLD2
using Printf

function load_and_analyze(L_values::Vector{Int})
    datadir = joinpath(@__DIR__, "data")
    
    all_results = Dict{Int, Any}()
    all_b_coeffs = Dict{Int, Dict{Symbol, Float64}}()
    
    for L in L_values
        filename = joinpath(datadir, "coherent_info_scar_L$(L).jld2")
        if isfile(filename)
            @load filename L p_values results b_coeffs
            all_results[L] = (p_values=p_values, I_c=results)
            all_b_coeffs[L] = b_coeffs
            println("Loaded L=$L")
        else
            println("Warning: $filename not found")
        end
    end
    
    # Print summary
    println("\n" * "="^60)
    println("b(L) Coefficient Summary")
    println("="^60)
    println("L\tb(L,X)\t\tb(L,Y)\t\tb(L,Z)")
    for L in sort(collect(keys(all_b_coeffs)))
        b = all_b_coeffs[L]
        @printf("%d\t%.6f\t%.6f\t%.6f\n", L, b[:X], b[:Y], b[:Z])
    end
    
    # Estimate effective scaling dimension
    if length(all_b_coeffs) >= 2
        println("\nScaling Analysis:")
        L_vals = sort(collect(keys(all_b_coeffs)))
        for α in [:X, :Y, :Z]
            b_vals = [all_b_coeffs[L][α] for L in L_vals]
            # Fit b(L) ~ L^(1-2Δ)
            if all(b_vals .> 0)
                log_L = log.(L_vals)
                log_b = log.(b_vals)
                # Linear regression
                n = length(log_L)
                slope = (n * sum(log_L .* log_b) - sum(log_L) * sum(log_b)) / 
                        (n * sum(log_L.^2) - sum(log_L)^2)
                Δ_eff = (1 - slope) / 2
                @printf("  %s-dephasing: Δ_eff ≈ %.3f (threshold if Δ > 0.5)\n", α, Δ_eff)
            end
        end
    end
    
    return all_results, all_b_coeffs
end

if abspath(PROGRAM_FILE) == @__FILE__
    load_and_analyze([4, 6, 8, 10])
end
```

- [ ] **Step 2: Commit**

```bash
git add exm/scar_qec/plot_qec_results.jl
git commit -m "feat(qec): add analysis and plotting script

Loads QEC scan results, computes scaling dimension from b(L)."
```

---

## Task 4: Run Full Computation

- [ ] **Step 1: Execute main scan**

```bash
cd /Users/cycling/Documents/projects/PXPConstrained
julia --project=. exm/scar_qec/run_qec_scan.jl
```

- [ ] **Step 2: Generate analysis**

```bash
julia --project=. exm/scar_qec/plot_qec_results.jl
```

- [ ] **Step 3: Commit results**

```bash
git add exm/scar_qec/data/
git commit -m "data(qec): add coherent information scan results

System sizes L=4,6,8,10
Dephasing types: X, Y, Z"
```

---

## Summary

| Task | Description | Files |
|------|-------------|-------|
| 1 | QEC Core Module | `src/QEC.jl`, update `src/PXPConstrained.jl` |
| 2 | Main Scan Driver | `exm/scar_qec/run_qec_scan.jl` |
| 3 | Analysis/Plotting | `exm/scar_qec/plot_qec_results.jl` |
| 4 | Run Computation | Execute scripts, save data |

**Key Integration Points:**
- Uses existing `PXP_basis()` for constrained basis
- Uses existing `ee()` pattern for entropy computation
- Uses existing `rdm_PXP()` approach for partial traces
- Follows existing module structure and documentation style
