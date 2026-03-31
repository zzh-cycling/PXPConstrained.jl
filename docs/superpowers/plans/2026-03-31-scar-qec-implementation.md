# Scar QEC Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement numerical computation of coherent information for scar vs thermal encodings under local dephasing, demonstrating that quantum many-body scars form approximate quantum error correcting codes.

**Architecture:** Build on existing PXP exact diagonalization infrastructure. Add coherent information computation module, dephasing channel implementation, and scaling analysis. Output data for plotting threshold behavior and scaling collapse.

**Tech Stack:** Julia, JLD for data storage, existing PXP Hamiltonian code

---

## File Structure

```
fig_plot_code/
├── scar_qec/
│   ├── coherent_information.jl    # Core I_c computation
│   ├── dephasing_channel.jl       # Noise channel implementation  
│   ├── scar_thermal_encoding.jl   # Code space preparation
│   ├── knill_laflamme.jl          # KL condition & b(L) computation
│   ├── run_qec_scan.jl            # Parameter sweep driver
│   └── plot_qec_results.jl        # Visualization
data/
├── scar_qec/
│   ├── coherent_info_scar_L{N}.jld
│   ├── coherent_info_thermal_L{N}.jld
│   └── scaling_collapse_data.jld
```

---

## Task 1: Dephasing Channel Implementation

**Files:**
- Create: `fig_plot_code/scar_qec/dephasing_channel.jl`

- [ ] **Step 1: Create module with single-site dephasing**

```julia
# fig_plot_code/scar_qec/dephasing_channel.jl
module DephasingChannel

using LinearAlgebra
using SparseArrays

export apply_single_site_dephasing, apply_uniform_dephasing, pauli_matrices

"""
Return Pauli matrices as sparse matrices
"""
function pauli_matrices()
    σx = sparse([0.0+0im 1.0; 1.0 0.0])
    σy = sparse([0.0+0im -1im; 1im 0.0])
    σz = sparse([1.0+0im 0.0; 0.0 -1.0])
    I2 = sparse([1.0+0im 0.0; 0.0 1.0])
    return σx, σy, σz, I2
end

"""
Apply single-site dephasing channel to density matrix.

N_p(ρ) = (1 - p/2)ρ + (p/2) σ_α ρ σ_α

Args:
    ρ: density matrix (can be for full system or subsystem)
    p: dephasing strength ∈ [0,1]
    σ: Pauli operator for dephasing (already in full Hilbert space)
"""
function apply_single_site_dephasing(ρ::AbstractMatrix, p::Real, σ::AbstractMatrix)
    return (1 - p/2) * ρ + (p/2) * σ * ρ * σ'
end

"""
Build single-site Pauli operator in full Hilbert space.

Args:
    pauli: single-site Pauli matrix (2x2)
    site: site index (1-based)
    L: total number of sites
    dim_local: local Hilbert space dimension (default 2)
"""
function single_site_operator(pauli::AbstractMatrix, site::Int, L::Int; dim_local::Int=2)
    ops = [i == site ? pauli : sparse(I, dim_local, dim_local) for i in 1:L]
    return reduce(kron, ops)
end

"""
Apply uniform dephasing to all sites.

N_{p,α}(ρ) = ⊗_j [(1-p/2)ρ + (p/2) σ_α^[j] ρ σ_α^[j]]

For efficiency, we apply site by site sequentially.

Args:
    ρ: density matrix
    p: dephasing strength
    α: Pauli type (:X, :Y, or :Z)
    L: number of sites
"""
function apply_uniform_dephasing(ρ::AbstractMatrix, p::Real, α::Symbol, L::Int)
    σx, σy, σz, _ = pauli_matrices()
    
    pauli = if α == :X
        σx
    elseif α == :Y
        σy
    elseif α == :Z
        σz
    else
        error("Unknown Pauli type: $α. Use :X, :Y, or :Z")
    end
    
    ρ_out = copy(ρ)
    for j in 1:L
        σj = single_site_operator(pauli, j, L)
        ρ_out = apply_single_site_dephasing(ρ_out, p, σj)
    end
    
    return ρ_out
end

end # module
```

- [ ] **Step 2: Test dephasing channel**

Create test file and verify:

```julia
# Test in REPL
include("fig_plot_code/scar_qec/dephasing_channel.jl")
using .DephasingChannel
using LinearAlgebra

# Test: p=0 should leave state unchanged
L = 2
ρ = [1.0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0] .+ 0im  # |00⟩⟨00|
ρ_out = apply_uniform_dephasing(ρ, 0.0, :Z, L)
@assert ρ ≈ ρ_out "p=0 should leave state unchanged"

# Test: p=1 Z-dephasing should fully dephase off-diagonals in Z basis
ρ_superpos = 0.5 * ones(ComplexF64, 4, 4)  # maximally mixed in computational basis
ρ_out = apply_uniform_dephasing(ρ_superpos, 1.0, :Z, L)
# After full Z-dephasing, off-diagonal elements should vanish
@assert maximum(abs.(ρ_out - Diagonal(diag(ρ_out)))) < 1e-10 "Full Z-dephasing should diagonalize"

println("✓ Dephasing channel tests passed")
```

- [ ] **Step 3: Commit**

```bash
git add fig_plot_code/scar_qec/dephasing_channel.jl
git commit -m "feat(qec): add dephasing channel implementation

Implements single-site and uniform dephasing channels for X, Y, Z noise.
Follows Sang-Zou convention: N_p(ρ) = (1-p/2)ρ + (p/2)σρσ"
```

---

## Task 2: Coherent Information Computation

**Files:**
- Create: `fig_plot_code/scar_qec/coherent_information.jl`

- [ ] **Step 1: Implement von Neumann entropy**

```julia
# fig_plot_code/scar_qec/coherent_information.jl
module CoherentInformation

using LinearAlgebra

export von_neumann_entropy, coherent_information, partial_trace

"""
Compute von Neumann entropy S(ρ) = -Tr(ρ log ρ)

Args:
    ρ: density matrix
    tol: eigenvalue tolerance (eigenvalues below this treated as 0)
"""
function von_neumann_entropy(ρ::AbstractMatrix; tol::Real=1e-12)
    # Ensure Hermitian
    ρ_herm = (ρ + ρ') / 2
    
    # Get eigenvalues
    λs = eigvals(Hermitian(Matrix(ρ_herm)))
    
    # Filter small/negative eigenvalues
    λs_pos = filter(λ -> λ > tol, real.(λs))
    
    if isempty(λs_pos)
        return 0.0
    end
    
    # S = -Σ λ log(λ)
    return -sum(λ -> λ * log(λ), λs_pos)
end

"""
Partial trace over subsystem B, keeping subsystem A.

For ρ_AB with dimensions (d_A * d_B, d_A * d_B),
returns ρ_A with dimensions (d_A, d_A).

Assumes tensor product ordering: A ⊗ B
"""
function partial_trace(ρ_AB::AbstractMatrix, d_A::Int, d_B::Int)
    @assert size(ρ_AB, 1) == d_A * d_B "Dimension mismatch"
    
    ρ_A = zeros(ComplexF64, d_A, d_A)
    
    for i in 1:d_A, j in 1:d_A
        for k in 1:d_B
            # Index in full space: (i-1)*d_B + k for row, (j-1)*d_B + k for col
            ρ_A[i, j] += ρ_AB[(i-1)*d_B + k, (j-1)*d_B + k]
        end
    end
    
    return ρ_A
end

"""
Compute coherent information I_c(R⟩Q) = S(ρ_Q) - S(ρ_RQ)

For a reference-code state |ψ_RQ⟩, after applying noise channel to Q:
- ρ_RQ = N(|ψ_RQ⟩⟨ψ_RQ|)  
- ρ_Q = Tr_R(ρ_RQ)

Args:
    ρ_RQ: joint density matrix after noise (d_R * d_Q, d_R * d_Q)
    d_R: dimension of reference system
    d_Q: dimension of code (physical) system
"""
function coherent_information(ρ_RQ::AbstractMatrix, d_R::Int, d_Q::Int)
    # Partial trace over R to get ρ_Q
    # Note: we need trace over R, so swap the order in partial_trace
    # ρ_RQ is in R⊗Q basis, we want to trace out R
    
    # Reshape and trace
    ρ_Q = partial_trace_first(ρ_RQ, d_R, d_Q)
    
    S_Q = von_neumann_entropy(ρ_Q)
    S_RQ = von_neumann_entropy(ρ_RQ)
    
    return S_Q - S_RQ
end

"""
Partial trace over the FIRST subsystem (R), keeping second (Q).
"""
function partial_trace_first(ρ_RQ::AbstractMatrix, d_R::Int, d_Q::Int)
    @assert size(ρ_RQ, 1) == d_R * d_Q "Dimension mismatch"
    
    ρ_Q = zeros(ComplexF64, d_Q, d_Q)
    
    for i in 1:d_Q, j in 1:d_Q
        for k in 1:d_R
            # In R⊗Q ordering: index = (r-1)*d_Q + q
            row_idx = (k-1)*d_Q + i
            col_idx = (k-1)*d_Q + j
            ρ_Q[i, j] += ρ_RQ[row_idx, col_idx]
        end
    end
    
    return ρ_Q
end

end # module
```

- [ ] **Step 2: Test coherent information**

```julia
include("fig_plot_code/scar_qec/coherent_information.jl")
using .CoherentInformation
using LinearAlgebra

# Test 1: Pure maximally entangled state should have I_c = log(d_R)
d_R, d_Q = 2, 2
# |Φ⟩ = (|00⟩ + |11⟩)/√2
ψ = [1, 0, 0, 1] / sqrt(2)
ρ_RQ = ψ * ψ'

I_c = coherent_information(ρ_RQ, d_R, d_Q)
@assert abs(I_c - log(2)) < 1e-10 "Pure Bell state should have I_c = log(2)"

# Test 2: Product state should have I_c = 0
ψ_prod = kron([1, 0], [1, 0])  # |0⟩|0⟩
ρ_prod = ψ_prod * ψ_prod'
I_c_prod = coherent_information(ρ_prod, d_R, d_Q)
@assert abs(I_c_prod) < 1e-10 "Product state should have I_c = 0"

# Test 3: Completely mixed state should have I_c ≤ 0
ρ_mixed = Matrix{ComplexF64}(I(4)) / 4
I_c_mixed = coherent_information(ρ_mixed, d_R, d_Q)
@assert I_c_mixed ≤ 1e-10 "Mixed state should have I_c ≤ 0"

println("✓ Coherent information tests passed")
```

- [ ] **Step 3: Commit**

```bash
git add fig_plot_code/scar_qec/coherent_information.jl
git commit -m "feat(qec): add coherent information computation

Implements von Neumann entropy and I_c = S(Q) - S(RQ).
Follows Sang-Zou definition for QEC threshold analysis."
```

---

## Task 3: Scar and Thermal State Encoding

**Files:**
- Create: `fig_plot_code/scar_qec/scar_thermal_encoding.jl`
- Reference: existing FSA code in `fig_plot_code/FSA_scar_thermal_ee_overlap.jl`

- [ ] **Step 1: Implement PXP constrained Hilbert space**

```julia
# fig_plot_code/scar_qec/scar_thermal_encoding.jl
module ScarThermalEncoding

using LinearAlgebra
using SparseArrays

export generate_pxp_basis, neel_state_z2, neel_state_z2prime
export prepare_scar_encoding, prepare_thermal_encoding
export reference_code_state

"""
Generate PXP constrained Hilbert space basis (no adjacent excitations).
Returns vector of valid basis states as integers.
"""
function generate_pxp_basis(L::Int)
    basis = Int[]
    for state in 0:(2^L - 1)
        valid = true
        for i in 0:(L-1)
            # Check adjacent bits (with PBC)
            bit_i = (state >> i) & 1
            bit_next = (state >> ((i+1) % L)) & 1
            if bit_i == 1 && bit_next == 1
                valid = false
                break
            end
        end
        if valid
            push!(basis, state)
        end
    end
    return basis
end

"""
Get index of Néel state |Z2⟩ = |101010...⟩ in PXP basis
"""
function neel_state_z2(L::Int, basis::Vector{Int})
    # |Z2⟩ = |101010...⟩ = alternating 1,0,1,0,...
    z2 = sum(1 << i for i in 0:2:(L-1))
    idx = findfirst(==(z2), basis)
    return idx
end

"""
Get index of Néel state |Z2'⟩ = |010101...⟩ in PXP basis
"""
function neel_state_z2prime(L::Int, basis::Vector{Int})
    # |Z2'⟩ = |010101...⟩ = alternating 0,1,0,1,...
    z2p = sum(1 << i for i in 1:2:(L-1))
    idx = findfirst(==(z2p), basis)
    return idx
end

"""
Prepare scar encoding: reference-code entangled state

|ψ_RQ⟩ = (1/√2)(|0⟩_R ⊗ |Z2⟩_Q + |1⟩_R ⊗ |Z2'⟩_Q)

Returns:
    ψ_RQ: state vector in R⊗Q space
    d_R: reference dimension (2)
    d_Q: code space dimension (PXP Hilbert space)
"""
function prepare_scar_encoding(L::Int)
    basis = generate_pxp_basis(L)
    d_Q = length(basis)
    d_R = 2
    
    # Get Néel state indices
    idx_z2 = neel_state_z2(L, basis)
    idx_z2p = neel_state_z2prime(L, basis)
    
    @assert idx_z2 !== nothing "Z2 state not found in basis"
    @assert idx_z2p !== nothing "Z2' state not found in basis"
    
    # Construct entangled state in R⊗Q space
    # |0⟩_R = [1,0], |1⟩_R = [0,1]
    # Full dimension: d_R * d_Q
    
    ψ_RQ = zeros(ComplexF64, d_R * d_Q)
    
    # |0⟩_R ⊗ |Z2⟩_Q: index = 0*d_Q + idx_z2 = idx_z2
    ψ_RQ[idx_z2] = 1/sqrt(2)
    
    # |1⟩_R ⊗ |Z2'⟩_Q: index = 1*d_Q + idx_z2p = d_Q + idx_z2p
    ψ_RQ[d_Q + idx_z2p] = 1/sqrt(2)
    
    return ψ_RQ, d_R, d_Q, basis
end

"""
Prepare thermal encoding: reference-code entangled state with thermal eigenstates

|ψ_RQ^th⟩ = (1/√2)(|0⟩_R ⊗ |th1⟩_Q + |1⟩_R ⊗ |th2⟩_Q)

Args:
    L: system size
    thermal_states: 2 x d_Q matrix where columns are thermal eigenstates at E≈0
"""
function prepare_thermal_encoding(L::Int, thermal_state1::Vector, thermal_state2::Vector)
    d_Q = length(thermal_state1)
    d_R = 2
    
    # Construct entangled state
    ψ_RQ = zeros(ComplexF64, d_R * d_Q)
    
    # |0⟩_R ⊗ |th1⟩_Q
    ψ_RQ[1:d_Q] = thermal_state1 / sqrt(2)
    
    # |1⟩_R ⊗ |th2⟩_Q  
    ψ_RQ[d_Q+1:2*d_Q] = thermal_state2 / sqrt(2)
    
    return ψ_RQ, d_R, d_Q
end

"""
Create density matrix from pure state
"""
function reference_code_state(ψ_RQ::Vector)
    return ψ_RQ * ψ_RQ'
end

end # module
```

- [ ] **Step 2: Test encoding**

```julia
include("fig_plot_code/scar_qec/scar_thermal_encoding.jl")
using .ScarThermalEncoding
using LinearAlgebra

# Test PXP basis generation
L = 6
basis = generate_pxp_basis(L)
println("L=$L: PXP basis size = $(length(basis)) (should be Fibonacci-like)")

# Test Néel states
idx_z2 = neel_state_z2(L, basis)
idx_z2p = neel_state_z2prime(L, basis)
println("Z2 index: $idx_z2, Z2' index: $idx_z2p")
@assert idx_z2 !== nothing && idx_z2p !== nothing

# Test scar encoding
ψ_RQ, d_R, d_Q, basis = prepare_scar_encoding(L)
@assert length(ψ_RQ) == d_R * d_Q
@assert abs(norm(ψ_RQ) - 1.0) < 1e-10 "State should be normalized"

ρ_RQ = reference_code_state(ψ_RQ)
@assert abs(tr(ρ_RQ) - 1.0) < 1e-10 "Density matrix should have trace 1"

println("✓ Scar encoding tests passed")
```

- [ ] **Step 3: Commit**

```bash
git add fig_plot_code/scar_qec/scar_thermal_encoding.jl
git commit -m "feat(qec): add scar and thermal state encoding

Implements PXP basis generation and Néel state encoding.
Prepares reference-code entangled states for I_c computation."
```

---

## Task 4: Knill-Laflamme Condition Analysis

**Files:**
- Create: `fig_plot_code/scar_qec/knill_laflamme.jl`

- [ ] **Step 1: Implement KL matrix elements and b(L) coefficient**

```julia
# fig_plot_code/scar_qec/knill_laflamme.jl
module KnillLaflamme

using LinearAlgebra
using SparseArrays

export compute_kl_matrix_elements, compute_b_coefficient
export analyze_kl_scaling

include("dephasing_channel.jl")
using .DephasingChannel: pauli_matrices, single_site_operator

"""
Compute Knill-Laflamme matrix elements for code states.

⟨ψ_i|E|ψ_j⟩ for error operators E = σ_α^[j]

Args:
    code_states: vector of code state vectors [|ψ_1⟩, |ψ_2⟩, ...]
    L: system size
    α: Pauli type (:X, :Y, :Z)
    
Returns:
    kl_matrix: array of size (n_states, n_states, L) where
               kl_matrix[i,j,site] = ⟨ψ_i|σ_α^[site]|ψ_j⟩
"""
function compute_kl_matrix_elements(code_states::Vector{<:Vector}, L::Int, α::Symbol)
    n_states = length(code_states)
    σx, σy, σz, _ = pauli_matrices()
    
    pauli = if α == :X
        σx
    elseif α == :Y
        σy
    elseif α == :Z
        σz
    else
        error("Unknown Pauli type: $α")
    end
    
    kl_matrix = zeros(ComplexF64, n_states, n_states, L)
    
    for site in 1:L
        σ_site = single_site_operator(pauli, site, L)
        for i in 1:n_states, j in 1:n_states
            kl_matrix[i, j, site] = code_states[i]' * σ_site * code_states[j]
        end
    end
    
    return kl_matrix
end

"""
Compute the b(L) coefficient from Sang-Zou.

b(L) = (1/D²) Σ_j [D·Tr(σ_j P σ_j P) - Tr(σ_j P)·Tr(σ_j P)]

where P is the code projector and D = dim(code space).

This determines threshold: b(L) → 0 implies threshold exists.
Scaling b(L) ∝ L^{1-2Δ} gives effective scaling dimension.

Args:
    code_states: orthonormal code state vectors
    L: system size
    α: Pauli type
"""
function compute_b_coefficient(code_states::Vector{<:Vector}, L::Int, α::Symbol)
    D = length(code_states)
    kl = compute_kl_matrix_elements(code_states, L, α)
    
    b = 0.0
    
    for j in 1:L
        # Tr(σ_j P σ_j P) = Σ_{i,k} |⟨ψ_i|σ_j|ψ_k⟩|²
        tr_σPσP = sum(abs2(kl[i, k, j]) for i in 1:D, k in 1:D)
        
        # Tr(σ_j P) = Σ_i ⟨ψ_i|σ_j|ψ_i⟩
        tr_σP = sum(kl[i, i, j] for i in 1:D)
        
        b += D * tr_σPσP - abs2(tr_σP)
    end
    
    b /= D^2
    return real(b)
end

"""
Compute ε_{ij} error terms for approximate KL conditions.

ε_{ij} = ⟨ψ_i|E|ψ_j⟩ - C_E δ_{ij}

where C_E = (1/D) Σ_i ⟨ψ_i|E|ψ_i⟩
"""
function compute_kl_errors(code_states::Vector{<:Vector}, L::Int, α::Symbol)
    D = length(code_states)
    kl = compute_kl_matrix_elements(code_states, L, α)
    
    errors = zeros(Float64, D, D, L)
    
    for j in 1:L
        # C_E for this site
        C_E = sum(kl[i, i, j] for i in 1:D) / D
        
        for i in 1:D, k in 1:D
            target = i == k ? C_E : 0.0
            errors[i, k, j] = abs(kl[i, k, j] - target)
        end
    end
    
    return errors
end

"""
Analyze KL condition scaling with system size.

Returns max |ε_{ij}| for different L values.
"""
function analyze_kl_scaling(L_values::Vector{Int}, α::Symbol, prepare_code_fn::Function)
    max_errors = Float64[]
    b_values = Float64[]
    
    for L in L_values
        code_states, _ = prepare_code_fn(L)
        
        errors = compute_kl_errors(code_states, L, α)
        push!(max_errors, maximum(errors))
        
        b = compute_b_coefficient(code_states, L, α)
        push!(b_values, b)
    end
    
    return max_errors, b_values
end

end # module
```

- [ ] **Step 2: Test KL computation**

```julia
include("fig_plot_code/scar_qec/knill_laflamme.jl")
include("fig_plot_code/scar_qec/scar_thermal_encoding.jl")
using .KnillLaflamme
using .ScarThermalEncoding

# Prepare scar code states for L=10
L = 10
basis = generate_pxp_basis(L)
d_Q = length(basis)

# Create code states as vectors in PXP basis
idx_z2 = neel_state_z2(L, basis)
idx_z2p = neel_state_z2prime(L, basis)

ψ_z2 = zeros(ComplexF64, d_Q)
ψ_z2[idx_z2] = 1.0

ψ_z2p = zeros(ComplexF64, d_Q)
ψ_z2p[idx_z2p] = 1.0

code_states = [ψ_z2, ψ_z2p]

# Compute b coefficient
b_Z = compute_b_coefficient(code_states, L, :Z)
b_X = compute_b_coefficient(code_states, L, :X)
println("L=$L: b(Z) = $b_Z, b(X) = $b_X")

# KL matrix elements
kl = compute_kl_matrix_elements(code_states, L, :Z)
println("KL diagonal (site 1): ", kl[1,1,1], ", ", kl[2,2,1])
println("KL off-diagonal (site 1): ", kl[1,2,1])

println("✓ Knill-Laflamme tests passed")
```

- [ ] **Step 3: Commit**

```bash
git add fig_plot_code/scar_qec/knill_laflamme.jl
git commit -m "feat(qec): add Knill-Laflamme condition analysis

Computes KL matrix elements and b(L) coefficient.
b(L) scaling determines effective scar 'scaling dimension'."
```

---

## Task 5: Main QEC Scan Driver

**Files:**
- Create: `fig_plot_code/scar_qec/run_qec_scan.jl`

- [ ] **Step 1: Implement parameter sweep**

```julia
# fig_plot_code/scar_qec/run_qec_scan.jl

using JLD
using LinearAlgebra
using Printf

# Include modules
include("dephasing_channel.jl")
include("coherent_information.jl")
include("scar_thermal_encoding.jl")
include("knill_laflamme.jl")

using .DephasingChannel
using .CoherentInformation
using .ScarThermalEncoding
using .KnillLaflamme

"""
Run coherent information scan for scar encoding.

Args:
    L: system size
    p_values: array of dephasing strengths
    α: Pauli type (:X, :Y, :Z)
"""
function run_scar_Ic_scan(L::Int, p_values::Vector{Float64}, α::Symbol)
    println("Running scar I_c scan for L=$L, α=$α")
    
    # Prepare encoding
    ψ_RQ, d_R, d_Q, basis = prepare_scar_encoding(L)
    ρ_RQ_init = reference_code_state(ψ_RQ)
    
    Ic_values = Float64[]
    
    for (i, p) in enumerate(p_values)
        # Apply dephasing to Q subsystem only
        # ρ_RQ is in R⊗Q space, need to apply dephasing to Q part
        ρ_RQ_noisy = apply_dephasing_to_Q(ρ_RQ_init, p, α, L, d_R, d_Q)
        
        # Compute coherent information
        Ic = coherent_information(ρ_RQ_noisy, d_R, d_Q)
        push!(Ic_values, Ic)
        
        if i % 10 == 0
            @printf("  p=%.2f: I_c = %.4f\n", p, Ic)
        end
    end
    
    return Ic_values
end

"""
Apply dephasing channel to Q subsystem of ρ_RQ.

ρ_RQ is (d_R * d_Q) x (d_R * d_Q) matrix in R⊗Q ordering.
We need to apply N_p,α to the Q subsystem.
"""
function apply_dephasing_to_Q(ρ_RQ::AbstractMatrix, p::Real, α::Symbol, L::Int, d_R::Int, d_Q::Int)
    σx, σy, σz, I2 = pauli_matrices()
    
    pauli = if α == :X
        σx
    elseif α == :Y
        σy
    elseif α == :Z
        σz
    else
        error("Unknown Pauli type")
    end
    
    ρ_out = copy(ρ_RQ)
    
    # For each site in Q, apply dephasing
    # The operator on R⊗Q is I_R ⊗ σ_j^Q
    I_R = sparse(I, d_R, d_R)
    
    for j in 1:L
        # Build σ_j in Q-space (PXP constrained)
        σ_j_Q = build_pauli_in_pxp_basis(pauli, j, L)
        
        # Full operator: I_R ⊗ σ_j_Q
        σ_j_full = kron(I_R, σ_j_Q)
        
        # Apply dephasing
        ρ_out = (1 - p/2) * ρ_out + (p/2) * σ_j_full * ρ_out * σ_j_full'
    end
    
    return ρ_out
end

"""
Build Pauli operator in PXP constrained basis.
"""
function build_pauli_in_pxp_basis(pauli_2x2::AbstractMatrix, site::Int, L::Int)
    basis = generate_pxp_basis(L)
    d = length(basis)
    
    σ = zeros(ComplexF64, d, d)
    
    for (i, state_i) in enumerate(basis)
        for (j, state_j) in enumerate(basis)
            # Compute ⟨state_i|σ_site|state_j⟩
            # σ acts on site, flipping/phasing the bit
            
            # Get bits at site
            bit_i = (state_i >> (site-1)) & 1
            bit_j = (state_j >> (site-1)) & 1
            
            # Check if states differ only at site
            mask = ~(1 << (site-1))
            if (state_i & mask) == (state_j & mask)
                # Matrix element from 2x2 Pauli
                σ[i, j] = pauli_2x2[bit_i+1, bit_j+1]
            end
        end
    end
    
    return sparse(σ)
end

"""
Main driver: run full parameter scan
"""
function main()
    # Parameters
    L_values = [10, 12, 14, 16, 18]
    p_values = collect(0.0:0.02:1.0)
    α_values = [:X, :Y, :Z]
    
    # Create output directory
    output_dir = "data/scar_qec"
    mkpath(output_dir)
    
    results = Dict()
    
    for L in L_values
        println("\n" * "="^50)
        println("System size L = $L")
        println("="^50)
        
        results[L] = Dict()
        
        for α in α_values
            println("\nDephasing type: $α")
            
            Ic_scar = run_scar_Ic_scan(L, p_values, α)
            results[L][α] = Ic_scar
            
            # Save intermediate results
            save_file = joinpath(output_dir, "coherent_info_scar_L$(L)_$(α).jld")
            save(save_file, "p_values", p_values, "Ic_values", Ic_scar, "L", L, "alpha", string(α))
            println("Saved to $save_file")
        end
        
        # Also compute b(L) coefficients
        basis = generate_pxp_basis(L)
        d_Q = length(basis)
        idx_z2 = neel_state_z2(L, basis)
        idx_z2p = neel_state_z2prime(L, basis)
        
        ψ_z2 = zeros(ComplexF64, d_Q); ψ_z2[idx_z2] = 1.0
        ψ_z2p = zeros(ComplexF64, d_Q); ψ_z2p[idx_z2p] = 1.0
        code_states = [ψ_z2, ψ_z2p]
        
        for α in α_values
            b = compute_b_coefficient(code_states, L, α)
            println("b($α) = $b")
            results[L][Symbol("b_", α)] = b
        end
    end
    
    # Save all results
    save(joinpath(output_dir, "all_results.jld"), "results", results, 
         "L_values", L_values, "p_values", p_values)
    
    println("\n✓ All scans complete!")
    return results
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
```

- [ ] **Step 2: Test on small system**

```julia
include("fig_plot_code/scar_qec/run_qec_scan.jl")

# Test on L=10
L = 10
p_values = [0.0, 0.1, 0.2, 0.5, 1.0]
Ic = run_scar_Ic_scan(L, p_values, :Z)
println("I_c values: ", Ic)

# Check: I_c(p=0) should be log(2)
@assert abs(Ic[1] - log(2)) < 0.1 "I_c(p=0) should be ≈ log(2)"
println("✓ QEC scan test passed")
```

- [ ] **Step 3: Commit**

```bash
git add fig_plot_code/scar_qec/run_qec_scan.jl
git commit -m "feat(qec): add main parameter scan driver

Implements coherent information scan over dephasing strength.
Computes I_c(p) for scar encoding under X, Y, Z dephasing."
```

---

## Task 6: Plotting and Analysis

**Files:**
- Create: `fig_plot_code/scar_qec/plot_qec_results.jl`

- [ ] **Step 1: Implement plotting functions**

```julia
# fig_plot_code/scar_qec/plot_qec_results.jl

using Plots
using JLD
using LaTeXStrings
using LsqFit

# Plot settings
default(
    tickfont=font(10),
    guidefont=font(14),
    legendfont=font(10),
    grid=false,
    widen=true
)

"""
Plot I_c vs p for multiple system sizes
"""
function plot_Ic_vs_p(data_dir::String, α::Symbol)
    fig = plot(xlabel=L"p", ylabel=L"I_c", title="$(α)-dephasing")
    
    colors = [:blue, :red, :green, :orange, :purple]
    
    files = filter(f -> occursin("coherent_info_scar", f) && occursin("_$(α).jld", f), readdir(data_dir))
    sort!(files)
    
    for (i, file) in enumerate(files)
        data = load(joinpath(data_dir, file))
        p_values = data["p_values"]
        Ic_values = data["Ic_values"]
        L = data["L"]
        
        plot!(p_values, Ic_values, label="L=$L", color=colors[mod1(i, length(colors))], 
              lw=2, marker=:circle, markersize=3)
    end
    
    # Reference lines
    hline!([log(2)], label=L"\log 2", linestyle=:dash, color=:black, lw=1)
    hline!([0], label=nothing, linestyle=:dot, color=:gray, lw=1)
    
    return fig
end

"""
Scaling collapse plot: I_c vs p*L^ν
"""
function plot_scaling_collapse(data_dir::String, α::Symbol, ν::Real)
    fig = plot(xlabel=L"p \cdot L^\nu", ylabel=L"I_c", 
               title="Scaling collapse (ν=$ν)")
    
    colors = [:blue, :red, :green, :orange, :purple]
    
    files = filter(f -> occursin("coherent_info_scar", f) && occursin("_$(α).jld", f), readdir(data_dir))
    sort!(files)
    
    for (i, file) in enumerate(files)
        data = load(joinpath(data_dir, file))
        p_values = data["p_values"]
        Ic_values = data["Ic_values"]
        L = data["L"]
        
        # Scaled x-axis
        x_scaled = p_values .* L^ν
        
        plot!(x_scaled, Ic_values, label="L=$L", color=colors[mod1(i, length(colors))],
              lw=2, marker=:circle, markersize=3)
    end
    
    return fig
end

"""
Plot b(L) coefficient scaling
"""
function plot_b_scaling(L_values::Vector{Int}, b_values::Vector{Float64}, α::Symbol)
    fig = plot(L_values, b_values, xlabel=L"L", ylabel=L"b(L)",
               title="$(α)-dephasing", marker=:circle, lw=2, label=nothing)
    
    # Fit power law: b(L) = A * L^(1-2Δ)
    @. model(x, p) = p[1] * x^p[2]
    fit = curve_fit(model, Float64.(L_values), b_values, [1.0, 0.5])
    
    A, exponent = fit.param
    Δ_eff = (1 - exponent) / 2
    
    L_fit = range(minimum(L_values), maximum(L_values), length=100)
    plot!(L_fit, model(L_fit, fit.param), linestyle=:dash, 
          label=latexstring("\\propto L^{$(round(exponent, digits=2))}, \\Delta_{\\rm eff}=$(round(Δ_eff, digits=2))"))
    
    return fig, Δ_eff
end

"""
Generate all figures
"""
function generate_all_figures(data_dir::String="data/scar_qec")
    # Load results
    all_data = load(joinpath(data_dir, "all_results.jld"))
    results = all_data["results"]
    L_values = all_data["L_values"]
    
    figs_dir = "figs/scar_qec"
    mkpath(figs_dir)
    
    # 1. I_c vs p plots
    for α in [:X, :Y, :Z]
        fig = plot_Ic_vs_p(data_dir, α)
        savefig(fig, joinpath(figs_dir, "Ic_vs_p_$(α).pdf"))
        println("Saved Ic_vs_p_$(α).pdf")
    end
    
    # 2. b(L) scaling plots
    for α in [:X, :Y, :Z]
        b_values = [results[L][Symbol("b_", α)] for L in L_values]
        fig, Δ_eff = plot_b_scaling(L_values, b_values, α)
        savefig(fig, joinpath(figs_dir, "b_scaling_$(α).pdf"))
        println("Saved b_scaling_$(α).pdf, Δ_eff = $Δ_eff")
    end
    
    # 3. Scaling collapse (try different ν values)
    for α in [:X, :Y, :Z]
        for ν in [-1.0, -0.5, 0.0, 0.5, 1.0]
            fig = plot_scaling_collapse(data_dir, α, ν)
            savefig(fig, joinpath(figs_dir, "scaling_collapse_$(α)_nu$(ν).pdf"))
        end
    end
    
    println("\n✓ All figures generated in $figs_dir")
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    generate_all_figures()
end
```

- [ ] **Step 2: Test plotting**

```julia
include("fig_plot_code/scar_qec/plot_qec_results.jl")

# Test with existing data (if available) or create mock data
data_dir = "data/scar_qec"
if isdir(data_dir) && !isempty(readdir(data_dir))
    fig = plot_Ic_vs_p(data_dir, :Z)
    display(fig)
    println("✓ Plotting test passed")
else
    println("No data yet - run main scan first")
end
```

- [ ] **Step 3: Commit**

```bash
git add fig_plot_code/scar_qec/plot_qec_results.jl
git commit -m "feat(qec): add plotting and analysis functions

Implements I_c vs p plots, scaling collapse, and b(L) scaling analysis.
Extracts effective scar 'scaling dimension' from b(L) fit."
```

---

## Task 7: Run Full Computation

- [ ] **Step 1: Execute main scan**

```bash
cd /Users/cycling/Documents/projects/RydbergErgotropy
julia fig_plot_code/scar_qec/run_qec_scan.jl
```

Expected output: Data files in `data/scar_qec/` for L=10,12,14,16,18

- [ ] **Step 2: Generate figures**

```bash
julia fig_plot_code/scar_qec/plot_qec_results.jl
```

Expected output: PDF figures in `figs/scar_qec/`

- [ ] **Step 3: Commit results**

```bash
git add data/scar_qec/ figs/scar_qec/
git commit -m "data(qec): add coherent information scan results

System sizes L=10,12,14,16,18
Dephasing types: X, Y, Z
Includes I_c(p) data and b(L) coefficients"
```

---

## Task 8: Documentation

**Files:**
- Update: `docs/superpowers/specs/2026-03-30-scar-qec-design.md`

- [ ] **Step 1: Add results section to spec**

Add to the spec document:

```markdown
## Numerical Results

### System Sizes Computed
- L = 10, 12, 14, 16, 18

### Key Findings

1. **Coherent information threshold**
   - X-dephasing: p_c = [value]
   - Y-dephasing: p_c = [value]  
   - Z-dephasing: p_c = [value]

2. **Effective scaling dimension**
   - From b(L) fit: Δ_eff = [value]
   - Threshold condition Δ > 1/2: [satisfied/not satisfied]

3. **Comparison with thermal encoding**
   - [Results]
```

- [ ] **Step 2: Commit documentation**

```bash
git add docs/superpowers/specs/2026-03-30-scar-qec-design.md
git commit -m "docs(qec): update spec with numerical results"
```

---

## Summary

| Task | Description | Key Output |
|------|-------------|------------|
| 1 | Dephasing channel | `dephasing_channel.jl` |
| 2 | Coherent information | `coherent_information.jl` |
| 3 | Scar/thermal encoding | `scar_thermal_encoding.jl` |
| 4 | Knill-Laflamme analysis | `knill_laflamme.jl` |
| 5 | Main scan driver | `run_qec_scan.jl` |
| 6 | Plotting | `plot_qec_results.jl` |
| 7 | Run computation | Data & figures |
| 8 | Documentation | Updated spec |
