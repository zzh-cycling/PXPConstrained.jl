"""
    QEC.jl

Quantum Error Correction analysis for PXP scar states.
Implements coherent information, dephasing channels, and Knill-Laflamme analysis
for demonstrating that quantum many-body scars form approximate quantum error 
correcting codes.

Works in the **constrained Hilbert space** (Fibonacci basis) where no two adjacent
sites can be simultaneously excited. Uses reference qubit technique similar to
FibonacciChain.jl for computing coherent information.

Based on:
- Brandão et al., "Quantum Error Correcting Codes in Eigenstates of Translation-Invariant Spin Chains"
- Sang & Zou, "Approximate quantum error correcting codes from conformal field theory"
"""

#=============================================================================
# Extended Basis with Reference Qubit
=============================================================================#

"""
    build_extended_basis(basis::Vector{T}) where {N, T <: BitStr{N}}

Construct the extended basis R⊗Q that prepends 1 reference qubit to the 
constrained PXP basis.

Reference bit occupies the highest-order position (leftmost) in binary representation.

# Arguments
- `basis::Vector{T}`: Sorted PXP basis (constrained Hilbert space)

# Returns
- `Vector{BitStr{N+1}}`: Extended basis with reference qubit prefix

# Example
```julia
basis = PXP_basis(6, true)
ext_basis = build_extended_basis(basis)  # |r⟩ ⊗ |q⟩ for r ∈ {0,1}
```
"""
function build_extended_basis(basis::Vector{T}) where {N, T <: BitStr{N}}
    newT = BitStr{N+1, Int}
    ext_basis = newT[]
    for r in 0:1
        for b in basis
            push!(ext_basis, join(BitStr{1,Int}(r), b))
        end
    end
    return sort(ext_basis)
end

#=============================================================================
# Dephasing Channel in Constrained Space
=============================================================================#

"""
    apply_Z_dephasing_constrained(ρ::AbstractMatrix, p::Real, L::Int, 
                                  basis::Vector{T}) where {T <: BitStr}

Apply uniform Z-dephasing channel to density matrix in constrained basis.

Z-dephasing preserves the constrained subspace since σ_Z is diagonal.
The channel is: N_p(ρ) = ∏_j [(1-p/2)ρ + (p/2) Z_j ρ Z_j]

# Arguments
- `ρ::AbstractMatrix`: Density matrix in constrained basis
- `p::Real`: Dephasing strength ∈ [0,1]
- `L::Int`: Number of sites
- `basis::Vector{T}`: PXP basis states (BitStr)

# Returns
- `Matrix{ComplexF64}`: Noisy density matrix in constrained basis
"""
function apply_Z_dephasing_constrained(ρ::AbstractMatrix, p::Real, L::Int, 
                                       basis::Vector{T}) where {N, T <: BitStr{N}}
    @assert 0 ≤ p ≤ 1 "Dephasing strength p must be in [0,1], got $p"
    @assert size(ρ, 1) == length(basis) "Dimension mismatch"
    
    d = length(basis)
    ρ_out = Matrix{ComplexF64}(ρ)
    
    for j in 1:L
        # Build diagonal Z_j operator in constrained basis
        # Z_j|b⟩ = (-1)^{b_j}|b⟩ where b_j is the j-th bit
        z_diag = [readbit(Int(b), N - j + 1) == 1 ? -1.0 : 1.0 for b in basis]
        
        # Apply dephasing: (1-p/2)ρ + (p/2) Z_j ρ Z_j
        for i in 1:d, k in 1:d
            ρ_out[i, k] = (1 - p/2) * ρ_out[i, k] + (p/2) * z_diag[i] * ρ_out[i, k] * z_diag[k]
        end
    end
    return ρ_out
end

"""
    apply_Z_dephasing_extended(ρ_RQ::AbstractMatrix, p::Real, L::Int,
                               ext_basis::Vector{T}) where {T <: BitStr}

Apply Z-dephasing to the Q (code) subsystem of a joint R⊗Q density matrix.

Acts as I_R ⊗ N_p where N_p is Z-dephasing on Q.

# Arguments
- `ρ_RQ::AbstractMatrix`: Joint density matrix in extended R⊗Q basis
- `p::Real`: Dephasing strength ∈ [0,1]
- `L::Int`: Number of sites in Q
- `ext_basis::Vector{T}`: Extended basis (R⊗Q)

# Returns
- `Matrix{ComplexF64}`: Noisy joint density matrix
"""
function apply_Z_dephasing_extended(ρ_RQ::AbstractMatrix, p::Real, L::Int,
                                    ext_basis::Vector{T}) where {N, T <: BitStr{N}}
    @assert 0 ≤ p ≤ 1 "Dephasing strength p must be in [0,1], got $p"
    @assert size(ρ_RQ, 1) == length(ext_basis) "Dimension mismatch"
    
    d = length(ext_basis)
    ρ_out = Matrix{ComplexF64}(ρ_RQ)
    
    # N = L + 1 (L sites for Q, 1 for R)
    # Reference qubit is at position N (leftmost), Q sites are at positions 1 to L
    for j in 1:L
        # Z_j acts on site j of Q (positions 1 to L in extended basis)
        z_diag = [readbit(Int(b), L - j + 1) == 1 ? -1.0 : 1.0 for b in ext_basis]
        
        for i in 1:d, k in 1:d
            ρ_out[i, k] = (1 - p/2) * ρ_out[i, k] + (p/2) * z_diag[i] * ρ_out[i, k] * z_diag[k]
        end
    end
    return ρ_out
end

#=============================================================================
# Coherent Information Computation
=============================================================================#

"""
    partial_trace_R(ρ_RQ::AbstractMatrix, d_Q::Int)

Partial trace over reference system R (first subsystem), keeping Q.

For the extended basis structure |r⟩⊗|q⟩ with r ∈ {0,1}, traces out R.

# Arguments
- `ρ_RQ::AbstractMatrix`: Joint density matrix in extended basis (2*d_Q, 2*d_Q)
- `d_Q::Int`: Dimension of code system (constrained basis size)

# Returns
- `Matrix{ComplexF64}`: Reduced density matrix ρ_Q (d_Q, d_Q)
"""
function partial_trace_R(ρ_RQ::AbstractMatrix, d_Q::Int)
    d_R = 2
    @assert size(ρ_RQ, 1) == d_R * d_Q "Dimension mismatch: expected $(d_R * d_Q), got $(size(ρ_RQ, 1))"
    
    # Extended basis is ordered as [|0⟩⊗|q_1⟩, |0⟩⊗|q_2⟩, ..., |1⟩⊗|q_1⟩, |1⟩⊗|q_2⟩, ...]
    # Tr_R(ρ_RQ)[i,j] = Σ_r ρ_RQ[r*d_Q + i, r*d_Q + j]
    ρ_Q = zeros(ComplexF64, d_Q, d_Q)
    for r in 0:(d_R-1)
        offset = r * d_Q
        ρ_Q .+= ρ_RQ[offset+1:offset+d_Q, offset+1:offset+d_Q]
    end
    return ρ_Q
end

"""
    von_neumann_entropy(ρ::AbstractMatrix; tol::Real=1e-12)

Compute von Neumann entropy S(ρ) = -Tr(ρ log ρ).

# Arguments
- `ρ::AbstractMatrix`: Density matrix (must be Hermitian)
- `tol::Real`: Eigenvalue tolerance (values below this treated as 0)

# Returns
- `Float64`: Von Neumann entropy
"""
function von_neumann_entropy(ρ::AbstractMatrix; tol::Real=1e-12)
    ρ_herm = (ρ + ρ') / 2  # Ensure Hermitian
    λs = eigvals(Hermitian(Matrix(ρ_herm)))
    λs_pos = filter(λ -> real(λ) > tol, λs)
    isempty(λs_pos) && return 0.0
    return -sum(λ -> real(λ) * log(real(λ)), λs_pos)
end

"""
    coherent_information_constrained(ρ_RQ::AbstractMatrix, d_Q::Int)

Compute coherent information I_c(R⟩Q) = S(ρ_Q) - S(ρ_RQ) for constrained basis.

# Arguments
- `ρ_RQ::AbstractMatrix`: Joint density matrix in extended basis
- `d_Q::Int`: Dimension of code system (constrained basis size)

# Returns
- `Float64`: Coherent information

# Example
```julia
basis = PXP_basis(6, true)
ψ_RQ, ext_basis = prepare_scar_encoding_constrained(6, basis)
ρ_RQ = ψ_RQ * ψ_RQ'
I_c = coherent_information_constrained(ρ_RQ, length(basis))  # Should be log(2)
```
"""
function coherent_information_constrained(ρ_RQ::AbstractMatrix, d_Q::Int)
    ρ_Q = partial_trace_R(ρ_RQ, d_Q)
    S_Q = von_neumann_entropy(ρ_Q)
    S_RQ = von_neumann_entropy(ρ_RQ)
    return S_Q - S_RQ
end

#=============================================================================
# Scar and Thermal State Encoding in Constrained Basis
=============================================================================#

"""
    neel_state_bitstr(L::Int)

Get Néel state BitStr representations.

- |Z2⟩ = |101010...⟩ (excitation on odd sites, counting from right)
- |Z2'⟩ = |010101...⟩ (excitation on even sites)

# Arguments
- `L::Int`: System size (must be even for PBC)

# Returns
- `Tuple{BitStr{L}, BitStr{L}}`: (|Z2⟩, |Z2'⟩) as BitStr
"""
function neel_state_bitstr(L::Int)
    T = BitStr{L, Int}
    # |Z2⟩ = |101010...⟩ has 1 at odd positions (1,3,5,...) counting from right starting at 1
    # In 0-indexed bit positions: 0,2,4,... → sum(1 << i for i in 0:2:(L-1))
    # But BitBasis counts from right, so bit position 0 is rightmost
    # |101010⟩ in binary means bit 1 at positions 1,3,5 (0-indexed)
    z2 = T(sum(1 << i for i in 1:2:(L-1)))      # |101010...⟩ = bits at odd 0-indexed positions
    z2p = T(sum(1 << i for i in 0:2:(L-1)))     # |010101...⟩ = bits at even 0-indexed positions
    return z2, z2p
end

"""
    prepare_scar_encoding_constrained(L::Int, basis::Vector{T}; pbc::Bool=true) where {T <: BitStr}

Prepare reference-code entangled state for scar encoding in constrained basis.

Constructs the state:
|ψ_RQ⟩ = (1/√2)(|0⟩_R ⊗ |Z2⟩_Q + |1⟩_R ⊗ |Z2'⟩_Q)

Works in the constrained Hilbert space using the extended basis R⊗Q.

# Arguments
- `L::Int`: System size
- `basis::Vector{T}`: PXP constrained basis
- `pbc::Bool=true`: Boundary conditions (must match basis)

# Returns
- `ψ_RQ::Vector{ComplexF64}`: Reference-code entangled state in extended basis
- `ext_basis::Vector{BitStr{L+1}}`: Extended basis (R⊗Q)

# Example
```julia
L = 6
basis = PXP_basis(L, true)
ψ_RQ, ext_basis = prepare_scar_encoding_constrained(L, basis)
ρ_RQ = ψ_RQ * ψ_RQ'
I_c = coherent_information_constrained(ρ_RQ, length(basis))  # Should be log(2)
```
"""
function prepare_scar_encoding_constrained(L::Int, basis::Vector{T}; pbc::Bool=true) where {N, T <: BitStr{N}}
    @assert N == L "Basis BitStr size must match L"
    
    # Get Néel states as BitStr
    z2, z2p = neel_state_bitstr(L)
    
    # Find indices in constrained basis
    idx_z2 = searchsortedfirst(basis, z2)
    idx_z2p = searchsortedfirst(basis, z2p)
    @assert basis[idx_z2] == z2 "Z2 state not in constrained basis"
    @assert basis[idx_z2p] == z2p "Z2' state not in constrained basis"
    
    # Build extended basis R⊗Q
    ext_basis = build_extended_basis(basis)
    d_Q = length(basis)
    
    # |0⟩_R ⊗ |Z2⟩_Q corresponds to index idx_z2 in the first d_Q entries
    # |1⟩_R ⊗ |Z2'⟩_Q corresponds to index d_Q + idx_z2p
    ψ_RQ = zeros(ComplexF64, 2 * d_Q)
    ψ_RQ[idx_z2] = 1/sqrt(2)           # |0⟩_R ⊗ |Z2⟩_Q
    ψ_RQ[d_Q + idx_z2p] = 1/sqrt(2)    # |1⟩_R ⊗ |Z2'⟩_Q
    
    return ψ_RQ, ext_basis
end

"""
    prepare_thermal_encoding_constrained(L::Int, th1::Vector, th2::Vector, 
                                         basis::Vector{T}) where {T <: BitStr}

Prepare reference-code entangled state for thermal encoding in constrained basis.

Constructs the state:
|ψ_RQ^th⟩ = (1/√2)(|0⟩_R ⊗ |th1⟩_Q + |1⟩_R ⊗ |th2⟩_Q)

where |th1⟩ and |th2⟩ are orthogonal thermal eigenstates at E≈0, 
expressed in the constrained basis.

# Arguments
- `L::Int`: System size
- `th1::Vector`: First thermal eigenstate (in constrained basis)
- `th2::Vector`: Second thermal eigenstate (orthogonal to th1)
- `basis::Vector{T}`: PXP constrained basis

# Returns
- `ψ_RQ::Vector{ComplexF64}`: Reference-code entangled state in extended basis
- `ext_basis::Vector{BitStr{L+1}}`: Extended basis (R⊗Q)
"""
function prepare_thermal_encoding_constrained(L::Int, th1::Vector, th2::Vector, 
                                              basis::Vector{T}) where {N, T <: BitStr{N}}
    @assert N == L "Basis BitStr size must match L"
    d_Q = length(basis)
    @assert length(th1) == d_Q "th1 must be in constrained basis"
    @assert length(th2) == d_Q "th2 must be in constrained basis"
    
    ext_basis = build_extended_basis(basis)
    
    ψ_RQ = zeros(ComplexF64, 2 * d_Q)
    ψ_RQ[1:d_Q] = th1 / sqrt(2)
    ψ_RQ[d_Q+1:2*d_Q] = th2 / sqrt(2)
    
    return ψ_RQ, ext_basis
end

#=============================================================================
# Knill-Laflamme Condition Analysis
=============================================================================#

"""
    knill_laflamme_coefficient_Z(L::Int, basis::Vector{T}) where {T <: BitStr}

Compute the b(L) coefficient for Z-dephasing errors in constrained basis.

Following Sang & Zou (arXiv:2406.09555), the coefficient is:
b(L) = (1/D²) Σ_j [D·Tr(Z_j P Z_j P) - Tr(Z_j P)²]

where:
- P = |Z2⟩⟨Z2| + |Z2'⟩⟨Z2'| is the code projector
- D = 2 is the code dimension
- The sum is over all sites j

For Z operators, this simplifies since Z_j|Z2⟩ = ±|Z2⟩ depending on site parity.

# Arguments
- `L::Int`: System size
- `basis::Vector{T}`: PXP constrained basis

# Returns
- `Float64`: The b(L) coefficient

# Example
```julia
basis = PXP_basis(10, true)
b = knill_laflamme_coefficient_Z(10, basis)
println("b(L=10, Z) = \$b")
```
"""
function knill_laflamme_coefficient_Z(L::Int, basis::Vector{T}) where {N, T <: BitStr{N}}
    @assert N == L "Basis BitStr size must match L"
    
    z2, z2p = neel_state_bitstr(L)
    idx_z2 = searchsortedfirst(basis, z2)
    idx_z2p = searchsortedfirst(basis, z2p)
    
    D = 2  # Code dimension
    b_sum = 0.0
    
    for j in 1:L
        # Z_j eigenvalue on |Z2⟩ and |Z2'⟩
        # |Z2⟩ = |101010...⟩ has 1 at odd positions (counting from right, 1-indexed)
        # Z_j|b⟩ = (-1)^{b_j}|b⟩
        z2_val = readbit(Int(z2), L - j + 1) == 1 ? -1.0 : 1.0
        z2p_val = readbit(Int(z2p), L - j + 1) == 1 ? -1.0 : 1.0
        
        # ⟨Z2|Z_j|Z2⟩ = z2_val, ⟨Z2'|Z_j|Z2'⟩ = z2p_val
        # Off-diagonal ⟨Z2|Z_j|Z2'⟩ = 0 since Z is diagonal
        
        # Tr(Z_j P Z_j P) = |⟨Z2|Z_j|Z2⟩|² + |⟨Z2'|Z_j|Z2'⟩|² = 1 + 1 = 2
        tr_ZPZ_P = z2_val^2 + z2p_val^2  # Always 2
        
        # Tr(Z_j P) = ⟨Z2|Z_j|Z2⟩ + ⟨Z2'|Z_j|Z2'⟩
        tr_ZP = z2_val + z2p_val
        
        b_sum += D * tr_ZPZ_P - abs2(tr_ZP)
    end
    
    return real(b_sum) / D^2
end

"""
    reference_code_state(ψ_RQ::Vector)

Create density matrix from pure reference-code state.

# Arguments
- `ψ_RQ::Vector`: Pure state vector

# Returns
- `Matrix{ComplexF64}`: Density matrix ρ_RQ = |ψ_RQ⟩⟨ψ_RQ|
"""
reference_code_state(ψ_RQ::Vector) = ψ_RQ * ψ_RQ'
