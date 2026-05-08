using PXPConstrained
using BitBasis
using LinearAlgebra

L = 11
H = PXP_Ham(L)
energy, states = eigen(H)

psi0 = zeros(length(PXP_basis(L)))  # all zero state
psi0[end] = 1  # Example initial state
times = collect(0:0.01:5.0)
stlis = wf_time_evolution(psi0, times, energy, states);
iso = iso_full2cons(L, true)
full_stlis = [iso*st for st in stlis]

# F = ⟨W^{†}V^†(t) W V(t)⟩ = ⟨ψ| W^{†}U V^†U^† W U^† V U|ψ⟩
# U = exp(-iHt) = S diag(exp(-iE*t)) S†, H = S diag(E) S†
# V(t) = S diag(Ṽ) S†, S 
# 1. |a⟩ = V(t)|ψ⟩ = S diag(exp(-iE*t)) S† |ψ⟩
# 2. |b⟩ = W|a⟩ = W S diag(exp(-iE*t)) S† |ψ⟩
# 3. |c⟩ = V†(t)|b⟩ = S diag(exp(iE*t)) S† W S diag(exp(-iE*t)) S† |ψ⟩
# 4. F = ⟨ψ| W† |c⟩ = ⟨W ψ|c⟩* = dot(W ψ, c)

"""
    OTOC(W::Matrix{ET}, V::Matrix{ET}, psi::Vector{ET}, times::Vector{Float64}, energy::Vector{Float64}, states::Matrix{Float64}) where {ET}

Calculate the out-of-time-ordered correlator (OTOC) using exact diagonalization.

Computes F(t) = ⟨ψ| W† U†(t) V† U(t) W U†(t) V U(t) |ψ⟩
where U(t) = exp(-iHt) is the time evolution operator.
Using the eigendecomposition H = S·diag(energy)·S†, we have U(t) = S·diag(exp(-i·energy·t))·S†.

# Arguments
- `W::Matrix{ET}`: First local operator in the PXP basis
- `V::Matrix{ET}`: Second local operator in the PXP basis
- `psi::Vector{ET}`: Initial quantum state
- `times::Vector{Float64}`: Time points to evaluate
- `energy::Vector{Float64}`: Eigenvalues of the Hamiltonian
- `states::Matrix{Float64}`: Eigenvectors of the Hamiltonian (columns)

# Returns
- `Vector{ComplexF64}`: OTOC values at each time point

# Example
```julia
N = 10
H = PXP_Ham(N)
energy, states = eigen(H)
psi = states[:, 1]  # ground state

# Build local Pauli-Z at site 1 in the PXP basis
basis = PXP_basis(N)
W = zeros(ComplexF64, length(basis), length(basis))
V = zeros(ComplexF64, length(basis), length(basis))
for (idx, str) in enumerate(basis)
    W[idx, idx] = str[N] == 0 ? 1.0 : -1.0  # Z at site 1 (rightmost bit)
    V[idx, idx] = str[N-1] == 0 ? 1.0 : -1.0  # Z at site 2
end

times = collect(0:0.1:10)
F = OTOC(W, V, psi, times, energy, states)
```
"""
function OTOC(W::Matrix{T1}, V::Matrix{T2}, psi::Vector{T3}, times::Vector{Float64}, energy::Vector{Float64}, states::Matrix{Float64}) where {T1, T2, T3}
    # H = states * diagm(energy) * states'
    # U(t) = states * diagm(exp.(-1im * energy * t)) * states'
    # U†(t) = states * diagm(exp.(1im * energy * t)) * states'
    
    # Precompute W† and V†
    Wdagger = W'
    Vdagger = V'
    
    # Transform operators to energy eigenbasis: Õ = S† · O · S
    # This makes time evolution trivial: U(t) = diag(exp(-i·energy·t))
    W_eig = states' * W * states
    V_eig = states' * V * states
    Wdagger_eig = states' * Wdagger * states
    Vdagger_eig = states' * Vdagger * states
    
    # Transform initial state to energy eigenbasis
    psi_eig = states' * psi
    
    otoc_values = zeros(ComplexF64, length(times))
    
    for (i, t) in enumerate(times)
        # U(t) in eigenbasis is diagonal
        exp_factors = exp.(-1im * t * energy)
        exp_factors_dag = exp.(1im * t * energy)
        
        # U(t) |ψ⟩ in eigenbasis
        psi_t = psi_eig .* exp_factors
        
        # V · U(t) |ψ⟩ in eigenbasis
        psi_v = V_eig * psi_t
        
        # U†(t) · V · U(t) |ψ⟩ in eigenbasis
        psi_v_ut = psi_v .* exp_factors_dag
        
        # W · U†(t) · V · U(t) |ψ⟩ in eigenbasis
        psi_w_v_ut = W_eig * psi_v_ut
        
        # U(t) · W · U†(t) · V · U(t) |ψ⟩ in eigenbasis
        psi_u_w_v_ut = psi_w_v_ut .* exp_factors
        
        # V† · U(t) · W · U†(t) · V · U(t) |ψ⟩ in eigenbasis
        psi_v_u_w_v_ut = Vdagger_eig * psi_u_w_v_ut
        
        # U†(t) · V† · U(t) · W · U†(t) · V · U(t) |ψ⟩ in eigenbasis
        psi_ut_v_u_w_v_ut = psi_v_u_w_v_ut .* exp_factors_dag
        
        # W† · U†(t) · V† · U(t) · W · U†(t) · V · U(t) |ψ⟩ in eigenbasis
        psi_w_ut_v_u_w_v_ut = Wdagger_eig * psi_ut_v_u_w_v_ut
        
        # Finally compute ⟨ψ| · (full operator chain) |ψ⟩
        otoc_values[i] = dot(psi_eig, psi_w_ut_v_u_w_v_ut)
    end
    
    return otoc_values
end

"""
    apply_operator_map_full(full_basis::Vector{T}, op_map::Function, state::Vector{ET}) where {N, T <: BitStr{N}, ET}

Apply an operator represented as a mapping function to a state vector in the full Hilbert space.

The operator map takes a basis state and returns a list of (output_state, amplitude) pairs.
Unlike the constrained version, this does NOT check if the output is in the basis — 
for full basis, all outputs are guaranteed to be in the basis.

# Arguments
- `full_basis::Vector{T}`: Full basis states (0 to 2^N-1, sorted)
- `op_map::Function`: Function that maps a basis state to Vector{Tuple{T, ET}} of (output_state, amplitude)
- `state::Vector{ET}`: Input state vector

# Returns
- `Vector{ET}`: Output state after applying the operator
"""
function apply_operator_map_full(full_basis::Vector{T}, op_map::Function, state::Vector{ET}) where {N, T <: BitStr{N}, ET}
    result = zeros(ET, length(full_basis))
    for (idx, str) in enumerate(full_basis)
        outputs = op_map(str)
        for (out_str, amp) in outputs
            # For full basis, out_str.buf + 1 is the index directly
            j = out_str.buf + 1
            result[j] += amp * state[idx]
        end
    end
    return result
end

"""
    Z_map_full(::Type{T}, i::Int64) where {N, T <: BitStr{N}}
    Z_map_full(N::Int64, i::Int64)

Create a Pauli-Z operator at site i as a mapping function for the full Hilbert space.

Z_i |...n_i...⟩ = (1 - 2*n_i) |...n_i...⟩ = ±1 * |...n_i...⟩
where n_i ∈ {0,1} is the occupation at site i.

Site index i counts from the left (1-based), consistent with the physical convention.

# Arguments
- `T::Type{BitStr{N}}` or `N::Int64`: System size specification
- `i::Int64`: Site index (1-based, counting from left)

# Returns
- `Function`: A function that maps a basis state `str` to Vector{Tuple{T, Float64}} pairs
"""
function Z_map_full(::Type{T}, i::Int64) where {N, T <: BitStr{N}}
    bit_idx = N + 1 - i
    
    function op_map(str::T)
        amp = str[bit_idx] == 0 ? 1.0 : -1.0
        return [(str, amp)]
    end
    
    return op_map
end
Z_map_full(N::Int64, i::Int64) = Z_map_full(BitStr{N, Int}, i)

"""
    X_map_full(::Type{T}, i::Int64) where {N, T <: BitStr{N}}
    X_map_full(N::Int64, i::Int64)

Create a Pauli-X operator at site i as a mapping function for the full Hilbert space.

X_i |...n_i...⟩ = flip the bit at site i.

Site index i counts from the left (1-based), consistent with the physical convention.
The implementation follows the style: `fl = bmask(T, N); flip(state, fl >> (i-1))`.

# Arguments
- `T::Type{BitStr{N}}` or `N::Int64`: System size specification
- `i::Int64`: Site index (1-based, counting from left)

# Returns
- `Function`: A function that maps a basis state `str` to Vector{Tuple{T, Float64}} pairs
"""
function X_map_full(::Type{T}, i::Int64) where {N, T <: BitStr{N}}
    fl = bmask(T, N)
    mask = fl >> (i - 1)
    
    function op_map(str::T)
        flipped = flip(str, mask)
        return [(flipped, 1.0)]
    end
    
    return op_map
end
X_map_full(N::Int64, i::Int64) = X_map_full(BitStr{N, Int}, i)

"""
    OTOC_full(W_map::Function, V_map::Function, psi_cons::Vector{ET}, times::Vector{Float64}, 
              energy::Vector{Float64}, states::Matrix{Float64}, iso::Matrix{T2}) where {ET, T2}

Calculate the out-of-time-ordered correlator (OTOC) in the full Hilbert space.

The state evolves under the PXP Hamiltonian (constrained space), but W and V act 
in the full Hilbert space. The isometry `iso` maps constrained → full space.

Computes F(t) = ⟨ψ_full| W† U†(t) V† U(t) W U†(t) V U(t) |ψ_full⟩
where |ψ_full⟩ = iso * |ψ_cons⟩ and U(t) = exp(-i H_PXP t).

# Arguments
- `W_map::Function`: Mapping function for operator W in full space, see `X_map_full` or `Z_map_full`
- `V_map::Function`: Mapping function for operator V in full space
- `psi_cons::Vector{ET}`: Initial quantum state in constrained PXP basis
- `times::Vector{Float64}`: Time points to evaluate
- `energy::Vector{Float64}`: Eigenvalues of the PXP Hamiltonian
- `states::Matrix{Float64}`: Eigenvectors of the PXP Hamiltonian (columns)
- `iso::Matrix`: Isometry mapping constrained space to full space (from iso_full2cons)

# Returns
- `Vector{ComplexF64}`: OTOC values at each time point

# Example
```julia
N = 10
H = PXP_Ham(N)
energy, states = eigen(H)
psi_cons = states[:, 1]  # ground state in constrained space
iso = iso_full2cons(N, true)

# W = X at site 1, V = Z at site 2 (in full space)
W = X_map_full(N, 1)
V = Z_map_full(N, 2)

times = collect(0:0.1:10)
F = OTOC_full(W, V, psi_cons, times, energy, states, iso)
```
"""
function OTOC_full(W_map::Function, V_map::Function, psi_cons::Vector{ET}, times::Vector{Float64}, 
                   energy::Vector{Float64}, states::Matrix{Float64}, iso::Matrix{T2}) where {ET, T2}
    
    # Transform initial state to energy eigenbasis
    psi_eig = states' * psi_cons
    
    # Isometry: constrained → full space
    # iso' : full → constrained
    dim_cons = length(energy)
    dim_full = size(iso, 1)
    full_basis = BitStr{L, Int}.(0:(dim_full-1))
    
    # Precompute operator matrix elements in the energy eigenbasis
    # W_{mn} = ⟨m| iso' · W · iso |n⟩ where |m⟩, |n⟩ are energy eigenstates
    # iso' · W · iso acts in constrained space, but W acts in full space
    
    W_eig = zeros(ComplexF64, dim_cons, dim_cons)
    V_eig = zeros(ComplexF64, dim_cons, dim_cons)
    
    for n in 1:dim_cons
        # |n⟩ in constrained basis
        state_n = states[:, n]
        
        # Map to full space, apply W, map back to constrained space
        state_n_full = iso * state_n
        W_state_full = apply_operator_map_full(full_basis, W_map, state_n_full)
        W_state_cons = iso' * W_state_full
        
        V_state_full = apply_operator_map_full(full_basis, V_map, state_n_full)
        V_state_cons = iso' * V_state_full
        
        # ⟨m|W|n⟩ = states[:, m]' * W_state_cons
        for m in 1:dim_cons
            W_eig[m, n] = dot(states[:, m], W_state_cons)
            V_eig[m, n] = dot(states[:, m], V_state_cons)
        end
    end
    
    # For Hermitian operators like X, Z: W† = W, V† = V
    Wdagger_eig = W_eig'
    Vdagger_eig = V_eig'
    
    otoc_values = zeros(ComplexF64, length(times))
    
    for (i, t) in enumerate(times)
        # U(t) in eigenbasis is diagonal
        exp_factors = exp.(-1im * t * energy)
        exp_factors_dag = exp.(1im * t * energy)
        
        # U(t) |ψ⟩ in eigenbasis
        psi_t = psi_eig .* exp_factors
        
        # V · U(t) |ψ⟩ in eigenbasis
        psi_v = V_eig * psi_t
        
        # U†(t) · V · U(t) |ψ⟩ in eigenbasis
        psi_v_ut = psi_v .* exp_factors_dag
        
        # W · U†(t) · V · U(t) |ψ⟩ in eigenbasis
        psi_w_v_ut = W_eig * psi_v_ut
        
        # U(t) · W · U†(t) · V · U(t) |ψ⟩ in eigenbasis
        psi_u_w_v_ut = psi_w_v_ut .* exp_factors
        
        # V† · U(t) · W · U†(t) · V · U(t) |ψ⟩ in eigenbasis
        psi_v_u_w_v_ut = Vdagger_eig * psi_u_w_v_ut
        
        # U†(t) · V† · U(t) · W · U†(t) · V · U(t) |ψ⟩ in eigenbasis
        psi_ut_v_u_w_v_ut = psi_v_u_w_v_ut .* exp_factors_dag
        
        # W† · U†(t) · V† · U(t) · W · U†(t) · V · U(t) |ψ⟩ in eigenbasis
        psi_w_ut_v_u_w_v_ut = Wdagger_eig * psi_ut_v_u_w_v_ut
        
        # Finally compute ⟨ψ| · (full operator chain) |ψ⟩
        otoc_values[i] = dot(psi_eig, psi_w_ut_v_u_w_v_ut)
    end
    
    return real.(otoc_values)
end

# W = X at site 1, V = Z at site 6 (in full space)

Wlis = [X_map_full(L, i) for i in vcat(1:5, 7:11)]
V = Z_map_full(L, 6)

otoc_matrix = zeros(L-1, length(times))
for i in 1:L-1
    W = Wlis[i]
    otoc_values = OTOC_full(W, V, psi0, times, energy, states, iso)
    otoc_matrix[i, :] = otoc_values
end
println("OTOC at t=0: ", otoc_matrix[1, 1])
println("OTOC at t=1: ", otoc_matrix[1, 101])
println("OTOC at t=5: ", otoc_matrix[1, end])

heatmap(1:L-1, times, otoc_matrix'; xlabel="Time", ylabel="Site i", title="OTOC F(t) for W=X_i and V=Z_6", colorbar_title="F(t)")