# Tests for QEC.jl - Quantum Error Correction module
using Test
using PXPConstrained
using LinearAlgebra
using BitBasis

@testset "build_extended_basis" begin
    L = 6
    basis = PXP_basis(L, true)
    ext_basis = build_extended_basis(basis)
    
    # Extended basis should have 2× the size (R = {0,1})
    @test length(ext_basis) == 2 * length(basis)
    @test issorted(ext_basis)
    
    # First half |0⟩_R, second half |1⟩_R
    d_Q = length(basis)
    for (i, b) in enumerate(ext_basis)
        r = readbit(Int(b), L + 1)
        @test (i <= d_Q) ? (r == 0) : (r == 1)
    end
end

@testset "neel_state_bitstr" begin
    # L=6: |Z2⟩ = |101010⟩ = 42, |Z2'⟩ = |010101⟩ = 21
    z2, z2p = neel_state_bitstr(6)
    @test Int(z2) == 0b101010
    @test Int(z2p) == 0b010101
    
    # Both should be valid PXP states
    basis = PXP_basis(6, true)
    @test z2 in basis
    @test z2p in basis
    
    # L=8: |Z2⟩ = |10101010⟩ = 170
    z2_8, z2p_8 = neel_state_bitstr(8)
    @test Int(z2_8) == 0b10101010
    @test Int(z2p_8) == 0b01010101
end

@testset "prepare_scar_encoding_constrained" begin
    L = 6
    basis = PXP_basis(L, true)
    ψ_RQ, ext_basis = prepare_scar_encoding_constrained(L, basis)
    d_Q = length(basis)
    
    @test norm(ψ_RQ) ≈ 1.0
    @test length(ψ_RQ) == 2 * d_Q
    
    # Exactly 2 non-zero entries at 1/√2
    non_zero = findall(x -> abs(x) > 1e-10, ψ_RQ)
    @test length(non_zero) == 2
    @test all(abs.(ψ_RQ[non_zero]) .≈ 1/sqrt(2))
    
    # Verify structure: |0⟩⊗|Z2⟩ and |1⟩⊗|Z2'⟩
    z2, z2p = neel_state_bitstr(L)
    idx_z2 = searchsortedfirst(basis, z2)
    idx_z2p = searchsortedfirst(basis, z2p)
    @test abs(ψ_RQ[idx_z2]) ≈ 1/sqrt(2)
    @test abs(ψ_RQ[d_Q + idx_z2p]) ≈ 1/sqrt(2)
end

@testset "partial_trace_R" begin
    d_Q = 4
    
    # Product state |0⟩_R ⊗ |0⟩_Q
    ψ = zeros(ComplexF64, 2 * d_Q); ψ[1] = 1.0
    ρ_Q = partial_trace_R(ψ * ψ', d_Q)
    @test ρ_Q[1,1] ≈ 1.0
    @test sum(abs.(ρ_Q)) ≈ 1.0
    
    # Maximally entangled: (|00⟩ + |11⟩)/√2
    ψ_ent = zeros(ComplexF64, 2 * d_Q)
    ψ_ent[1] = 1/sqrt(2); ψ_ent[d_Q + 2] = 1/sqrt(2)
    ρ_Q_ent = partial_trace_R(ψ_ent * ψ_ent', d_Q)
    @test ρ_Q_ent[1,1] ≈ 0.5
    @test ρ_Q_ent[2,2] ≈ 0.5
    @test tr(ρ_Q_ent) ≈ 1.0
end

@testset "von_neumann_entropy" begin
    # Pure state: S = 0
    ψ = [1.0, 0.0, 0.0, 0.0]
    @test von_neumann_entropy(ψ * ψ') ≈ 0.0 atol=1e-10
    
    # Maximally mixed 2×2: S = log(2)
    @test von_neumann_entropy([0.5 0; 0 0.5]) ≈ log(2) atol=1e-10
    
    # Maximally mixed 4×4: S = log(4)
    @test von_neumann_entropy(I(4) / 4) ≈ log(4) atol=1e-10
end

@testset "coherent_information_constrained" begin
    L = 6
    basis = PXP_basis(L, true)
    d_Q = length(basis)
    
    # Pure maximally entangled: I_c = log(2)
    ψ_RQ, _ = prepare_scar_encoding_constrained(L, basis)
    @test coherent_information_constrained(ψ_RQ * ψ_RQ', d_Q) ≈ log(2) atol=1e-10
    
    # Product state: I_c = 0
    ψ_prod = zeros(ComplexF64, 2 * d_Q); ψ_prod[1] = 1.0
    @test coherent_information_constrained(ψ_prod * ψ_prod', d_Q) ≈ 0.0 atol=1e-10
end

@testset "apply_Z_dephasing" begin
    L = 6
    basis = PXP_basis(L, true)
    d_Q = length(basis)
    
    z2, _ = neel_state_bitstr(L)
    idx_z2 = searchsortedfirst(basis, z2)
    ψ = zeros(ComplexF64, d_Q); ψ[idx_z2] = 1.0
    ρ = ψ * ψ'
    
    # p=0: unchanged
    @test apply_Z_dephasing_constrained(ρ, 0.0, L, basis) ≈ ρ
    
    # Trace and Hermiticity preserved
    ρ_p05 = apply_Z_dephasing_constrained(ρ, 0.5, L, basis)
    @test tr(ρ_p05) ≈ 1.0 atol=1e-10
    @test ρ_p05 ≈ ρ_p05'
end

@testset "apply_Z_dephasing_extended" begin
    L = 6
    basis = PXP_basis(L, true)
    d_Q = length(basis)
    
    ψ_RQ, ext_basis = prepare_scar_encoding_constrained(L, basis)
    ρ_RQ = ψ_RQ * ψ_RQ'
    
    # p=0: unchanged
    @test apply_Z_dephasing_extended(ρ_RQ, 0.0, L, ext_basis) ≈ ρ_RQ
    
    # Trace preserved, I_c decreases
    ρ_p05 = apply_Z_dephasing_extended(ρ_RQ, 0.5, L, ext_basis)
    @test tr(ρ_p05) ≈ 1.0 atol=1e-10
    @test coherent_information_constrained(ρ_p05, d_Q) < coherent_information_constrained(ρ_RQ, d_Q)
end

@testset "knill_laflamme_coefficient_Z" begin
    # For Z-dephasing on Néel states: b(L) = L
    # Because Tr(Z_j P) = 0 (opposite parities) and Tr(Z_j P Z_j P) = 2
    # So b = (1/4) * Σ_j [2*2 - 0] = L
    for L in [6, 8, 10]
        basis = PXP_basis(L, true)
        @test knill_laflamme_coefficient_Z(L, basis) ≈ L atol=1e-10
    end
end

@testset "coherent_information_vs_noise" begin
    L = 6
    basis = PXP_basis(L, true)
    d_Q = length(basis)
    
    ψ_RQ, ext_basis = prepare_scar_encoding_constrained(L, basis)
    ρ_RQ_0 = ψ_RQ * ψ_RQ'
    
    p_values = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    I_c_prev = Inf
    
    for p in p_values
        ρ_noisy = apply_Z_dephasing_extended(ρ_RQ_0, p, L, ext_basis)
        I_c = coherent_information_constrained(ρ_noisy, d_Q)
        @test I_c <= I_c_prev + 1e-10  # Monotonically decreasing
        I_c_prev = I_c
    end
end
