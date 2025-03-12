using Test
using PXPConstrained
using BitBasis
using LinearAlgebra

#writing test to get the correct sparse Hamiltonian
@testset "ishermitian for PXP_MSS_Ham_sparse" begin
    for L in 6:20
        @test ishermitian(PXP_MSS_Ham_sparse(L, 0))
        @test ishermitian(PXP_MSS_Ham_sparse(L, div(L, 2)))
    end
end

@testset "PXP_MSS_Ham_sparse" begin
    @test eigvals(Matrix(PXP_MSS_Ham_sparse(16, 0))) ≈ eigvals(PXP_MSS_Ham(16, 0))
    @test eigvals(Matrix(PXP_MSS_Ham_sparse(16, 8))) ≈ eigvals(PXP_MSS_Ham(16, 8))
end

@testset "time_evolution_mss" begin
    state=zeros(455)
    state[end]=1 # Noted that here the "Z2 state" is the superposition of Z2 and Z2tilde. 
    tlis=collect(0:0.1:20)
    H_MSS=PXP_MSS_Ham_sparse(20, 0)
    wflis=wf_time_evolution_mss(20, 0, state, tlis)
    Elis=similar(tlis)
    for i in 1:length(tlis)
        Elis[i]=real(wflis[i]'*H_MSS*wflis[i])
    end
    @test norm.(wflis) ≈ ones(length(wflis)) #normalized
    @test isapprox(Elis, zeros(length(Elis)), atol=1e-10)
    # energy conservation
end

@testset "ergotropy in MSS" begin
    state=zeros(455)
    state[end]=1
    @test ergotropy_PXP_MSS_state(20, 10, state) ≈ 0
end