using Test
using PXPConstrained
using LinearAlgebra

@testset "ishermitian for PXP_Ham_sparse" begin
    for L in 6:14
        @test ishermitian(PXP_Ham_sparse(L))
    end
end

@testset "PXP_Ham_sparse" begin
    @test eigvals(Matrix(PXP_Ham_sparse(12))) ≈ eigvals(PXP_Ham(12))
end

@testset "ishermitian for PXP_K_Ham_sparse" begin
    for L in 6:18
        @test ishermitian(PXP_K_Ham_sparse(L, 0))
        @test ishermitian(PXP_K_Ham_sparse(L, div(L, 2)))
    end
end

@testset "PXP_MSS_Ham_sparse" begin
    @test eigvals(Matrix(PXP_K_Ham_sparse(14, 0))) ≈ eigvals(PXP_K_Ham(14, 0))
    @test eigvals(Matrix(PXP_K_Ham_sparse(14, 7))) ≈ eigvals(PXP_K_Ham(14, 7))
end

@testset "ishermitian for PXP_MSS_Ham_sparse" begin
    for L in 6:20
        @test ishermitian(PXP_MSS_Ham_sparse(L, 0))
        @test ishermitian(PXP_MSS_Ham_sparse(L, div(L, 2)))
    end
end

@testset "PXP_MSS_Ham_sparse" begin
    @test eigvals(Matrix(PXP_MSS_Ham_sparse(16, 0))) ≈ eigvals(PXP_MSS_Ham(16, 0))
    @test eigvals(Matrix(PXP_MSS_Ham_sparse(16, 8))) ≈ eigvals(PXP_MSS_Ham(16, 8))
    @test PXP_K_Ham_sparse(8,4) ≈ PXP_MSS_Ham_sparse(8,4)
end


