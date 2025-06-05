using Test
using PXPConstrained
using LinearAlgebra

@testset "FSA" begin
    N = 12
    # map from total basis to FSA basis, isometry which should satisfy u'*u=I, belike:
    u = iso_total2FSA(N)
    @test size(u) == (322, 13)
    @test u'*u ≈ I(13)
    P = u*u'
    @test isapprox(P*P, P, atol=1e-10)
end

@testset "ishermitian" begin
    for L in 6:14
        @test ishermitian(PXP_Ham_sparse(L))
    end
    for L in 8:16
        @test ishermitian(PXP_K_Ham_sparse(L, 0))
        @test ishermitian(PXP_K_Ham_sparse(L, div(L, 2)))
    end
    for L in 10:18
        @test ishermitian(PXP_MSS_Ham_sparse(L, 0))
        @test ishermitian(PXP_MSS_Ham_sparse(L, div(L, 2)))
    end
end

@testset "sparse hamiltonian" begin
    @test eigvals(Matrix(PXP_Ham_sparse(12))) ≈ eigvals(PXP_Ham(12))
    @test eigvals(Matrix(PXP_K_Ham_sparse(14, 0))) ≈ eigvals(PXP_K_Ham(14, 0))
    @test eigvals(Matrix(PXP_K_Ham_sparse(14, 7))) ≈ eigvals(PXP_K_Ham(14, 7))
    @test eigvals(Matrix(PXP_MSS_Ham_sparse(16, 0))) ≈ eigvals(PXP_MSS_Ham(16, 0))
    @test eigvals(Matrix(PXP_MSS_Ham_sparse(16, 8))) ≈ eigvals(PXP_MSS_Ham(16, 8))
    @test PXP_K_Ham_sparse(8,4) ≈ PXP_MSS_Ham_sparse(8,4,-1)

    vals= eigvals(Matrix(PXP_Ham_sparse(12)))
    @test isapprox(vals.+reverse(vals), zeros(length(vals)), atol=1e-10)
    vals= eigvals(Matrix(PXP_K_Ham_sparse(14, 0)))
    @test isapprox(vals.+reverse(vals), zeros(length(vals)), atol=1e-10)
    vals= eigvals(Matrix(PXP_K_Ham_sparse(14, 7)))
    @test isapprox(vals.+reverse(vals), zeros(length(vals)), atol=1e-10)
    vals= eigvals(Matrix(PXP_MSS_Ham_sparse(16, 0)))
    @test isapprox(vals.+reverse(vals), zeros(length(vals)), atol=1e-10)
    vals= eigvals(Matrix(PXP_MSS_Ham_sparse(16, 8)))
    @test isapprox(vals.+reverse(vals), zeros(length(vals)), atol=1e-10)
end

@testset "sparse iso" begin
    u= iso_total2K_sparse(12, 6)
    P=u*u'
    @test isapprox(P*P, P, atol=1e-10)
    @test u'*u ≈ I(size(u, 2))
    @test Matrix(u) ≈ iso_total2K(12, 6)

    u= iso_total2K_sparse(12, 0)
    P=u*u'
    @test isapprox(P*P, P, atol=1e-10)
    @test u'*u ≈ I(size(u, 2))
    @test Matrix(u) ≈ iso_total2K(12, 0)
    @test u[1,1] ≈ 1.0
    @test u[end,end] ≈ 1/√2
    @test u[233,end] ≈ 1/√2

    u = iso_K2MSS_sparse(12, 6)
    P=u*u'
    @test isapprox(P*P, P, atol=1e-10)
    @test u'*u ≈ I(size(u, 2))
    @test Matrix(u) ≈ iso_K2MSS(12, 6)

    u = iso_K2MSS_sparse(12, 0)
    P=u*u'
    @test isapprox(P*P, P, atol=1e-10)
    @test u'*u ≈ I(size(u, 2))
    @test Matrix(u) ≈ iso_K2MSS(12, 0)
    non_zero_lis = reshape(mapslices(col -> count(x -> x != 0, col), u, dims=1), 26)
    @test non_zero_lis == PXP_MSS_basis(12, 0)[3]

    u = iso_total2MSS_sparse(14, 7)
    P=u*u'
    @test isapprox(P*P, P, atol=1e-10)
    @test u'*u ≈ I(size(u, 2))
    @test Matrix(u) ≈ iso_total2MSS(14, 7)

    u = iso_total2MSS_sparse(14, 0)
    P=u*u'
    @test isapprox(P*P, P, atol=1e-10)
    @test u'*u ≈ I(size(u, 2))
    @test Matrix(u) ≈ iso_total2MSS(14, 0)
    @test u[1,1] ≈ 1.0
    @test u[end,end] ≈ 1/√2
    
    u = iso_total2MSS_sparse(14, 0, -1)
    P=u*u'
    @test isapprox(P*P, P, atol=1e-10)
    @test u'*u ≈ I(size(u, 2))
    @test Matrix(u) ≈ iso_total2MSS(14, 0, -1)
end


