using Test
using PXPConstrained, BitBasis    
using LinearAlgebra

@testset "pxp functions" begin
    N=12
    state=BitStr{N}(0)
    res = PXPConstrained.Fibonacci_chain_OBC(BitStr{N})
    # constrained space dims scales as Fibonacci sequence
    @test length(res) == 377
    @test BitStr{N}(0) in res
    @test BitStr{N, Int}(0b100101010101) in res

    basis = PXP_basis(BitStr{N, Int})
    @test length(basis) == 322

    # PBC and OBC basis have different dimensions
    basis_obc = PXP_basis(BitStr{N, Int}, false)
    @test length(basis_obc) == 377

    res = actingH_PXP(BitStr{N, Int}, state)
    @test length(res) == 12

    res_obc = actingH_PXP(BitStr{N, Int}, state,false)
    @test length(res_obc) == 12

    H = PXP_Ham(BitStr{N, Int})
    @test size(H) == (322, 322)

    hk = PXP_K_Ham(BitStr{N, Int}, 0)
    @test size(hk) == (31, 31)

    # map from total basis to FSA basis, isometry which should satisfy u'*u=I, belike:
    res = iso_total2FSA(BitStr{N, Int})
    @test size(res) == (322, 13)
    @test res'*res ≈ I(13)

    # Dims of PBC's rdm should = Dims of OBC
    rdm = rdm_PXP(BitStr{N, Int}, BitStr{div(N,2), Int},collect(1:div(N,2)), ones(322))
    @test size(rdm) == (21, 21)
    @test length(PXP_basis(BitStr{6, Int}, false)) == 21

    map_idx = iso_total2K(BitStr{12, Int},0)
    @test size(map_idx) == (322, 31)
    @test map_idx'*map_idx ≈ I(31)
    u = iso_total2K(14, 7)
    P = u*u'
    @test isapprox(P*P, P, atol=1e-10)
    @test isapprox(u'*u, I(size(u, 2)), atol=1e-8)

    u = iso_total2MSS(14, 7)
    P = u*u'
    @test isapprox(P*P, P, atol=1e-10)
    @test isapprox(u'*u, I(size(u, 2)), atol=1e-10)

    u = iso_total2MSS(14, 0)
    P = u*u'
    @test isapprox(P*P, P, atol=1e-10)
    @test isapprox(u'*u, I(size(u, 2)), atol=1e-10)

    u = iso_total2K(14, 0)
    P = u*u'
    @test isapprox(P*P, P, atol=1e-10)
    @test isapprox(u'*u, I(size(u, 2)), atol=1e-10)

    u = iso_K2MSS(14, 0)
    P = u*u'
    @test isapprox(P*P, P, atol=1e-10)
    @test isapprox(u'*u, I(size(u, 2)), atol=1e-10)

    rdm_K = rdm_PXP_K(BitStr{24, Int}, BitStr{12, Int}, collect(1:12), ones(4341),0)
    @test size(rdm_K) == (377, 377)
    @test length(PXP_basis(BitStr{12, Int}, false)) == 377

    # It means that translate N times is the same as identity, and inversion matrix squared is the identity
    T=translation_matrix(BitStr{N, Int})
    @test isapprox(T^N, I(size(T)[1]), atol=1e-6)
    Inv=inversion_matrix(BitStr{N, Int})

    Inv=inversion_matrix(BitStr{N, Int})
    @test isapprox(Inv^2, I(size(Inv)[1]))

    MSS_basis = PXP_MSS_basis(BitStr{N, Int}, 0)[1]
    @test length(MSS_basis) == 26
    @test length(PXP_MSS_basis(BitStr{8, Int}, 0)[1]) == length(PXP_K_basis(BitStr{8, Int}, 0)[1])

    # The energy spectrum is symmetric about zero
    H_MSS = PXP_MSS_Ham(BitStr{N, Int}, 0)
    MSS_vals, MSS_vecs = eigen(H_MSS)
    @test isapprox(reverse(MSS_vals) .+ MSS_vals, zeros(26), atol=1e-6)

    # The MSS dims + MSSminv dims = K dims
    @test length(PXP_MSS_basis(BitStr{N, Int}, 0)[1]) + length(PXP_MSS_basis(BitStr{N, Int}, 0, -1)[1]) == length(PXP_K_basis(BitStr{N, Int}, 0)[1])
    H_MSSminv = PXP_MSS_Ham(BitStr{N, Int}, 0, -1)
    MSSminv_vals, MSSminv_vecs = eigen(H_MSSminv)
    @test isapprox(MSSminv_vals[3], 0.0, atol=1e-6)

    map_K2MSS = iso_K2MSS(BitStr{N, Int},0)
    map_K2MSSminv = iso_K2MSS(BitStr{N, Int},0, -1)
    @test size(map_K2MSS) == (31, 26)
    @test map_K2MSS'*map_K2MSS ≈ I(26)
    @test size(map_K2MSSminv) == (31, 5)
    @test map_K2MSSminv'*map_K2MSSminv ≈ I(5)

    map_total2MSS = iso_total2MSS(BitStr{N, Int},0)
    @test size(map_total2MSS) == (322, 26)
    @test map_total2MSS'*map_total2MSS ≈ I(26)

    MSS_vec=MSS_vecs[:,12]
    rdm_MSS = rdm_PXP_MSS(BitStr{N, Int}, BitStr{div(N,2), Int}, collect(1:6), MSS_vec,0)
    @test size(rdm_MSS) == (21, 21) == (length(PXP_basis(BitStr{6, Int}, false)) ,length(PXP_basis(BitStr{6, Int}, false)))

    # Fit the scar's central charge, may change depends on the machine and basic linear algebra package.
    splitlis = collect(1:N-1)
    EE_lis=zeros(length(splitlis))
    for m in eachindex(EE_lis)
        subrho=rdm_PXP_MSS(BitStr{N, Int}, BitStr{length(1:splitlis[m]), Int}, collect(1:splitlis[m]), MSS_vec, 0)
        EE_lis[m]=ee(subrho)
    end

    cent, _= fitpage_curve(EE_lis; mincut=1)
    @test isapprox(cent, 0.587, atol=1e-2)
end
