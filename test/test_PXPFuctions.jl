using Test
using PXPConstrained, BitBasis    
using LinearAlgebra

@testset "pxp functions" begin
    N=12
    state=BitStr{N}(0)
    res = PXPConstrained.Fibonacci_chain_OBC(BitStr{N})
    # constrained space dims scales as Fibonacci sequence
    @test length(res) == 377
    @test BitStr{N}(0)==res[1]
    @test BitStr{N, Int}(0b100101010101) in res

    basis = PXP_basis(N)
    @test length(basis) == 322
    @test BitStr{N}(0)==basis[1]
    @test basis[end].buf==2730
    @test basis[233].buf==2730>>1

    # PBC and OBC basis have different dimensions
    basis_obc = PXP_basis(N, false)
    @test length(basis_obc) == 377
    @test basis_obc[end].buf==2730

    res = actingH_PXP(BitStr{N, Int}, state)
    @test length(res) == 12

    res_obc = actingH_PXP(BitStr{N, Int}, state,false)
    @test length(res_obc) == 12

    H = PXP_Ham(N)
    @test size(H) == (322, 322)
    @test ishermitian(H)
    @test H[1, 1]==H[end,end] ≈ 0.0

    hk = PXP_K_Ham(N, 0)
    @test size(hk) == (31, 31)
    @test ishermitian(hk)
    @test hk[1, 1] == hk[end, end]≈ 0.0

    hmss=PXP_MSS_Ham(N, 0)
    @test size(hmss) == (26, 26)
    @test ishermitian(hmss)
    @test hmss[1, 1] == hmss[end, end] ≈ 0.0

    # map from total basis to FSA basis, isometry which should satisfy u'*u=I, belike:
    u = iso_total2FSA(N)
    @test size(u) == (322, 13)
    @test u'*u ≈ I(13)
    P = u*u'
    @test isapprox(P*P, P, atol=1e-10)

    # Dims of PBC's rdm should = Dims of OBC
    rdm = rdm_PXP(N, [collect(1:div(N,2))], ones(322))
    @test ishermitian(rdm)
    @test isapprox(tr(rdm), 322.0, atol=1e-6)
    @test size(rdm) == (21, 21)
    @test length(PXP_basis(6, false)) == 21

    map_idx = iso_total2K(BitStr{12, Int},0)
    @test size(map_idx) == (322, 31)
    @test map_idx'*map_idx ≈ I(31)
    u = iso_total2K(14, 7)
    P = u*u'
    @test isapprox(P*P, P, atol=1e-10)
    @test isapprox(u'*u, I(size(u, 2)), atol=1e-8)
    @test map_idx[1, 1] == 1.0
    @test map_idx[end, end] ≈ 1/√2
    @test map_idx[233, end] ≈ 1/√2

    iso = iso_K2MSS(12, 0)
    @test size(iso) == (31, 26)
    @test iso'*iso ≈ I(26)
    P=iso*iso'
    @test isapprox(P*P, P, atol=1e-10)
    non_zero_lis = reshape(mapslices(col -> count(x -> x != 0, col), iso, dims=1), 26)
    @test non_zero_lis == PXP_MSS_basis(12, 0)[3]
    u = iso_K2MSS(14, 0)
    P = u*u'
    @test isapprox(P*P, P, atol=1e-10)
    @test isapprox(u'*u, I(size(u, 2)), atol=1e-10)

    u = iso_total2MSS(12, 0)
    P = u*u'
    @test isapprox(P*P, P, atol=1e-10)
    @test isapprox(u'*u, I(size(u, 2)), atol=1e-10)
    @test u[1,1] == 1.0
    @test u[end, end] ≈ 1/√2
    @test u[233, end] ≈ 1/√2

    u = iso_total2MSS(14, 7)
    P = u*u'
    @test isapprox(P*P, P, atol=1e-10)
    @test isapprox(u'*u, I(size(u, 2)), atol=1e-10)

    u = iso_total2K(14, 0)
    P = u*u'
    @test isapprox(P*P, P, atol=1e-10)
    @test isapprox(u'*u, I(size(u, 2)), atol=1e-10)
    @test u[1,1] == 1.0
    @test u[end, end] ≈ 1/√2

    rdm_K = rdm_PXP_K(24, [collect(1:12)], ones(4341),0)
    @test size(rdm_K) == (377, 377)
    @test length(PXP_basis(12, false)) == 377

    # It means that translate N times is the same as identity, and inversion matrix squared is the identity
    T=translation_matrix(N)
    @test isapprox(T^N, I(size(T)[1]), atol=1e-6)
    Inv=inversion_matrix(N)
    @test isapprox(Inv^2, I(size(Inv)[1]))

    MSS_basis = PXP_MSS_basis(N, 0)[1]
    @test length(MSS_basis) == 26
    @test length(PXP_MSS_basis(8, 0)[1]) == length(PXP_K_basis(8, 0)[1])

    # The energy spectrum is symmetric about zero
    H_MSS = PXP_MSS_Ham(N, 0)
    MSS_vals, MSS_vecs = eigen(H_MSS)
    @test isapprox(reverse(MSS_vals) .+ MSS_vals, zeros(26), atol=1e-6)

    # The MSS dims + MSSminv dims = K dims
    @test length(PXP_MSS_basis(N, 0)[1]) + length(PXP_MSS_basis(N, 0, -1)[1]) == length(PXP_K_basis(N, 0)[1])
    H_MSSminv = PXP_MSS_Ham(N, 0, -1)
    MSSminv_vals, MSSminv_vecs = eigen(H_MSSminv)
    @test isapprox(MSSminv_vals[3], 0.0, atol=1e-6)

    map_K2MSS = iso_K2MSS(N, 0)
    map_K2MSSminv = iso_K2MSS(N, 0, -1)
    @test size(map_K2MSS) == (31, 26)
    @test map_K2MSS'*map_K2MSS ≈ I(26)
    @test size(map_K2MSSminv) == (31, 5)
    @test map_K2MSSminv'*map_K2MSSminv ≈ I(5)

    map_total2MSS = iso_total2MSS(N,0)
    @test size(map_total2MSS) == (322, 26)
    @test map_total2MSS'*map_total2MSS ≈ I(26)

    MSS_vec=MSS_vecs[:,12]
    rdm_MSS = rdm_PXP_MSS(N, [collect(1:6)], MSS_vec,0)
    @test size(rdm_MSS) == (21, 21) == (length(PXP_basis(6, false)) ,length(PXP_basis(6, false)))

    # Fit the scar's central charge, may change depends on the machine and basic linear algebra package.
    splitlis = collect(1:N-1)
    EE_lis=zeros(length(splitlis))
    for m in eachindex(EE_lis)
        subrho=rdm_PXP_MSS(N, [collect(1:splitlis[m])], MSS_vec, 0)
        EE_lis[m]=ee(subrho)
    end

    cent, _= fitCCEntEntScal(EE_lis; mincut=1, pbc=true)
    @test isapprox(cent, 1.924, atol=1e-2)
end
