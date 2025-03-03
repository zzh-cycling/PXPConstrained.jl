using Test
using PXPConstrained, BitBasis    
    
@testset "pxp functions" begin
    N=12
    state=BitStr{N}(0)
    res = PXPConstrained.Fibonacci_chain_OBC(BitStr{N})
    @test length(res) == 377
    @test BitStr{N}(0) in res
    @test BitStr{N, Int}(0b100101010101) in res

    basis = PXP_basis(BitStr{N, Int})
    @test length(basis) == 322

    res = actingH_PXP(BitStr{N, Int}, state)
    @test length(res) == 12

    h = PXP_Ham(BitStr{N, Int})
    @test size(h) == (322, 322)

    hk = PXP_K_Ham(BitStr{N, Int}, 0)
    @test size(hk) == (31, 31)

    res = iso_total2FSA(BitStr{N, Int})
    @test size(res) == (322, 13)

    rdm = rdm_PXP(BitStr{N, Int}, collect(1:6), ones(322))
    @test size(rdm) == (21, 21)

    map_idx = iso_total2K(BitStr{12, Int},0)
    @test size(map_idx) == (322, 31)

    rdm_K = rdm_PXP_K(BitStr{24, Int}, collect(1:12), ones(4341),0)
    @test size(rdm_K) == (377, 377)

    T=translation_matrix(BitStr{N, Int})
    @test isapprox(T^N, I(size(T)[1]), atol=1e-6)
    Inv=inversion_matrix(BitStr{N, Int})

    Inv=inversion_matrix(BitStr{N, Int})
    @test isapprox(Inv^2, I(size(Inv)[1]))

    # map_idx_MSS=iso_total2MSS(BitStr{12, Int},0)
end