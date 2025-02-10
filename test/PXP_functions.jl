using Test
using PXPConstrained, BitBasis

@testset "Fibonacci chain" begin
    N=12
    state=BitStr{N}(0)
    res = PXPConstrained.Fibonacci_chain_OBC(BitStr{N})
    @test length(res) == 377
    @test BitStr{N}(0) in res
    @test BitStr{N, Int}(0b100101010101) in res

    basis = PXPConstrained.PXP_basis(BitStr{N, Int})
    @test length(basis) == 322

    res = PXPConstrained.actingH_PXP(BitStr{N, Int}, state)
    @test length(res) == 12

    h = PXPConstrained.PXP_Ham(BitStr{N, Int})
    @test size(h) == (322, 322)

    hk = PXPConstrained.PXP_K_Ham(BitStr{N, Int}, 0)
    @test size(hk) == (31, 31)

    res = PXPConstrained.iso_total2FSA(BitStr{N, Int})
    @test size(res) == (322, 13)

    rdm = PXPConstrained.rdm_PXP(BitStr{N, Int}, collect(1:6), ones(322))
    @test size(rdm) == (21, 21)
end

