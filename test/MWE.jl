using Test
using PXPConstrainted, BitBasis

using Profile, BenchmarkTools

@testset "Fibonacci chain" begin
    N=12
    state=BitStr{N}(0)
    res = PXPConstrainted.Fibonacci_chain_OBC(BitStr{N})
    @test length(res) == 377
    @test BitStr{N}(0) in res
    @test BitStr{N, Int}(0b100101010101) in res

    basis = PXPConstrainted.PXP_basis(BitStr{N, Int})
    @test length(basis) == 322

    res = PXPConstrainted.actingH_PXP(BitStr{N, Int}, state)
    @test length(res) == 12

    h = PXPConstrainted.PXP_Ham(BitStr{N, Int})
    @test size(h) == (322, 322)

    hk = PXPConstrainted.PXP_K_Ham(BitStr{N, Int}, 0)
    @test size(hk) == (31, 31)

    res = PXPConstrainted.iso_total2FSA(BitStr{N, Int})
    @test size(res) == (322, 13)

    rdm = PXPConstrainted.rdm_PXP(BitStr{N, Int}, collect(1:6), ones(322))
    @test size(rdm) == (21, 21)
end

