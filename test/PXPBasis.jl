using Test
using PXPConstrained, BitBasis    
using LinearAlgebra

@testset "get_representative" begin
    # In case overflow, bit"100000000001110011000000000011100100" > typemax(Int64)
    PXPConstrained.get_representative(bit"100000000001110011000000000011100100") == (bit"000000000011100100100000000001110011", 18)
end

@testset "pxp ham" begin
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

    
end


@testset "process_join" begin
    # [[000 ₍₂₎, 001 ₍₂₎, 010 ₍₂₎, 100 ₍₂₎, 101 ₍₂₎], [000 ₍₂₎, 001 ₍₂₎, 010 ₍₂₎, 100 ₍₂₎, 101 ₍₂₎]]
    lis1 = BitStr{2}[0, 1, 2]
    lis2 = BitStr{3}[0, 1, 2, 4, 5]
    res = PXPConstrained.process_join(lis1, lis2) 
    @test res == vec([join(l2, l1) for l1 in lis1, l2 in lis2])

    # joint_pxp_basis
    res = PXPConstrained.joint_pxp_basis([2, 3])
    @test res == vec([join(l2, l1) for l1 in lis1, l2 in lis2])

    # move_subsystem
    res = PXPConstrained.move_subsystem(BitStr{5, Int}, BitStr{3, Int}(0b101), [1, 2, 5])
    @test res == BitStr{5}(0b10001)

    # takeenviron
    bs, mask = BitStr{5}(0b11001), BitStr{5}(0b10001)
    env = PXPConstrained.takeenviron(bs, mask)
    sys = PXPConstrained.takesystem(bs, mask)
    @test env == BitStr{5}(0b01000)
    @test sys == BitStr{5}(0b10001)
end

@testset "connected components" begin
    v = [1, 2, 4, 5, 7]
    @test PXPConstrained.connected_components(v) == [[1, 2], [4, 5], [7]]
    @test PXPConstrained.connected_components([1,2,3,7,8,9]) == [[1, 2, 3], [7, 8, 9]]
end

@testset "pxp rdm" begin
    N = 12
    # Dims of PBC's rdm should = Dims of OBC
    rdm = rdm_PXP(N, collect(1:div(N,2)), ones(322))
    @test ishermitian(rdm)
    @test isapprox(tr(rdm), 322.0, atol=1e-6)
    @test size(rdm) == (21, 21)
    @test length(PXP_basis(6, false)) == 21
end