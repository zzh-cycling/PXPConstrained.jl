using Test
using PXPConstrained, BitBasis
using ITensors
using LinearAlgebra

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

@testset "scar separate" begin
    N=12
    T = BitStr{N, Int}
    proj=proj_FSA(T)
    @test isapprox(proj, I(size(proj)[1]), atol=1e-6)

    proj=proj_FSA2total(T)
    @test size(proj) == (322, 322)

    H = PXP_Ham(T)
    energy, states= eigen(H)
    scar, total_states = sep_scar_FSA(T, energy, states)

    Ob=randn((322, 322))
    proj=proj_Ob(energy, states, Ob)
    @test size(proj)==(22,22)

    scar, thermal = sep_scar_FSA(T, energy, states)
    total_st=hcat(scar, thermal)
    @test isapprox(total_st'*total_st, I(22))

    scar=gene_scar(N)
    scar1=storage(scar)

    basis= PXP_basis(T)
    basis_int = [i.buf for i in basis]
    scar1=scar1[basis_int.+1]
    @test isapprox(norm(H*scar1), 0, atol=1e-12)

    exact_scar,exact_scar_prime,thermal_ensemble=sep_scar_exact(T, energy, states)
    total_st = hcat(exact_scar, exact_scar_prime, thermal_ensemble)
    @test isapprox(total_st'*total_st, I(322))
end

@testset "Observables" begin
    scar_indexlis16=[1, 2, 9, 27, 82, 202, 408, 728, 1075, 1480, 1800, 2006, 2126, 2181, 2199, 2206, 2207];
    @test isapprox(EE(ones((100,100))),-460.51701859880916)

    N=12
    T = BitStr{N, Int}
    H = PXP_Ham(T)
    energy, states= eigen(H)
    state=states[:,1]
    splitlis=Vector(1:N-1)
    @test isapprox(fitCCEntEntScal(EE_PXP_state(T,splitlis,state); mincut=1, pbc=true)[1], 0.19684135629232746)

    A, B, C = collect(1:3), collect(4:6), collect(7:9)
    scar, thermal = sep_scar_FSA(T, energy, states)
    @test isapprox(Tri_mutual_information(T, scar, (A, B, C))/log(2), 0.85760514, atol=1e-6)
    @test isapprox(Mutual_information(T, scar, (A, C)), 0.97392703, atol=1e-6)

    z2state=zeros(322)
    z2state[end]=1
    @test isapprox(QFI(domain_wall(T), z2state),0)

    Ob=diagm(domain_wall(T))
    @test isapprox(QFI(Ob, z2state),0)
    @test QFI(domain_wall(BitStr{4, Int}), ones(7)) == 48.0

    P=particlenumber(BitStr{4, Int})
    @test isapprox(QFI(particlenumber(T), z2state),0)
    @test diag(P)==[0.0, 1.0, 1.0, 1.0, 2.0, 1.0, 2.0]
    
    vals, vecs= eigen(PXP_Ham(BitStr{16, Int}))

    @test isapprox(vals[728], -1.34002, atol=1e-6)
    @test isapprox(vals[408], -2.6701435, atol=1e-6)
    @test isapprox(vals[202], -3.97957011, atol=1e-6)
    @test isapprox(vals[82], -5.255931, atol=1e-6)

    vec_728 = vecs[:, 728]
    vec_408 = vecs[:, 408]
    vec_202 = vecs[:, 202]
    vec_82 = vecs[:, 82]

    qfi_728 = QFI(domain_wall(BitStr{16, Int}), vec_728) / 16
    qfi_408 = QFI(domain_wall(BitStr{16, Int}), vec_408) / 16
    qfi_202 = QFI(domain_wall(BitStr{16, Int}), vec_202) / 16
    qfi_82 = QFI(domain_wall(BitStr{16, Int}), vec_82) / 16
    
    @test isapprox(qfi_728, 7.899781, atol=1e-6)
    @test isapprox(qfi_408, 7.6411086, atol=1e-6)
    @test isapprox(qfi_202, 7.1952766, atol=1e-6)
    @test isapprox(qfi_82, 6.55122, atol=1e-6)
end