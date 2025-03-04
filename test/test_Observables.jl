using Test
using PXPConstrained, BitBasis 
using LinearAlgebra

@testset "Observables" begin
    scar_indexlis16=[1, 2, 9, 27, 82, 202, 408, 728, 1075, 1480, 1800, 2006, 2126, 2181, 2199, 2206, 2207];
    @test isapprox(ee(ones((100,100))),-460.51701859880916)
    
    
    N=12
    T = BitStr{N, Int}
    H = PXP_Ham(T)
    energy, states= eigen(H)
    state=states[:,1]
    splitlis=Vector(1:N-1)
    @test isapprox(fitCCEntEntScal(ee_PXP_state(T,splitlis,state); mincut=1, pbc=true)[1], 0.19684135629232746)
    @test isapprox(fitCCEntEntScal(ee_PXP_idx(T,splitlis,1); mincut=1, pbc=true)[1], 0.19684135629232746)

    cent_cc, _ = ee_PXP_scaling_fig(T, states[:,1], "CC")
    cent_page, _ = ee_PXP_scaling_fig(T, states[:,3], "Page")
    @test isapprox(cent_cc, 0.19684135629232746, atol=1e-3) 
    @test cent_page > 0

    A, B, C = collect(1:3), collect(4:6), collect(7:9)
    scar, thermal = sep_scar_FSA(T, energy, states)
    @test isapprox(tri_mutual_information(T, scar, (A, B, C))/log(2), 0.85760514, atol=1e-6)
    @test isapprox(mutual_information(T, scar, (A, C)), 0.97392703, atol=1e-6)

    z2state=zeros(322)
    z2state[end]=1
    @test isapprox(qfi(domain_wall(T), z2state),0)

    Ob=diagm(domain_wall(T))
    @test isapprox(qfi(Ob, z2state),0)
    @test qfi(domain_wall(BitStr{4, Int}), ones(7)) == 48.0

    P=particlenumber(BitStr{4, Int})
    @test isapprox(qfi(particlenumber(T), z2state),0)
    @test diag(P)==[0.0, 1.0, 1.0, 1.0, 2.0, 1.0, 2.0]
    
    vals, vecs = eigen(on_siten(T,1))
    @test sum(vals) == 89.0

    vals, vecs= eigen(PXP_Ham(BitStr{16, Int}))

    @test isapprox(vals[728], -1.34002, atol=1e-6)
    @test isapprox(vals[408], -2.6701435, atol=1e-6)
    @test isapprox(vals[202], -3.97957011, atol=1e-6)
    @test isapprox(vals[82], -5.255931, atol=1e-6)

    vec_728 = vecs[:, 728]
    vec_408 = vecs[:, 408]
    vec_202 = vecs[:, 202]
    vec_82 = vecs[:, 82]

    qfi_728 = qfi(domain_wall(BitStr{16, Int}), vec_728) / 16
    qfi_408 = qfi(domain_wall(BitStr{16, Int}), vec_408) / 16
    qfi_202 = qfi(domain_wall(BitStr{16, Int}), vec_202) / 16
    qfi_82 = qfi(domain_wall(BitStr{16, Int}), vec_82) / 16
    
    @test isapprox(qfi_728, 7.899781, atol=1e-6) # via Pappalardi et al. 2018, PHYSICAL REVIEW LETTERS 129, 020601 (2022)
    @test isapprox(qfi_408, 7.6411086, atol=1e-6)
    @test isapprox(qfi_202, 7.1952766, atol=1e-6)
    @test isapprox(qfi_82, 6.55122, atol=1e-6)

    GS_energy, subenergy, passive_energy = ergotropy_PXP_state(BitStr{16, Int}, 8, vecs[:,1074])
    W = (GS_energy - passive_energy)/16
    @test isapprox(W, 0.16915174104683417, atol=1e-6)

    GS_energy, subenergy, passive_energy = ergotropy_PXP_idx(BitStr{16, Int}, 8, 1074)
    W = (GS_energy - passive_energy)/16
    @test isapprox(W, 0.16915174104683417, atol=1e-6)

    GS_energy, subenergy, passive_energy = ergotropy_PXP_idx_OBC(BitStr{16, Int}, 8, 1300)
    W = (GS_energy - passive_energy)/16
    @test isapprox(W, 0.16076624524701522, atol=1e-6)
end

