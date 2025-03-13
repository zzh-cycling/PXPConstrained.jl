using Test
using PXPConstrained, BitBasis 
using LinearAlgebra

@testset "Observables" begin
    # This is the natural order index for scar state in PXP model's eigen
    scar_indexlis16=[1, 2, 9, 27, 82, 202, 408, 728, 1075, 1480, 1800, 2006, 2126, 2181, 2199, 2206, 2207];
    @test isapprox(ee(ones((100,100))),-460.51701859880916)
    
    
    N=12
    T = BitStr{N, Int}
    H = PXP_Ham(T)
    energy, states= eigen(H)
    state=states[:,1]
    splitlis=Vector(1:N-1)

    # Fit the central charge of 1st scar
    @test isapprox(fitCCEntEntScal(ee_PXP_state(N,splitlis,state); mincut=1, pbc=true)[1], 0.19684135629232746)
    @test isapprox(fitCCEntEntScal(ee_PXP_idx(N,splitlis,1); mincut=1, pbc=true)[1], 0.19684135629232746)

    cent_cc, _ = ee_PXP_scaling_fig(N, states[:,1], "CC")
    cent_page, _ = ee_PXP_scaling_fig(N, states[:,3], "Page")
    @test isapprox(cent_cc, 0.19684135629232746, atol=1e-3) 
    @test cent_page > 0

    # calculate the FSA tmi and mi.
    A, B, C = collect(1:3), collect(4:6), collect(7:9)
    scar, thermal = sep_scar_FSA(N, energy, states)
    @test isapprox(tri_mutual_information(N, (A, B, C), scar)/log(2), 0.85760514, atol=1e-6)
    @test isapprox(mutual_information(N, (A, C), scar), 0.97392703, atol=1e-6)

    # The qfi_dw of Z2 state should be 0
    z2state=zeros(322)
    z2state[end]=1
    @test isapprox(qfi(domain_wall(N), z2state),0)

    # The qfi_dw of Z2 state should be 0(use matrix form function), and qfi_dw of 1 state should be 48
    Ob=diagm(domain_wall(N))
    @test isapprox(qfi(Ob, z2state),0)
    @test qfi(domain_wall(4), ones(7)) == 48.0

    # The qfi_particlenumber of Z2 state should be 0, and particlenumber is [0.0, 1.0, 1.0, 1.0, 2.0, 1.0, 2.0]
    P=particlenumber(4)
    @test isapprox(qfi(particlenumber(N), z2state),0)
    @test diag(P)==[0.0, 1.0, 1.0, 1.0, 2.0, 1.0, 2.0]
    
    vals, vecs = eigen(on_siten(N, 1))
    @test sum(vals) == 89.0

    vals, vecs= eigen(PXP_Ham(BitStr{16, Int}))

    # According to the FSA prediction, theire energy should be equal distant.
    @test isapprox(vals[728], -1.34002, atol=1e-6)
    @test isapprox(vals[408], -2.6701435, atol=1e-6)
    @test isapprox(vals[202], -3.97957011, atol=1e-6)
    @test isapprox(vals[82], -5.255931, atol=1e-6)

    vec_728 = vecs[:, 728]
    vec_408 = vecs[:, 408]
    vec_202 = vecs[:, 202]
    vec_82 = vecs[:, 82]

    qfi_728 = qfi(domain_wall(16), vec_728) / 16
    qfi_408 = qfi(domain_wall(16), vec_408) / 16
    qfi_202 = qfi(domain_wall(16), vec_202) / 16
    qfi_82 = qfi(domain_wall(16), vec_82) / 16
    
    @test isapprox(qfi_728, 7.899781, atol=1e-6) # qfi value via Pappalardi et al. 2018, PHYSICAL REVIEW LETTERS 129, 020601 (2022)
    @test isapprox(qfi_408, 7.6411086, atol=1e-6)
    @test isapprox(qfi_202, 7.1952766, atol=1e-6)
    @test isapprox(qfi_82, 6.55122, atol=1e-6)

    # calculate their ergotropy.
    GS_energy, subenergy, passive_energy = ergotropy_PXP_state(16, 8, vecs[:,1074])
    W = (GS_energy - passive_energy)/16
    @test isapprox(W, 0.16915174104683417, atol=1e-6)

    GS_energy, subenergy, passive_energy = ergotropy_PXP_idx(16, 8, 1074)
    W = (GS_energy - passive_energy)/16
    @test isapprox(W, 0.16915174104683417, atol=1e-6)

    GS_energy, subenergy, passive_energy = ergotropy_PXP_idx_OBC(16, 8, 1200)
    W = (GS_energy - passive_energy)/16
    @test isapprox(W,0.19446000922247064, atol=1e-3)
end

@testset "ergotropy in MSS" begin
    state=zeros(455)
    state[end]=1
    @test ergotropy_PXP_MSS_state(20, 10, state) ≈ 0
    total_state=iso_total2MSS(20, state)
end