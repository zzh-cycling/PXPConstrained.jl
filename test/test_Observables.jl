using Test
using PXPConstrained, BitBasis 
using LinearAlgebra
using Yao: arrayreg, density_matrix, von_neumann_entropy, expect, matblock

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

    GS_energy, subenergy, passive_energy = ergotropy_PXP_idx(16, 8, 1200, false)
    W = (GS_energy - passive_energy)/16
    @test isapprox(W,0.19446000922247064, atol=1e-3)
end

@testset "tmi,mi" begin
    # calculate the FSA tmi and mi.
    N=12
    energy, states= eigen(PXP_Ham(N))
    A, B, C = collect(1:3), collect(4:6), collect(7:9)
    scar, thermal = sep_scar_FSA(N, energy, states)
    scar_total=zeros(2^12)
    scar_total[[i.buf+1 for i in PXP_basis(12)]]=scar
    @test norm(scar_total) ≈ 1
    
    reg = arrayreg(scar_total)
    ρ_A = density_matrix(reg, A)
    ρ_B = density_matrix(reg, B)
    ρ_C = density_matrix(reg, C)
    ρ_AB = density_matrix(reg, vcat(A, B))
    ρ_BC = density_matrix(reg, vcat(B, C))
    ρ_AC = density_matrix(reg, vcat(A, C))
    ρ_ABC = density_matrix(reg, vcat(A, B, C))
    
    valsA=filter(x -> !isapprox(x, 0.0; atol=1e-10), eigvals(ρ_A.state))
    valsB=filter(x -> !isapprox(x, 0.0; atol=1e-10), eigvals(ρ_B.state))
    valsC=filter(x -> !isapprox(x, 0.0; atol=1e-10), eigvals(ρ_C.state))
    valsAB=filter(x -> !isapprox(x, 0.0; atol=1e-10), eigvals(ρ_AB.state))
    valsBC=filter(x -> !isapprox(x, 0.0; atol=1e-10), eigvals(ρ_BC.state))
    valsAC=filter(x -> !isapprox(x, 0.0; atol=1e-10), eigvals(ρ_AC.state))
    valsABC=filter(x -> !isapprox(x, 0.0; atol=1e-10), eigvals(ρ_ABC.state))
    @test valsA ≈ valsB ≈ valsC  ≈ valsABC
    @test eigvals(rdm_PXP(N, [A], scar)) ≈ valsA
    @test eigvals(rdm_PXP(N, [B], scar)) ≈ valsB
    @test eigvals(rdm_PXP(N, [C], scar)) ≈ valsC
    @test eigvals(rdm_PXP(N, [A, B], scar)) ≈ valsAB
    @test_broken eigvals(rdm_PXP(N, [B, C], scar)) ≈ valsBC
    @test_broken eigvals(rdm_PXP(N, [A, C], scar)) ≈ valsAC
    @test_broken eigvals(rdm_PXP(N, [A, B, C], scar)) ≈ valsABC
    
    S_A= von_neumann_entropy(ρ_A)
    S_B= von_neumann_entropy(ρ_B)
    S_C= von_neumann_entropy(ρ_C)
    S_AB= von_neumann_entropy(ρ_AB)
    S_BC= von_neumann_entropy(ρ_BC)
    S_AC= von_neumann_entropy(ρ_AC)
    S_ABC= von_neumann_entropy(ρ_ABC)
    # Calculate the mutual information
    I_ABC = S_A + S_B + S_C - S_AB - S_BC - S_AC + S_ABC
    
    # Calculate the mutual information
    I_AB = S_A + S_C - S_AC
    @test isapprox(I_AB, 0.97392703, atol=1e-6)
    @test isapprox(I_ABC/log(2), 0.85760514, atol=1e-6)
    @test isapprox(tri_mutual_information(N, (A, B, C), scar)/log(2), 0.85760514, atol=1e-6)
    @test isapprox(mutual_information(N, (A, C), scar), 0.97392703, atol=1e-6)
end

@testset "ergotropy in MSS" begin
    state=zeros(455)
    state[end]=1
    GS_energy, subenergy, passive_energy= ergotropy_PXP_MSS_state(20, 10, state)
    W = (GS_energy - passive_energy)/20
    @test isapprox(W, 0.28935025927242286, atol=1e-6)
    total_state=iso_total2MSS(20, 0)*state
    GS_energy1, subenergy1, passive_energy1= ergotropy_PXP_state(20, 10, total_state)
    W = (GS_energy1 - passive_energy1)/20
    @test isapprox(W, 0.28935025927242286, atol=1e-6)
    @test isapprox(GS_energy, GS_energy1, atol=1e-6)
    @test isapprox(passive_energy, passive_energy1, atol=1e-6)
    @test isapprox(subenergy, subenergy1, atol=1e-6)
end

