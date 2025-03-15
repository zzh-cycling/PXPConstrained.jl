using Test
using PXPConstrained
using LinearAlgebra

@testset "time_evolution" begin
    # Test the time evolution of the wavefunction, it should be normalized under evolution, and the energy should be conserved.
    psi0=zeros(322)
    psi0[end]=1
    N = 12
    H = PXP_Ham(BitStr{N, Int})
    @test size(H) == (322, 322)
    energy, states = eigen(H)
    times = collect(range(0, stop=10, length=100))
    wflis = wf_time_evolution(psi0, times, energy, states)
    @test length(wflis) == 100
    @test isapprox(norm.(wflis), ones(100), atol = 1e-10)
    @test isapprox([norm(wf'*H*wf) for wf in wflis], zeros(100), atol = 1e-10)
end

@testset "time_evolution_mss" begin
    state=zeros(455)
    state[end]=1 # Noted that here the "Z2 state" is the superposition of Z2 and Z2tilde. 
    tlis=collect(0:0.1:20)
    H_MSS=PXP_MSS_Ham_sparse(20, 0)
    wflis=wf_time_evolution_mss(20, 0, state, tlis)
    Elis=similar(tlis)
    for i in 1:length(tlis)
        Elis[i]=norm(wflis[i]'*H_MSS*wflis[i])
    end
    @test norm.(wflis) ≈ ones(length(wflis)) #normalized
    @test isapprox(Elis, zeros(length(Elis)), atol=1e-10)
    # energy conservation
end

@testset "rotated_psi_state" begin
    # Test the rotated state by the on-site rotation exp(i θ/2 Y)
    N = 12
    # for θ in 0.0:0.1:π
    #     rotated_state = rotated_psi_state_mss(N, θ)
    #     @test isapprox(norm(rotated_state), 1.0, atol=1e-10)
    # end
end

@test "rotated_psi_state_mss" begin
    # Test the rotated state by the on-site rotation exp(i θ/2 Y)
    N = 12
    # for θ in 0.0:0.1:π
    #     rotated_state = rotated_psi_state_mss(N, θ)
    #     @test isapprox(norm(rotated_state), 1.0, atol=1e-10)
    # end
    # @test rotated_psi_state_mss(12, 0) ≈ rotated_psi_state_mss(12, π)
end