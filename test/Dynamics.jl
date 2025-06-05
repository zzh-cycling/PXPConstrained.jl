using Test
using PXPConstrained
using LinearAlgebra

# TEST need to test for 4N and 4N+2
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

    psi0=zeros(843)
    psi0[end]=1
    N = 14
    H = PXP_Ham(BitStr{N, Int})
    @test size(H) == (843, 843)
    energy, states = eigen(H)
    times = collect(range(0, stop=10, length=100))
    wflis = wf_time_evolution(psi0, times, energy, states)
    @test length(wflis) == 100
    @test isapprox(norm.(wflis), ones(100), atol = 1e-10)
    @test isapprox([norm(wf'*H*wf) for wf in wflis], zeros(100), atol = 1e-10)
end

@testset "time_evolution_sparse" begin
    N =20
    state=zeros(455)
    state[end]=1 # Noted that here the "Z2 state" is the superposition of Z2 and Z2tilde. 
    tlis=collect(0:0.1:20)
    H_MSS=PXP_MSS_Ham_sparse(N, 0)
    wflis=wf_time_evolution_sparse(N, 0, state, tlis)
    Elis=similar(tlis)
    for i in 1:length(tlis)
        Elis[i]=norm(wflis[i]'*H_MSS*wflis[i])
    end
    @test norm.(wflis) ≈ ones(length(wflis)) #normalized
    @test isapprox(Elis, zeros(length(Elis)), atol=1e-10)
    # energy conservation

    N = 22
    state=zeros(1022)
    state[end]=1 # Noted that here the "Z2 state" is the superposition of Z2 and Z2tilde.
    tlis=collect(0:0.1:20)
    H_MSS=PXP_MSS_Ham_sparse(N, 0)
    wflis=wf_time_evolution_sparse(N, 0, state, tlis)
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
    for θ in 0.0:0.1:π
        rotated_state = rotated_psi_state(N, θ)
        @test isapprox(norm(rotated_state), 1.0, atol=1e-10)
    end
    @test isapprox(norm(rotated_psi_state(12, π)), 1.0, atol=1e-10)
    @test isapprox(rotated_psi_state(12, 0)[end], rotated_psi_state(12, π)[233], atol=1e-10)
    N =14 
    for θ in 0.0:0.1:π
        rotated_state = rotated_psi_state(N, θ)
        @test isapprox(norm(rotated_state), 1.0, atol=1e-10)
    end
    @test isapprox(norm(rotated_psi_state(14, 0)), 1.0, atol=1e-10)
    @test isapprox(rotated_psi_state(14, 0)[end], abs(rotated_psi_state(14, π)[610]), atol=1e-10)
end

@testset "rotated_psi_state_mss" begin
    # Test the rotated state by the on-site rotation exp(i θ/2 Y)
    N = 12
    for θ in 0.0:0.1:π
        rotated_state = rotated_psi_state_mss(N, 0, θ)
        @test isapprox(norm(rotated_state), 1.0, atol=1e-10)
    end
    @test norm(rotated_psi_state_mss(12,0, 0)) ≈ 1.0 atol = 1e-10
    @test rotated_psi_state_mss(12,0, 0) ≈ rotated_psi_state_mss(12,0, π) atol = 1e-10
    N = 14
    for θ in 0.0:0.1:π
        rotated_state = rotated_psi_state_mss(N, 0, θ)
        @test isapprox(norm(rotated_state), 1.0, atol=1e-10)
    end
    @test norm(rotated_psi_state_mss(14,0, 0)) ≈ 1.0 atol = 1e-10
    @test abs.(rotated_psi_state_mss(14,0, 0)) ≈ abs.(rotated_psi_state_mss(14,0, π)) atol = 1e-10

end

# @testset "rotated_psi_state_mss_minus" begin
#     # Test the rotated state by the on-site rotation exp(i θ/2 Y)
#     N = 12
#     for θ in 0.0:0.1:π
#         rotated_state = rotated_psi_state_mss(N, 0, θ)
#         @test isapprox(norm(rotated_state), 1.0, atol=1e-10)
#     end
#     @test norm(rotated_psi_state_mss(12,0, 0)) ≈ 1.0 atol = 1e-10
#     @test rotated_psi_state_mss(12,0, 0) ≈ rotated_psi_state_mss(12,0, π) atol = 1e-10
#     N = 14
#     for θ in 0.0:0.1:π
#         rotated_state = rotated_psi_state_mss(N, 0, θ)
#         @test isapprox(norm(rotated_state), 1.0, atol=1e-10)
#     end
#     @test norm(rotated_psi_state_mss(14,0, 0)) ≈ 1.0 atol = 1e-10
#     @test abs.(rotated_psi_state_mss(14,0, 0)) ≈ abs.(rotated_psi_state_mss(14,0, π)) atol = 1e-10

# end

#even_zeros, even_ones, odd_zeros, odd_ones
@testset "count_zeros_and_ones" begin
    # Z2tilde state
    @test PXPConstrained.count_zeros_and_ones(BitStr{12, Int}(0b010101010101)) == (6, 0, 0, 6)
    # Z2 state：10101010  
    @test PXPConstrained.count_zeros_and_ones(BitStr{12, Int}(0b101010101010)) == (0, 6, 6, 0)
    @test PXPConstrained.count_zeros_and_ones(BitStr{12, Int}(0b111111111111)) == (0, 6, 0, 6)
    @test PXPConstrained.count_zeros_and_ones(BitStr{12, Int}(0b000000000000)) == (6, 0, 6, 0)
end

#10101010
@testset "Z2_overlap" begin
    @test PXPConstrained.Z2_overlap((0, 6, 6, 0), 0) ≈ 1.0
    @test PXPConstrained.Z2_overlap((0, 6, 6, 0), π) ≈ 0.0 atol = 1e-20
end

#01010101
@testset "Z2tilde_overlap" begin
    @test PXPConstrained.Z2tilde_overlap((6, 0, 0, 6), 0) ≈ 1.0
    @test PXPConstrained.Z2tilde_overlap((6, 0, 0, 6), π) ≈ 0.0 atol = 1e-20
end
