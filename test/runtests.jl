using Test
using PXPConstrained, BitBasis
using ITensors
using LinearAlgebra

@testset "pxp functions" begin
    include("./test_PXPFuctions.jl")
end

@testset "scar separate" begin
    include("./test_ScarSeparate.jl")
end

@testset "Observables" begin
    include("./test_Observables.jl")
end

@testset "PXP_MSS_Ham_sparse" begin
    include("./test_PXPSparse.jl")
end

@testset "Dynamics" begin
    include("./test_Dynamics.jl")
end 