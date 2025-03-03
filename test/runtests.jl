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