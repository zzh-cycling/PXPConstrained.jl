using Test
using PXPConstrained, BitBasis
using ITensors
using LinearAlgebra
@testset "pxp basis" begin
    include("./PXPBasis.jl")
    include("./PXPSymmetry.jl")
end

@testset "scar separate" begin
    include("./ScarSeparate.jl")
end

@testset "Observables" begin
    include("./Observables.jl")
end

@testset "PXP_MSS_Ham_sparse" begin
    include("./PXPSparse.jl")
end

@testset "Dynamics" begin
    include("./Dynamics.jl")
end 