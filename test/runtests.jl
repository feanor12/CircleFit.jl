using Test
using StatsBase
using CircleFit
using Statistics

@testset "Analytic Fits" begin
    @testset "Pratt" begin
        include("pratt.jl")
    end
    @testset "Taubin" begin
        include("taubin.jl")
    end
    @testset "Kasa" begin
        include("kasa.jl")
    end
    @testset "GRAF" begin
        include("graf.jl")
    end
end
@testset "StatsModel" begin
    include("Circle.jl")
end