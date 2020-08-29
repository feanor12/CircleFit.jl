using Test
using CircleFit

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
end