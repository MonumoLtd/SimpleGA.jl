#Test suite for Geometric Algebra

using GeometricAlgebra
using LinearAlgebra
using Test

include("testfuncs.jl")

@testset "GA Tests" begin
    #! format: off
    @testset "GA(2, 0)" begin include("test20.jl") end
    @testset "GA(3, 0)" begin include("test30.jl") end
    @testset "STA" begin include("testSTA.jl") end
    @testset "GA(4, 0)" begin include("test40.jl") end
    @testset "GA(3, 1)" begin include("test31.jl") end
    @testset "PGA" begin include("testPGA.jl") end
    @testset "CGA" begin include("testCGA.jl") end
    @testset "GA(3, 3)" begin include("test33.jl") end
    @testset "GA(2, 4)" begin include("test24.jl") end
    @testset "GA(4, 4)" begin include("test44.jl") end
    @testset "GA(32, 32)" begin include("test3232.jl") end
    #! format: on
end
