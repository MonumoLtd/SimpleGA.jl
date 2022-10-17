using GADraft
using Test

include("common.jl")

@testset "GADraft.jl" begin

    @testset "core" begin
        @test !Odd == Even
        @test !Even == Odd
    end

    @testset "ga20" begin include("ga20.jl") end
    @testset "ga30" begin include("ga30.jl") end
end
