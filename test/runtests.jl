using Test
using DataFrames

include(joinpath(@__DIR__, "..", "src", "pd39", "PD39.jl"))
using .PD39

@testset "PD39 preregistration scaffold" begin
    @test CANDIDATE_SG_BUSES == [30, 32, 33, 34, 35, 36, 37, 38]
    @test SLACK_BUS == 31
    @test length(candidate_portfolios()) == 256
    @test length(uncertainty_scenarios()) == 9
    @test length(unique(getfield.(uncertainty_scenarios(), :id))) == 9
    candidates = candidate_table()
    @test nrow(candidates) == 8
    @test !(SLACK_BUS in candidates.bus)
end
