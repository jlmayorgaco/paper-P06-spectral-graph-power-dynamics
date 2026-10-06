include(joinpath(@__DIR__, "..", "..", "src", "bnd_graph", "BNDGraph.jl"))

const REPOSITORY_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
BNDGraph.Reporting.run_experiment(REPOSITORY_ROOT)
