using Test
using LinearAlgebra
using Random

include(joinpath(@__DIR__, "..", "..", "src", "bnd_graph", "BNDGraph.jl"))

const BG = BNDGraph
include("test_graph_basis.jl")
include("test_energy_identity.jl")
include("test_schur_identity.jl")
include("test_coupling_scaling.jl")
include("test_pole_predictor.jl")
include("test_expA_adapter.jl")
