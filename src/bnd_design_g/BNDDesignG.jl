module BNDDesignG

using LinearAlgebra
using CSV
using DataFrames
using SHA
using TOML
using Statistics
using Printf

# The design implementation depends only on audited analytical sources.
# PowerDynamics is deliberately absent from this module.
include(joinpath(@__DIR__, "..", "bnd_design_e", "PhysicalData.jl"))
include(joinpath(@__DIR__, "..", "bnd_design_e", "CollectiveModel.jl"))
include(joinpath(@__DIR__, "..", "bnd_design_e", "ClosureSpectrum.jl"))
include("Endpoint.jl")
include("Robustness.jl")
include("Transients.jl")
include("GraphAblation.jl")

export PhysicalData, CollectiveModel, ClosureSpectrum, file_sha256,
       gauge_vector, zero_structure, spectrum_with_gauge, pll_reference,
       validate_sensitivities, gain_correction, enumerate_surrogate,
       low_rank_boundary, modal_residue_bound, robustness_certificate, transient_analysis,
       graph_features

file_sha256(path) = bytes2hex(sha256(read(path)))

end
