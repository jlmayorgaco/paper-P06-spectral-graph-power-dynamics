using CSV, DataFrames, TOML, SHA
const EXP_ROOT=normpath(joinpath(@__DIR__,"..",".."))
const EXP_OUT=joinpath(EXP_ROOT,"reports","experiment_P","P5")
const CANDIDATE=joinpath(EXP_OUT,"Z_P_NOMINAL_FINAL.toml")
bytes2hex(sha256(read(CANDIDATE)))==strip(read(CANDIDATE*".sha256",String)) ||
    error("frozen P5 candidate SHA check failed")
include(joinpath(EXP_ROOT,"src","bnd_design_p","ExpP.jl"))
using .ExpP
include(joinpath(EXP_ROOT,"src","bnd_opt_expP","LinearPulse.jl"))
using .LinearPulse
c=TOML.parsefile(CANDIDATE)
ctx=ExpP.PDExactDesignN.design_context(EXP_ROOT)
frames=CSV.read(joinpath(EXP_OUT,"TABLE_P5_tds_trajectories.csv"),DataFrame)
comparison=LinearPulse.compare_frozen_pulses(ctx,Float64.(c["rho"]),
    Float64.(c["Kp"]),Float64.(c["Ki"]),frames)
CSV.write(joinpath(EXP_OUT,"TABLE_P5_linear_nonlinear_comparison.csv"),comparison)
println("P5_LINEAR_TDS_COMPARISON_ROWS=",nrow(comparison))
println("P5_LINEAR_TDS_MAX_FREQUENCY_REL_ERROR=",maximum(comparison.frequency_relative_error))
println("P5_LINEAR_TDS_MAX_ROCOF_REL_ERROR=",maximum(comparison.rocof_relative_error))
println("P5_LINEAR_TDS_MATCH_PASS=",all(comparison.match_pass))
