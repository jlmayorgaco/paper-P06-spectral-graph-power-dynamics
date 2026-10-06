include("ReducedDAE.jl")
using .ReducedDAE, LinearAlgebra, TOML, CSV, DataFrames
const R=ReducedDAE
BLAS.set_num_threads(1);ctx=R.N.design_context(R.ROOT)
d=TOML.parsefile(joinpath(R.ROOT,"reports","codesign_validation_20261001","candidate_improved.toml"))
m=R.model(ctx,d["rho"],d["Kp"],d["Ki"];bus=8,delta=100.,dc_convention=:physical_supply)
r=R.simulate(m;horizon=5.,dt=.0025,tol=1e-10,wall_limit=120.)
met=R.metrics(m,r;dt=.0025,monitor_buses=collect(1:39),savepath=joinpath(R.OUT,"physical_counterexample.csv"))
CSV.write(joinpath(R.OUT,"physical_counterexample_metrics.csv"),DataFrame([met]))
println("PHYSICAL_COUNTEREXAMPLE ",met);flush(stdout)
