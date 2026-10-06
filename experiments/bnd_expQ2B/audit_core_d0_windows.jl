using CSV,DataFrames,TOML,SHA
const ROOT=normpath(joinpath(@__DIR__,"..",".."));const CORE=joinpath(ROOT,"reports","experiment_Q2B","CORE")
include(joinpath(ROOT,"src","bnd_design_p","ExpP.jl"));include(joinpath(ROOT,"src","bnd_expQ","LinearSecurity.jl"));include(joinpath(ROOT,"src","bnd_expQ2B","FiniteWindow.jl"))
const N=ExpP.PDExactDesignN
d=TOML.parsefile(joinpath(CORE,"Z_Q2B_CORE_d000.toml"));ctx=N.design_context(ROOT)
rho=Float64.(d["rho"]);kp=Float64.(d["Kp"]);ki=Float64.(d["Ki"])
m=N.descriptor(ctx,rho,kp,ki)
out=FiniteWindow.design_metrics(ctx,m,rho;load_bus=16,disturbance_MW=100.0,
    windows=(.2,.5,1.,2.),dt_s=.01,horizon_s=60.,gauge_vector=N.gauge_vector)
CSV.write(joinpath(CORE,"TABLE_Q2B_core_d0_rocof_at100.csv"),out.metrics)
println("CORE_D0_AT100 windows=",join([string(x.window_s,":F=",x.F_peak_Hz,",R=",x.R_peak_Hz_s) for x in out.metrics],";"))
