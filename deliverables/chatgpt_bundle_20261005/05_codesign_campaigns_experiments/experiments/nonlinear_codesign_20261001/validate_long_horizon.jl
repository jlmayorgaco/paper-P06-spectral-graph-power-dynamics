include("ReducedDAE.jl")
using .ReducedDAE, LinearAlgebra, CSV, DataFrames, TOML, SHA
const R=ReducedDAE
BLAS.set_num_threads(1)
ctx=R.N.design_context(R.ROOT)
candidate=abspath(ARGS[1]);d=TOML.parsefile(candidate)
dest=joinpath(R.OUT,length(ARGS)>1 ? ARGS[2] : "long_horizon");mkpath(dest)
rows=NamedTuple[]
for bus in (8,16,29),delta in (100.,-100.)
    println("LONG_START ",bus," ",delta);flush(stdout)
    m=R.model(ctx,d["rho"],d["Kp"],d["Ki"];bus,delta,dc_convention=Symbol(d["dc_convention"]))
    s=R.simulate(m;horizon=600.,dt=.01,tol=1e-9,wall_limit=180.,maxiters=200000,rotating_frame=true)
    met=R.metrics(m,s;dt=.01,window=.5,monitor_buses=collect(1:39))
    push!(rows,(;bus,delta,met...))
    CSV.write(joinpath(dest,"events.csv"),DataFrame(rows))
    println("LONG_RESULT ",last(rows));flush(stdout)
end
open(joinpath(dest,"provenance.toml"),"w") do io
    TOML.print(io,Dict("candidate_sha256"=>bytes2hex(sha256(read(candidate))),
        "horizon_s"=>600.,"scope"=>"Extended finite-horizon validation in the exact moving frame; not an infinite-horizon invariant-set certificate."))
end
