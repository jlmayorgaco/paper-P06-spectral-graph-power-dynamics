include("ReducedDAE.jl")
using .ReducedDAE, LinearAlgebra, TOML, CSV, DataFrames, SHA
const R=ReducedDAE
BLAS.set_num_threads(1)
path=abspath(ARGS[1]);d=TOML.parsefile(path);ctx=R.N.design_context(R.ROOT)
m=R.model(ctx,d["rho"],d["Kp"],d["Ki"];bus=29,delta=-100.,dc_convention=Symbol(d["dc_convention"]))
println("WEAK_TAIL_START");flush(stdout)
s=R.simulate(m;horizon=10000.,dt=.1,tol=1e-9,wall_limit=240.,maxiters=1000000,rotating_frame=true)
met=R.metrics(m,s;dt=.1,window=.5,monitor_buses=collect(1:39))
dest=joinpath(R.OUT,"weak_tail_final");mkpath(dest)
CSV.write(joinpath(dest,"events.csv"),DataFrame([(;bus=29,delta=-100.,met...)]))
open(joinpath(dest,"provenance.toml"),"w") do io
    TOML.print(io,Dict("candidate_sha256"=>bytes2hex(sha256(read(path))),"horizon_s"=>10000.,
        "sample_step_s"=>.1,"scope"=>"Supplemental long-tail numerical audit for the weakest post-disturbance equilibrium. Coarser output mesh; early peaks are checked separately at 0.005/0.01 s. No invariant-region certificate."))
end
println("WEAK_TAIL_RESULT ",met);flush(stdout)
