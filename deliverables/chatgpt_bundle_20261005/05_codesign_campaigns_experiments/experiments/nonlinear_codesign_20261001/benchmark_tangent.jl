include("ReducedDAE.jl")
using .ReducedDAE, LinearAlgebra, TOML
const R=ReducedDAE
BLAS.set_num_threads(1);ctx=R.N.design_context(R.ROOT);mkpath(R.OUT)
m=R.model(ctx,fill(.8,10),fill(R.N.K0P,10),fill(R.N.K0I,10);bus=29,delta=-100.)
elapsed=@elapsed trajectory=R.tangent_simulate(m;horizon=60.,dt=.05,monitor_buses=collect(1:39))
peaks=R.tangent_peaks(trajectory)
println("BENCHMARK wall=",elapsed," peaks=",[(p.name,p.peak,p.bus,p.time) for p in peaks]," V=",(trajectory.vmin,trajectory.vmax)," limiter=",trajectory.limiter_fraction);flush(stdout)
adaptive=R.simulate(m;horizon=60.,dt=.005,tol=1e-9,wall_limit=100.)
met=R.metrics(m,adaptive;dt=.005,monitor_buses=collect(1:39))
println("REFERENCE ",met);flush(stdout)
open(joinpath(R.OUT,"tangent_discretization.toml"),"w") do io
    TOML.print(io,Dict("step_s"=>.05,"horizon_s"=>60.,"runtime_s"=>elapsed,
        "monitored_buses"=>collect(1:39),"F_tangent"=>peaks[1].peak,"R_tangent"=>peaks[2].peak,
        "F_adaptive"=>met.Fpeak_Hz,"R_adaptive"=>met.Rpeak_Hz_s,"reference_complete"=>adaptive.ok,
        "F_absolute_error"=>abs(peaks[1].peak-met.Fpeak_Hz),"R_absolute_error"=>abs(peaks[2].peak-met.Rpeak_Hz_s)))
end
