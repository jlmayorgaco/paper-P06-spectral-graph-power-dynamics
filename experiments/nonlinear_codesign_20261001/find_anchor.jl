include("ReducedDAE.jl")
using .ReducedDAE, LinearAlgebra, CSV, DataFrames, TOML
const R=ReducedDAE
BLAS.set_num_threads(1);mkpath(R.OUT)
ctx=R.N.design_context(R.ROOT);rows=NamedTuple[]
kp=fill(R.N.K0P,10);ki=fill(R.N.K0I,10)
for share in (.8,.9,.93)
    rho=fill(share,10)
    sp=R.N.spectrum(ctx,rho,kp,ki)
    for bus in (8,16,29),delta in (100.,-100.)
        println("ANCHOR_START rho=",share," bus=",bus," delta=",delta," alpha=",sp.alpha);flush(stdout)
        m=R.model(ctx,rho,kp,ki;bus,delta)
        r=R.simulate(m;horizon=60.,dt=.01,tol=1e-8,wall_limit=100.)
        met=R.metrics(m,r;dt=.01,savepath=joinpath(R.OUT,"anchor_$(share)_$(bus)_$(Int(delta)).csv"))
        row=(;rho=share,bus,delta,alpha=sp.alpha,met...);push!(rows,row)
        CSV.write(joinpath(R.OUT,"anchor_cases.csv"),DataFrame(rows))
        println("ANCHOR_RESULT ",row);flush(stdout)
        (!r.ok || met.Fpeak_Hz>.5 || met.Rpeak_Hz_s>.5) && break
    end
end
