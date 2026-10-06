include("ReducedDAE.jl")
using .ReducedDAE, LinearAlgebra, CSV, DataFrames, TOML
const R=ReducedDAE
BLAS.set_num_threads(1)
ctx=R.N.design_context(R.ROOT)
rows=NamedTuple[]
for (bus,delta) in [(8,0.),(8,100.),(8,-100.),(16,100.),(16,-100.),(29,100.),(29,-100.)]
    net=R.network(ctx,bus,delta)
    H=Symmetric(-(net.reduced+net.reduced')/2)
    gamma=eigmin(H)
    passive_sigma=minimum(svdvals(net.Y[1:58,1:58]))
    push!(rows,(;bus,delta,port_accretivity_lower_bound=gamma,passive_sigma,
        inverse_port_norm_upper_bound=gamma>0 ? 1/gamma : Inf))
end
machines=[ctx.net.sg[b].op.parameters for b in 30:39]
conditions=all(p.rs==0 && p.xdpp==p.xqpp for p in machines)
CSV.write(joinpath(R.OUT,"algebraic_regularity.csv"),DataFrame(rows))
open(joinpath(R.OUT,"algebraic_regularity.toml"),"w") do io
    TOML.print(io,Dict("zero_stator_resistance_equal_subtransient_reactances"=>conditions,
        "all_declared_cases_positive_port_accretivity"=>all(r.port_accretivity_lower_bound>0 for r in rows),
        "scope"=>"State-independent port regularity conditional on nonzero SG speed and fixed Z-load admittance. Numerical eigenvalue evaluation, not interval arithmetic; no claim for other load models or a continuous disturbance family."))
end
show(stdout,"text/plain",DataFrame(rows));println()
