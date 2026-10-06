using CSV, DataFrames, TOML, LinearAlgebra
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const NEV=include(joinpath(ROOT,"experiments","delay_dressed_replacement_frontier_20261002","DDE_NEV.jl"))
const D=NEV.D;const R=D.R
function main()
    d=TOML.parsefile(joinpath(@__DIR__,"seed_uniform_875.toml"));ctx=R.N.design_context(ROOT)
    L=D.linearization(ctx,Float64.(d["rho"]),Float64.(d["Kp"]),Float64.(d["Ki"]))
    rows=NamedTuple[]
    for tau_ms in (20.0,40.0)
        roots=NEV.track_roots(L,fill(tau_ms/1000,length(L.Ai));real_cut=-0.25,
            continuation_steps=6,max_roots=12,maxiter=20)
        for q in roots
            push!(rows,(;tau_ms,root_id=q.root_id,zero_delay_real=real(q.initial),
                real=real(q.s),imag=imag(q.s),frequency_Hz=abs(imag(q.s))/(2pi),
                nev_residual=q.residual,converged=q.converged,iterations=q.iterations,
                right_of_margin=real(q.s)>-0.05,claim="SUPPORTED_LOCAL; root coverage not established"))
        end
    end
    CSV.write(joinpath(@__DIR__,"TABLE_Q03_LOCAL_ROOT_TRACKING.csv"),DataFrame(rows))
    println(DataFrame(rows))
end
main()
