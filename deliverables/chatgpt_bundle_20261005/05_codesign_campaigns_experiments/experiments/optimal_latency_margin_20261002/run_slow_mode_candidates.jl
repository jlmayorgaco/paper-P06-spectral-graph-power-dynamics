using CSV, DataFrames, LinearAlgebra, TOML

module DDE
include(joinpath(@__DIR__,"..","delay_dressed_replacement_frontier_20261002","DDE_NEV.jl"))
end

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const MEGA=joinpath(ROOT,"experiments","analytical_delay_codesign_mega_20261002")
BLAS.set_num_threads(1)

function main()
    ctx=DDE.DDENEV.D.R.N.design_context(ROOT)
    grid=(0.0,5.0,10.0,15.0,20.0,25.0,30.0,32.0,34.0,36.0,38.0,40.0)
    rows=NamedTuple[]
    for (design_id,path) in (("seed_875",joinpath(MEGA,"seed_uniform_875.toml")),
                             ("best_zero_delay_88455",joinpath(MEGA,"M1_ZERO_DELAY_DESIGN.toml")))
        d=TOML.parsefile(path)
        L=DDE.DDENEV.D.linearization(ctx,Float64.(d["rho"]),Float64.(d["Kp"]),Float64.(d["Ki"]))
        eig=eigen(ComplexF64.(L.A));positive=findall(imag.(eig.values).>0)
        idx=positive[argmax(real.(eig.values[positive]))]
        s=eig.values[idx];v=eig.vectors[:,idx]
        tau_previous_ms=0.0
        for tau_ms in grid
            points=tau_ms==tau_previous_ms ? Float64[] :
                collect(range(tau_previous_ms,tau_ms,length=max(2,ceil(Int,(tau_ms-tau_previous_ms)/2)+1)))[2:end]
            converged=true;residual=0.0
            for point in points
                rr=DDE.DDENEV.newton_root(L,s,fill(point/1000,length(L.Ai)),v;
                    tol=2e-10,maxiter=25,max_step=2.0)
                s,v,residual=rr.s,rr.v,rr.residual
                if !rr.converged
                    converged=false;break
                end
            end
            push!(rows,(;design_id,tau_ms,root_real=real(s),root_imag=imag(s),
                frequency_hz=imag(s)/(2pi),residual,converged,
                family_id="slow_ode_$(idx)",status=converged ? "SUPPORTED_LOCAL" : "TRACKING_FAILED"))
            CSV.write(joinpath(@__DIR__,"T01_SLOW_MODE_CANDIDATES.csv"),DataFrame(rows))
            println("SLOW_MODE ",last(rows));flush(stdout)
            converged || break
            tau_previous_ms=tau_ms
        end
    end
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
