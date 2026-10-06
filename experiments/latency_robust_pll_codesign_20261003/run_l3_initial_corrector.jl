using CSV, DataFrames, LinearAlgebra, TOML

module Scan
include(joinpath(@__DIR__, "run_l2_family_scan.jl"))
end

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const PREV=joinpath(ROOT,"experiments","optimal_latency_margin_20261002")
BLAS.set_num_threads(1)

function classify(L,tau_ms)
    q=Scan.Roots.MegaOracle.trace_count(L,fill(tau_ms/1000,10);gamma=-0.05,step=10.0)
    good=abs(real(q.estimate)-q.nearest)<0.01 && abs(imag(q.estimate))<0.01 &&
        abs(q.phase_estimate-q.nearest)<0.01 && q.largest_phase_step<pi/2 &&
        q.quadrature_error_estimate<0.01
    (;status=good ? (q.nearest==0 ? "SAFE" : "UNSAFE") : "INDETERMINATE",
      count=good ? q.nearest : missing,trace_count_real=real(q.estimate),
      phase_count=q.phase_estimate,quadrature_error=q.quadrature_error_estimate,
      min_sigma=q.min_sigma_small,evaluations=q.evaluations)
end

function main()
    ctx=Scan.Roots.MegaOracle.DC.R.N.design_context(ROOT)
    prior=CSV.read(joinpath(PREV,"T08_FAST_MODE_PROVENANCE.csv"),DataFrame)
    catalogue=prior[abs.(prior.tau_ms.-40.0).<1e-8,:]
    seeds=ComplexF64.(catalogue.root_real.+im.*catalogue.root_imag)
    candidates=CSV.read(joinpath(@__DIR__,"L3_INITIAL_PREDICTORS.csv"),DataFrame)
    summaries=NamedTuple[];rootrows=NamedTuple[];counts=NamedTuple[]
    for item in eachrow(candidates)
        id=String(item.candidate_id)
        d=TOML.parsefile(joinpath(@__DIR__,"designs",id*".toml"))
        rho,kp,ki=Float64.(d["rho"]),Float64.(d["Kp"]),Float64.(d["Ki"])
        L=Scan.Roots.MegaOracle.DC.linearization(ctx,rho,kp,ki)
        m=Scan.Roots.reduced_model(L)
        found=NamedTuple[]
        for (sid,s0) in enumerate(seeds)
            r40=Scan.Roots.refine(m,s0,fill(0.040,10);tol=1e-10)
            r40.converged || continue
            any(abs(r40.s-x.s40)<1e-3 for x in found) && continue
            cross=Scan.root_crossing(m,r40.s)
            cross===nothing && continue
            push!(found,(;s40=r40.s,seed_id=sid,tau_ms=cross.tau_ms,s=cross.s,
                frequency_hz=imag(cross.s)/(2pi),residual=cross.residual))
        end
        sort!(found,by=x->x.tau_ms)
        for (rank,r) in enumerate(found)
            push!(rootrows,(;candidate_id=id,rank,seed_id=r.seed_id,
                local_crossing_ms=r.tau_ms,critical_real=real(r.s),
                critical_imag=imag(r.s),critical_frequency_hz=r.frequency_hz,
                small_factor_residual=r.residual,status="LOCAL_EXACT_CHARACTERISTIC_ROOT"))
        end
        CSV.write(joinpath(@__DIR__,"L3_INITIAL_CORRECTOR_ROOTS.csv"),DataFrame(rootrows))
        if isempty(found)
            println("L3_NO_ROOTS ",id);flush(stdout)
            continue
        end
        c=first(found);alpha0=maximum(real,eigvals(L.A))
        qbelow="INDETERMINATE";qabove="INDETERMINATE"
        for tau_ms in (30.0,c.tau_ms-0.03,c.tau_ms+0.03)
            println("L3_COUNT_START ",id," tau_ms=",tau_ms);flush(stdout)
            try
                q=classify(L,tau_ms)
                push!(counts,merge((;candidate_id=id,tau_ms),q))
                if abs(tau_ms-(c.tau_ms-0.03))<1e-8;qbelow=q.status;end
                if abs(tau_ms-(c.tau_ms+0.03))<1e-8;qabove=q.status;end
            catch err
                push!(counts,(;candidate_id=id,tau_ms,status="INDETERMINATE",count=missing,
                    trace_count_real=NaN,phase_count=NaN,quadrature_error=NaN,
                    min_sigma=NaN,evaluations=missing))
                println("L3_COUNT_ERROR ",id," ",sprint(showerror,err));flush(stdout)
            end
            CSV.write(joinpath(@__DIR__,"L3_INITIAL_CORRECTOR_COUNTS.csv"),DataFrame(counts))
        end
        status=alpha0<=-0.05 && qbelow=="SAFE" && qabove=="UNSAFE" ?
            "SPECTRAL_CORRECTOR_PASS_ZERO_DELAY_EVENTS_PENDING" : "SPECTRAL_CORRECTOR_FAIL_OR_INDETERMINATE"
        push!(summaries,(;candidate_id=id,parent=item.parent,method=item.method,
            predicted_gain_ms=item.predicted_gain_ms,alpha_zero_delay=alpha0,
            first_local_crossing_ms=c.tau_ms,critical_frequency_hz=c.frequency_hz,
            second_local_crossing_ms=length(found)>1 ? found[2].tau_ms : NaN,
            fast_roots_discovered=length(found),safe_below_status=qbelow,unsafe_above_status=qabove,
            all_five_events="PENDING",status))
        CSV.write(joinpath(@__DIR__,"L3_INITIAL_CORRECTOR.csv"),DataFrame(summaries))
        println("L3_CORRECTOR ",last(summaries));flush(stdout)
    end
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
