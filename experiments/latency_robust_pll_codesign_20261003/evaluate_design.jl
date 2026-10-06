using CSV, DataFrames, LinearAlgebra, TOML

module Scan
include(joinpath(@__DIR__,"run_l2_family_scan.jl"))
end

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const PREV=joinpath(ROOT,"experiments","optimal_latency_margin_20261002")
BLAS.set_num_threads(1)

function classify(L,tau_ms)
    r=Scan.Roots.MegaOracle.trace_count(L,fill(tau_ms/1000,10);gamma=-0.05,step=10.0)
    agreed=abs(real(r.estimate)-r.nearest)<0.01 && abs(imag(r.estimate))<0.01 &&
        abs(r.phase_estimate-r.nearest)<0.01 && r.largest_phase_step<pi/2 &&
        r.quadrature_error_estimate<0.01
    (;status=agreed ? (r.nearest==0 ? "SAFE" : "UNSAFE") : "INDETERMINATE",
      root_count=agreed ? r.nearest : missing,trace_real=real(r.estimate),
      phase_count=r.phase_estimate,quadrature_error=r.quadrature_error_estimate,
      min_sigma=r.min_sigma_small,evaluations=r.evaluations)
end

function main()
    length(ARGS)==1 || error("usage: evaluate_design.jl design.toml")
    design_path=abspath(ARGS[1]);id=splitext(basename(design_path))[1]
    out=joinpath(@__DIR__,"evaluations",id);mkpath(out)
    d=TOML.parsefile(design_path)
    rho,kp,ki=Float64.(d["rho"]),Float64.(d["Kp"]),Float64.(d["Ki"])
    ctx=Scan.Roots.MegaOracle.DC.R.N.design_context(ROOT)
    L=Scan.Roots.MegaOracle.DC.linearization(ctx,rho,kp,ki)
    m=Scan.Roots.reduced_model(L)
    old=CSV.read(joinpath(PREV,"T08_FAST_MODE_PROVENANCE.csv"),DataFrame)
    prior=old[abs.(old.tau_ms.-40.0).<1e-8,:]
    seeds=ComplexF64.(prior.root_real.+im.*prior.root_imag)
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
    isempty(found) && error("no fast root crossings discovered")
    rootrows=[(;candidate_id=id,rank=j,seed_id=x.seed_id,local_crossing_ms=x.tau_ms,
        critical_real=real(x.s),critical_imag=imag(x.s),frequency_hz=x.frequency_hz,
        small_factor_residual=x.residual,status="LOCAL_ROOT_REFINED_EXACT_CHARACTERISTIC")
        for (j,x) in enumerate(found)]
    CSV.write(joinpath(out,"ROOTS.csv"),DataFrame(rootrows))
    c=first(found);counts=NamedTuple[]
    for tau_ms in (20.0,30.0,c.tau_ms-0.03,c.tau_ms+0.03)
        println("EVALUATE_COUNT_START ",id," ",tau_ms);flush(stdout)
        try
            r=classify(L,tau_ms)
            push!(counts,merge((;candidate_id=id,tau_ms),r))
        catch err
            push!(counts,(;candidate_id=id,tau_ms,status="INDETERMINATE",
                root_count=missing,trace_real=NaN,phase_count=NaN,
                quadrature_error=NaN,min_sigma=NaN,evaluations=missing))
            println("EVALUATE_COUNT_ERROR ",sprint(showerror,err));flush(stdout)
        end
        CSV.write(joinpath(out,"ROOT_COUNTS.csv"),DataFrame(counts))
    end
    alpha0=maximum(real,eigvals(L.A))
    pass=alpha0<=-0.05 && counts[1].status=="SAFE" &&
        counts[2].status=="SAFE" && counts[3].status=="SAFE" && counts[4].status=="UNSAFE"
    summary=(;candidate_id=id,alpha_zero_delay_s_inv=alpha0,
        local_tau_crit_ms=c.tau_ms,critical_frequency_hz=c.frequency_hz,
        critical_root_real=real(c.s),critical_root_imag=imag(c.s),
        second_local_crossing_ms=length(found)>1 ? found[2].tau_ms : NaN,
        third_local_crossing_ms=length(found)>2 ? found[3].tau_ms : NaN,
        fast_roots_discovered=length(found),complete_contour_bracket_pass=pass,
        zero_delay_events="PENDING",status=pass ?
            "SPECTRAL_CORRECTOR_PASS_EVENTS_PENDING" : "SPECTRAL_CORRECTOR_FAIL_OR_INDETERMINATE")
    CSV.write(joinpath(out,"SPECTRAL.csv"),DataFrame([summary]))
    println("EVALUATE_DONE ",summary);flush(stdout)
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
