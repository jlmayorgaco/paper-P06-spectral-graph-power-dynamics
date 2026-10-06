using CSV, DataFrames, LinearAlgebra, TOML

module Scan
include(joinpath(@__DIR__,"run_l2_family_scan.jl"))
end

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
BLAS.set_num_threads(1)

function main()
    length(ARGS) in (1,2,3) || error("usage: discover_boundary_roots.jl design.toml [fine] [high]")
    path=abspath(ARGS[1]);id=splitext(basename(path))[1]
    out=joinpath(@__DIR__,"evaluations",id);mkpath(out)
    d=TOML.parsefile(path)
    rho,kp,ki=Float64.(d["rho"]),Float64.(d["Kp"]),Float64.(d["Ki"])
    ctx=Scan.Roots.MegaOracle.DC.R.N.design_context(ROOT)
    L=Scan.Roots.MegaOracle.DC.linearization(ctx,rho,kp,ki)
    model=Scan.Roots.reduced_model(L)
    counts=CSV.read(joinpath(out,"ROOT_COUNTS.csv"),DataFrame)
    high=length(ARGS)==3 && ARGS[3]=="high"
    unsafe=counts[counts.status.=="UNSAFE",:]
    tau_ms=isempty(unsafe) ? NaN : (high ? maximum(unsafe.tau_ms) : minimum(unsafe.tau_ms))
    expected=isempty(unsafe) ? 0 : only(unsafe[unsafe.tau_ms.==tau_ms,:].root_count)
    if high && isfile(joinpath(out,"ADDITIONAL_COUNTS.csv"))
        extra=CSV.read(joinpath(out,"ADDITIONAL_COUNTS.csv"),DataFrame)
        extra_unsafe=extra[extra.status.=="UNSAFE",:]
        if !isempty(extra_unsafe) && (!isfinite(tau_ms) || maximum(extra_unsafe.tau_ms)>tau_ms)
            tau_ms=maximum(extra_unsafe.tau_ms)
            expected=only(extra_unsafe[extra_unsafe.tau_ms.==tau_ms,:].root_count)
        end
    end
    isfinite(tau_ms) || error("need known UNSAFE contour in ROOT_COUNTS or ADDITIONAL_COUNTS")
    tau=fill(tau_ms/1000,10)
    starts=Tuple{Float64,ComplexF64}[]
    fine=length(ARGS)==2 && ARGS[2]=="fine"
    grid=collect(20.0:(fine ? 0.025 : 0.125):48.0)
    real_lines=fine ? (-0.05,-0.04,-0.03,-0.02,0.0,0.05,0.2,0.5,1.0,2.0) :
        (-0.05,0.0,0.2,0.5,1.0,2.0)
    for re in real_lines
        vals=Float64[]
        for omega in grid
            F=Scan.Roots.small_factor(model,re+im*omega,tau)[1]
            push!(vals,minimum(svdvals(F)))
        end
        for j in 2:length(grid)-1
            vals[j]<vals[j-1] && vals[j]<vals[j+1] && vals[j]<1.0 &&
                push!(starts,(vals[j],re+im*grid[j]))
        end
    end
    suffix=high ? "_HIGH" : ""
    CSV.write(joinpath(out,"DISCOVERY_SEEDS"*suffix*".csv"),DataFrame(
        [(;candidate_id=id,tau_ms,score,initial_real=real(s),initial_imag=imag(s)) for (score,s) in starts]))
    sort!(starts,by=first)
    roots=NamedTuple[]
    for (score,s0) in first(starts,min(length(starts),250))
        r=Scan.Roots.refine(model,s0,tau;tol=1e-10)
        r.converged && real(r.s)>-0.05 && imag(r.s)>0 || continue
        any(abs(r.s-z.s)<1e-5 for z in roots) && continue
        D=Scan.Roots.MegaOracle.DC.delta_matrix(L,r.s,tau)
        v=svd(D).V[:,end]
        residual=norm(D*v)/max(1,norm(D)*norm(v))
        residual<1e-7 || continue
        push!(roots,(;s=r.s,score,small_residual=r.sv,full_residual=residual))
    end
    sort!(roots,by=x->real(x.s),rev=true)
    rows=NamedTuple[]
    for (j,r) in enumerate(roots)
        crossing=Scan.root_crossing(model,r.s;initial=tau_ms)
        push!(rows,(;candidate_id=id,root_id=j,unsafe_tau_ms=tau_ms,
            root_real_at_unsafe=real(r.s),root_imag_at_unsafe=imag(r.s),
            frequency_at_unsafe_hz=imag(r.s)/(2pi),
            small_factor_residual=r.small_residual,full_characteristic_residual=r.full_residual,
            local_crossing_ms=crossing===nothing ? NaN : crossing.tau_ms,
            frequency_at_crossing_hz=crossing===nothing ? NaN : imag(crossing.s)/(2pi),
            crossing_root_residual=crossing===nothing ? NaN : crossing.residual,
            status=crossing===nothing ? "UNSAFE_ROOT_FOUND_CROSSING_UNRESOLVED" :
                "FULL_CHARACTERISTIC_ROOT_AND_LOCAL_CROSSING"))
    end
    CSV.write(joinpath(out,"DISCOVERED_ROOTS"*suffix*".csv"),DataFrame(rows))
    observed=2*length(roots)
    status=observed==expected && all(isfinite(r.local_crossing_ms) for r in rows) ?
        "COMPLETE_NUMERICAL_ROOT_COVERAGE_AT_UNSAFE_CONTOUR" : "INCOMPLETE_ROOT_COVERAGE"
    CSV.write(joinpath(out,"ROOT_COVERAGE"*suffix*".csv"),DataFrame([(;candidate_id=id,tau_ms,
        contour_root_count=expected,refined_conjugate_root_count=observed,
        distinct_positive_imag_roots=length(roots),status)]))
    println("DISCOVER_DONE ",id," expected=",expected," found=",observed," status=",status)
    for row in rows;println("DISCOVER_ROOT ",row);end
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
