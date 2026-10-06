using CSV, DataFrames, LinearAlgebra, TOML

module MegaOracle
include(joinpath(@__DIR__,"..","analytical_delay_codesign_mega_20261002","m3_a_trace_integral.jl"))
end
module FastRoots
include(joinpath(@__DIR__,"run_fast_root_locus.jl"))
end

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const MEGA=joinpath(ROOT,"experiments","analytical_delay_codesign_mega_20261002")
const TRACE=joinpath(@__DIR__,"T02_BISECTION_TRACE.csv")
const OUT=joinpath(@__DIR__,"T02_FIXED_DESIGN_DELAY_MARGINS.csv")
BLAS.set_num_threads(1)

function classify(L,tau_ms)
    if tau_ms==0
        alpha=maximum(real,eigvals(L.A))
        return (;status=alpha<=-0.05 ? "SAFE" : "UNSAFE",count=count(z->real(z)>-0.05,eigvals(L.A)),
            phase_count=NaN,trace_count=NaN,min_sigma=NaN)
    end
    r=MegaOracle.trace_count(L,fill(tau_ms/1000,length(L.Ai));gamma=-0.05,step=20.0)
    agree=abs(real(r.estimate)-r.nearest)<0.01 &&
        abs(imag(r.estimate))<0.01 &&
        abs(r.phase_estimate-r.nearest)<0.01 &&
        r.largest_phase_step<pi/2 && r.quadrature_error_estimate<0.01
    (;status=agree ? (r.nearest==0 ? "SAFE" : "UNSAFE") : "INDETERMINATE",
      count=agree ? r.nearest : missing,phase_count=r.phase_estimate,
      trace_count=real(r.estimate),min_sigma=r.min_sigma_small)
end

function critical_candidates(L,model,design_id,tau_ms,locus)
    sample=locus[locus.design_id.==design_id,:]
    sample=sample[sample.status.=="LOCAL_EXACT_CHARACTERISTIC_CONTINUATION",:]
    isempty(sample) && return (NaN,NaN,0,NaN)
    roots=NamedTuple[]
    for group in groupby(sample,:family_id)
        j=argmin(abs.(group.tau_ms .-tau_ms))
        seed=group[j,:]
        r=FastRoots.refine(model,seed.root_real+im*seed.root_imag,
                           fill(tau_ms/1000,length(L.Ai)))
        r.converged && push!(roots,(;family_id=Int(seed.family_id),s=r.s,residual=r.sv))
    end
    isempty(roots) && return (NaN,NaN,0,NaN)
    best=roots[argmax(real.(getproperty.(roots,:s)))]
    (real(best.s),imag(best.s)/(2pi),best.family_id,best.residual)
end

function main()
    grid=CSV.read(joinpath(@__DIR__,"T01_DELAY_TRANSITION_REPRODUCTION.csv"),DataFrame)
    locus=CSV.read(joinpath(@__DIR__,"T08_FAST_MODE_PROVENANCE.csv"),DataFrame)
    ctx=MegaOracle.DC.R.N.design_context(ROOT)
    trace=NamedTuple[];margins=NamedTuple[]
    for (design_id,path) in (("seed_875",joinpath(MEGA,"seed_uniform_875.toml")),
                             ("best_zero_delay_88455",joinpath(MEGA,"M1_ZERO_DELAY_DESIGN.toml")))
        d=TOML.parsefile(path)
        L=MegaOracle.DC.linearization(ctx,Float64.(d["rho"]),Float64.(d["Kp"]),Float64.(d["Ki"]))
        model=FastRoots.reduced_model(L)
        g=sort(grid[grid.design_id.==design_id,:],:tau_ms)
        findfirst(g.status.=="INDETERMINATE")===nothing || error("T01 contains indeterminate count")
        firstbad=findfirst(g.status.=="UNSAFE")
        firstbad===nothing && error("no unsafe bracket for $(design_id)")
        firstbad>1 || error("first grid point already unsafe")
        lower=Float64(g.tau_ms[firstbad-1]);upper=Float64(g.tau_ms[firstbad])
        for tau_ms in (0.0,5.0,10.0,15.0)
            println("T02_EARLY_CHECK ",design_id," ",tau_ms);flush(stdout)
            q=classify(L,tau_ms)
            push!(trace,(;design_id,stage="EARLY_INTERVAL_CHECK",iteration=0,tau_ms,
                status=q.status,roots_violating_margin=q.count,phase_count=q.phase_count,
                trace_count=q.trace_count,min_sigma_small=q.min_sigma,
                candidate_rightmost_real=NaN,candidate_rightmost_frequency_hz=NaN,
                candidate_family_id=0,candidate_root_residual=NaN,
                current_safe_lower_ms=lower,current_unsafe_upper_ms=upper))
            CSV.write(TRACE,DataFrame(trace))
            q.status=="SAFE" || error("early-delay check not safe; cannot call grid boundary first transition")
        end
        iteration=0
        while upper-lower>0.1
            iteration+=1;tau_ms=(lower+upper)/2
            println("T02_BISECT ",design_id," iteration=",iteration," tau_ms=",tau_ms);flush(stdout)
            q=classify(L,tau_ms)
            rr,ff,fam,res=critical_candidates(L,model,design_id,tau_ms,locus)
            if q.status=="SAFE";lower=tau_ms
            elseif q.status=="UNSAFE";upper=tau_ms
            else;error("indeterminate contour at $(tau_ms) ms")
            end
            push!(trace,(;design_id,stage="BISECTION",iteration,tau_ms,
                status=q.status,roots_violating_margin=q.count,phase_count=q.phase_count,
                trace_count=q.trace_count,min_sigma_small=q.min_sigma,
                candidate_rightmost_real=rr,candidate_rightmost_frequency_hz=ff,
                candidate_family_id=fam,candidate_root_residual=res,
                current_safe_lower_ms=lower,current_unsafe_upper_ms=upper))
            CSV.write(TRACE,DataFrame(trace))
        end
        tau_ms=(lower+upper)/2
        rr,ff,fam,res=critical_candidates(L,model,design_id,tau_ms,locus)
        push!(margins,(;design_id,replacement_percent=design_id=="seed_875" ? 87.5 : Float64(d["GFL_percent"]),
            safe_lower_ms=lower,unsafe_upper_ms=upper,bracket_width_ms=upper-lower,
            tau_critical_estimate_ms=tau_ms,critical_root_real_candidate=rr,
            critical_frequency_candidate_hz=ff,critical_family_id=fam,
            candidate_root_residual=res,
            status="NUMERICAL_FIRST_OBSERVED_TRANSITION_NOT_GLOBAL_DELAY_CERTIFICATE"))
        CSV.write(OUT,DataFrame(margins))
        println("T02_MARGIN ",last(margins));flush(stdout)
    end
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
