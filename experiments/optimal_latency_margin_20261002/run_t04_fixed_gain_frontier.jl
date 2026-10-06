using CSV, DataFrames, LinearAlgebra, TOML

module MegaOracle
include(joinpath(@__DIR__,"..","analytical_delay_codesign_mega_20261002","m3_a_trace_integral.jl"))
end

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const OUT=joinpath(@__DIR__,"T04_FIXED_GAIN_RHO_TAU_FRONTIER.csv")
const TRACE=joinpath(@__DIR__,"T04_FIXED_GAIN_COUNT_TRACE.csv")
BLAS.set_num_threads(1)

function count_margin(L,tau_ms)
    if tau_ms==0
        vals=eigvals(L.A);n=count(z->real(z)>-0.05,vals)
        return (;status=n==0 ? "SAFE" : "UNSAFE",count=n,
            phase_count=NaN,trace_count=NaN,min_sigma_small=NaN)
    end
    r=MegaOracle.trace_count(L,fill(tau_ms/1000,length(L.Ai));gamma=-0.05,step=20.0)
    okay=abs(real(r.estimate)-r.nearest)<0.01 &&
        abs(imag(r.estimate))<0.01 &&
        abs(r.phase_estimate-r.nearest)<0.01 &&
        r.largest_phase_step<pi/2 && r.quadrature_error_estimate<0.01
    (;status=okay ? (r.nearest==0 ? "SAFE" : "UNSAFE") : "INDETERMINATE",
      count=okay ? r.nearest : missing,phase_count=r.phase_estimate,
      trace_count=real(r.estimate),min_sigma_small=r.min_sigma_small)
end

function main()
    frozen=joinpath(@__DIR__,"fixed_gain_designs")
    ctx=MegaOracle.DC.R.N.design_context(ROOT)
    output=NamedTuple[];logs=NamedTuple[]
    for path in sort(readdir(frozen;join=true))
        endswith(path,".toml") || continue
        d=TOML.parsefile(path)
        design_id=d["label"];percent=Float64(d["target_GFL_percent"])
        L=MegaOracle.DC.linearization(ctx,Float64.(d["rho"]),Float64.(d["Kp"]),Float64.(d["Ki"]))
        earlysafe=true
        for tau_ms in (0.0,10.0,20.0,30.0,40.0)
            println("T04_COARSE ",design_id," tau_ms=",tau_ms);flush(stdout)
            q=count_margin(L,tau_ms)
            push!(logs,(;design_id,replacement_percent=percent,stage="COARSE",iteration=0,
                tau_ms,status=q.status,roots_violating_margin=q.count,
                phase_count=q.phase_count,trace_count=q.trace_count,
                min_sigma_small=q.min_sigma_small))
            CSV.write(TRACE,DataFrame(logs))
            if tau_ms<=30 && q.status!="SAFE";earlysafe=false;end
        end
        recent=logs[end-4:end]
        if !earlysafe || recent[end].status!="UNSAFE"
            push!(output,(;design_id,replacement_percent=percent,
                safe_lower_ms=NaN,unsafe_upper_ms=NaN,bracket_width_ms=NaN,
                tau_critical_estimate_ms=NaN,upper_endpoint_roots_violating_margin=missing,
                status="INDETERMINATE_NO_30_TO_40_BRACKET"))
            CSV.write(OUT,DataFrame(output));continue
        end
        lower,upper=30.0,40.0;iteration=0;lastcount=recent[end].roots_violating_margin
        while upper-lower>0.1
            iteration+=1;tau_ms=(lower+upper)/2
            println("T04_BISECT ",design_id," k=",iteration," tau_ms=",tau_ms);flush(stdout)
            q=count_margin(L,tau_ms)
            push!(logs,(;design_id,replacement_percent=percent,stage="BISECTION",iteration,
                tau_ms,status=q.status,roots_violating_margin=q.count,
                phase_count=q.phase_count,trace_count=q.trace_count,
                min_sigma_small=q.min_sigma_small))
            CSV.write(TRACE,DataFrame(logs))
            if q.status=="SAFE";lower=tau_ms
            elseif q.status=="UNSAFE";upper=tau_ms;lastcount=q.count
            else;error("T04 indeterminate count for $(design_id) at $(tau_ms) ms")
            end
        end
        push!(output,(;design_id,replacement_percent=percent,
            safe_lower_ms=lower,unsafe_upper_ms=upper,bracket_width_ms=upper-lower,
            tau_critical_estimate_ms=(lower+upper)/2,
            upper_endpoint_roots_violating_margin=lastcount,
            status="NUMERICAL_FIRST_OBSERVED_MARGIN_CROSSING_FIXED_SEED_GAINS"))
        CSV.write(OUT,DataFrame(output))
        println("T04_RESULT ",last(output));flush(stdout)
    end
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
