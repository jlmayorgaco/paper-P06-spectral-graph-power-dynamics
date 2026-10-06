using CSV, DataFrames, LinearAlgebra, TOML

module FastRoots
include(joinpath(@__DIR__,"run_fast_root_locus.jl"))
end

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
BLAS.set_num_threads(1)

function root_at(model,seedtable,family_id,tau_ms)
    rows=seedtable[seedtable.family_id.==family_id,:]
    j=argmin(abs.(rows.tau_ms .-tau_ms))
    z=rows[j,:]
    FastRoots.refine(model,z.root_real+im*z.root_imag,
        fill(tau_ms/1000,size(model.right,2));tol=1e-11)
end

function main()
    frontier=CSV.read(joinpath(@__DIR__,"T04_FIXED_GAIN_RHO_TAU_FRONTIER.csv"),DataFrame)
    locus=CSV.read(joinpath(@__DIR__,"T08_FAST_MODE_PROVENANCE.csv"),DataFrame)
    seedtable=locus[locus.design_id.=="seed_875",:]
    ctx=FastRoots.MegaOracle.DC.R.N.design_context(ROOT)
    output=NamedTuple[]
    for item in eachrow(frontier)
        path=joinpath(@__DIR__,"fixed_gain_designs",String(item.design_id)*".toml")
        d=TOML.parsefile(path)
        L=FastRoots.MegaOracle.DC.linearization(ctx,Float64.(d["rho"]),Float64.(d["Kp"]),Float64.(d["Ki"]))
        model=FastRoots.reduced_model(L)
        candidates=NamedTuple[]
        for family_id in sort(unique(seedtable.family_id))
            lo=Float64(item.safe_lower_ms);hi=Float64(item.unsafe_upper_ms)
            rlo=root_at(model,seedtable,family_id,lo)
            rhi=root_at(model,seedtable,family_id,hi)
            (rlo.converged && rhi.converged && real(rlo.s)<-0.05<real(rhi.s)) || continue
            iterations=0
            while hi-lo>1e-4
                iterations+=1;mid=(lo+hi)/2
                r=root_at(model,seedtable,family_id,mid)
                r.converged || error("fixed-gain branch failed at $mid")
                if real(r.s)<-0.05;lo=mid;else;hi=mid;end
            end
            t=(lo+hi)/2;r=root_at(model,seedtable,family_id,t)
            push!(candidates,(;family_id,critical_delay_ms=t,lower_ms=lo,upper_ms=hi,
                critical_root_real=real(r.s),critical_root_imag=imag(r.s),
                critical_frequency_hz=imag(r.s)/(2pi),root_residual=r.sv,
                d_real_lambda_d_uniform_delay=real(r.root_tau_sensitivity),iterations))
        end
        if isempty(candidates)
            push!(output,(;design_id=String(item.design_id),replacement_percent=item.replacement_percent,
                fixed_gain_tau_crit_ms=NaN,critical_family_id=missing,
                critical_frequency_hz=NaN,critical_root_real=NaN,root_residual=NaN,
                contour_safe_lower_ms=item.safe_lower_ms,contour_unsafe_upper_ms=item.unsafe_upper_ms,
                status="NO_TRACKED_BRANCH_BRACKETS_CONTOUR"))
        else
            c=candidates[argmin(getproperty.(candidates,:critical_delay_ms))]
            push!(output,(;design_id=String(item.design_id),replacement_percent=item.replacement_percent,
                fixed_gain_tau_crit_ms=c.critical_delay_ms,critical_family_id=c.family_id,
                critical_frequency_hz=c.critical_frequency_hz,critical_root_real=c.critical_root_real,
                root_residual=c.root_residual,
                contour_safe_lower_ms=item.safe_lower_ms,contour_unsafe_upper_ms=item.unsafe_upper_ms,
                status="LOCAL_ROOT_CROSSING_WITH_CONTOUR_BRACKET"))
        end
        CSV.write(joinpath(@__DIR__,"T04_PRECISE_CROSSINGS.csv"),DataFrame(output))
        println("T04_PRECISE ",last(output));flush(stdout)
    end
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
