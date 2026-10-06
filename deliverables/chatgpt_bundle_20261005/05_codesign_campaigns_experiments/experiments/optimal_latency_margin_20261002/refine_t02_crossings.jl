using CSV, DataFrames, LinearAlgebra, TOML

module FastRoots
include(joinpath(@__DIR__,"run_fast_root_locus.jl"))
end

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const MEGA=joinpath(ROOT,"experiments","analytical_delay_codesign_mega_20261002")
BLAS.set_num_threads(1)

function at_delay(model,series,tau_ms)
    j=argmin(abs.(series.tau_ms .-tau_ms))
    z=series[j,:]
    r=FastRoots.refine(model,z.root_real+im*z.root_imag,
        fill(tau_ms/1000,size(model.right,2));tol=1e-11)
    r.converged || error("local root refinement failed at $tau_ms ms")
    r
end

function main()
    margins=CSV.read(joinpath(@__DIR__,"T02_FIXED_DESIGN_DELAY_MARGINS.csv"),DataFrame)
    locus=CSV.read(joinpath(@__DIR__,"T08_FAST_MODE_PROVENANCE.csv"),DataFrame)
    ctx=FastRoots.MegaOracle.DC.R.N.design_context(ROOT)
    rows=NamedTuple[]
    for row in eachrow(margins)
        path=row.design_id=="seed_875" ? joinpath(MEGA,"seed_uniform_875.toml") :
             joinpath(MEGA,"M1_ZERO_DELAY_DESIGN.toml")
        d=TOML.parsefile(path)
        L=FastRoots.MegaOracle.DC.linearization(ctx,Float64.(d["rho"]),Float64.(d["Kp"]),Float64.(d["Ki"]))
        model=FastRoots.reduced_model(L)
        series=locus[(locus.design_id.==row.design_id).&(locus.family_id.==row.critical_family_id),:]
        lo=Float64(row.safe_lower_ms);hi=Float64(row.unsafe_upper_ms)
        rlo=at_delay(model,series,lo);rhi=at_delay(model,series,hi)
        real(rlo.s)<-0.05<real(rhi.s) || error("active branch does not bracket margin for $(row.design_id)")
        iterations=0
        while hi-lo>1e-4
            iterations+=1
            mid=(lo+hi)/2
            r=at_delay(model,series,mid)
            if real(r.s)<-0.05;lo=mid;else;hi=mid;end
        end
        center=(lo+hi)/2;r=at_delay(model,series,center)
        push!(rows,(;design_id=row.design_id,replacement_percent=row.replacement_percent,
            original_contour_safe_lower_ms=row.safe_lower_ms,
            original_contour_unsafe_upper_ms=row.unsafe_upper_ms,
            local_root_safe_lower_ms=lo,local_root_unsafe_upper_ms=hi,
            local_root_crossing_ms=center,local_bracket_width_ms=hi-lo,
            critical_root_real=real(r.s),critical_root_imag=imag(r.s),
            critical_frequency_hz=imag(r.s)/(2pi),critical_family_id=row.critical_family_id,
            small_factor_root_residual=r.sv,
            d_real_lambda_d_uniform_delay_s_inv_per_s=real(r.root_tau_sensitivity),
            iterations,status="LOCAL_ROOT_CROSSING_WITH_CONTOUR_BRACKET"))
        CSV.write(joinpath(@__DIR__,"T02_PRECISE_CROSSINGS.csv"),DataFrame(rows))
        println("PRECISE_CROSSING ",last(rows));flush(stdout)
    end
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
