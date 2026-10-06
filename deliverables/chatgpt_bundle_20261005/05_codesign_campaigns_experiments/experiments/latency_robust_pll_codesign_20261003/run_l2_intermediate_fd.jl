using CSV, DataFrames, LinearAlgebra, TOML

module FD
include(joinpath(@__DIR__,"run_l2_gradient_fd.jl"))
end

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
BLAS.set_num_threads(1)

function main()
    ctx=FD.Roots.MegaOracle.DC.R.N.design_context(ROOT)
    z=TOML.parsefile(joinpath(@__DIR__,"designs","Z_zero_delay_tuned.toml"))
    n=TOML.parsefile(joinpath(@__DIR__,"designs","N_nominal.toml"))
    rho=Float64.(z["rho"])
    scan=CSV.read(joinpath(@__DIR__,"L2_FAMILY_SCAN.csv"),DataFrame)
    gradients=CSV.read(joinpath(@__DIR__,"L2_MULTIMODE_GRADIENTS.csv"),DataFrame)
    tests=Dict(0.5=>[(:Kp,35),(:Kp,36),(:Ki,35),(:Ki,36)],
               0.9=>[(:Kp,30),(:Ki,30),(:Kp,35),(:Kp,36)])
    out=NamedTuple[]
    for eta in (0.5,0.9)
        kp=exp.((1-eta).*log.(Float64.(z["Kp"])).+eta.*log.(Float64.(n["Kp"])))
        ki=exp.((1-eta).*log.(Float64.(z["Ki"])).+eta.*log.(Float64.(n["Ki"])))
        branches=scan[(abs.(scan.eta.-eta).<1e-8).&(scan.rank_by_local_crossing.<=2),:]
        sort!(branches,:rank_by_local_crossing)
        nrow(branches)==2 || error("need two branches at eta=$eta")
        for (kind,bus) in tests[eta]
            i=bus-29;k=kind==:Kp ? kp[i] : ki[i];h=0.001*k
            crossings=Dict{String,Vector{Float64}}()
            for (side,sign) in (("plus",1),("minus",-1))
                kp1=copy(kp);ki1=copy(ki)
                if kind==:Kp;kp1[i]+=sign*h;else;ki1[i]+=sign*h;end
                branch_crossings=Float64[]
                for b in eachrow(branches)
                    s=ComplexF64(b.critical_real+2pi*im*b.frequency_hz)
                    r=FD.root_crossing(ctx,rho,kp1,ki1,s,Float64(b.local_crossing_ms))
                    push!(branch_crossings,r===nothing ? NaN : r.tau_ms)
                end
                crossings[side]=branch_crossings
            end
            gp=only(eachrow(gradients[(abs.(gradients.eta.-eta).<1e-8).&
                (gradients.mode_rank.==1).&(gradients.bus.==bus),:]))
            analytic=kind==:Kp ? gp.d_margin_ms_d_Kp : gp.d_margin_ms_d_Ki
            p,m=crossings["plus"],crossings["minus"]
            valid=all(isfinite,vcat(p,m))
            fdbranch=valid ? (p[1]-m[1])/(2h) : NaN
            fdenvelope=valid ? (minimum(p)-minimum(m))/(2h) : NaN
            relbranch=valid ? abs(fdbranch-analytic)/max(abs(analytic),abs(fdbranch),1e-12) : NaN
            relenv=valid ? abs(fdenvelope-analytic)/max(abs(analytic),abs(fdenvelope),1e-12) : NaN
            switched=valid && (argmin(p)!=argmin(m))
            status=!valid ? "INDETERMINATE_ROOT_REFINEMENT" : switched ?
                "NONSMOOTH_ENVELOPE_MODE_SWITCH" : relenv<0.05 ?
                "LOCAL_ENVELOPE_GRADIENT_VALIDATED" : "LOCAL_ENVELOPE_GRADIENT_FAIL"
            push!(out,(;eta,bus,gain_kind=String(kind),gain_value=k,perturbation=h,
                analytic_active_branch_gradient_ms_per_gain=analytic,
                centered_FD_same_branch_gradient_ms_per_gain=fdbranch,
                centered_FD_two_mode_envelope_gradient_ms_per_gain=fdenvelope,
                branch_relative_error=relbranch,envelope_relative_error=relenv,
                plus_first_mode_crossing_ms=p[1],plus_second_mode_crossing_ms=p[2],
                minus_first_mode_crossing_ms=m[1],minus_second_mode_crossing_ms=m[2],
                plus_active_mode=argmin(p),minus_active_mode=argmin(m),
                mode_switched_across_FD=switched,status,
                scope="TWO_CATALOGUED_FAST_MODES;FULL_CONTOUR_REQUIRED_FOR_GLOBAL_ACTIVE_LABEL"))
            CSV.write(joinpath(@__DIR__,"L2_INTERMEDIATE_GRADIENT_VALIDATION.csv"),DataFrame(out))
            println("L2_INTERMEDIATE_FD ",last(out));flush(stdout)
        end
    end
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
