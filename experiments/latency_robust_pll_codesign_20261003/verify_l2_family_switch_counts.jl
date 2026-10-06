using CSV, DataFrames, LinearAlgebra, TOML

module Oracle
include(joinpath(@__DIR__,"..","analytical_delay_codesign_mega_20261002","m3_a_trace_integral.jl"))
end

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
BLAS.set_num_threads(1)

function main()
    ctx=Oracle.DC.R.N.design_context(ROOT)
    z=TOML.parsefile(joinpath(@__DIR__,"designs","Z_zero_delay_tuned.toml"))
    n=TOML.parsefile(joinpath(@__DIR__,"designs","N_nominal.toml"))
    active=CSV.read(joinpath(@__DIR__,"L2_ACTIVE_FAMILY_SCAN.csv"),DataFrame)
    out=NamedTuple[]
    for eta in (0.5,0.8,0.9)
        kp=exp.((1-eta).*log.(Float64.(z["Kp"])).+eta.*log.(Float64.(n["Kp"])))
        ki=exp.((1-eta).*log.(Float64.(z["Ki"])).+eta.*log.(Float64.(n["Ki"])))
        L=Oracle.DC.linearization(ctx,Float64.(z["rho"]),kp,ki)
        row=only(eachrow(active[abs.(active.eta.-eta).<1e-8,:]))
        c=Float64(row.local_earliest_crossing_ms)
        for (side,offset) in (("below",-0.03),("above",0.03))
            tau_ms=c+offset
            println("L2_SWITCH_COUNT_START ",eta," ",side," ",tau_ms);flush(stdout)
            try
                r=Oracle.trace_count(L,fill(tau_ms/1000,10);gamma=-0.05,step=10.0)
                good=abs(real(r.estimate)-r.nearest)<0.01 && abs(imag(r.estimate))<0.01 &&
                    abs(r.phase_estimate-r.nearest)<0.01 && r.largest_phase_step<pi/2 &&
                    r.quadrature_error_estimate<0.01
                push!(out,(;eta,side,tau_ms,status=good ? (r.nearest==0 ? "SAFE" : "UNSAFE") : "INDETERMINATE",
                    root_count=good ? r.nearest : missing,trace_count_real=real(r.estimate),
                    phase_count=r.phase_estimate,min_sigma=r.min_sigma_small,
                    quadrature_error=r.quadrature_error_estimate,error=""))
            catch err
                push!(out,(;eta,side,tau_ms,status="INDETERMINATE",root_count=missing,
                    trace_count_real=NaN,phase_count=NaN,min_sigma=NaN,quadrature_error=NaN,
                    error=sprint(showerror,err)))
            end
            CSV.write(joinpath(@__DIR__,"L2_FAMILY_SWITCH_COUNTS.csv"),DataFrame(out))
            println("L2_SWITCH_COUNT_DONE ",last(out));flush(stdout)
        end
    end
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
