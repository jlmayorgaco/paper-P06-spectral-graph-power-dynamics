using CSV, DataFrames, LinearAlgebra, TOML

module MegaOracle
include(joinpath(@__DIR__,"..","analytical_delay_codesign_mega_20261002","m3_a_trace_integral.jl"))
end

const ROOT = normpath(joinpath(@__DIR__,"..",".."))
const MEGA = joinpath(ROOT,"experiments","analytical_delay_codesign_mega_20261002")
const GRID_MS = (20.0,25.0,30.0,32.0,34.0,36.0,38.0,40.0)
const OUT = joinpath(@__DIR__,"T01_DELAY_TRANSITION_REPRODUCTION.csv")
BLAS.set_num_threads(1)

function count_one(L,tau_ms)
    tau=fill(tau_ms/1000,length(L.Ai))
    r=MegaOracle.trace_count(L,tau;gamma=-0.05,step=20.0)
    agreed=abs(real(r.estimate)-r.nearest)<0.01 &&
        abs(imag(r.estimate))<0.01 &&
        abs(r.phase_estimate-r.nearest)<0.01 &&
        r.largest_phase_step<pi/2 && r.quadrature_error_estimate<0.01
    (;status=agreed ? (r.nearest==0 ? "SAFE" : "UNSAFE") : "INDETERMINATE",
      roots_violating_margin=agreed ? r.nearest : missing,
      trace_count_real=real(r.estimate),phase_count=r.phase_estimate,
      largest_phase_step_rad=r.largest_phase_step,
      integral_error_estimate=r.quadrature_error_estimate,
      min_sigma_small=r.min_sigma_small,
      contour_radius_rad_s=r.radius,
      contour_evaluations=r.evaluations)
end

function main()
    ctx=MegaOracle.DC.R.N.design_context(ROOT)
    rows=NamedTuple[]
    designs=(("seed_875",joinpath(MEGA,"seed_uniform_875.toml")),
             ("best_zero_delay_88455",joinpath(MEGA,"M1_ZERO_DELAY_DESIGN.toml")))
    for (design_id,path) in designs
        d=TOML.parsefile(path)
        rho,kp,ki=Float64.(d["rho"]),Float64.(d["Kp"]),Float64.(d["Ki"])
        L=MegaOracle.DC.linearization(ctx,rho,kp,ki)
        replacement_percent=design_id=="seed_875" ? 87.5 : Float64(d["GFL_percent"])
        for tau_ms in GRID_MS
            println("T01_START ",design_id," tau_ms=",tau_ms);flush(stdout)
            try
                q=count_one(L,tau_ms)
                push!(rows,merge((;design_id,replacement_percent,tau_ms,
                    margin_real_s_inv=-0.05,source="NEW_TARGETED_FLOAT64_CONTOUR"),q))
            catch err
                push!(rows,(;design_id,replacement_percent,tau_ms,margin_real_s_inv=-0.05,
                    source="NEW_TARGETED_FLOAT64_CONTOUR",status="INDETERMINATE",
                    roots_violating_margin=missing,trace_count_real=NaN,phase_count=NaN,
                    largest_phase_step_rad=NaN,integral_error_estimate=NaN,
                    min_sigma_small=NaN,contour_radius_rad_s=NaN,contour_evaluations=missing))
                println("T01_ERROR ",sprint(showerror,err));flush(stdout)
            end
            CSV.write(OUT,DataFrame(rows))
            println("T01_DONE ",last(rows));flush(stdout)
        end
    end
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
