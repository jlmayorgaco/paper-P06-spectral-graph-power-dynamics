using CSV, DataFrames, LinearAlgebra, TOML

module Oracle
include(joinpath(@__DIR__, "..", "analytical_delay_codesign_mega_20261002", "m3_a_trace_integral.jl"))
end

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const OUT = joinpath(@__DIR__, "T04_ENDPOINT_BOUNDARY_VERIFICATION_0P01MS.csv")
BLAS.set_num_threads(1)

function main()
    ctx = Oracle.DC.R.N.design_context(ROOT)
    crossings = CSV.read(joinpath(@__DIR__, "T04_PRECISE_CROSSINGS.csv"), DataFrame)
    rows = NamedTuple[]
    for j in (1, nrow(crossings))
        c = crossings[j, :]
        d = TOML.parsefile(joinpath(@__DIR__, "fixed_gain_designs", "rho_" * lpad(j == 1 ? "0" : "4", 2, '0') * ".toml"))
        L = Oracle.DC.linearization(ctx, Float64.(d["rho"]), Float64.(d["Kp"]), Float64.(d["Ki"]))
        for (side, offset) in (("below", -0.01), ("above", 0.01))
            tau_ms = c.fixed_gain_tau_crit_ms + offset
            println("T04_ENDPOINT_START ", c.design_id, " ", side, " ", tau_ms); flush(stdout)
            try
                r = Oracle.trace_count(L, fill(tau_ms / 1000, length(L.Ai)); gamma=-0.05, step=10.0)
                agreed = abs(real(r.estimate)-r.nearest) < 0.01 &&
                    abs(imag(r.estimate)) < 0.01 &&
                    abs(r.phase_estimate-r.nearest) < 0.01 &&
                    r.largest_phase_step < pi/2 && r.quadrature_error_estimate < 0.01
                status = agreed ? (r.nearest == 0 ? "SAFE" : "UNSAFE") : "INDETERMINATE"
                push!(rows, (; design_id=c.design_id, replacement_percent=c.replacement_percent,
                    side, tau_ms, status, roots_violating_margin=agreed ? r.nearest : missing,
                    trace_count_real=real(r.estimate), phase_count=r.phase_estimate,
                    min_sigma_small=r.min_sigma_small, max_phase_step=r.largest_phase_step,
                    quadrature_error_estimate=r.quadrature_error_estimate,
                    contour_evaluations=r.evaluations, error=""))
            catch err
                push!(rows, (; design_id=c.design_id, replacement_percent=c.replacement_percent,
                    side, tau_ms, status="INDETERMINATE", roots_violating_margin=missing,
                    trace_count_real=NaN, phase_count=NaN, min_sigma_small=NaN,
                    max_phase_step=NaN, quadrature_error_estimate=NaN,
                    contour_evaluations=missing, error=sprint(showerror, err)))
            end
            CSV.write(OUT, DataFrame(rows))
            println("T04_ENDPOINT_DONE ", last(rows)); flush(stdout)
        end
    end
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
