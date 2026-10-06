"""Freeze the saturated-rho analytic candidates before any new PowerDynamics call.

The admissible replacement fraction is rho in [0,1]. If the exact frozen
full-state matrix at rho=1 and nominal PLL settings meets every non-gauge
spectral constraint, the objective upper bound is attained. This gives a
global certificate without solving an unnecessary pole-boundary polynomial.

This script reads only frozen ExpC matrices and the copied IEEE-39 data. It
does not import or call PowerDynamics.
"""

using CSV, DataFrames, LinearAlgebra, SHA, Dates

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const REPORT = joinpath(ROOT, "reports", "experiment_D")
const C_MATRICES = joinpath(ROOT, "reports", "experiment_C", "matrices")
const C_TABLES = joinpath(ROOT, "reports", "experiment_C", "tables")
const SIGMA_REQ = 0.05
const PLL_SCALE_BOUNDS = (0.9, 1.1)
const KP0 = 5.0 * 2pi
const KI0 = KP0^2 / 4
const BUS_MATRIX = Dict(30 => "C30_Ared.csv", 33 => "C33_regenerated_Ared.csv",
                        35 => "C35_Ared.csv", 37 => "C37_Ared.csv")

file_sha(path) = bytes2hex(sha256(read(path)))
is_gauge(z) = abs(z) <= 1e-8

function main()
    busdata = CSV.read(joinpath(REPORT, "inputs", "bus.csv"), DataFrame)
    tables = NamedTuple[]
    candidates = NamedTuple[]
    input_hashes = Dict{String,String}()
    all_pass = true

    for bus in (33, 30, 35, 37)
        matrix_name = BUS_MATRIX[bus]
        matrix_path = joinpath(C_MATRICES, matrix_name)
        data_path = joinpath(REPORT, "inputs", "bus.csv")
        input_hashes["reports/experiment_C/matrices/$(matrix_name)"] = file_sha(matrix_path)
        matrix_df = CSV.read(matrix_path, DataFrame)
        A = Matrix{Float64}(matrix_df[:, Not(:row_index)])
        size(A, 1) == size(A, 2) || error("Ared for bus $bus is not square")
        λ = ComplexF64.(eigvals(A))
        all(isfinite, real.(λ)) && all(isfinite, imag.(λ)) || error("nonfinite spectrum at bus $bus")
        nongauge = filter(z -> !is_gauge(z), λ)
        isempty(nongauge) && error("no non-gauge poles at bus $bus")
        rightmost = nongauge[argmax(real.(nongauge))]
        alpha = real(rightmost)
        mode_ok = alpha <= -SIGMA_REQ + 1e-10
        all_pass &= mode_ok

        busrow = busdata[findfirst(==(bus), Int.(busdata.bus)), :]
        p0_mw = 100.0 * Float64(busrow.P)
        mw = p0_mw
        kprow = 1.0 * KP0
        kirow = 1.0 * KI0
        beta = 1.0
        first(PLL_SCALE_BOUNDS) <= beta <= last(PLL_SCALE_BOUNDS) ||
            error("nominal PLL scale is outside its domain")
        push!(candidates, (
            bus=bus, rho_star=1.0, P0_MW=p0_mw, MW_star=mw,
            Kp_star=kprow, Ki_star=kirow, pll_scale=beta,
            rightmost_pole_real=real(rightmost), rightmost_pole_imag=imag(rightmost),
            alpha=alpha, required_margin=SIGMA_REQ,
            non_gauge_mode_count=length(nongauge), total_finite_mode_count=length(λ),
            all_modes_pass=mode_ok, active_constraint="rho_upper_bound",
            matrix_source="reports/experiment_C/matrices/$(matrix_name)"))

        for (idx, z) in enumerate(λ)
            gauge = is_gauge(z)
            pass = gauge || real(z) <= -SIGMA_REQ + 1e-10
            push!(tables, (bus=bus, mode=idx, lambda_real=real(z), lambda_imag=imag(z),
                damping_ratio=(abs(z) > 1e-8 && abs(imag(z)) > 1e-12 ? -real(z)/abs(z) : missing),
                is_gauge=gauge, pass=pass, required_alpha=-SIGMA_REQ,
                source="reports/experiment_C/matrices/$(matrix_name)"))
        end
    end

    all_pass || error("at least one rho=1 endpoint fails the declared all-mode margin")
    input_hashes["reports/experiment_D/inputs/bus.csv"] = file_sha(joinpath(REPORT, "inputs", "bus.csv"))
    input_hashes["reports/experiment_C/RESULTS_EXP_C.json"] = file_sha(joinpath(ROOT, "reports", "experiment_C", "RESULTS_EXP_C.json"))
    input_hashes["experiments/bnd_expD/freeze_endpoint_candidates.jl"] = file_sha(@__FILE__)

    primary = only(filter(c -> c.bus == 33, candidates))
    cert = (
        type="SATURATED_PHYSICAL_BOUND",
        objective="P0_MW * rho",
        rho_interval=[0.0, 1.0],
        upper_bound_MW=primary.P0_MW,
        achieved_MW=primary.MW_star,
        optimality_gap_MW=0.0,
        witness_rho=primary.rho_star,
        endpoint_all_modes_pass=primary.all_modes_pass,
        argument="For every admissible point, P0_MW*rho <= P0_MW. The exact rho=1 endpoint is feasible, so this upper bound is attained globally; no interior or pole-boundary branch can improve it.",
        gain_tie_break="minimize squared PLL-scale distance to nominal; beta=1 is the unique zero-distance choice.",
        endpoint_structure="At rho=1 the SG has zero share and is removed from the physical realization; its dormant local states are not counted as network poles." )

    candidate = (
        experiment="BND_EXP_D",
        candidate_frozen=true,
        status="FROZEN",
        frozen_at_utc=Dates.format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SS.sssZ"),
        primary_bus=33,
        objective="maximize active SG dispatch replaced by GFL, P0_MW*rho",
        domain=(rho=[0.0, 1.0], pll_scale=collect(PLL_SCALE_BOUNDS),
            PLL_Kp="beta * 5*2pi rad/s", PLL_Ki="(beta*5*2pi)^2/4 rad/s^2",
            sigma_req=SIGMA_REQ, gauge_tolerance=1e-8),
        primary_candidate=primary,
        cross_bus_candidates=candidates,
        global_certificate=cert,
        analysis_method=(name="exact frozen full-state endpoint spectrum plus saturated-bound proof",
            powerdynamics_calls=0, optimizer_calls=0,
            limitation="The proof certifies the rho objective by its physical upper bound. It does not enumerate nonbinding pole-boundary branches because the upper bound is attained."),
        input_sha256=input_hashes)

    outpath = joinpath(REPORT, "ANALYTIC_CANDIDATE.json")
    open(outpath, "w") do io
        write(io, BNDDesign.AnalyticDeviceModel.json_string(candidate))
    end
    open(joinpath(REPORT, "ANALYTIC_CANDIDATE.sha256"), "w") do io
        write(io, file_sha(outpath), "\n")
    end
    CSV.write(joinpath(REPORT, "tables", "TABLE_D10_multimode_check.csv"), DataFrame(tables))
    CSV.write(joinpath(REPORT, "tables", "TABLE_D09_analytic_optimum.csv"),
        DataFrame([merge(c, (omega_star=abs(c.rightmost_pole_imag),
            boundary_real=missing, boundary_imag=missing)) for c in candidates]))

    println("ANALYTIC_CANDIDATE_FROZEN: YES")
    println("PRIMARY_BUS=33; rho_star=1.0; MW_star=$(primary.MW_star)")
    println("Kp_star=$(primary.Kp_star); Ki_star=$(primary.Ki_star)")
    println("rightmost_pole=$(primary.rightmost_pole_real) + $(primary.rightmost_pole_imag)im")
    println("all_modes_pass=$(primary.all_modes_pass); required_margin=$(SIGMA_REQ)")
    println("candidate_sha256=$(file_sha(outpath))")
    println("POWERDYNAMICS_IMPORTED: NO")
end

include(joinpath(ROOT, "src", "bnd_design", "BNDDesign.jl"))
main()
