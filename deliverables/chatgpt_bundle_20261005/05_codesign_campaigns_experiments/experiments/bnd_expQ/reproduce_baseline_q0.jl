using CSV, DataFrames, LinearAlgebra, SHA, TOML
using NetworkDynamics, PowerDynamics, OrdinaryDiffEqRosenbrock, SciMLBase

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const OUT = joinpath(ROOT, "reports", "experiment_Q", "Q0")
mkpath(OUT)

include(joinpath(ROOT, "src", "bnd_design_p", "ExpP.jl"))
using .ExpP
include(joinpath(ROOT, "src", "bnd_design_p", "AlgebraicFrontiers.jl"))
using .AlgebraicFrontiers
include(joinpath(ROOT, "src", "bnd_model_expN", "PDReferenceN.jl"))
include(joinpath(ROOT, "src", "bnd_model_expN", "PDExactDesignN.jl"))
const PD39 = PDReferenceN.PD39

function event_network(nw, bus, deltaP_MW; start=1.0, base_mva=100.0)
    vertices, edges = PD39.PD39Model.copy_network_components(nw)
    defaults = get_defaults_dict(vertices[bus])
    pkeys = [s for s in keys(defaults) if occursin("Pset", string(s))]
    qkeys = [s for s in keys(defaults) if occursin("Qset", string(s))]
    length(pkeys) == 1 && length(qkeys) == 1 ||
        error("load P/Q setpoints absent at bus $bus")
    pkey, qkey = only(pkeys), only(qkeys)
    p0, q0 = defaults[pkey], defaults[qkey]
    delta_p_pu = deltaP_MW / base_mva
    affect = (u, p, ctx) -> begin
        ctx.t >= start && (p[pkey] = p0 - delta_p_pu)
        p[qkey] = q0
    end
    callback = PresetTimeComponentCallback([start],
        ComponentAffect(affect, (), (pkey, qkey)))
    set_callback!(vertices[bus], callback)
    net = Network(vertices, edges)
    set_jac_prototype!(net)
    net, (; p0, q0, delta_p_pu)
end

function main()
    nroot = joinpath(ROOT, "reports", "experiment_N")
    proot = joinpath(ROOT, "reports", "experiment_P")
    toml = joinpath(nroot, "Z_N_NOMINAL_FINAL.toml")
    p5toml = joinpath(proot, "P5", "Z_P_NOMINAL_FINAL.toml")
    nsha = bytes2hex(sha256(read(toml)))
    psha = bytes2hex(sha256(read(p5toml)))
    nsha == strip(read(toml * ".sha256", String)) || error("ExpN freeze mismatch")
    psha == strip(read(p5toml * ".sha256", String)) || error("ExpP candidate freeze mismatch")
    nsha == "3915a8f57da0f552779b11464572f14a86f0b648b532be77124f668eeff9f263" ||
        error("unexpected ExpN candidate identity")
    psha == "dbe727713acccec9b0802619f4b2e10763a0a88b9799ba24a74411ac88ff7046" ||
        error("unexpected ExpP candidate identity")

    c = TOML.parsefile(toml)
    rho, kp, ki = Float64.(c["rho"]), Float64.(c["Kp"]), Float64.(c["Ki"])
    epsv = 1 .- rho
    # The ExpP algebraic identity accepts its own module-scoped context type.
    ctx = ExpP.PDExactDesignN.design_context(ROOT)
    an = ExpP.PDExactDesignN.spectrum(ctx, rho, kp, ki)
    retained = dot(ctx.power, epsv)
    abs(retained - 1.124344477653483) < 1e-8 || error("ExpN retained MW changed")
    abs(an.alpha - (-0.0500000011050941)) < 1e-8 || error("ExpN analytic alpha changed")

    # Independently rebuild the mixed PowerDynamics architecture and equilibrium.
    base = PDReferenceN.frozen_baseline()
    nw = PDReferenceN.build_architecture(base, rho, kp, ki)
    state = PDReferenceN.trim_state(nw, base, rho, kp, ki)
    trim = PDReferenceN.residual_audit(nw, state)
    power = PDReferenceN.direct_power_audit(state, base, rho)
    max_p = maximum(Float64.(power.max_P_error_pu))
    max_q = maximum(Float64.(power.max_Q_error_pu))
    pdlin = linearize_network(state)
    pdpoles = jacobian_eigenvals(pdlin)
    gauge_idx = argmin(abs.(pdpoles))
    pdalpha = maximum(real.(pdpoles[[i for i in eachindex(pdpoles) if i != gauge_idx]]))
    abs(pdalpha - (-0.050000000138500054)) < 1e-7 ||
        error("ExpN PD alpha did not reproduce")

    # Recompute the ExpP bus-38 algebraic root from the exact retained-SG update.
    roots = AlgebraicFrontiers.retention_roots(ctx, kp, ki, 38)
    candidate_eps = epsv[38 - 29]
    rootrow = roots.roots[argmin(abs.(roots.roots.epsilon .- candidate_eps)), :]
    root_mw = Float64(rootrow.retained_SG_MW)
    abs(root_mw - retained) < 1e-6 || error("ExpP conditional root did not recover ExpN")

    # Re-evaluate the pointwise robustness witness stored in ExpP from its frozen CSV.
    robust = CSV.read(joinpath(proot, "P4", "TABLE_P4_resolvent_observations.csv"), DataFrame)
    nominal = robust[robust.candidate .== "ExpN nominal", :]
    nrow(nominal) == 1 || error("nominal ExpP beta witness missing")
    beta_upper = Float64(only(nominal.beta_upper_from_observation))
    beta_req = Float64(only(nominal.beta_req))
    beta_upper < beta_req || error("ExpN pointwise beta witness no longer rejects incumbent")

    # Verify separate SG/GFL component dispatch, proportional Q share, and ZIP loads.
    contract = power[!, [:bus, :rho, :SG_P_MW, :SG_Q_Mvar, :GFL_P_MW, :GFL_Q_Mvar,
                         :ZIP_P_MW, :ZIP_Q_Mvar, :target_SG_P_MW, :target_SG_Q_Mvar,
                         :target_GFL_P_MW, :target_GFL_Q_Mvar, :max_P_error_pu,
                         :max_Q_error_pu]]
    CSV.write(joinpath(OUT, "TABLE_Q00_component_contract.csv"), contract)
    max_contract_error = max(maximum(contract.max_P_error_pu), maximum(contract.max_Q_error_pu))
    for bus in (31, 39)
        row = only(eachrow(contract[contract.bus .== bus, :]))
        base_row = base.rows[bus]
        expected_load_p = 100.0 * Float64(state[VIndex(bus, :ZIPLoad₊P)])
        expected_load_q = 100.0 * Float64(state[VIndex(bus, :ZIPLoad₊Q)])
        abs(row.ZIP_P_MW - expected_load_p) < 1e-10 || error("bus $bus load P was conflated")
        abs(row.ZIP_Q_Mvar - expected_load_q) < 1e-10 || error("bus $bus load Q was conflated")
        abs(row.target_SG_P_MW + row.target_GFL_P_MW - base_row.P) < 1e-9 ||
            error("bus $bus P share contract changed")
        abs(row.target_SG_Q_Mvar + row.target_GFL_Q_Mvar - base_row.Q) < 1e-9 ||
            error("bus $bus Q share contract changed")
    end

    # Independently reproduce ExpP's declared 100 MW sustained bus-16 event,
    # writing all new evidence only below reports/experiment_Q.
    scenario = TOML.parsefile(joinpath(ROOT, "experiments", "bnd_expG",
        "configs", "DESIGN_SCENARIO_FROZEN.toml"))
    bus = Int(scenario["event_bus"])
    deltaP = Float64(scenario["deltaP_MW"])
    active, profile = event_network(nw, bus, deltaP;
        base_mva=Float64(scenario["base_MVA"]))
    dt, horizon = 0.01, 60.0
    prob = SciMLBase.ODEProblem(active, state, (0.0, horizon))
    sol = SciMLBase.solve(prob, OrdinaryDiffEqRosenbrock.Rodas5P();
        callback=get_callbacks(active), initializealg=SciMLBase.NoInit(),
        saveat=dt, abstol=1e-9, reltol=1e-9)
    SciMLBase.successful_retcode(sol.retcode) || error("Q0 baseline TDS failed: $(sol.retcode)")
    times = collect(0.0:dt:horizon)
    omega0 = Float64(state[VIndex(38, :ctrld_gen₊machine₊ω)])
    frames = NamedTuple[]
    freq = Float64[]
    for t in times
        s = NetworkDynamics.NWState(sol, t)
        f = 60 * (Float64(s[VIndex(38, :ctrld_gen₊machine₊ω)]) - omega0)
        push!(freq, f)
        push!(frames, (; time_s=t, frequency_deviation_Hz=f))
    end
    rocof = vcat(0.0, diff(freq) ./ dt)
    after = findall(times .>= 1.0)
    peakf = maximum(abs.(freq[after]))
    peakr = maximum(abs.(rocof[after]))
    threshold = max(0.02 * peakf, 1e-9)
    settling = NaN
    for j in after
        if maximum(abs.(freq[j:end])) <= threshold
            settling = times[j] - 1.0
            break
        end
    end
    CSV.write(joinpath(OUT, "TABLE_Q00_baseline_event_trajectory.csv"),
              DataFrame(time_s=times, frequency_deviation_Hz=freq, RoCoF_Hz_s=rocof))
    p5event = CSV.read(joinpath(proot, "P5", "TABLE_P5_declared_event_validation.csv"), DataFrame)
    abs(peakf - only(p5event.peak_COI_frequency_Hz)) < 0.03 ||
        error("Q0 reproduced frequency event differs from ExpP")
    abs(peakr - only(p5event.peak_RoCoF_Hz_s)) < 0.01 ||
        error("Q0 reproduced RoCoF event differs from ExpP")
    isfinite(settling) == isfinite(only(p5event.settling_time_s)) ||
        error("Q0 reproduced settling decision differs from ExpP")

    rows = DataFrame(
        check=["expN_candidate_sha", "expP_candidate_sha", "analytic_alpha",
               "PD_alpha", "complete_PD_poles", "trim_residual", "max_P_error_pu",
               "max_Q_error_pu", "bus38_algebraic_root_retained_MW",
               "nominal_beta_pointwise_upper", "beta_req", "event_peak_frequency_Hz",
               "event_peak_RoCoF_Hz_s", "event_settling_s", "max_contract_error_pu"],
        value=[nsha, psha, an.alpha, pdalpha, length(pdpoles)-1, trim.maximum,
               max_p, max_q, root_mw, beta_upper, beta_req, peakf, peakr,
               settling, max_contract_error],
        pass=[true, true, abs(an.alpha + 0.0500000011050941)<1e-8,
              abs(pdalpha + 0.050000000138500054)<1e-7, length(pdpoles)-1==101,
              trim.maximum<1e-10, max_p<1e-9, max_q<1e-9,
              abs(root_mw-retained)<1e-6, beta_upper<beta_req, true,
              abs(peakf-39.33340125578227)<0.03,
              abs(peakr-1.072208443371192)<0.01, !isfinite(settling),
              max_contract_error<1e-9])
    all(rows.pass) || error("Q0 baseline gate failed: $(rows[.!rows.pass,:])")
    CSV.write(joinpath(OUT, "TABLE_Q00_baseline_reproduction.csv"), rows)
    result = Dict(
        "stage"=>"Q0", "status"=>"PASS_BASELINE",
        "model_sha"=>"e2f104608f1eb1705f0ad2c07764e7beab0a5df4deb4d4ce4f359b1e79947f0a",
        "ExpN_candidate_sha256"=>nsha, "ExpP_candidate_sha256"=>psha,
        "alpha_analytic"=>an.alpha, "alpha_PD"=>pdalpha, "physical_poles"=>length(pdpoles)-1,
        "trim_residual"=>trim.maximum, "max_P_error_pu"=>max_p, "max_Q_error_pu"=>max_q,
        "retained_SG_MW"=>retained, "GFL_MW"=>sum(ctx.power)-retained,
        "bus38_root_retained_MW"=>root_mw, "pointwise_beta_upper"=>beta_upper,
        "beta_req"=>beta_req, "step_bus"=>bus, "step_MW"=>deltaP,
        "step_peak_frequency_Hz"=>peakf, "step_peak_RoCoF_Hz_s"=>peakr,
        "step_settling_s"=>settling, "step_simulation_s"=>horizon,
        "PD_rebuilds"=>1, "PD_event_simulations"=>1)
    ExpP.write_json(joinpath(OUT, "Q0_RESULTS.json"), result)
    write(joinpath(OUT, "STAGE_SUMMARY.md"), """
    # Q0 — Frozen ExpN/ExpP baseline reproduction

    Status: **PASS_BASELINE**. Exact ExpN and ExpP candidate hashes were verified before model calls. No input candidate, frozen model, Project, or Manifest was modified.

    - Rebuilt the ExpP mixed PD architecture and trimmed the declared P/Q shares. Maximum trim residual: `$(trim.maximum)`; maximum component P/Q contract errors: `$(max_p)` / `$(max_q)` pu.
    - Recomputed the full physical spectrum: analytic alpha `$(an.alpha)` s⁻¹, PD alpha `$(pdalpha)` s⁻¹, 101 physical poles.
    - Recomputed the ExpP bus-38 conditional algebraic root: `$(root_mw)` MW versus incumbent `$(retained)` MW.
    - Rechecked the nominal frequency-zero robustness witness: beta upper `$(beta_upper)` < required `$(beta_req)`.
    - Independently simulated the sustained 100 MW bus-16 load step for 60 s using PD. Peak frequency deviation `$(peakf)` Hz; peak sampled RoCoF `$(peakr)` Hz/s; settling `$(settling)` s. Q/P loads remain separate at buses 31 and 39.
    - Runtime records are in `Q0_RESULTS.json`; event TDS used one independent PD build and one 60 s simulation.
    """)
    println("EXP_Q_Q0_STATUS=PASS_BASELINE alpha_an=", an.alpha,
        " alpha_PD=", pdalpha, " event_df=", peakf, " event_RoCoF=", peakr,
        " settle=", settling)
end

main()
