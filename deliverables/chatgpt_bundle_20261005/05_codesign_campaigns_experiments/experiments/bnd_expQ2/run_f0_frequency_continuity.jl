using CSV, DataFrames, LinearAlgebra, TOML
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const OUT=joinpath(ROOT,"reports","experiment_Q2","F0")
mkpath(OUT)
include(joinpath(ROOT,"src","bnd_design_p","ExpP.jl"))
include(joinpath(ROOT,"src","bnd_expQ","LinearSecurity.jl"))
include(joinpath(ROOT,"src","bnd_expQ","FrequencyMetrics.jl"))
include(joinpath(ROOT,"src","bnd_expQ2","GridFrequency.jl"))
using .ExpP, .LinearSecurity, .FrequencyMetrics, .GridFrequency
const N=ExpP.PDExactDesignN
const MODEL_SHA="e2f104608f1eb1705f0ad2c07764e7beab0a5df4deb4d4ce4f359b1e79947f0a"
const RHO_TEST=(0.0,1e-8,1e-7,1e-6,1e-5,1e-4)
const WINDOWS=(0,2,5,10,20)
const DT=0.01
const HORIZON=60.0

function metrics_dict(rows)
    d=Dict{String,Float64}()
    for row in rows
        e=String(row.estimator)
        d["$(e):F_peak_Hz"]=Float64(row.F_peak_Hz)
        d["$(e):R_peak_Hz_s"]=Float64(row.R_peak_Hz_s)
        d["$(e):F_inf_Hz"]=Float64(row.F_inf_Hz)
        d["$(e):phase_jump_peak_rad"]=Float64(row.phase_jump_peak_rad)
    end
    d
end

function main()
    modelsha=ExpP.verify_expN_freeze(ROOT)
    modelsha==MODEL_SHA || error("ExpN frozen model SHA mismatch")
    ctx=N.design_context(ROOT)
    kp=fill(N.K0P,10);ki=fill(N.K0I,10);rho0=zeros(10)
    times=collect(-0.5:DT:HORIZON)
    t0=time()
    m0=N.descriptor(ctx,rho0,kp,ki)
    s0=grid_step_signals(ctx,m0,rho0,times;load_bus=16,disturbance_MW=1.0,
        gauge_vector=N.gauge_vector)
    rows0=grid_step_metrics(s0;half_windows=WINDOWS,dt_s=DT,event_time_s=0.0)
    md0=metrics_dict(rows0)
    records=NamedTuple[];continuity=NamedTuple[]
    function savecase(bus,rho_val,m,sig,met)
        for x in met
            push!(records,(;bus, rho_GFL=rho_val,epsilon_SG=1-rho_val,
                estimator=x.estimator,half_window=x.half_window,
                filter_support_s=x.window_support_s,derivative_minus3dB_Hz=x.derivative_minus3dB_Hz,
                F_peak_Hz=x.F_peak_Hz,R_peak_Hz_s=x.R_peak_Hz_s,F_inf_Hz=x.F_inf_Hz,
                phase_jump_peak_rad=x.phase_jump_peak_rad,alpha_per_s=sig.alpha,
                condition_A=sig.condition_A,eigenvector_condition=sig.eigenvector_condition,
                gauge_residual=sig.reduction_gauge_residual))
        end
    end
    for bus in 30:39
        savecase(bus,0.0,m0,s0,rows0)
        i=bus-29
        for r in RHO_TEST[2:end]
            rho=copy(rho0);rho[i]=r
            m=N.descriptor(ctx,rho,kp,ki)
            sig=grid_step_signals(ctx,m,rho,times;load_bus=16,disturbance_MW=1.0,
                gauge_vector=N.gauge_vector)
            met=grid_step_metrics(sig;half_windows=WINDOWS,dt_s=DT,event_time_s=0.0)
            savecase(bus,r,m,sig,met)
            md=metrics_dict(met)
            append!(continuity,continuity_rows(bus,r,md0,md;condition=sig.condition_A,
                ref_condition=s0.condition_A,rho_value=r,alpha=sig.alpha))
        end
    end
    tbl=DataFrame(records);ct=DataFrame(continuity)
    CSV.write(joinpath(OUT,"TABLE_Q2_F0_frequency_metric_continuity.csv"),tbl)
    CSV.write(joinpath(OUT,"TABLE_Q2_F0_continuity_errors.csv"),ct)
    atedge=ct[ct.rho_GFL.==1e-8,:]
    pass=all(atedge.pass_within_tolerance)
    status=pass ? "PASS_CONTINUOUS" : "FAIL_DISCONTINUOUS"
    bandwidth=DataFrame(estimator=[x.estimator for x in rows0],
        half_window=[x.half_window for x in rows0],support_s=[x.window_support_s for x in rows0],
        minus3dB_Hz=[x.derivative_minus3dB_Hz for x in rows0],
        F_peak_1MW_Hz=[x.F_peak_Hz for x in rows0],R_peak_1MW_Hz_s=[x.R_peak_Hz_s for x in rows0],
        F_inf_1MW_Hz=[x.F_inf_Hz for x in rows0],phase_jump_rad=[x.phase_jump_peak_rad for x in rows0])
    CSV.write(joinpath(OUT,"TABLE_Q2_F0_bandwidth_sensitivity_allSG.csv"),bandwidth)
    result=(;stage="F0",status,grid_frequency_metric_status=status,model_sha=modelsha,
        inherited_measurement="ExpG inertia-weighted retained-SG COI; ExpP retained-SG COI; no architecture-invariant bus-voltage frequency output or frozen measurement bandwidth was found",
        primary_grid_frequency="frequency from unwrapped bus-voltage phase at buses 30:39, same passive measurement in every architecture",
        sampling_s=DT,phase_differentiators=collect(WINDOWS),disturbance_MW=1.0,
        event_bus=16,horizon_s=HORIZON,architecture_comparisons=length(RHO_TEST)*10,
        continuity_rows=nrow(ct),edge_rho=1e-8,
        maximum_edge_error=maximum(atedge.absolute_error),maximum_tolerance=maximum(atedge.tolerance_condition_aware),
        maximum_error_by_rho=combine(groupby(ct,:rho_GFL),:absolute_error=>maximum=>:max_abs_error),
        unfiltered_phase_jump_max_rad=maximum(abs.(s0.phase_jump_rad)),
        unfiltered_step_rocof_impulse=maximum(abs.(s0.phase_jump_rad))>0,
        rocof_limit_status="METRIC_DEPENDENT: inherited 0.5 Hz/s limits were applied to SG COI and no frozen bus-measurement bandwidth was found",
        elapsed_s=time()-t0,PowerDynamics_calls=0)
    open(joinpath(OUT,"F0_RESULTS.toml"),"w") do io
        TOML.print(io,Dict("stage"=>result.stage,"status"=>result.status,
            "model_sha"=>result.model_sha,"inherited_measurement"=>result.inherited_measurement,
            "primary_grid_frequency"=>result.primary_grid_frequency,"sampling_s"=>result.sampling_s,
            "phase_differentiators"=>result.phase_differentiators,"disturbance_MW"=>result.disturbance_MW,
            "event_bus"=>result.event_bus,"horizon_s"=>result.horizon_s,
            "architecture_comparisons"=>result.architecture_comparisons,"continuity_rows"=>result.continuity_rows,
            "edge_rho"=>result.edge_rho,"maximum_edge_error"=>result.maximum_edge_error,
            "maximum_tolerance"=>result.maximum_tolerance,"unfiltered_phase_jump_max_rad"=>result.unfiltered_phase_jump_max_rad,
            "unfiltered_step_rocof_impulse"=>result.unfiltered_step_rocof_impulse,
            "rocof_limit_status"=>result.rocof_limit_status,"elapsed_s"=>result.elapsed_s,"PowerDynamics_calls"=>0))
    end
    write(joinpath(OUT,"STAGE_SUMMARY.md"),"""
# F0 — Architecture-invariant voltage-frequency metric

**Gate: $status.** The audit finds that ExpG used an inertia-weighted retained-SG COI output and ExpP used retained-SG COI at the retained machine; neither freezes an architecture-invariant bus-frequency sensor or filter bandwidth. ExpQ's Savitzky–Golay bus derivative was a diagnostic, not a frozen primary requirement.

For design measurements, this stage uses the same passive quantity for every architecture: unwrapped busbar-voltage phase at generator buses 30–39, differentiated at a fixed 100 Hz sample rate. The 1 MW bus-16 event is evaluated for all ten buses at `rho=0, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4`; all-SG is the reference. A sampled unfiltered finite difference is retained along with cubic Savitzky–Golay phase differentiators with 0.04, 0.10, 0.20 and 0.40 s supports. The filter outputs are a bandwidth sensitivity study, not a post-hoc choice of the passing metric.

At `rho=0`, the algebraic bus-voltage phase has a jump of up to $(maximum(abs.(s0.phase_jump_rad))) rad for the 1 MW step. Its ideal unfiltered derivative contains an impulse; consequently the continuous-time unfiltered RoCoF supremum is unbounded. Sampled and filtered RoCoF values are finite but depend on estimator bandwidth. The inherited 0.5 Hz/s design limit was previously applied to SG COI and has no frozen bus-measurement bandwidth; its application to bus-frequency RoCoF is therefore **METRIC_DEPENDENT**.

Continuity uses condition-aware tolerance `max(1e-8 × scale, 100 eps × max(cond(A),cond(A_all-SG)) × scale)`. The rho=1e-8 endpoint maximum absolute discrepancy is $(maximum(atedge.absolute_error)) versus tolerance $(maximum(atedge.tolerance_condition_aware)); pass is $(pass). The metric is continuous at the tested GFL-insertion boundary only to this finite-horizon, sampled linear test. See the detailed tables for each output/window, rho level, and conditioning.

No PowerDynamics call or optimization was used. Evaluation cost: $(nrow(ct)) per-output continuity comparisons over $(length(RHO_TEST)*10) architectures; elapsed $(result.elapsed_s) s. This gate allows proceeding only if the endpoint continuity checks pass. It does not resolve a defensible primary RoCoF bandwidth.
""")
    println("F0 status=$status rows=$(nrow(ct)) max_edge_error=$(maximum(atedge.absolute_error)) tol=$(maximum(atedge.tolerance_condition_aware)) jump=$(maximum(abs.(s0.phase_jump_rad))) elapsed=$(result.elapsed_s)")
end
main()
