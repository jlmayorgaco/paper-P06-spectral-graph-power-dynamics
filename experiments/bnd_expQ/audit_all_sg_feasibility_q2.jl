using CSV, DataFrames, LinearAlgebra, SHA, TOML
using NetworkDynamics, PowerDynamics, OrdinaryDiffEqRosenbrock, SciMLBase

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const OUT=joinpath(ROOT,"reports","experiment_Q","Q2_REFERENCE")
mkpath(OUT)
include(joinpath(ROOT,"src","bnd_design_p","ExpP.jl"))
include(joinpath(ROOT,"src","bnd_model_expN","PDReferenceN.jl"))
include(joinpath(ROOT,"src","bnd_opt_expP","RobustAudit.jl"))
include(joinpath(ROOT,"src","bnd_expQ","LinearSecurity.jl"))
include(joinpath(ROOT,"src","bnd_expQ","FrequencyMetrics.jl"))
using .ExpP, .PDReferenceN, .RobustAudit, .LinearSecurity, .FrequencyMetrics
const N=ExpP.PDExactDesignN

function add_event(nw,bus,p0,q0,amp;start=1.0,base_mva=100.0)
    vertices,edges=PDReferenceN.PD39.PD39Model.copy_network_components(nw)
    defaults=get_defaults_dict(vertices[bus])
    pk=only([x for x in keys(defaults) if occursin("Pset",string(x))])
    qk=only([x for x in keys(defaults) if occursin("Qset",string(x))])
    affect=(u,p,ctx)->begin
        p[pk]=ctx.t>=start ? p0-amp/base_mva : p0
        p[qk]=q0
    end
    cb=PresetTimeComponentCallback([start],ComponentAffect(affect,(),(pk,qk)))
    set_callback!(vertices[bus],cb)
    active=Network(vertices,edges);set_jac_prototype!(active)
    active
end

function main()
    req=1.6991206999182038e-6; sigma=0.05; deltaP=100.0
    cand_path=joinpath(OUT,"Z_Q_ALL_SG_SECURE_REFERENCE.toml")
    cand_sha=bytes2hex(sha256(read(cand_path)))
    cand_sha==strip(read(cand_path*".sha256",String)) || error("all-SG reference candidate hash mismatch")
    cand=TOML.parsefile(cand_path)
    cand["model_sha"]==ExpP.verify_expN_freeze(ROOT) || error("all-SG reference model SHA mismatch")
    ctx=N.design_context(ROOT); rho=Float64.(cand["rho"]);kp=Float64.(cand["Kp"]);ki=Float64.(cand["Ki"])
    t0=time();sp=N.spectrum(ctx,rho,kp,ki)
    as=Matrix{Float64}(transpose(sp.quotient)*sp.model.Ared*sp.quotient)+sigma*I
    beta0=RobustAudit.sigma_min_at(as,0.0)
    # Bounded-real Riccati/LMI certificate, evaluated in Float64. This is a
    # numerical certificate; the run does not use directed-rounding arithmetic.
    care=RobustAudit.hinf_upper_certificate(as,1/req)
    model=reduced_step_model(ctx,sp.model,rho,16;gauge_vector=N.gauge_vector)
    metrics=step_metrics(model;disturbance_MW=deltaP,horizon_s=60.0,dt_s=0.02)
    robust_pass=care.certified
    println("ALLSG alpha=$(sp.alpha) beta0=$(beta0) beta_lower=$(care.beta_lower) CARE=$(care.status) Fpeak=$(metrics.frequency_peak_MW) Rpeak=$(metrics.rocof_peak_MW)");flush(stdout)
    pre=(;retained_SG_MW=sum(ctx.power),alpha=sp.alpha,physical_poles=length(sp.lambda),
        beta_req=req,beta_pointwise_omega0_upper=beta0,beta_observed_upper=NaN,
        beta_numeric_lower=care.certified ? care.beta_lower : NaN,beta_numeric_status=care.status,
        care_relative_residual=care.care_relative_residual,care_max_eigenvalue=care.care_max_eigenvalue,
        min_P_eigenvalue=care.min_P_eigenvalue,closed_loop_abscissa=care.closed_loop_abscissa,
        beta_formal_certified=false,linear_F_inf_Hz=metrics.steady_frequency_MW,
        linear_F_peak_Hz=metrics.frequency_peak_MW,linear_R_peak_Hz_s=metrics.rocof_peak_MW)
    ExpP.write_json(joinpath(OUT,"Q2_ALL_SG_PRE_TDS.json"),Dict(string(k)=>v for (k,v) in pairs(pre)))

    base=PDReferenceN.frozen_baseline()
    nw=PDReferenceN.build_architecture(base,rho,kp,ki)
    state=PDReferenceN.trim_state(nw,base,rho,kp,ki)
    powers=PDReferenceN.direct_power_audit(state,base,rho)
    residual=PDReferenceN.residual_audit(nw,state)
    bus=16;load=CSV.read(joinpath(ROOT,"reports","experiment_D","inputs","load.csv"),DataFrame)
    lr=only(eachrow(load[load.bus.==bus,:]))
    active=add_event(nw,bus,Float64(lr.Pset),Float64(lr.Qset),deltaP)
    started=time()
    sol=SciMLBase.solve(SciMLBase.ODEProblem(active,deepcopy(state),(0.0,60.0)),
        OrdinaryDiffEqRosenbrock.Rodas5P();callback=get_callbacks(active),
        initializealg=SciMLBase.NoInit(),saveat=0.01,abstol=1e-9,reltol=1e-9)
    SciMLBase.successful_retcode(sol.retcode)||error("allSG reference TDS failed: $(sol.retcode)")
    tt=collect(0.0:0.01:60.0); f=fill(NaN,length(tt),10); volts=Matrix{ComplexF64}(undef,length(tt),10)
    for (k,t) in enumerate(tt)
        s=NetworkDynamics.NWState(sol,t)
        for b in 30:39
            j=b-29; key=b==39 ? :machine₊ω : :ctrld_gen₊machine₊ω
            f[k,j]=sg_frequency_deviation(Float64(s[VIndex(b,key)]))
            ur=Float64(s[VIndex(b,:busbar₊u_r)]);ui=Float64(s[VIndex(b,:busbar₊u_i)])
            volts[k,j]=complex(ur,ui)
        end
    end
    rocof=similar(f)
    for j in 1:10;rocof[:,j].=sg_polynomial_derivative(f[:,j],0.01;half_window=5,degree=3);end
    activeidx=findall(tt.>=1.0); fpeak=maximum(abs.(f[activeidx,:])); rpeak=maximum(abs.(rocof[activeidx,:]))
    fb=hcat([bus_frequency_deviation(volts[:,j],0.01;half_window=5,degree=3) for j in 1:10]...)
    buspeak=maximum(abs.(fb[activeidx,:]))
    td=DataFrame(time_s=tt)
    for j in 1:10;td[!,Symbol("bus$(29+j)_SG_df_Hz")]=f[:,j];td[!,Symbol("bus$(29+j)_SG_RoCoF_Hz_s")]=rocof[:,j];td[!,Symbol("bus$(29+j)_bus_df_Hz")]=fb[:,j];end
    CSV.write(joinpath(OUT,"TABLE_Q2_all_sg_reference_trajectory.csv"),td)
    pqerr=max(maximum(powers.max_P_error_pu),maximum(powers.max_Q_error_pu))
    summary=(;candidate_sha256=cand_sha,model_sha=ExpP.verify_expN_freeze(ROOT),retained_SG_MW=sum(ctx.power),converted_GFL_MW=0.0,
        alpha=sp.alpha,physical_poles=length(sp.lambda),sigma_req=sigma,beta_req=req,
        beta_pointwise_omega0_upper=beta0,beta_observed_upper=NaN,
        beta_observed_omega=NaN,beta_numeric_lower=care.certified ? care.beta_lower : NaN,
        beta_numeric_status=care.status,beta_formal_certified=false,
        beta_CARE_relative_residual=care.care_relative_residual,
        beta_CARE_LMI_max_eigenvalue=care.care_max_eigenvalue,
        beta_CARE_min_P_eigenvalue=care.min_P_eigenvalue,
        beta_CARE_closed_loop_abscissa=care.closed_loop_abscissa,beta_pass=robust_pass,
        linear_F_inf=metrics.steady_frequency_MW,linear_F_peak=metrics.frequency_peak_MW,
        linear_R_peak=metrics.rocof_peak_MW,PD_F_peak=fpeak,PD_R_peak=rpeak,
        PD_bus_frequency_peak=buspeak,frequency_limit_Hz=0.5,rocof_limit_Hz_s=0.5,
        frequency_pass=fpeak<=0.5,rocof_pass=rpeak<=0.5,bus_frequency_diagnostic_peak=buspeak,
        trim_residual=residual.maximum,max_PQ_error_pu=pqerr,solver_status=string(sol.retcode),
        tds_elapsed_s=time()-started,total_elapsed_s=time()-t0,
        reference_feasible=robust_pass && fpeak<=0.5 && rpeak<=0.5)
    ExpP.write_json(joinpath(OUT,"Q2_ALL_SG_REFERENCE.json"),Dict(string(k)=>v for (k,v) in pairs(summary)))
    CSV.write(joinpath(OUT,"TABLE_Q2_all_sg_reference_metrics.csv"),DataFrame([summary]))
    write(joinpath(OUT,"STAGE_SUMMARY.md"),"""
# Q2 feasibility reference — all SG

This is an upper-bound reference, not an optimized Q design. It uses the frozen ExpN model, controller bounds and the exact +100 MW bus-16 sustained event with Q unchanged.

- Complete physical spectrum: alpha=$(sp.alpha) s⁻¹ across $(length(sp.lambda)) physical poles.
- Robustness: beta requirement $(req); omega=0 pointwise upper witness $(beta0). The Bounded Real CARE/LMI result is $(care.status), beta lower bound $(care.certified ? care.beta_lower : NaN), CARE relative residual $(care.care_relative_residual), minimum P eigenvalue $(care.min_P_eigenvalue), maximum CARE/LMI residual eigenvalue $(care.care_max_eigenvalue), and closed-loop abscissa $(care.closed_loop_abscissa). It uses Float64 arithmetic without directed rounding, so it is numerically certified, not formally validated.
- Analytic step model: F_inf=$(metrics.steady_frequency_MW) Hz, F_peak=$(metrics.frequency_peak_MW) Hz, R_peak=$(metrics.rocof_peak_MW) Hz/s. Independent PD TDS over 60 s: local SG frequency peak $(fpeak) Hz, SG RoCoF $(rpeak) Hz/s, bus-frequency diagnostic peak $(buspeak) Hz.
- Trim residual $(residual.maximum), P/Q maximum error $(pqerr) pu. TDS wall time $(time()-started) s.
- Reference feasible under declared constraints: $(summary.reference_feasible). No optimization, candidate freeze, or globality claim is made.
""")
end

main()
