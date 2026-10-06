using CSV, DataFrames, LinearAlgebra, SHA, TOML
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const OUT=joinpath(ROOT,"reports","experiment_Q","Q2_REFERENCE")
include(joinpath(ROOT,"src","bnd_design_p","ExpP.jl"))
include(joinpath(ROOT,"src","bnd_opt_expP","RobustAudit.jl"))
include(joinpath(ROOT,"src","bnd_expQ","LinearSecurity.jl"))
using .ExpP, .RobustAudit, .LinearSecurity
const N=ExpP.PDExactDesignN

function main()
    src=joinpath(ROOT,"reports","experiment_G","Z_G_FINAL.toml")
    srcsha=bytes2hex(sha256(read(src)))
    srcsha==strip(read(src*".sha256",String)) || error("frozen ExpG seed hash mismatch")
    old=TOML.parsefile(src); gens=old["generator"]
    sort!(gens;by=x->Int(x["bus"]))
    Int[x["bus"] for x in gens]==collect(30:39) || error("ExpG seed bus inventory mismatch")
    rho=Float64[x["rho"] for x in gens];kp=Float64[x["Kp"] for x in gens];ki=Float64[x["Ki"] for x in gens]
    ctx=N.design_context(ROOT); all(ctx.kpmin .<= kp .<= ctx.kpmax) || error("Kp outside frozen ExpK box")
    all(ctx.kimin .<= ki .<= ctx.kimax) || error("Ki outside frozen ExpK box")
    sigma=0.05; beta_req=1.6991206999182038e-6
    sp=N.spectrum(ctx,rho,kp,ki)
    model=reduced_step_model(ctx,sp.model,rho,16;gauge_vector=N.gauge_vector)
    met=step_metrics(model;disturbance_MW=100.0,horizon_s=60.0,dt_s=0.02)
    as=Matrix{Float64}(transpose(sp.quotient)*sp.model.Ared*sp.quotient)+sigma*I
    care=RobustAudit.hinf_upper_certificate(as,1/beta_req)
    retained=dot(ctx.power,1 .- rho)
    analytic_pass=sp.alpha<=-sigma && care.certified &&
        met.frequency_peak_MW<=0.5 && met.rocof_peak_MW<=0.5
    row=(;source="frozen ExpG coordinate seed, re-evaluated with ExpN PD-exact model",
        source_sha256=srcsha,model_sha=ExpP.verify_expN_freeze(ROOT),retained_SG_MW=retained,
        converted_GFL_MW=sum(ctx.power)-retained,rho=join(rho,";"),Kp=join(kp,";"),Ki=join(ki,";"),
        alpha=sp.alpha,physical_poles=length(sp.lambda),beta_lower=care.certified ? care.beta_lower : NaN,
        beta_status=care.status,CARE_relative_residual=care.care_relative_residual,
        CARE_LMI_max_eigenvalue=care.care_max_eigenvalue,CARE_min_P_eigenvalue=care.min_P_eigenvalue,
        F_inf_Hz=met.steady_frequency_MW,F_peak_Hz=met.frequency_peak_MW,
        R_peak_Hz_s=met.rocof_peak_MW,F_peak_active_bus=Int(met.table.bus[argmax(met.table.F_peak_unit_Hz)]),
        R_peak_active_bus=Int(met.table.bus[argmax(met.table.R_peak_unit_Hz_s)]),
        analytic_constraints_pass=analytic_pass,PD_called=false)
    CSV.write(joinpath(OUT,"TABLE_Q2_expG_seed_exact_model.csv"),DataFrame([row]))
    if analytic_pass
        cand=Dict{String,Any}("candidate_kind"=>"re-trimmed ExpG seed; not Q optimum",
            "candidate_frozen_before_independent_PD_validation"=>true,
            "source_ExpG_candidate_sha256"=>srcsha,"model_sha"=>ExpP.verify_expN_freeze(ROOT),
            "rho"=>rho,"Kp"=>kp,"Ki"=>ki,"retained_SG_MW"=>retained,
            "converted_GFL_MW"=>sum(ctx.power)-retained,"sigma_req_per_s"=>sigma,
            "beta_req_normalized"=>beta_req,"event_bus"=>16,"event_delta_P_MW"=>100.0,
            "frequency_limit_Hz"=>0.5,"rocof_limit_Hz_s"=>0.5)
        dest=joinpath(OUT,"Z_Q_EXP_G_RETRIMMED_REFERENCE.toml")
        open(dest,"w") do io; TOML.print(io,cand;sorted=true); end
        sha=bytes2hex(sha256(read(dest)));write(dest*".sha256",sha*"\n")
    end
    ExpP.write_json(joinpath(OUT,"Q2_EXPG_SEED_ANALYTIC.json"),Dict(
        "source_sha256"=>srcsha,"analytic_constraints_pass"=>analytic_pass,
        "retained_SG_MW"=>retained,"alpha"=>sp.alpha,"beta_numeric_lower"=>(care.certified ? care.beta_lower : nothing),
        "beta_status"=>care.status,"F_inf_Hz"=>met.steady_frequency_MW,
        "F_peak_Hz"=>met.frequency_peak_MW,"R_peak_Hz_s"=>met.rocof_peak_MW,
        "candidate_frozen"=>analytic_pass,"PD_called"=>false))
    println("EXPG_SEED retained=$(retained) alpha=$(sp.alpha) beta=$(care.certified ? care.beta_lower : NaN) Fpeak=$(met.frequency_peak_MW) Rpeak=$(met.rocof_peak_MW) pass=$(analytic_pass)")
end

main()
