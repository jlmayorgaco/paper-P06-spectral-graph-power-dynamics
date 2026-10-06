using CSV, DataFrames, LinearAlgebra, TOML
include(joinpath(@__DIR__, "..", "..", "src", "bnd_design_p", "ExpP.jl"))
using .ExpP
include(joinpath(@__DIR__, "..", "..", "src", "bnd_opt_expP", "RobustAudit.jl"))
using .RobustAudit

function main()
    started=time()
    root=ExpP.ROOT
    ctx=ExpP.PDExactDesignN.design_context(root)
    candidate,_=ExpP.read_candidate(root)
    p4=read(joinpath(root,"reports","experiment_P","P4","P4_RESULTS.json"),String)
    m=match(r"\"conditional_point_epsilon\"\s*:\s*([0-9.eE+-]+)",p4)
    m===nothing && error("P4 conditional point not recorded")
    epsilon=parse(Float64,m.captures[1])
    kp=Float64.(candidate["Kp"]);ki=Float64.(candidate["Ki"])
    mb=match(r"\"beta_req\"\s*:\s*([0-9.eE+-]+)",p4)
    mb===nothing && error("P4 beta requirement not recorded")
    beta=parse(Float64,mb.captures[1])
    rho=ones(10);rho[38-29]=1-epsilon
    sp=ExpP.design_spectrum(ctx,rho,kp,ki)
    sh=RobustAudit.shifted_quotient(ctx,38,epsilon,kp,ki)
    cert=RobustAudit.certify_radius(sh.As,beta;max_nodes=50000,max_depth=64)
    row=(;candidate="P4 conditional scalar point",epsilon,beta_req=beta,
        sigma_req=0.05,alpha=sp.alpha,physical_poles=length(sp.lambda),
        nominal_spectrum_pass=sp.alpha<=-0.05+1e-10,
        interval_status=cert.status,certified=cert.certified,
        beta_lower_cert=cert.beta_lower_cert,beta_upper_observed=cert.beta_upper_observed,
        gamma_lower_observed=cert.gamma_lower_observed,
        omega_witness=cert.omega_witness,nodes=cert.nodes,max_depth=cert.max_depth,
        tail_lower=cert.tail_lower,max_nodes=50000,elapsed_s=time()-started)
    out=joinpath(root,"reports","experiment_P","P4","TABLE_P4_interval_certificate_conditional.csv")
    CSV.write(out,DataFrame([row]))
    println(row)
end

main()
