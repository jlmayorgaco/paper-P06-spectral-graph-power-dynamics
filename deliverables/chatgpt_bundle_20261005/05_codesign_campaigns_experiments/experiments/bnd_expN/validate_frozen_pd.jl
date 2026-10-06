using SHA, TOML
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const OUT=joinpath(ROOT,"reports","experiment_N")
const CANDIDATE=joinpath(OUT,"Z_N_NOMINAL_FINAL.toml")
isfile(CANDIDATE) && isfile(CANDIDATE*".sha256") || error("frozen candidate missing")
bytes2hex(sha256(read(CANDIDATE)))==strip(read(CANDIDATE*".sha256",String)) ||
    error("candidate SHA mismatch")
candidate=TOML.parsefile(CANDIDATE)
modeltext=read(joinpath(OUT,"MODEL_FREEZE.json"),String)
m=match(r"\"MODEL_SHA\":\s*\"([0-9a-f]+)\"",modeltext)
m!==nothing && candidate["model_sha"]==m.captures[1] || error("model SHA mismatch")

# Import PowerDynamics only after the analytical candidate hash is verified.
using CSV, DataFrames, DelimitedFiles, LinearAlgebra, NetworkDynamics, PowerDynamics
include(joinpath(ROOT,"src","bnd_model_expN","PDReferenceN.jl"))
include(joinpath(ROOT,"src","bnd_model_expN","PDExactDesignN.jl"))

function save_matrix(path,x)
    writedlm(path,Matrix(x),',')
end

function main()
    rho=Float64.(candidate["rho"])
    kp=Float64.(candidate["Kp"])
    ki=Float64.(candidate["Ki"])
    base=PDReferenceN.frozen_baseline()
    nw=PDReferenceN.build_architecture(base,rho,kp,ki)
    state=PDReferenceN.trim_state(nw,base,rho,kp,ki)
    residue=PDReferenceN.residual_audit(nw,state)
    power=PDReferenceN.direct_power_audit(state,base,rho)
    CSV.write(joinpath(OUT,"TABLE_N13_powerdynamics_component_power.csv"),power)
    maxp=maximum(Float64.(power.max_P_error_pu))
    maxq=maximum(Float64.(power.max_Q_error_pu))
    bounds=all(power.bounds_status .== "PASS")
    pd=linearize_network(state)
    λpd=jacobian_eigenvals(pd)
    jg=argmin(abs.(λpd))
    λphys=λpd[[i for i in eachindex(λpd) if i!=jg]]
    αpd=maximum(real.(λphys))
    analytic=PDExactDesignN.spectrum(PDExactDesignN.design_context(ROOT),rho,kp,ki)
    path=joinpath(OUT,"postfreeze_pd");mkpath(path)
    n=size(pd.A,1)
    mass=pd.M isa UniformScaling ? Matrix{Float64}(I,n,n) : Matrix(pd.M)
    save_matrix(joinpath(path,"PD_M.csv"),mass)
    save_matrix(joinpath(path,"PD_A.csv"),pd.A)
    CSV.write(joinpath(path,"PD_state_map.csv"),DataFrame(
        index=1:length(pd.sym),state_name=string.(pd.sym),differential=diag(mass).==1))
    CSV.write(joinpath(path,"PD_poles.csv"),DataFrame(real=real.(λpd),imag=imag.(λpd)))
    save_matrix(joinpath(path,"AN_Ared.csv"),analytic.model.Ared)
    CSV.write(joinpath(path,"AN_state_map.csv"),analytic.model.state_inventory)
    CSV.write(joinpath(path,"AN_poles.csv"),DataFrame(
        real=real.(analytic.lambda),imag=imag.(analytic.lambda)))
    row=(;candidate_sha=strip(read(CANDIDATE*".sha256",String)),
        trim_residual=residue.maximum,max_P_error_pu=maxp,max_Q_error_pu=maxq,
        physical_bounds_pass=bounds,PD_gauge_real=real(λpd[jg]),
        PD_gauge_imag=imag(λpd[jg]),PD_finite_count=length(λphys),
        AN_finite_count=length(analytic.lambda),PD_alpha=αpd,
        AN_alpha=analytic.alpha,alpha_error=abs(αpd-analytic.alpha),
        PD_margin_pass=αpd<=-0.05,AN_margin_pass=analytic.alpha<=-0.05)
    CSV.write(joinpath(OUT,"TABLE_N13_powerdynamics_raw.csv"),DataFrame([row]))
    println("PD_POSTFREEZE alpha=",αpd," analytic=",analytic.alpha,
        " difference=",abs(αpd-analytic.alpha)," residual=",residue.maximum,
        " bounds=",bounds);flush(stdout)
end

main()
