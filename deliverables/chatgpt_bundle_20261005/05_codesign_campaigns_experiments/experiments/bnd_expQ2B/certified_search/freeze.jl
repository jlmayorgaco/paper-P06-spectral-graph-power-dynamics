include("LocalOracle.jl")
using .LocalOracle, LinearAlgebra, TOML, CSV, DataFrames, SHA, Dates
const L=LocalOracle
Base.include(L,joinpath(@__DIR__,"MultiPeak.jl"))
include("SQP.jl")
include("adjacent_exact_limits.jl")
cert=TOML.parsefile(joinpath(L.OUT,"LOCAL_CERTIFICATE_FEASIBLE.toml"))
br=TOML.parsefile(joinpath(L.OUT,"FULL_BAND_CERTIFICATE.toml"))
ta=TOML.parsefile(joinpath(L.OUT,"ALL_TIME_CERTIFICATE.toml"))
cert["local_KKT_certificate"]&&br["strict_negative_LMI"]&&ta["pass"] || error("Pre-freeze gates not passed")
d=TOML.parsefile(joinpath(L.OUT,"FEASIBLE_NUMERICAL_REPRESENTATIVE.toml"))
a=L.architecture(d["support"],d);x=L.encode(a.support,d);v=L.multi_evaluate(a,x;derivatives=true)
limits=adjacent_exact_limits(d)
CSV.write(joinpath(L.OUT,"FROZEN_ANALYTICAL_POLES.csv"),DataFrame(real=real.(v.lam),imag=imag.(v.lam)))
CSV.write(joinpath(L.OUT,"FROZEN_STATE_MAP.csv"),v.m.state_inventory)
CSV.write(joinpath(L.OUT,"FROZEN_PHYSICAL_A.csv"),DataFrame(v.A,:auto))
ts=collect(0.:.01:60.);df=DataFrame(event_time_s=ts)
for k in 1:10,kind in (:F,:R)
    df[!,Symbol("$(kind)_bus$(k+29)")]=[L.output(v,t,kind,k) for t in ts]
end
CSV.write(joinpath(L.OUT,"FROZEN_LINEAR_TRAJECTORIES.csv"),df)
au=audit(a,x,v);mu=zeros(length(v.g))
for (id,m) in zip(au.ids,au.mu);id>0 && (mu[id]=m);end
authorities=NamedTuple[]
for (j,bus) in enumerate(a.support)
    p=L.CTX.power[bus-29]
    # Objective is J/1000 and constraints are normalized. Therefore the
    # physical multiplier is 1000*mu; the resulting weighted authority is 1.
    contributions=-1000mu.*v.Jac[:,j]/p
    push!(authorities,(;bus,modal=contributions[1],robust=contributions[2],steady=contributions[3],
        peak_frequency=sum(contributions[findall(startswith("F_"),v.names)]),
        RoCoF=sum(contributions[findall(startswith("R_"),v.names)]),total=sum(contributions)))
end
CSV.write(joinpath(L.OUT,"MARGINAL_SECURITY_VALUES.csv"),DataFrame(authorities))
candidate=merge(d,Dict("local_certified"=>true,"certificate_type"=>"NUMERICAL_FIXED_ARCHITECTURE_KKT_SOSC",
    "epsilon"=>1 .-d["rho"],"total_original_SG_MW"=>sum(L.CTX.power),
    "converted_GFL_MW"=>sum(L.CTX.power)-v.J,"GFL_fraction"=>1-v.J/sum(L.CTX.power),
    "beta_requirement"=>L.BETA,"sigma_requirement"=>.05,"frequency_limit_Hz"=>.5,
    "RoCoF_limit_Hz_s"=>.5,"measurement_window_s"=>.5,"disturbance_bus"=>16,
    "disturbance_MW"=>100.,"event"=>"Sustained initialized ZIP Pset change; Qset unchanged",
    "primary_frequency"=>"Ten generator-bus voltage phases: causal phase difference / (2pi T); RoCoF second phase difference / (2pi T^2)",
    "frozen_at_UTC"=>string(now(UTC)),"model_sha"=>"e2f104608f1eb1705f0ad2c07764e7beab0a5df4deb4d4ce4f359b1e79947f0a",
    "global_certified"=>false,"PowerDynamics_validation_at_freeze"=>"NOT_RUN"))
path=joinpath(L.OUT,"Z_LOCAL_SECURE_FINAL.toml");isfile(path)&&error("Immutable candidate already exists")
open(path,"w") do io;TOML.print(io,candidate);end
write(path*".sha256",bytes2hex(sha256(read(path)))*"\n")
println("FROZEN ",path," SHA ",strip(read(path*".sha256",String)));flush(stdout)
