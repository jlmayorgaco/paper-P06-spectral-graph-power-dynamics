include("LocalOracle.jl")
using .LocalOracle, LinearAlgebra, TOML, CSV, DataFrames, SHA
const L=LocalOracle
Base.include(L,joinpath(@__DIR__,"MultiPeak.jl"))
include("SQP.jl")
include("Certificate.jl")
include("GlobalAudit.jl")
d=TOML.parsefile(joinpath(L.OUT,"MULTIPEAK_SOC_FINAL.toml"))
a=L.architecture(d["support"],d);x=L.encode(a.support,d)
cert=local_certificate(a,x)
cert["local_KKT_certificate"] || error("Local certificate failed")
safe=feasible_roundoff_representative(a,x)
certsafe=local_certificate(a,safe.x;label="LOCAL_CERTIFICATE_FEASIBLE")
certsafe["local_KKT_certificate"] || error("Numerical representative certificate failed")
