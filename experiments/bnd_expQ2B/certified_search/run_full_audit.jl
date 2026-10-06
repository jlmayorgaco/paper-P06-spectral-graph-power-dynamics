include("LocalOracle.jl")
using .LocalOracle, LinearAlgebra, TOML, CSV, DataFrames, SHA
const L=LocalOracle
Base.include(L,joinpath(@__DIR__,"MultiPeak.jl"))
include("SQP.jl")
include("GlobalAudit.jl")
include("RiccatiCertificate.jl")
d=TOML.parsefile(joinpath(L.OUT,"FEASIBLE_NUMERICAL_REPRESENTATIVE.toml"))
a=L.architecture(d["support"],d);x=L.encode(a.support,d)
v=L.multi_evaluate(a,x;derivatives=false)
println("FULL_AUDIT_START J=$(v.J)");flush(stdout)
CSV.write(joinpath(L.OUT,"ANALYTICAL_A.csv"),DataFrame(v.A,:auto))
br=riccati_certificate(v.A)
open(joinpath(L.OUT,"FULL_BAND_CERTIFICATE.toml"),"w") do io;TOML.print(io,br);end
println(br);flush(stdout)
ta=time_full_audit(v)
open(joinpath(L.OUT,"ALL_TIME_CERTIFICATE.toml"),"w") do io;TOML.print(io,Dict(string(k)=>value for (k,value) in pairs(ta)));end
println(ta);flush(stdout)
for line in eachline(stdin)
    strip(line)=="EXIT" && break
    try;include_string(Main,line,"interactive_audit");catch err;showerror(stdout,err,catch_backtrace());println();end
    println("COMMAND_DONE");flush(stdout)
end
