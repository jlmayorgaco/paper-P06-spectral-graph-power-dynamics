ENV["BND_SEARCH_SUBDIRECTORY"]=length(ARGS)>1 ? ARGS[2] : "stationarity"
ENV["BND_MODAL_BUNDLE_SIZE"]=length(ARGS)>2 ? ARGS[3] : "1"
include("search_oracle.jl")
d=TOML.parsefile(abspath(ARGS[1]))
p=vcat(Float64.(d["rho"]),log.(Float64.(d["Kp"])),log.(Float64.(d["Ki"])))
out=evaluate(p;grad=true,gradient_horizon=60.,gradient_dt=.025)
println("FINAL_GRADIENT_COMPLETE constraints=",out[1:NBUNDLE+5length(CASES)]);flush(stdout)
