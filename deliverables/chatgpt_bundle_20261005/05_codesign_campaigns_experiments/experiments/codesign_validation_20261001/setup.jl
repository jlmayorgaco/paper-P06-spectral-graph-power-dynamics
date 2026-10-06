using LinearAlgebra, CSV, DataFrames, TOML, SHA
BLAS.set_num_threads(1)
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const OUT=joinpath(ROOT,"reports","codesign_validation_20261001")
mkpath(OUT)
include(joinpath(ROOT,"src","bnd_design_p","ExpP.jl"))
include(joinpath(ROOT,"src","bnd_expQ","LinearSecurity.jl"))
include(joinpath(@__DIR__,"FiniteWindow.jl"))
include(joinpath(@__DIR__,"CoreDesign.jl"))
const N=ExpP.PDExactDesignN
const CTX=N.design_context(ROOT)
const ORIGINAL=TOML.parsefile(joinpath(@__DIR__,"frozen","original_candidate.toml"))
function save_toml(name,d)
    open(joinpath(OUT,name),"w") do io;TOML.print(io,d);end
end
function save_candidate(name,e,kp,ki; extra=Dict{String,Any}())
    d=Dict{String,Any}("epsilon"=>e,"rho"=>1 .-e,"Kp"=>kp,"Ki"=>ki,
      "support"=>collect(30:39)[e.>0],"retained_SG_MW"=>dot(CTX.power,e),
      "converted_GFL_MW"=>dot(CTX.power,1 .-e),"model_sha"=>ORIGINAL["model_sha"])
    merge!(d,extra); save_toml(name,d)
    open(joinpath(OUT,name*".sha256"),"w") do io
      print(io,bytes2hex(sha256(read(joinpath(OUT,name)))))
    end
end
