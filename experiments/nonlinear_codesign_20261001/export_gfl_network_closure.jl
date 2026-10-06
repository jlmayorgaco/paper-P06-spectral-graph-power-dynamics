include("ReducedDAE.jl")
using .ReducedDAE, LinearAlgebra, TOML, SHA
const R=ReducedDAE
const OUT=joinpath(R.OUT,"full_gfl_sector_contract")
ctx=R.N.design_context(R.ROOT)
candidate=joinpath(R.OUT,"candidate_final_physical.toml")
dd=TOML.parsefile(candidate)
m=R.model(ctx,Float64.(dd["rho"]),Float64.(dd["Kp"]),Float64.(dd["Ki"]);dc_convention=:physical_supply)
G=copy(m.net.reduced)
for i in 1:10
    D,c=R.norton(m.x0[m.sgidx[i]],m.sp[i])
    ix=2i-1:2i
    G[ix,ix].+=(1-m.rho[i]).*D
end
Z=-inv(G)
out=Dict("status"=>"EXACT_AFFINE_CURRENT_SLICE_OF_NETWORK_DAE",
         "rho"=>m.rho,"theta0"=>[R.N.trim_gfl(ctx,b).x[3] for b in 30:39],
         "buses"=>collect(30:39),"G"=>[collect(r) for r in eachrow(G)],
         "current_to_voltage"=>[collect(r) for r in eachrow(Z)],
         "candidate_sha256"=>bytes2hex(sha256(read(candidate))),
         "model_source_sha256"=>bytes2hex(sha256(read(joinpath(@__DIR__,"ReducedDAE.jl")))),
         "script_sha256"=>bytes2hex(sha256(read(@__FILE__))),
         "scope"=>"SG states fixed at nominal; GFL phase increments zero; currents vary within local certificates. A closure failure is not an instability proof.")
open(joinpath(OUT,"network_closure.toml"),"w") do io;TOML.print(io,out);end
println("EXPORTED network-current slice; inverse residual ",norm(G*Z+I,Inf))
