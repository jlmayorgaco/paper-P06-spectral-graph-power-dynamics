include("ReducedDAE.jl")
using .ReducedDAE, LinearAlgebra, TOML, SHA, Random
const R=ReducedDAE
const OUT=joinpath(R.OUT,"full_gfl_sector_contract")
mkpath(OUT)
ctx=R.N.design_context(R.ROOT)
devices=Dict[]
for bus in 30:39
    op=R.N.trim_gfl(ctx,bus)
    p=op.pars;V0=norm(op.u);id0=op.target_p/V0;iq0=-op.target_q/V0
    vid0=p.cc_ki*op.x[2];viq0=p.cc_ki*op.x[1]
    residual=norm(R.gfl_rhs(op.x,op.u,p,R.N.K0P,R.N.K0I,:physical_supply),Inf)
    @assert residual<1e-8
    push!(devices,Dict("bus"=>bus,"V0"=>V0,"id0"=>id0,"iq0"=>iq0,
        "vid0"=>vid0,"viq0"=>viq0,"trim_residual"=>residual,
        "parameters"=>Dict(String(k)=>v for (k,v) in pairs(p)),
        "Kp_nominal"=>R.N.K0P,"Ki_nominal"=>R.N.K0I,
        "x0"=>op.x,"u0"=>op.u))
end
files=[joinpath(@__DIR__,"ReducedDAE.jl"),joinpath(R.ROOT,"src","bnd_model_expN","PDExactDesignN.jl")]
out=Dict("status"=>"FROZEN_FULL_GFL_MODELS_EXPORTED",
         "devices"=>devices,
         "source_hashes"=>Dict(relpath(f,R.ROOT)=>bytes2hex(sha256(read(f))) for f in files),
         "script_sha256"=>bytes2hex(sha256(read(@__FILE__))))
open(joinpath(OUT,"models.toml"),"w") do io;TOML.print(io,out);end
println("EXPORTED ",length(devices)," models; maximum trim residual ",maximum(d["trim_residual"] for d in devices))
