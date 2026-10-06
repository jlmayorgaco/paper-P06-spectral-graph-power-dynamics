using CSV, DataFrames, DelimitedFiles, LinearAlgebra, TOML, SHA
include(joinpath(@__DIR__,"..","..","src","bnd_design_k","BNDDesignK.jl"))
using .BNDDesignK

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const CASE=length(ARGS)==1 ? ARGS[1] : error("pass case")
const OUT=joinpath(ROOT,"reports","experiment_M","matrices",CASE)
mkpath(OUT)

function frozen(path)
    bytes2hex(sha256(read(path)))==strip(read(path*".sha256",String)) || error("candidate hash mismatch")
    TOML.parsefile(path)
end

ctx=BNDDesignK.design_context(ROOT)
rho=zeros(10);kp=copy(ctx.kp0);ki=copy(ctx.ki0)
if CASE=="ExpG_candidate"
    c=frozen(joinpath(ROOT,"reports","experiment_G","Z_G_FINAL.toml"))
    rows=sort(c["generator"];by=x->Int(x["bus"]))
    rho=Float64.([r["rho"] for r in rows])
    kp=Float64.([r["Kp"] for r in rows])
    ki=Float64.([r["Ki"] for r in rows])
elseif CASE=="ExpK_nominal"
    c=frozen(joinpath(ROOT,"reports","experiment_K","Z_K_NOMINAL_FINAL.toml"))
    rho=Float64.(c["rho"]);kp=Float64.(c["Kp"]);ki=Float64.(c["Ki"])
elseif CASE!="all_SG"
    error("unknown case")
end

cm=BNDDesignK.BNDDesignG.CollectiveModel
mod=cm.mixed_jacobian(ctx.net,rho,kp,ki)
nx=mod.n_dynamic;ny=mod.n_algebraic
E=zeros(nx+ny,nx+ny);E[1:nx,1:nx].=Matrix{Float64}(I,nx,nx)
A=[mod.A mod.B;mod.C mod.Gy]
for (name,m) in (("AN_M",E),("AN_A",A),("AN_Ared",mod.Ared),
    ("AN_A_local",mod.A),("AN_B_port",mod.B),("AN_C_port",mod.C),
    ("AN_D_port",mod.D),("AN_Ystatic",ctx.net.y_static),("AN_Gy",mod.Gy))
    writedlm(joinpath(OUT,name*".csv"),m,',')
end

state_rows=NamedTuple[]
for row in eachrow(mod.state_map)
    names=row.kind=="SG" ? (Int(row.bus)==39 ? cm.AnalyticSG.uncontrolled_state_names() :
        cm.AnalyticSG.state_names()) : cm.AnalyticGFLPLL.state_names()
    x=row.kind=="SG" ? ctx.net.sg[Int(row.bus)].op.x : ctx.net.gfl[Int(row.bus)].op.x
    length(names)==Int(row.last)-Int(row.first)+1 || error("local state count mismatch")
    for (k,name) in enumerate(names)
        push!(state_rows,(index=Int(row.first)+k-1,bus=Int(row.bus),kind=String(row.kind),
            state_name=String(name),differential=true,equilibrium_value=Float64(x[k])))
    end
end
for bus in 1:39, coord in ("u_r","u_i")
    push!(state_rows,(index=nx+2bus-(coord=="u_r" ? 1 : 0),bus,kind="network",
        state_name="busbar_"*coord,differential=false,
        equilibrium_value=coord=="u_r" ? real(ctx.net.voltage[bus]) : imag(ctx.net.voltage[bus])))
end
CSV.write(joinpath(OUT,"AN_state_map.csv"),DataFrame(state_rows))
ev=BNDDesignK.evaluate(ctx,1 .-rho,kp,ki)
CSV.write(joinpath(OUT,"AN_poles.csv"),DataFrame(real=real.(ev.lambda),imag=imag.(ev.lambda)))
println("AN_CASE_DONE ",CASE," n=",nx+ny," dynamic=",nx," algebraic=",ny,
    " alpha=",ev.alpha)
