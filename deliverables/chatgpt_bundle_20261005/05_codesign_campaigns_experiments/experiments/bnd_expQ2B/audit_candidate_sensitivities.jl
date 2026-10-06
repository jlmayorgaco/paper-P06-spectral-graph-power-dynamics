using CSV,DataFrames,LinearAlgebra,TOML
const ROOT=normpath(joinpath(@__DIR__,"..",".."));const CORE=joinpath(ROOT,"reports","experiment_Q2B","CORE")
include(joinpath(ROOT,"src","bnd_design_p","ExpP.jl"));include(joinpath(ROOT,"src","bnd_expQ","LinearSecurity.jl"));include(joinpath(ROOT,"src","bnd_expQ2B","FiniteWindow.jl"));include(joinpath(ROOT,"src","bnd_expQ2B","CoreDesign.jl"))
const N=ExpP.PDExactDesignN;const CD=CoreDesign;const BREQ=1.6991206999182038e-6
d=TOML.parsefile(joinpath(CORE,"Z_Q2B_SECURE_T05_FINAL.toml"));ctx=N.design_context(ROOT);s=Int.(d["support"])
eps=Float64.(d["epsilon"]);kp=Float64.(d["Kp"]);ki=Float64.(d["Ki"]);y=CD.encode(ctx,eps,kp,ki)
ev=CD.design_eval(ctx,s,eps,kp,ki;disturbance_MW=100.0,rocof_window=.5,dt_s=.05,horizon_s=30.0)
der=CD.constraint_jacobian(ctx,s,y,ev;rocof_window=.5,disturbance_MW=100.0,
    dt_s=.05,horizon_s=30.0,fd_step=2e-4,active_tol=.8)
P=Float64.(ctx.power);rows=NamedTuple[]
for (j,bus) in enumerate(s)
    dalpha=der.J[1,j]*.01
    dbeta=-der.J[2,j]*BREQ
    dfinf=der.J[3,j]*.5
    dfpeak=der.J[4,j]*.5
    drocof=der.J[5,j]*.5
    push!(rows,(;bus,Pgen0_MW=P[bus-29],
        d_alpha_d_epsilon_s_inv=dalpha,
        a_alpha_per_MW=-dalpha/P[bus-29],
        d_beta0_d_epsilon=dbeta,a_beta0_per_MW=dbeta/(BREQ*P[bus-29]),
        dFinf_d_epsilon_Hz=dfinf,a_Finf_per_MW=-dfinf/P[bus-29],
        dFpeak_d_epsilon_Hz=dfpeak,a_Fpeak_per_MW=-dfpeak/P[bus-29],
        dR05_d_epsilon_Hz_s=drocof,a_R05_per_MW=-drocof/P[bus-29],
        derivative_status="LOCAL_FINITE_DIFFERENCE; beta uses active omega=0; KKT multipliers unavailable"))
end
CSV.write(joinpath(CORE,"TABLE_Q2B_candidate_sensitivities.csv"),DataFrame(rows))
println("SENSITIVITY rows=",length(rows)," constraints=",ev.g,
    " active_gradient_rows=",der.active_rows," KKT_not_certified=true")
