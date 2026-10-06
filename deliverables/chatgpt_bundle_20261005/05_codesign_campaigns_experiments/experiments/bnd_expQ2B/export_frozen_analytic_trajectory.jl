using CSV,DataFrames,LinearAlgebra,SHA,TOML
const ROOT=normpath(joinpath(@__DIR__,"..",".."));const CORE=joinpath(ROOT,"reports","experiment_Q2B","CORE")
const F=joinpath(CORE,"Z_Q2B_SECURE_T05_FINAL.toml");const EXPECTED="f53250383cdb906a2d0cb3130ca9c1e72c63b472717c5c82aea1d5838ac785ce"
bytes2hex(sha256(read(F)))==EXPECTED||error("candidate SHA mismatch")
include(joinpath(ROOT,"src","bnd_design_p","ExpP.jl"));include(joinpath(ROOT,"src","bnd_expQ","LinearSecurity.jl"));include(joinpath(ROOT,"src","bnd_expQ2B","FiniteWindow.jl"))
const N=ExpP.PDExactDesignN
c=TOML.parsefile(F);rho=Float64.(c["rho"]);eps=Float64.(c["epsilon"]);kp=Float64.(c["Kp"]);ki=Float64.(c["Ki"]);ctx=N.design_context(ROOT)
m=N.descriptor(ctx,rho,kp,ki);dt=.01;times=collect(-1.0:dt:60.0)
sig=FiniteWindow.phase_step(ctx,m,rho,times;load_bus=16,disturbance_MW=100.0,gauge_vector=N.gauge_vector)
T=.5;lag=Int(round(T/dt));zero=findfirst(==(0.0),times);f=zeros(10,length(times));r=zeros(10,length(times))
for j in zero:length(times)
    jm=j-lag;jm2=j-2lag
    thm=jm>=zero ? sig.phase_rad[:,jm] : zeros(10)
    thm2=jm2>=zero ? sig.phase_rad[:,jm2] : zeros(10)
    f[:,j]=(sig.phase_rad[:,j]-thm)/(2pi*T)
    r[:,j]=(sig.phase_rad[:,j]-2thm+thm2)/(2pi*T^2)
end
df=DataFrame(event_time_s=times)
for k in 1:10
    bus=29+k;df[!,Symbol("phase_bus$(bus)_rad")]=collect(view(sig.phase_rad,k,:))
    df[!,Symbol("fGrid_bus$(bus)_T05_Hz")]=collect(view(f,k,:))
    df[!,Symbol("RGrid_bus$(bus)_T05_Hz_s")]=collect(view(r,k,:))
end
df[!,:max_abs_fGrid_T05_Hz]=[maximum(abs.(f[:,j])) for j in eachindex(times)]
df[!,:max_abs_RGrid_T05_Hz_s]=[maximum(abs.(r[:,j])) for j in eachindex(times)]
CSV.write(joinpath(CORE,"TABLE_Q2B_analytic_TDS_timeseries.csv"),df)
CSV.write(joinpath(CORE,"TABLE_Q2B_candidate_by_bus.csv"),DataFrame(
    bus=collect(30:39),Pgen0_MW=ctx.power,epsilon=eps,rho=rho,Kp=kp,Ki=ki,
    retained_SG_MW=ctx.power.*eps,GFL_MW=ctx.power.*rho))
println("ANALYTIC_TRAJECTORY rows=",nrow(df)," total_SG_MW=",sum(ctx.power.*eps),
    " GFL_MW=",sum(ctx.power.*rho)," max_f=",maximum(df.max_abs_fGrid_T05_Hz),
    " max_R=",maximum(df.max_abs_RGrid_T05_Hz_s))
