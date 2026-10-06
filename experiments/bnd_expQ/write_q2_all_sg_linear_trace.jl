using CSV, DataFrames, LinearAlgebra, TOML
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const OUT=joinpath(ROOT,"reports","experiment_Q","Q2_REFERENCE")
include(joinpath(ROOT,"src","bnd_design_p","ExpP.jl"))
include(joinpath(ROOT,"src","bnd_expQ","LinearSecurity.jl"))
using .ExpP, .LinearSecurity
const N=ExpP.PDExactDesignN

function main()
    cand=TOML.parsefile(joinpath(OUT,"Z_Q_ALL_SG_SECURE_REFERENCE.toml"))
    ctx=N.design_context(ROOT)
    rho=Float64.(cand["rho"]); kp=Float64.(cand["Kp"]); ki=Float64.(cand["Ki"])
    sp=N.spectrum(ctx,rho,kp,ki)
    model=reduced_step_model(ctx,sp.model,rho,16;gauge_vector=N.gauge_vector)
    eig=eigen(model.A); z=eig.vectors\model.B; cv=model.C*eig.vectors
    pd=CSV.read(joinpath(OUT,"TABLE_Q2_all_sg_reference_trajectory.csv"),DataFrame)
    times=Float64.(pd.time_s)
    sgcols=[Symbol("bus$(b)_SG_df_Hz") for b in 30:39]
    linear=zeros(length(times)); nonlinear=zeros(length(times))
    for j in eachindex(times)
        tau=max(times[j]-1.0,0.0)
        factors=[abs(lam)<1e-13 ? tau : expm1(lam*tau)/lam for lam in eig.values]
        y=real.(cv*(factors.*z))*100.0
        linear[j]=maximum(abs.(y))
        nonlinear[j]=maximum(abs(Float64(pd[j,c])) for c in sgcols)
    end
    out=DataFrame(time_s=times,linear_max_abs_local_SG_Hz=linear,
        PD_max_abs_local_SG_Hz=nonlinear)
    CSV.write(joinpath(OUT,"TABLE_Q2_all_sg_linear_vs_PD.csv"),out)
    println("samples=$(length(times)) linear_peak=$(maximum(linear)) PD_peak=$(maximum(nonlinear)) max_abs_difference=$(maximum(abs.(linear.-nonlinear)))")
end
main()
