using CSV, DataFrames, TOML, Statistics
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const OUT=joinpath(ROOT,"reports","experiment_Q2","F0")
include(joinpath(ROOT,"src","bnd_design_p","ExpP.jl"))
include(joinpath(ROOT,"src","bnd_expQ","LinearSecurity.jl"))
include(joinpath(ROOT,"src","bnd_expQ","FrequencyMetrics.jl"))
include(joinpath(ROOT,"src","bnd_expQ2","GridFrequency.jl"))
using .ExpP, .LinearSecurity, .FrequencyMetrics, .GridFrequency
const N=ExpP.PDExactDesignN

function main()
    ctx=N.design_context(ROOT);rho=zeros(10);kp=fill(N.K0P,10);ki=fill(N.K0I,10)
    m=N.descriptor(ctx,rho,kp,ki)
    pd=CSV.read(joinpath(ROOT,"reports","experiment_Q","Q2_REFERENCE",
        "TABLE_Q2_all_sg_reference_trajectory.csv"),DataFrame)
    rel=Float64.(pd.time_s).-1.0
    sig=grid_step_signals(ctx,m,rho,rel;load_bus=16,disturbance_MW=100.0,
        gauge_vector=N.gauge_vector)
    dt=0.01; h=5;degree=3
    lin=similar(sig.phase_rad)
    for k in axes(lin,1)
        lin[k,:].=sg_polynomial_derivative(view(sig.phase_rad,k,:),dt;
            half_window=h,degree)./(2pi)
    end
    rows=NamedTuple[]
    for bus in 30:39
        col=Symbol("bus$(bus)_bus_df_Hz")
        y=Float64.(pd[!,col]);k=bus-29
        inds=findall((Float64.(pd.time_s).>=1.1).&(Float64.(pd.time_s).<=60.0))
        err=lin[k,inds].-y[inds]
        push!(rows,(;bus,disturbance_MW=100.0,window_half_samples=h,
            support_s=2h*dt,minus3dB_Hz=savgol_cutoff_hz(dt;half_window=h,degree),
            linear_peak_Hz=maximum(abs.(lin[k,inds])),PD_peak_Hz=maximum(abs.(y[inds])),
            max_abs_error_Hz=maximum(abs.(err)),rmse_Hz=sqrt(mean(abs2,err)),
            n_samples=length(inds),linear_F_inf_Hz=sig.F_inf_Hz[k]))
    end
    CSV.write(joinpath(OUT,"TABLE_Q2_F0_allSG_linear_PD_bus_frequency.csv"),DataFrame(rows))
    println("allSG 100MW bus-frequency max error=$(maximum(x.max_abs_error_Hz for x in rows)) Hz")
end
main()
