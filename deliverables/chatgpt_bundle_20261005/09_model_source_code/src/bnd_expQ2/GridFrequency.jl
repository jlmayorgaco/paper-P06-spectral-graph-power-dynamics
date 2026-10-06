module GridFrequency

using LinearAlgebra
import ..LinearSecurity
import ..FrequencyMetrics

export grid_step_signals, grid_step_metrics, continuity_rows

"""Architecture-invariant algebraic grid-voltage phase measurement model.

The reduced plant is the ExpN DAE closure. Outputs are bus-voltage phase
derivatives at buses 30:39. `f_regular` is the finite t>0 derivative for a
constant step. `phase_jump_rad` records the algebraic voltage-angle jump at
the step; ideal unfiltered RoCoF therefore contains an impulse when this is
nonzero. Sampled/filter metrics are diagnostics with their estimator stated.
"""
function grid_step_signals(ctx, m, rho, times_s; load_bus=16, disturbance_MW=1.0,
                           gauge_vector, system_base_mva=100.0)
    red=LinearSecurity.reduced_step_model(ctx,m,rho,load_bus;
        gauge_vector,system_base_mva)
    eig=eigen(red.A)
    λ=eig.values; V=eig.vectors; z=V\red.B
    CV=red.C_bus*V
    nbus=10; nt=length(times_s)
    f_regular=zeros(Float64,nbus,nt)
    phase=zeros(Float64,nbus,nt)
    yload=LinearSecurity.load_input_vector(ctx,m,load_bus;system_base_mva)
    dv_jump=-(m.Gy\yload).*disturbance_MW
    phase_jump=zeros(Float64,nbus)
    for bus in 30:39
        k=bus-29; v=ctx.net.voltage[bus]
        ar=zeros(1,size(m.Gy,1))
        ar[1,2bus-1]=-imag(v)/abs2(v)
        ar[1,2bus]= real(v)/abs2(v)
        phase_jump[k]=(ar*dv_jump)[1]
    end
    phase_jump ./= 2pi
    d=Float64(disturbance_MW)
    finf=real.(red.C_bus*(-(red.A\red.B)) .+ red.D_bus).*d
    for (j,t0) in enumerate(times_s)
        t=Float64(t0)
        if t < 0
            continue
        end
        stepfac=[abs(x)<1e-13 ? t : expm1(x*t)/x for x in λ]
        intfac=[abs(x)<1e-13 ? t^2/2 : (expm1(x*t)-x*t)/x^2 for x in λ]
        f_regular[:,j].=real.(CV*(stepfac.*z)).*d .+ red.D_bus.*d
        phase[:,j].=phase_jump .+ 2pi.*(real.(CV*(intfac.*z)).*d .+
            red.D_bus.*d.*t)
    end
    (;times_s=Float64.(times_s),phase_rad=phase,f_regular_Hz=f_regular,
      phase_jump_rad=phase_jump,A=red.A,B=red.B,C_bus=red.C_bus,D_bus=red.D_bus,
      poles=λ,alpha=maximum(real.(λ)),condition_A=cond(red.A),
      eigenvector_condition=cond(V),reduction_gauge_residual=red.gauge_residual,
      F_inf_Hz=finf)
end

function _sampled_derivatives(x,dt)
    n=size(x,1); nt=size(x,2)
    f=zeros(Float64,n,nt); r=zeros(Float64,n,nt)
    for k in 1:n
        f[k,2:end].=diff(view(x,k,:))./(2pi*dt)
        f[k,1]=f[k,2]
        r[k,2:end].=diff(view(f,k,:))./dt
        r[k,1]=r[k,2]
    end
    f,r
end

"""Peak/steady metrics at several fixed offline phase-differentiator windows.

`half_windows=0` means unfiltered sampled first/second differences. Positive
values apply the same cubic local-polynomial derivative to phase and then to
frequency (the latter for RoCoF). The step's phase jump remains present.
"""
function grid_step_metrics(signals; half_windows=(0,2,5,10,20), dt_s=0.01,
                           event_time_s=0.0)
    t=signals.times_s; active=findall(t .>= event_time_s)
    rows=NamedTuple[]
    for h in half_windows
        if h==0
            f,r=_sampled_derivatives(signals.phase_rad,dt_s)
            cutoff=NaN; filter_support=0.0; label="UNFILTERED_100HZ_FINITE_DIFFERENCE"
        else
            f=similar(signals.phase_rad);r=similar(signals.phase_rad)
            degree=min(3,2h)
            for k in axes(f,1)
                f[k,:].=FrequencyMetrics.sg_polynomial_derivative(
                    view(signals.phase_rad,k,:),dt_s;half_window=h,degree)./(2pi)
                r[k,:].=FrequencyMetrics.sg_polynomial_derivative(
                    view(f,k,:),dt_s;half_window=h,degree)
            end
            cutoff=FrequencyMetrics.savgol_cutoff_hz(dt_s;half_window=h,degree)
            filter_support=2h*dt_s;label="SAVGOL_HALF_WINDOW_$(h)"
        end
        fpeak=maximum(abs,view(f,:,active))
        rpeak=maximum(abs,view(r,:,active))
        fss=maximum(abs,signals.F_inf_Hz)
        push!(rows,(;estimator=label,half_window=h,window_support_s=filter_support,
            derivative_minus3dB_Hz=cutoff,F_peak_Hz=fpeak,R_peak_Hz_s=rpeak,
            F_inf_Hz=fss,phase_jump_peak_rad=maximum(abs,signals.phase_jump_rad),
            includes_event_phase_jump=true,unfiltered_RoCoF_impulse_if_jump=true))
    end
    rows
end

function continuity_rows(bus, rho, ref_metrics, trial_metrics; condition, ref_condition,
                         rho_value, alpha)
    rows=NamedTuple[]
    scale=max(1.0,maximum(abs,Float64.(values(ref_metrics))))
    tol=max(1e-8*scale,100eps(Float64)*max(condition,ref_condition)*scale)
    for key in keys(ref_metrics)
        ref=Float64(ref_metrics[key]); val=Float64(trial_metrics[key])
        err=abs(val-ref)
        push!(rows,(;bus,rho_GFL=rho_value,epsilon_SG=1-rho_value,
            estimator=String(key),metric_reference=ref,metric_rho=val,
            absolute_error=err,tolerance_condition_aware=tol,
            pass_within_tolerance=isfinite(err)&&err<=tol,
            alpha_per_s=alpha,condition_A=condition))
    end
    rows
end

end
