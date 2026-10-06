module FiniteWindow

using LinearAlgebra
import ..LinearSecurity

export phase_step, window_metrics, design_metrics

"""Return the phase deviation at generator-bus voltages for a bus-load step.

The phase signal is computed from the ExpN descriptor's exact algebraic
voltage jump and the full physical finite-pole response after the validated
rotational quotient. Values before the event are zero.
"""
function phase_step(ctx, m, rho, times_s; load_bus=16, disturbance_MW=100.0,
                    gauge_vector, system_base_mva=100.0)
    red=LinearSecurity.reduced_step_model(ctx,m,rho,load_bus;
        gauge_vector,system_base_mva)
    F=eigen(red.A); λ=F.values; V=F.vectors; z=V\red.B
    CV=red.C_bus*V
    yload=LinearSecurity.load_input_vector(ctx,m,load_bus;system_base_mva)
    dvjump=-(m.Gy\yload).*disturbance_MW
    phasejump=zeros(10)
    for bus in 30:39
        k=bus-29; v=ctx.net.voltage[bus]
        angrow=zeros(1,size(m.Gy,1))
        angrow[1,2bus-1]=-imag(v)/abs2(v)
        angrow[1,2bus]=real(v)/abs2(v)
        phasejump[k]=(angrow*dvjump)[1] # voltage-angle differential is in radians
    end
    ts=Float64.(times_s); phase=zeros(10,length(ts)); d=Float64(disturbance_MW)
    active=findall(ts .>= 0)
    if !isempty(active)
        ta=ts[active]
        intfac=similar(ComplexF64.(λ .* ta'))
        for j in eachindex(λ),k in eachindex(ta)
            x=λ[j]*ta[k]
            intfac[j,k]=abs(x)<1e-7 ? ta[k]^2/2 + λ[j]*ta[k]^3/6 :
                (expm1(x)-x)/λ[j]^2
        end
        integ=real.(CV*(intfac .* reshape(z,:,1))).*d
        phase[:,active].=phasejump .+ 2pi.*(integ .+ red.D_bus.*(d.*reshape(ta,1,:)))
    end
    (;times_s=ts,phase_rad=phase,phase_jump_rad=phasejump,
      F_inf_Hz=real.(red.C_bus*(-(red.A\red.B)).+red.D_bus).*d,
      poles=λ,alpha=maximum(real.(λ)),A=red.A,B=red.B,C_bus=red.C_bus,
      D_bus=red.D_bus,condition_A=cond(red.A),eigenvector_condition=cond(V),
      quotient=red)
end

"""Apply the causal phase-difference frequency and RoCoF estimators.

The evaluation grid must include t=0 and have `T/dt` integral. The routine
includes one-sided samples at moving-window break points in its peak scan.
"""
function window_metrics(sig; windows=(0.2,0.5,1.0,2.0), dt_s=0.02,
                        horizon_s=40.0, disturbance_MW=100.0)
    rows=NamedTuple[]; t=sig.times_s; n=length(t); dt=Float64(dt_s)
    zero=findfirst(==(0.0),t); zero===nothing && error("grid must contain t=0")
    d=Float64(disturbance_MW)
    for T0 in windows
        T=Float64(T0); lagf=T/dt; abs(lagf-round(lagf))<1e-8 ||
            error("window must align to dt: T=$T, dt=$dt")
        lag=Int(round(lagf)); f=zeros(10,n); r=zeros(10,n)
        for k in zero:n
            km=k-lag; km2=k-2lag
            θm=km>=zero ? sig.phase_rad[:,km] : zeros(10)
            θm2=km2>=zero ? sig.phase_rad[:,km2] : zeros(10)
            f[:,k].=(sig.phase_rad[:,k]-θm)./(2pi*T)
            r[:,k].=(sig.phase_rad[:,k]-2θm+θm2)./(2pi*T^2)
        end
        # Steady-state is a required candidate extremum, even if the finite
        # simulation horizon has not fully reached it.
        fs=maximum(abs.(sig.F_inf_Hz))
        fi=maximum(abs.(f[:,zero:end])); ir=argmax(abs.(vec(f[:,zero:end])))
        rb=maximum(abs.(r[:,zero:end])); rr=argmax(abs.(vec(r[:,zero:end])))
        fp=max(fi,fs)
        push!(rows,(;window_s=T,F_inf_Hz=fs,F_peak_Hz=fp,
          F_peak_sample_Hz=fi,F_peak_bus=29+mod1(ir,10),
          F_peak_time_s=t[zero+div(ir-1,10)],
          R_peak_Hz_s=rb,R_peak_bus=29+mod1(rr,10),
          R_peak_time_s=t[zero+div(rr-1,10)],dt_s=dt,
          horizon_s=Float64(horizon_s),disturbance_MW=d,
          estimator="CAUSAL_UNWRAPPED_PHASE_DIFFERENCE"))
    end
    rows
end

function design_metrics(ctx,m,rho;load_bus=16,disturbance_MW=100.0,
                        windows=(0.2,0.5,1.0,2.0),dt_s=0.02,
                        horizon_s=40.0,gauge_vector)
    Tmax=maximum(windows); dt=Float64(dt_s)
    nlag=Int(round(2Tmax/dt)); ntime=Int(round(horizon_s/dt))
    times=collect(-nlag:ntime).*dt
    sig=phase_step(ctx,m,rho,times;load_bus,disturbance_MW,gauge_vector)
    (;signals=sig,metrics=window_metrics(sig;windows,dt_s,horizon_s,disturbance_MW))
end

end
