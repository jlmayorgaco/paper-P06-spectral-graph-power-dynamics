module LinearPulse

using LinearAlgebra, CSV, DataFrames
using ..ExpP
const N=ExpP.PDExactDesignN

export compare_frozen_pulses

function _linear_step(ctx,rho,kp,ki,bus;duration=60.0,dt=0.1,base_mva=100.0,f0=60.0)
    sp=N.spectrum(ctx,rho,kp,ki);m=sp.model;Q=sp.quotient
    A=transpose(Q)*m.Ared*Q
    Cf=zeros(size(m.Ared,1));h_total=0.0
    for b in 30:39
        eps=1-rho[b-29];eps>0 || continue
        r=findfirst(x->Int(x.bus)==b && x.kind=="SG",eachrow(m.state_map))
        r===nothing && continue
        row=m.state_map[r,:]
        state_index=Int(row.last)-1
        par=ctx.net.sg[b].op.parameters
        h=eps*par.inertia*par.rating_mva
        Cf[state_index]=f0*h
        h_total+=h
    end
    h_total>0 || error("COI frequency output requested without retained SG")
    Cf./=h_total

    loads=CSV.read(joinpath(ctx.root,"reports","experiment_D","inputs","load.csv"),DataFrame)
    loadrow=only(eachrow(loads[loads.bus.==bus,:]))
    p0=Float64(loadrow.Pset);q0=Float64(loadrow.Qset)
    qratio=abs(q0/p0)
    V=ctx.net.voltage[bus]
    Iunit=conj((1+im*qratio)/(base_mva*V))
    iv=zeros(Float64,size(m.Gy,1));iv[2bus-1]=real(Iunit);iv[2bus]=imag(Iunit)
    Bd=m.B*(m.Gy\iv)
    Aq=transpose(Q)*m.Ared*Q
    Bq=transpose(Q)*Bd;Cq=transpose(Cf)*Q
    M=exp(Aq*dt)
    step_increment=Aq\((M-Matrix{Float64}(I,size(Aq,1),size(Aq,2)))*Bq)
    ts=collect(0.0:dt:duration);xs=zeros(Float64,length(Bq));out=zeros(Float64,length(ts))
    for k in 2:length(ts)
        xs=M*xs+step_increment
        out[k]=dot(Cq,xs)
    end
    (;times=ts,step_frequency_per_MW=out,load_p_MW=-p0*base_mva,
      load_q_Mvar=-q0*base_mva,q_to_p_ratio=qratio,
      alpha=sp.alpha,physical_poles=length(sp.lambda),input_complex_current_pu=Iunit)
end

"""Compare frozen PD nonlinear pulses with the exact-closure linear response.

The same 0.1 s pulse support `[1.0,1.1)` and P/Q load ratio are used. The
linear response is a difference of two step responses, sampled at TDS output
times; RoCoF is computed with the same 0.1 s backward difference.
"""
function compare_frozen_pulses(ctx,rho,kp,ki,tds_frames;
        buses=(8,16,29),fractions=(1e-5,2e-5),start=1.0,stop=1.1,
        duration=60.0,dt=0.1)
    rows=NamedTuple[]
    for bus in buses
        lin=_linear_step(ctx,rho,kp,ki,bus;duration,dt)
        for frac in fractions
            actual=tds_frames[(tds_frames.event_bus.==bus).&
                              (tds_frames.pulse_fraction.==frac),:]
            nrow(actual)>0 || error("missing nonlinear TDS trace bus=$bus fraction=$frac")
            sort!(actual,:time_s)
            times=Float64.(actual.time_s)
            ΔP=lin.load_p_MW*frac
            function step_at(t)
                t<0 && return 0.0
                idx=round(Int,t/dt)+1
                abs(times[clamp(idx,1,length(times))]-t)<dt/100 ||
                    error("linear and TDS sampling grids differ")
                lin.step_frequency_per_MW[clamp(idx,1,length(times))]
            end
            flinear=[ΔP*(step_at(t-start)-step_at(t-stop)) for t in times]
            rlinear=vcat(0.0,diff(flinear)./dt)
            fnl=Float64.(actual.COI_frequency_deviation_Hz)
            rnl=Float64.(actual.RoCoF_Hz_s)
            ferr=maximum(abs.(flinear-fnl))/max(maximum(abs.(fnl)),1e-12)
            rerr=maximum(abs.(rlinear-rnl))/max(maximum(abs.(rnl)),1e-12)
            push!(rows,(;event_bus=bus,pulse_fraction=frac,delta_P_MW=ΔP,
              P_load_MW=lin.load_p_MW,Q_load_Mvar=lin.load_q_Mvar,
              q_to_p_ratio=lin.q_to_p_ratio,pulse_start_s=start,pulse_stop_s=stop,
              peak_frequency_linear_Hz=maximum(abs.(flinear)),
              peak_frequency_nonlinear_Hz=maximum(abs.(fnl)),
              peak_rocof_linear_Hz_s=maximum(abs.(rlinear)),
              peak_rocof_nonlinear_Hz_s=maximum(abs.(rnl)),
              frequency_relative_error=ferr,rocof_relative_error=rerr,
              linear_alpha=lin.alpha,physical_poles=lin.physical_poles,
              match_pass=ferr<0.05 && rerr<0.05))
        end
    end
    DataFrame(rows)
end

end
