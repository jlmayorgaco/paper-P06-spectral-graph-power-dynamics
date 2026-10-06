using CSV, DataFrames, LinearAlgebra, Statistics, TOML
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const OUT=joinpath(ROOT,"reports","experiment_Q2","F2");mkpath(OUT)
include(joinpath(ROOT,"src","bnd_design_p","ExpP.jl"))
include(joinpath(ROOT,"src","bnd_expQ","LinearSecurity.jl"))
include(joinpath(ROOT,"src","bnd_expQ","FrequencyMetrics.jl"))
include(joinpath(ROOT,"src","bnd_expQ2","GridFrequency.jl"))
const N=ExpP.PDExactDesignN

function arg_response(q,lam,D,t,d)
    d*(D+real(sum(q[j]*(abs(lam[j])<1e-13 ? t : expm1(lam[j]*t)/lam[j])
        for j in eachindex(lam))))
end
function rocof_response(q,lam,t,d)
    d*real(sum(q[j]*exp(lam[j]*t) for j in eachindex(lam)))
end
function golden_max(f,a,b;it=70)
    phi=(sqrt(5)-1)/2;c=b-phi*(b-a);e=a+phi*(b-a);fc=f(c);fe=f(e)
    for _ in 1:it
        if fc>fe;b=e;e=c;fe=fc;c=b-phi*(b-a);fc=f(c)
        else;a=c;c=e;fc=fe;e=a+phi*(b-a);fe=f(e);end
    end
    x=(a+b)/2;(x,f(x))
end
function peak_response(f,grid)
    vals=f.(grid); ix=Int[]
    for k in 2:length(grid)-1
        abs(vals[k])>=abs(vals[k-1]) && abs(vals[k])>=abs(vals[k+1]) && push!(ix,k)
    end
    best=(abs(vals[1]),grid[1],vals[1])
    for k in ix
        t,v=golden_max(x->abs(f(x)),grid[k-1],grid[k+1])
        abs(v)>best[1] && (best=(abs(v),t,f(t)))
    end
    abs(vals[end])>best[1] && (best=(abs(vals[end]),grid[end],vals[end]))
    best
end
function main()
    ctx=ExpP.PDExactDesignN.design_context(ROOT)
    c=TOML.parsefile(joinpath(ROOT,"reports","experiment_N","Z_N_NOMINAL_FINAL.toml"))
    rho=Float64.(c["rho"]);kp=Float64.(c["Kp"]);ki=Float64.(c["Ki"])
    m=N.descriptor(ctx,rho,kp,ki)
    red=LinearSecurity.reduced_step_model(ctx,m,rho,16;gauge_vector=N.gauge_vector)
    F=eigen(red.A);lam=F.values;V=F.vectors;W=inv(V)
    d=100.0; dt=0.01;grid=collect(0.0:dt:60.0)
    Q=Matrix{ComplexF64}(undef,length(lam),10)
    for b in 1:10
        Q[:,b].=vec(red.C_bus[b,:]'*V).*vec(W*red.B)
    end
    metric_times=collect(-0.5:dt:60.0)
    metric_sig=GridFrequency.grid_step_signals(ctx,m,rho,metric_times;load_bus=16,
        disturbance_MW=d,gauge_vector=N.gauge_vector)
    bandwidth=GridFrequency.grid_step_metrics(metric_sig;half_windows=(0,2,5,10,20),
        dt_s=dt,event_time_s=0.0)
    CSV.write(joinpath(OUT,"TABLE_Q2_F2_ExpN_bandwidth_sensitivity.csv"),DataFrame(bandwidth))
    freq_peak=(0.0,NaN,NaN);rocof_peak=(0.0,NaN,NaN);steady_peak=(0.0,NaN,NaN)
    for b in 1:10
        q=Q[:,b];df=real(red.D_bus[b])
        fp=peak_response(t->arg_response(q,lam,df,t,d),grid)
        rp=peak_response(t->rocof_response(q,lam,t,d),grid)
        ss=d*real(df-sum(q./lam))
        fp[1]>freq_peak[1] && (freq_peak=(fp[1],29+b,fp[2]))
        rp[1]>rocof_peak[1] && (rocof_peak=(rp[1],29+b,rp[2]))
        abs(ss)>steady_peak[1] && (steady_peak=(abs(ss),29+b,ss))
    end
    checktimes=[0.0,0.01,0.1,1.0,5.0,60.0]
    direct=GridFrequency.grid_step_signals(ctx,m,rho,checktimes;load_bus=16,
        disturbance_MW=d,gauge_vector=N.gauge_vector)
    recerr=0.0;rocof0err=0.0
    for b in 1:10
        q=Q[:,b]
        rocof0err=max(rocof0err,abs(d*real(sum(q))-d*real((red.C_bus[b,:]'*red.B)[1])))
        for (k,t) in enumerate(checktimes)
            fmodal=arg_response(q,lam,real(red.D_bus[b]),t,d)
            recerr=max(recerr,abs(fmodal-direct.f_regular_Hz[b,k]))
        end
    end
    bR=Int(rocof_peak[2])-29;bF=Int(freq_peak[2])-29;bS=Int(steady_peak[2])-29
    rows=NamedTuple[]
    for j in eachindex(lam)
        lamj=lam[j];qR=Q[j,bR];qF=Q[j,bF];qS=Q[j,bS]
        damp=abs(lamj)>0 ? -real(lamj)/abs(lamj) : NaN
        cR=d*real(qR*exp(lamj*rocof_peak[3]))
        cF=d*real(qF*(abs(lamj)<1e-13 ? freq_peak[3] : expm1(lamj*freq_peak[3])/lamj))
        cS=d*real(-qS/lamj)
        push!(rows,(;pole_index=j,lambda_real_per_s=real(lamj),lambda_imag_per_s=imag(lamj),
            frequency_Hz=abs(imag(lamj))/(2pi),damping_ratio=damp,
            residue_rocof_magnitude_Hz_s_per_MW=abs(qR),
            residue_rocof_phase_deg=rad2deg(angle(qR)),
            residue_frequency_magnitude_Hz_per_MW=abs(qF),
            initial_rocof_contribution_Hz_s=cR,
            frequency_peak_contribution_Hz=cF,
            steady_offset_contribution_Hz=cS,
            is_spectral_abscissa_mode=j==argmax(real.(lam))))
    end
    CSV.write(joinpath(OUT,"TABLE_Q2_F2_modal_frequency_contributions.csv"),DataFrame(rows))
    result=Dict("stage"=>"F2","candidate"=>"ExpN nominal, bus 16, +100 MW sustained",
        "complete_physical_modes"=>length(lam),"alpha_s_inv"=>maximum(real.(lam)),
        "rightmost_pole_real_s_inv"=>real(lam[argmax(real.(lam))]),
        "grid_frequency_peak_Hz"=>freq_peak[1],"grid_frequency_peak_bus"=>freq_peak[2],
        "grid_frequency_peak_time_s_after_step"=>freq_peak[3],
        "smooth_grid_RoCoF_peak_Hz_s"=>rocof_peak[1],"RoCoF_peak_bus"=>rocof_peak[2],
        "RoCoF_peak_time_s_after_step"=>rocof_peak[3],"ideal_unfiltered_rocof_impulse"=>true,
        "phase_jump_max_rad_at_100MW"=>maximum(abs,direct.phase_jump_rad),
        "phase_jump_bus"=>29+argmax(abs.(direct.phase_jump_rad)),
        "phase_jump_all_buses_rad"=>direct.phase_jump_rad,
        "RoCoF_measurement_bandwidth_status"=>"METRIC_DEPENDENT",
        "measurement_bandwidth_sensitivity_rows"=>length(bandwidth),
        "steady_offset_peak_Hz"=>steady_peak[1],"steady_offset_bus"=>steady_peak[2],
        "feedthrough_Dbus_Hz_per_MW"=>maximum(abs,red.D_bus),
        "F_peak_modal_expansion_reconstruction_error_Hz"=>recerr,
        "eigenvector_condition"=>cond(V),"smooth_RoCoF_t0_max_abs_Hz_s"=>
            maximum(abs,d .* (red.C_bus*red.B)),
        "RoCoF_residue_t0_reconstruction_error_Hz_s"=>rocof0err,
        "modal_expansion_note"=>"finite-mode state response; direct algebraic feedthrough reported separately; ideal phase jump creates distributional unfiltered RoCoF")
    open(joinpath(OUT,"F2_RESULTS.toml"),"w") do io;TOML.print(io,result);end
    write(joinpath(OUT,"STAGE_SUMMARY.md"),"""
    # F2 — ExpN grid-frequency modal residues

    The full physical quotient spectrum of the frozen ExpN point was expanded at buses 30–39. Finite pole residues use q=(C r)(lᴴB)/(lᴴr); the direct algebraic feedthrough D is separately retained. The phase jump from the step produces a distributional component and is not represented by the finite-pole residue sum.

    - Smooth filtered-free grid frequency peak (state response plus D): $(freq_peak[1]) Hz at bus $(freq_peak[2]), t=$(freq_peak[3]) s after the step.
    - Smooth finite-state RoCoF peak: $(rocof_peak[1]) Hz/s at bus $(rocof_peak[2]), t=$(rocof_peak[3]) s; ideal unfiltered RoCoF also has an impulse because phase jumps.
    - Steady offset magnitude: $(steady_peak[1]) Hz at bus $(steady_peak[2]).
    - Active finite mode: $(lam[argmax(real.(lam))]). Residue condition number: $(cond(V)); six-time modal reconstruction error: $(recerr) Hz; zero-time RoCoF reconstruction error: $(rocof0err) Hz/s.
    - Contribution rows contain per-mode residue and contributions evaluated at each metric's independently determined maximum.
    """)
    println("F2_ALPHA=",maximum(real.(lam))," Fpeak=",freq_peak," smooth_Rpeak=",rocof_peak," Finf=",steady_peak)
end
main()
