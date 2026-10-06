module LinearSecurity

using LinearAlgebra, DataFrames, CSV

export reduced_step_model, frequency_outputs, load_input_vector,
       step_response_eigensystem, step_metrics, frequency_authority_from_governors

function _quotient(m, gauge_vector)
    g=gauge_vector(m)
    norm(m.Ared*g)/max(norm(m.Ared)*norm(g),eps()) < 1e-8 ||
        error("gauge vector failed the DAE reduction residual")
    Q=nullspace(reshape(g/norm(g),1,:))
    A=transpose(Q)*m.Ared*Q
    (;Q,A,gauge_residual=norm(m.Ared*g)/max(norm(m.Ared)*norm(g),eps()))
end

function frequency_outputs(m, rho; f0_hz=60.0)
    inv=m.state_inventory
    rows=NamedTuple[]
    C=zeros(Float64,0,size(m.Ared,1))
    function append_row(bus,kind,state_name,scale)
        mask=(Int.(inv.bus).==bus).&(String.(inv.kind).==kind).&
             (String.(inv.state_name).==state_name)
        idx=findall(mask)
        length(idx)==1 || error("state $state_name at $kind bus $bus is not unique")
        row=zeros(Float64,1,size(m.Ared,1))
        row[1,idx[1]]=scale
        C=vcat(C,row)
        push!(rows,(;bus,kind,state_name,scale))
        nothing
    end
    for bus in 30:39
        i=bus-29
        if rho[i] < 1
            append_row(bus,"SG","machine_omega",Float64(f0_hz))
        end
        if rho[i] > 0
            append_row(bus,"GFL","delta_omega_rad_s",1/(2pi))
        end
    end
    (;C,metadata=DataFrame(rows))
end

"""Realified load injection vector for a +1 MW constant-impedance P step."""
function load_input_vector(ctx,m,bus; system_base_mva=100.0)
    1 <= bus <= length(ctx.net.voltage) || throw(BoundsError(ctx.net.voltage,bus))
    v=ctx.net.voltage[bus]
    dy=-1/(Float64(system_base_mva)*abs2(v))
    y=zeros(Float64,size(m.Gy,1))
    y[2bus-1]=real(dy*v)
    y[2bus]=imag(dy*v)
    y
end

"""
Construct the exact fixed-architecture reduced linear load-to-frequency map.

The nodal equation is `G_y δv + C δx + δY_load v0 d = 0`; this routine
therefore uses the algebraic Schur complement rather than a fitted transfer
function. A unit input is one MW of increased constant-impedance load.
"""
function reduced_step_model(ctx,m,rho,load_bus; gauge_vector,
                             system_base_mva=100.0, f0_hz=60.0)
    q=_quotient(m,gauge_vector)
    yload=load_input_vector(ctx,m,load_bus;system_base_mva)
    Bfull= -m.B*(m.Gy\yload)
    Bq=transpose(q.Q)*Bfull
    outputs=frequency_outputs(m,rho;f0_hz)
    Cq=outputs.C*q.Q
    Cdot=-(m.Gy\(m.C*m.Ared))
    Ddot=-(m.Gy\(m.C*Bfull))

    # Bus-angle-rate diagnostic is measurement-only and has no feedback path.
    Cbus=zeros(Float64,10,size(m.Ared,1))
    Dbus=zeros(Float64,10)
    for bus in 30:39
        i=bus-29;v=ctx.net.voltage[bus]
        angrow=zeros(Float64,1,size(m.Gy,1))
        angrow[1,2bus-1]=-imag(v)/(abs2(v)*2pi)
        angrow[1,2bus]= real(v)/(abs2(v)*2pi)
        Cbus[i,:].=(angrow*Cdot)[:]
        Dbus[i]=(angrow*Ddot)[1]
    end
    Cbusq=Cbus*q.Q
    Aq=q.A
    condA=cond(Aq)
    (;A=Aq,B=Bq,C=Cq,output_metadata=outputs.metadata,
      C_bus=Cbusq,D_bus=Dbus,B_full=Bfull,
      A_full=m.Ared,B_algebraic=m.B,C_algebraic=m.C,Gy=m.Gy,
      condition_A=condA,gauge_residual=q.gauge_residual,
      load_bus=Int(load_bus),load_per_unit_MW=1/Float64(system_base_mva),
      system_base_mva=Float64(system_base_mva))
end

function step_response_eigensystem(model)
    F=eigen(model.A)
    V=F.vectors
    λ=F.values
    z=V\model.B
    (;lambda=λ,V,z,condition_V=cond(V))
end

function _step_factors(lambda,t)
    [abs(λ)<1e-13 ? t : expm1(λ*t)/λ for λ in lambda]
end

function _at(eig,C,t; derivative_order=0)
    factors = derivative_order==0 ? _step_factors(eig.lambda,t) :
        _exp_factors(eig.lambda,t) .* eig.lambda.^(derivative_order-1)
    C*eig.V*(factors.*eig.z)
end
_exp_factors(lambda,t)=exp.(lambda .* t)

function _bisect_zero(fun,a,b; tol=1e-9,maxiter=80)
    fa,fb=fun(a),fun(b)
    fa==0 && return a
    fb==0 && return b
    signbit(fa)==signbit(fb) && return NaN
    for _ in 1:maxiter
        c=(a+b)/2;fc=fun(c)
        if abs(fc)<1e-13 || abs(b-a)<tol
            return c
        elseif signbit(fc)==signbit(fa)
            a,fa=c,fc
        else
            b,fb=c,fc
        end
    end
    (a+b)/2
end

function _roots_for_channel(fun,grid)
    roots=Float64[]
    for k in 1:length(grid)-1
        a,b=grid[k],grid[k+1]
        fa,fb=fun(a),fun(b)
        if isfinite(fa) && isfinite(fb) && signbit(fa)!=signbit(fb)
            r=_bisect_zero(fun,a,b)
            isfinite(r) && (isempty(roots) || abs(r-last(roots))>1e-7) && push!(roots,r)
        end
    end
    roots
end

"""Continuous-time unit-step extrema using modal evaluations and root refinement.
The time grid brackets extrema; sign-changing roots are refined by bisection."""
function step_metrics(model; disturbance_MW=100.0,horizon_s=60.0,dt_s=0.01)
    e=step_response_eigensystem(model)
    channels=size(model.C,1)
    # Cache output modal residues once. Re-forming C*V at every time sample
    # turns a dense validation trace into an avoidable cubic-cost loop.
    CV=model.C*e.V
    tgrid=collect(0.0:dt_s:horizon_s)
    F=zeros(Float64,channels,length(tgrid))
    R=zeros(Float64,channels,length(tgrid))
    for (j,t) in enumerate(tgrid)
        F[:,j].=real.(CV*(_step_factors(e.lambda,t).*e.z))
        R[:,j].=real.(CV*(_exp_factors(e.lambda,t).*e.z))
    end
    output_rows=NamedTuple[]
    for ch in 1:channels
        cv=CV[ch,:]
        funit=t->real(sum(cv .* _step_factors(e.lambda,t) .* e.z))
        runit=t->real(sum(cv .* _exp_factors(e.lambda,t) .* e.z))
        rdot=t->real(sum(cv .* _exp_factors(e.lambda,t) .* e.lambda .* e.z))
        tf=_roots_for_channel(runit,tgrid)
        tr=_roots_for_channel(rdot,tgrid)
        fss=real((-model.C[ch:ch,:]*(model.A\model.B))[1])
        ftime=vcat([0.0,horizon_s],tf)
        fvals=[abs(funit(t)) for t in ftime]
        rtime=vcat([0.0,horizon_s],tr)
        rvals=[abs(runit(t)) for t in rtime]
        bestf=argmax(fvals);bestr=argmax(rvals)
        if abs(fss)>fvals[bestf]
            fpeak=abs(fss); fpeaktime=Inf
        else
            fpeak=fvals[bestf]; fpeaktime=ftime[bestf]
        end
        meta=model.output_metadata[ch,:]
        push!(output_rows,(;bus=Int(meta.bus),kind=String(meta.kind),
            F_inf_unit_Hz=fss,F_peak_unit_Hz=max(fpeak,maximum(abs.(F[ch,:]))),
            F_peak_time_s=fpeaktime,
            R_peak_unit_Hz_s=max(rvals[bestr],maximum(abs.(R[ch,:]))),
            R_peak_time_s=rtime[bestr]))
    end
    table=DataFrame(output_rows)
    (;table,lambda=e.lambda,eigenvector_condition=e.condition_V,
      frequency_peak_unit=maximum(table.F_peak_unit_Hz),
      rocof_peak_unit=maximum(table.R_peak_unit_Hz_s),
      steady_frequency_unit=maximum(abs.(table.F_inf_unit_Hz)),
      frequency_peak_MW=maximum(table.F_peak_unit_Hz)*disturbance_MW,
      rocof_peak_MW=maximum(table.R_peak_unit_Hz_s)*disturbance_MW,
      steady_frequency_MW=maximum(abs.(table.F_inf_unit_Hz))*disturbance_MW,
      grid_dt_s=dt_s,horizon_s=horizon_s)
end

"""TGOV1 zero-frequency governor droop authority, MW/Hz, per retained fraction."""
function frequency_authority_from_governors(root,original)
    gov=CSV.read(joinpath(root,"reports","experiment_D","inputs","gov.csv"),DataFrame)
    rows=NamedTuple[]
    for r in eachrow(original)
        bus=Int(r.bus)
        j=findfirst(==(bus),Int.(gov.bus))
        if j===nothing
            c=0.0; R=NaN; DT=NaN; status="NO_TGOV1"
        else
            gr=gov[j,:];R=Float64(gr.R);DT=Float64(gr.DT)
            c=Float64(r.Sn_original_MVA)*(1/R+DT)/60.0
            status="TGOV1_UNSATURATED"
        end
        push!(rows,(;bus,Pgen0_MW=Float64(r.P_gen_MW),
            Sn0_MVA=Float64(r.Sn_original_MVA),R,DT,
            c_MW_per_Hz_per_epsilon=c,
            authority_per_dispatch_MW_per_Hz=c/Float64(r.P_gen_MW),
            status))
    end
    DataFrame(rows)
end

end
