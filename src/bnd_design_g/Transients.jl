using LinearAlgebra
using Statistics

function bisection_root(f,a,b;tol=1e-12,maxiter=80)
    fa=f(a); fb=f(b)
    fa==0 && return a
    fb==0 && return b
    fa*fb>0 && return NaN
    for _ in 1:maxiter
        c=(a+b)/2; fc=f(c)
        abs(fc)<tol || abs(b-a)<tol && return c
        if fa*fc<=0; b=c; fb=fc else; a=c; fa=fc end
    end
    (a+b)/2
end

"Analytical modal and direct-matrix step responses for a retained-SG COI frequency."
function transient_analysis(net,model,eps,buses;event_bus=16,system_base_mva=100.0,
                            f0=60.0,rocof_limit=0.5,frequency_limit=0.5)
    h=zeros(length(buses)); states=Int[]; weights=Float64[]
    for (j,b) in enumerate(buses)
        p=net.sg[b].op.parameters
        h[j]=p.inertia*p.rating_mva
        row=findfirst(r->Int(r.bus)==b && r.kind=="SG",eachrow(model.state_map))
        if row!==nothing && eps[j]>0
            rr=model.state_map[row,:]
            idx=Int(rr.first)+(Int(rr.last)-Int(rr.first))-1
            push!(states,idx); push!(weights,eps[j]*h[j])
        end
    end
    totalH=sum(weights)
    totalH>0 || return (;status="BLOCKED_NO_RETAINED_SYNCHRONOUS_INERTIA",
        H_total_MVA_s=0.0,rocof_unit_gain=NaN,frequency_unit_gain=NaN)
    Cf=zeros(Float64,model.n_dynamic)
    for (i,w) in zip(states,weights); Cf[i]=f0*w/totalH; end
    V=net.voltage[event_bus]
    Iinj=conj((1/system_base_mva)/V) # one MW active injection; sign does not affect capacities
    iv=zeros(Float64,size(model.Gy,1)); iv[2event_bus-1]=real(Iinj); iv[2event_bus]=imag(Iinj)
    Bd=model.B*(model.Gy\iv)
    A=model.Ared; λ,R=eigen(A); L=inv(R)'
    q=ComplexF64[]
    for j in eachindex(λ)
        den=dot(L[:,j],R[:,j])
        push!(q,(dot(Cf,R[:,j])*dot(L[:,j],Bd))/den)
    end
    stable=maximum(real.(λ)) < -1e-9 || (abs(maximum(real.(λ)))<1e-8)
    decay=minimum(-real(z) for z in λ if real(z)<-1e-8)
    tmax=min(600.0,max(30.0,12/max(decay,1e-5)))
    maxomega=maximum(abs.(imag.(λ)),init=0.0)
    dt=min(0.02,2pi/(16max(maxomega,1.0)))
    ntime=min(250_001,max(5_001,ceil(Int,tmax/dt)+1))
    ts=collect(range(0.0,tmax,length=ntime))
    rocof_modal(t)=real(sum(q[j]*exp(λ[j]*t) for j in eachindex(λ)))
    freq_modal(t)=real(sum(abs(λ[j])<1e-10 ? q[j]*t :
        q[j]*expm1(λ[j]*t)/λ[j] for j in eachindex(λ)))
    Rm=rocof_modal.(ts); Fm=freq_modal.(ts)
    # Direct exp(A t) validation on a deterministic coarser mesh.
    check_ts=collect(range(0.0,tmax,length=41))
    directR=Float64[]; directF=Float64[]
    Aug=zeros(Float64,size(A,1)+1,size(A,2)+1)
    Aug[1:end-1,1:end-1].=A; Aug[1:end-1,end].=Bd
    for t in check_ts
        X=exp(A*t)*Bd
        push!(directR,dot(Cf,X))
        Z=exp(Aug*t)
        push!(directF,dot(Cf,Z[1:end-1,end]))
    end
    modalR=rocof_modal.(check_ts); modalF=freq_modal.(check_ts)
    rerr=maximum(abs.(directR-modalR))/max(maximum(abs.(directR)),1e-12)
    ferr=maximum(abs.(directF-modalF))/max(maximum(abs.(directF)),1e-12)
    # Deterministic root brackets for exact modal RoCoF extrema.
    roots=Float64[]
    for i in 1:length(ts)-1
        if Rm[i]*Rm[i+1]<0
            z=bisection_root(rocof_modal,ts[i],ts[i+1])
            isfinite(z) && push!(roots,z)
        end
    end
    candidates=unique(vcat(0.0,roots,tmax))
    fvals=freq_modal.(candidates)
    nadir_idx=argmin(fvals); peak_idx=argmax(fvals)
    qR=sum(abs.(q))
    qF=sum(abs(q[j])*2/max(abs(λ[j]),1e-12) for j in eachindex(λ) if abs(λ[j])>1e-10)
    peak_rocof=max(maximum(abs.(Rm)),abs(rocof_modal(tmax)))
    peak_frequency=max(abs(minimum(Fm)),abs(maximum(Fm)))
    qdom=argmax(abs.(q)); domλ=λ[qdom]
    pair=findall(j->abs(λ[j]-conj(domλ))<1e-7,eachindex(λ))
    pair=unique(vcat(qdom,pair)); tail=sum(abs(q[j]) for j in eachindex(q) if !(j in pair))
    initial=abs(rocof_modal(0.0))
    cap_init=2rocof_limit*totalH/f0
    cap_peak=rocof_limit/max(qR,1e-15)
    cap_freq=frequency_limit/max(qF,1e-15)
    residue=DataFrame(mode=collect(eachindex(λ)),lambda_real=real.(λ),lambda_imag=imag.(λ),
        frequency_hz=abs.(imag.(λ))/(2pi),
        damping_ratio=[abs(imag(z))<1e-12 ? NaN : -real(z)/abs(z) for z in λ],
        q_real=real.(q),q_imag=imag.(q),q_abs=abs.(q),
        rocof_participation=abs.(q)./max(qR,Base.eps(Float64)),
        frequency_excursion_bound=[abs(q[j])*2/max(abs(λ[j]),1e-12) for j in eachindex(λ)])
    (;status=stable ? "EVALUATED" : "UNSTABLE_ANALYTICAL_MODEL",
      A=A,Bd=Bd,Cf=Cf,eigenvalues=λ,residues=q,times=ts,rocof=Rm,
      frequency=Fm,unit_rocof_peak=peak_rocof,unit_initial_rocof=initial,
      unit_frequency_peak=peak_frequency,frequency_nadir=minimum(Fm),
      frequency_nadir_time=candidates[nadir_idx],frequency_peak_time=candidates[peak_idx],
      rocof_unit_gain=qR,frequency_unit_gain=qF,H_total_MVA_s=totalH,
      deltaP_max_initial_rocof_MW=cap_init,deltaP_max_peak_rocof_MW=cap_peak,
      deltaP_max_frequency_MW=cap_freq,
      deltaP_max_total_MW=min(cap_init,cap_peak,cap_freq),
      direct_modal_rocof_relative_error=rerr,direct_modal_frequency_relative_error=ferr,
      dominant_mode=qdom,dominant_pair=pair,dominant_pair_tail_bound=tail,
      dominant_pair_residue=sum(abs(q[j]) for j in pair),dominant_pair_justified=
          tail<=0.1*max(sum(abs(q[j]) for j in pair),Base.eps(Float64)),
      modal_residues=residue,frequency_limit_Hz=frequency_limit,
      rocof_limit_Hz_s=rocof_limit,checked_times=check_ts,
      checked_direct_rocof=directR,checked_modal_rocof=modalR,
      checked_direct_frequency=directF,checked_modal_frequency=modalF)
end
