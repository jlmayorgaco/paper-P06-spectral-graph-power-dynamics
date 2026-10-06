include("ReducedDAE.jl")
using .ReducedDAE, LinearAlgebra, TOML, SHA, Random, ForwardDiff, CSV, DataFrames
const R=ReducedDAE
const OUT=joinpath(R.OUT,"coupled_dq_contract")
mkpath(OUT); BLAS.set_num_threads(1)
ctx=R.N.design_context(R.ROOT)
candidate=joinpath(R.OUT,"candidate_final_physical.toml")
dd=TOML.parsefile(candidate)
m=R.model(ctx,Float64.(dd["rho"]),Float64.(dd["Kp"]),Float64.(dd["Ki"]);dc_convention=:physical_supply)
@assert all(0 .< m.rho .<1)
@assert all(p.rs==0 && p.xdpp==p.xqpp for p in m.sp)
xbase=copy(m.x0)
for ix in m.gfidx
    sn,cs=sincos(m.x0[ix[3]]);ii=m.x0[ix[6]];ir=m.x0[ix[7]]
    xbase[ix[6]]=-sn*ir+cs*ii;xbase[ix[7]]=cs*ir+sn*ii
end
ig=last(m.sgidx[end]); keep=setdiff(1:length(m.x0),[ig]); nx=length(keep)+78
tau_f=.1; tau_r=.1
vnom=R.voltage(m.x0,m;allbus=true); phase0=atan.(vnom[2:2:end],vnom[1:2:end])
function full_rhs(z,w)
    T=promote_type(eltype(z),eltype(w)); x=T.(xbase);x[keep].+=view(z,1:length(keep))
    G=Matrix{T}(m.net.Y);hh=T.(view(w,1:78))
    for i in 1:10
        vi=(2(i+29)-1):(2(i+29));D,c=R.norton(view(x,m.sgidx[i]),m.sp[i])
        G[vi,vi].+=(1-m.rho[i])*D;hh[vi].+=(1-m.rho[i])*c
        ix=m.gfidx[i];sn,cs=sincos(x[ix[3]])
        hh[vi].+=m.rho[i].*[cs*x[ix[7]]-sn*x[ix[6]],sn*x[ix[7]]+cs*x[ix[6]]]
    end
    v=-(G\hh); dx=zeros(T,length(x));nu=m.sp[end].omega_base*(x[ig-1]-m.sp[end].omega_frame)
    for i in 1:10
        vi=(2(i+29)-1):(2(i+29));u=view(v,vi);xi=m.sgidx[i];p=m.sp[i]
        dx[xi].=p.controlled ? R.SG.rhs(view(x,xi),u,p) : R.SG.rhs_uncontrolled(view(x,xi),u,p)
        dx[last(xi)]-=nu
        ix=m.gfidx[i];pp=m.gp[i]
        gq,gd,theta,om,xi,iq,id,vdi,vdc=x[ix];sn,cs=sincos(theta)
        ud=cs*u[1]+sn*u[2];uq=-sn*u[1]+cs*u[2]
        ed=(vdc-pp.Vdc)*pp.dc_kp+vdi-id;eq=pp.iset_q-iq
        vid=pp.cc_kp*ed+pp.cc_ki*gd;viq=pp.cc_kp*eq+pp.cc_ki*gq
        dx[ix].=[eq,ed,om-nu,(xi+m.kp[i]*uq-om)/pp.pll_tau,m.ki[i]*uq,
                 pp.omega_base/pp.Xf*(viq-uq-pp.Rf*iq)-(pp.omega_base*pp.omega_frame+om)*id,
                 pp.omega_base/pp.Xf*(vid-ud-pp.Rf*id)+(pp.omega_base*pp.omega_frame+om)*iq,
                 (vdc-pp.Vdc)*pp.dc_ki,(pp.Pdc+w[78+i]-vid*id-viq*iq)/(pp.Cdc*vdc)]
    end
    out=vcat(dx[keep],zeros(T,78));ff=zeros(T,39);rr=zeros(T,39)
    for b in 1:39
        ur,ui=v[2b-1:2b];cs=cos(phase0[b]);sn=sin(phase0[b])
        ph=atan((cs*ui-sn*ur)/(cs*ur+sn*ui))
        eta=z[length(keep)+b];chi=z[length(keep)+39+b]
        ff[b]=(ph-eta)/(2pi*tau_f);rr[b]=(ff[b]-chi)/tau_r
        out[length(keep)+b]=(ph-eta)/tau_f-nu
        out[length(keep)+39+b]=rr[b]
    end
    (;out,v,ff,rr)
end
function global_rhs(z,w)
    T=promote_type(eltype(z),eltype(w)); x=T.(m.x0);x[keep].+=view(z,1:length(keep))
    G=Matrix{T}(m.net.Y);hh=T.(view(w,1:78))
    for i in 1:10
        vi=(2(i+29)-1):(2(i+29));D,c=R.norton(view(x,m.sgidx[i]),m.sp[i])
        G[vi,vi].+=(1-m.rho[i])*D;hh[vi].+=(1-m.rho[i])*c
        ix=m.gfidx[i];hh[vi].+=m.rho[i].*x[ix[[7,6]]]
    end
    v=-(G\hh); dx=zeros(T,length(x));nu=m.sp[end].omega_base*(x[ig-1]-m.sp[end].omega_frame)
    for i in 1:10
        vi=(2(i+29)-1):(2(i+29));u=view(v,vi);xi=m.sgidx[i];p=m.sp[i]
        dx[xi].=p.controlled ? R.SG.rhs(view(x,xi),u,p) : R.SG.rhs_uncontrolled(view(x,xi),u,p)
        dx[last(xi)]-=nu
        ix=m.gfidx[i];pp=m.gp[i]
        dx[ix].=R.gfl_rhs(view(x,ix),u,pp,m.kp[i],m.ki[i],:physical_supply)
        dx[ix[3]]-=nu;dx[ix[6]]-=nu*x[ix[7]];dx[ix[7]]+=nu*x[ix[6]]
        dx[ix[9]]+=w[78+i]/(pp.Cdc*x[ix[9]])
    end
    out=vcat(dx[keep],zeros(T,78));ff=zeros(T,39);rr=zeros(T,39)
    for b in 1:39
        ur,ui=v[2b-1:2b];cs=cos(phase0[b]);sn=sin(phase0[b])
        ph=atan((cs*ui-sn*ur)/(cs*ur+sn*ui))
        eta=z[length(keep)+b];chi=z[length(keep)+39+b]
        ff[b]=(ph-eta)/(2pi*tau_f);rr[b]=(ff[b]-chi)/tau_r
        out[length(keep)+b]=(ph-eta)/tau_f-nu
        out[length(keep)+39+b]=rr[b]
    end
    (;out,v,ff,rr)
end
z=zeros(nx);w=zeros(88);zeroout=full_rhs(z,w)
println("EXPORT_BEGIN nx=",nx," trim=",norm(zeroout.out,Inf));flush(stdout)
A=ForwardDiff.jacobian(z->full_rhs(z,w).out,z)
B=ForwardDiff.jacobian(w->full_rhs(z,w).out,w)
println("JACOBIANS alpha=",maximum(real,eigvals(A)));flush(stdout)
SGPAR=[Dict(String(k)=>getfield(p,k) for k in fieldnames(typeof(p))) for p in m.sp]
GFPAR=[Dict(String(k)=>v for (k,v) in pairs(p)) for p in m.gp]
files=[candidate,joinpath(@__DIR__,"ReducedDAE.jl"),joinpath(R.ROOT,"src","bnd_design","AnalyticSG.jl"),
       joinpath(R.ROOT,"src","bnd_model_expN","PDExactDesignN.jl")]
data=Dict("status"=>"EXACT_NONLINEAR_FULL_NETWORK_MODEL_EXPORTED",
    "rho"=>m.rho,"Kp"=>m.kp,"Ki"=>m.ki,"power_MW"=>ctx.power,
    "Y"=>collect(eachrow(m.net.Y)),"x0"=>xbase,"sg_indices"=>collect.(m.sgidx),
    "gfl_indices"=>collect.(m.gfidx),"keep_indices"=>keep,"gauge_index"=>ig,
    "sg_parameters"=>SGPAR,"gfl_parameters"=>GFPAR,"phase0"=>phase0,
    "tau_f"=>tau_f,"tau_r"=>tau_r,"f0"=>zeroout.out,"v0"=>zeroout.v,
    "A"=>collect(eachrow(A)),"B"=>collect(eachrow(B)),
    "state_count"=>nx,"input_count"=>length(w),"modal_alpha"=>maximum(real,eigvals(A)),
    "source_hashes"=>Dict(relpath(f,R.ROOT)=>bytes2hex(sha256(read(f))) for f in files),
    "script_sha256"=>bytes2hex(sha256(read(@__FILE__))))
open(joinpath(OUT,"model.toml"),"w") do io; TOML.print(io,data);end
rng=MersenneTwister(20261002);rows=NamedTuple[];max_mapping_error=0.
for k in 1:30
    zz=1e-5.*randn(rng,nx);ww=1e-5.*randn(rng,88);r=full_rhs(zz,ww)
    xq=copy(xbase);xq[keep].+=zz[1:length(keep)];xg=copy(xq)
    for ix in m.gfidx
        sn,cs=sincos(xq[ix[3]]);iq=xq[ix[6]];id=xq[ix[7]]
        xg[ix[6]]=sn*id+cs*iq;xg[ix[7]]=cs*id-sn*iq
    end
    zr=vcat(xg[keep]-m.x0[keep],zz[length(keep)+1:end]);gr=global_rhs(zr,ww)
    fx=zeros(length(xbase));fx[keep].=gr.out[1:length(keep)]
    for ix in m.gfidx
        sn,cs=sincos(xq[ix[3]]);di=fx[ix[6]];dr=fx[ix[7]];rate=fx[ix[3]]
        fx[ix[6]]=-sn*dr+cs*di-xq[ix[7]]*rate
        fx[ix[7]]=cs*dr+sn*di+xq[ix[6]]*rate
    end
    mapped=vcat(fx[keep],gr.out[length(keep)+1:end])
    err=norm(mapped-r.out,Inf);@assert err<1e-6
    global max_mapping_error=max(max_mapping_error,err)
    push!(rows,(;case=k,z=zz,w=ww,f=r.out,v=r.v,frequency=r.ff,rocof=r.rr))
end
open(joinpath(OUT,"julia_probes.toml"),"w") do io
    TOML.print(io,Dict("probes"=>[Dict(String(k)=>v for (k,v) in pairs(r)) for r in rows],
       "max_global_coordinate_mapping_error"=>max_mapping_error,
       "model_sha256"=>bytes2hex(sha256(read(joinpath(OUT,"model.toml"))))))
end
println("EXPORTED_DQ_NETWORK; mapped global-coordinate error ",max_mapping_error)
