module PDExactDesignN

using CSV, DataFrames, ForwardDiff, LinearAlgebra, SHA
include(joinpath(@__DIR__, "..", "bnd_design_e", "CollectiveModel.jl"))
using .CollectiveModel

export DesignContext, design_context, trim_gfl, gfl_jacobians,
    descriptor, spectrum, port_transfer, gain_bounds, simple_mode_sensitivities

const SBASE = 100.0
const K0P = 2pi*5
const K0I = (2pi*5)^2/4

struct DesignContext
    root::String
    net::Any
    original::DataFrame
    power::Vector{Float64}
    kpmin::Vector{Float64}
    kpmax::Vector{Float64}
    kimin::Vector{Float64}
    kimax::Vector{Float64}
end

function design_context(root::AbstractString)
    path=joinpath(root,"reports","experiment_N","TABLE_N01_original_operating_point.csv")
    hashpath=path*".sha256"
    bytes2hex(sha256(read(path)))==strip(read(hashpath,String)) || error("N01 freeze hash mismatch")
    original=CSV.read(path,DataFrame)
    sort!(original,:bus)
    Int.(original.bus)==collect(30:39) || error("N01 bus set differs from generator set")
    net=CollectiveModel.frozen_network(root)
    for row in eachrow(original)
        bus=Int(row.bus)
        op=net.sg[bus].op
        v=net.voltage[bus]
        abs(real(v)-Float64(row.V_real_pu))<1e-9 || error("SG voltage mismatch at bus $bus")
        abs(imag(v)-Float64(row.V_imag_pu))<1e-9 || error("SG voltage mismatch at bus $bus")
        diagnosis=net.sg[bus].diagnostic
        abs(100*diagnosis.network_power.p-Float64(row.P_gen_MW))<1e-7 ||
            error("SG direct P mismatch at bus $bus")
        abs(100*diagnosis.network_power.q-Float64(row.Q_gen_Mvar))<1e-7 ||
            error("SG direct Q mismatch at bus $bus")
    end
    n=10
    DesignContext(String(root),net,original,Float64.(original.P_gen_MW),
        fill(0.25K0P,n),fill(4K0P,n),fill(0.25K0I,n),fill(4K0I,n))
end

gain_bounds(ctx::DesignContext)=(;kpmin=ctx.kpmin,kpmax=ctx.kpmax,
    kimin=ctx.kimin,kimax=ctx.kimax)

"""Full-dispatch per-unit GFL trim derived from installed PLL/CC1/DC equations."""
function trim_gfl(ctx::DesignContext,bus::Integer)
    row=ctx.original[findfirst(==(bus),Int.(ctx.original.bus)),:]
    ur=Float64(row.V_real_pu);ui=Float64(row.V_imag_pu)
    p=Float64(row.P_gen_MW)/SBASE;q=Float64(row.Q_gen_Mvar)/SBASE
    vm=hypot(ur,ui);vm>0 || error("zero voltage at bus $bus")
    theta=atan(ui,ur)
    ir=(p*ur+q*ui)/vm^2; ii=(p*ui-q*ur)/vm^2
    id=p/vm;iq=-q/vm
    rf=0.01;xf=0.03
    vid=vm+rf*id-xf*iq;viq=rf*iq+xf*id
    cki=(xf/(2pi*60))*(600*2pi)^2/4
    x=Float64[viq/cki,vid/cki,theta,0,0,ii,ir,id,2.5]
    pars=(;Rf=rf,Xf=xf,omega_base=2pi*60,omega_frame=1.0,
        pll_tau=1/(300*2pi),cc_kp=(xf/(2pi*60))*(600*2pi),
        cc_ki=cki,Cdc=1.25,Vdc=2.5,
        dc_kp=2.5*1.25*(5*2pi),dc_ki=2.5*1.25*(5*2pi)^2/4,
        iset_q=iq,Pdc=vid*id+viq*iq)
    (;x,u=Float64[ur,ui],pars,current=Float64[ir,ii],
        target_p=p,target_q=q)
end

function gfl_rhs(x,u,p,Kp,Ki)
    γq,γd,θ,Δω,Δωi,ifi,ifr,vdi,vdc=x
    ur,ui=u;c=cos(θ);s=sin(θ)
    id=c*ifr+s*ifi;iq=-s*ifr+c*ifi
    irefd=(p.Vdc-vdc)*p.dc_kp+vdi
    ed=irefd-id;eq=p.iset_q-iq
    vid=p.cc_kp*ed+p.cc_ki*γd
    viq=p.cc_kp*eq+p.cc_ki*γq
    vir=c*vid-s*viq;vii=s*vid+c*viq
    eang=-s*ur+c*ui
    pac=vid*id+viq*iq
    [eq,ed,Δω,(Δωi+Kp*eang-Δω)/p.pll_tau,Ki*eang,
     (p.omega_base/p.Xf)*(vii-ui-p.Rf*ifi-p.omega_frame*p.Xf*ifr),
     (p.omega_base/p.Xf)*(vir-ur-p.Rf*ifr+p.omega_frame*p.Xf*ifi),
     (p.Vdc-vdc)*p.dc_ki,(pac-p.Pdc)/(p.Cdc*vdc)]
end

function gfl_jacobians(op,Kp,Ki)
    x=op.x;u=op.u;p=op.pars
    fz=z->gfl_rhs(z,u,p,Kp,Ki)
    fu=v->gfl_rhs(x,v,p,Kp,Ki)
    A=ForwardDiff.jacobian(fz,x)
    B=ForwardDiff.jacobian(fu,u)
    C=zeros(2,9);C[1,7]=1;C[2,6]=1
    D=zeros(2,2)
    (;A,B,C,D)
end

"""Fixed-architecture physical port closure; never calls PowerDynamics."""
function descriptor(ctx::DesignContext,rho::AbstractVector,
                    kp::AbstractVector=fill(K0P,10),ki::AbstractVector=fill(K0I,10))
    length(rho)==length(kp)==length(ki)==10 || throw(DimensionMismatch("ten triples required"))
    all(0 .<= rho .<= 1) || throw(ArgumentError("rho outside [0,1]"))
    all(ctx.kpmin .<= kp .<= ctx.kpmax) || throw(ArgumentError("Kp outside ExpK box"))
    all(ctx.kimin .<= ki .<= ctx.kimax) || throw(ArgumentError("Ki outside ExpK box"))
    blocks=NamedTuple[]
    for bus in 30:39
        r=Float64(rho[bus-29])
        if r<1
            push!(blocks,(;bus,kind="SG",share=1-r,J=ctx.net.sg[bus].J,
                equilibrium=ctx.net.sg[bus].op.x))
        end
        if r>0
            op=trim_gfl(ctx,bus)
            J=gfl_jacobians(op,Float64(kp[bus-29]),Float64(ki[bus-29]))
            push!(blocks,(;bus,kind="GFL",share=r,J,equilibrium=op.x))
        end
    end
    nx=sum(size(b.J.A,1) for b in blocks)
    ny=size(ctx.net.y_static,1)
    A=zeros(nx,nx);B=zeros(nx,ny);C=zeros(ny,nx);D=zeros(ny,ny)
    map=NamedTuple[];inventory=NamedTuple[];off=0
    for b in blocks
        d=size(b.J.A,1);xi=(off+1):(off+d);yi=(2b.bus-1):(2b.bus)
        A[xi,xi].=b.J.A;B[xi,yi].=b.J.B
        C[yi,xi].+=b.share.*b.J.C
        D[yi,yi].+=b.share.*b.J.D
        push!(map,(;bus=b.bus,kind=b.kind,first=first(xi),last=last(xi)))
        names=if b.kind=="SG"
            b.bus==39 ? CollectiveModel.AnalyticSG.uncontrolled_state_names() :
                        CollectiveModel.AnalyticSG.state_names()
        else
            ["gamma_q","gamma_d","theta","delta_omega_rad_s",
             "delta_omega_i_rad_s","i_f_i","i_f_r","v_dc_i","v_dc_state"]
        end
        length(names)==d || error("physical state name count mismatch at bus $(b.bus)")
        for (j,name) in enumerate(names)
            push!(inventory,(;index=off+j,bus=b.bus,kind=b.kind,state_name=name))
        end
        off+=d
    end
    Gy=ctx.net.y_static+D
    Ared=A-B*(Gy\C)
    (;A,B,C,D,Gy,Ared,state_map=DataFrame(map),state_inventory=DataFrame(inventory),blocks,
        n_dynamic=nx,n_algebraic=ny)
end

function gauge_vector(m)
    g=zeros(size(m.Ared,1))
    for block in m.blocks
        row=only(eachrow(m.state_map[(m.state_map.bus .== block.bus) .&
            (m.state_map.kind .== block.kind),:]))
        first=Int(row.first);last=Int(row.last)
        if block.kind=="SG"
            g[last]=1
        else
            x=block.equilibrium
            g[first+2]=1
            g[first+5]=x[7]
            g[first+6]=-x[6]
        end
    end
    g
end

function spectrum(ctx::DesignContext,rho,kp=fill(K0P,10),ki=fill(K0I,10))
    m=descriptor(ctx,rho,kp,ki)
    g=gauge_vector(m)
    gauge_res=norm(m.Ared*g)/max(norm(m.Ared)*norm(g),eps())
    gauge_res<1e-8 || error("gauge residual $gauge_res")
    Q=nullspace(reshape(g/norm(g),1,:))
    Aq=transpose(Q)*m.Ared*Q
    λ=eigvals(Aq)
    (;model=m,lambda=λ,alpha=maximum(real.(λ)),gauge_residual=gauge_res,
        gauge=g,quotient=Q)
end

function port_transfer(J,s)
    J.D+J.C*((s*I(size(J.A,1))-J.A)\J.B)
end

"""Total simple-pole derivatives within a fixed physical architecture.

The declared per-unit internal trim and fixed bus voltage do not change with
rho, Kp or Ki. This is the explicit solution of the equilibrium IFT in these
coordinates: dx*/dz=0. The SG rating and GFL aggregate port factor enter C,D;
PLL gains enter the local GFL A,B through the phase detector derivatives.
"""
function simple_mode_sensitivities(ctx::DesignContext,rho,kp,ki;mode=:rightmost)
    s=spectrum(ctx,rho,kp,ki)
    m=s.model;Q=s.quotient;Aq=transpose(Q)*m.Ared*Q
    F=eigen(Aq)
    j=mode===:rightmost ? argmax(real.(F.values)) : Int(mode)
    λ=F.values[j];right=F.vectors[:,j]
    FL=eigen(adjoint(Aq))
    left=FL.vectors[:,argmin(abs.(FL.values .- conj(λ)))]
    denominator=dot(left,right)
    abs(denominator)>1e-12 || error("ill-conditioned left/right pole pairing")
    gyro=m.Gy\m.C
    gr=fill(ComplexF64(NaN),10)
    gp=fill(ComplexF64(NaN),10)
    gi=fill(ComplexF64(NaN),10)
    modal(dAz)=dot(left,transpose(Q)*dAz*Q*right)/denominator
    for bus in 30:39
        i=bus-29;yi=(2bus-1):(2bus)
        sgidx=findfirst(b->b.bus==bus && b.kind=="SG",m.blocks)
        gfidx=findfirst(b->b.bus==bus && b.kind=="GFL",m.blocks)
        if sgidx!==nothing && gfidx!==nothing
            sb=m.blocks[sgidx];fb=m.blocks[gfidx]
            sr=m.state_map[sgidx,:];fr=m.state_map[gfidx,:]
            sx=Int(sr.first):Int(sr.last);fx=Int(fr.first):Int(fr.last)
            dC=zeros(size(m.C));dD=zeros(size(m.D))
            dC[yi,sx].=-sb.J.C
            dC[yi,fx].=fb.J.C
            dD[yi,yi].=fb.J.D-sb.J.D
            gr[i]=modal(m.B*(m.Gy\(dD*gyro-dC)))
        end
        if gfidx!==nothing
            fb=m.blocks[gfidx];fr=m.state_map[gfidx,:]
            fx=Int(fr.first):Int(fr.last)
            op=trim_gfl(ctx,bus)
            θ=op.x[3];ur=op.u[1];ui=op.u[2]
            sinθ=sin(θ);cosθ=cos(θ)
            de_dθ=-cosθ*ur-sinθ*ui
            for (dest,row,scale) in ((gp,4,1/op.pars.pll_tau),(gi,5,1.0))
                dAr=zeros(size(m.A));dBr=zeros(size(m.B))
                dAr[fx[row],fx[3]]=scale*de_dθ
                dBr[fx[row],yi[1]]=scale*(-sinθ)
                dBr[fx[row],yi[2]]=scale*cosθ
                dest[i]=modal(dAr-dBr*gyro)
            end
        end
    end
    condition=norm(left)*norm(right)/abs(denominator)
    (;lambda=λ,rho=gr,Kp=gp,Ki=gi,condition,denominator,
        gauge_residual=s.gauge_residual)
end

end
