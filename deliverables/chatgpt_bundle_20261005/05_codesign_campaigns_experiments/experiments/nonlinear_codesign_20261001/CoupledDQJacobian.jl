"Exact local-device Jacobians plus full algebraic Schur derivative; no model truncation."
function dq_local_rhs(x,u,p,kp,ki,wdc)
    gq,gd,theta,om,xi,iq,id,vdi,vdc=x;sn,cs=sincos(theta)
    ud=cs*u[1]+sn*u[2];uq=-sn*u[1]+cs*u[2]
    ed=(vdc-p.Vdc)*p.dc_kp+vdi-id;eq=p.iset_q-iq
    vid=p.cc_kp*ed+p.cc_ki*gd;viq=p.cc_kp*eq+p.cc_ki*gq
    [eq,ed,om,(xi+kp*uq-om)/p.pll_tau,ki*uq,
     p.omega_base/p.Xf*(viq-uq-p.Rf*iq)-(p.omega_base*p.omega_frame+om)*id,
     p.omega_base/p.Xf*(vid-ud-p.Rf*id)+(p.omega_base*p.omega_frame+om)*iq,
     (vdc-p.Vdc)*p.dc_ki,(p.Pdc+wdc-vid*id-viq*iq)/(p.Cdc*vdc)]
end

function dq_augmented_jacobian(z,w)
    m=Q.m;R=Q.R;nraw=length(Q.xbase);nk=length(Q.keep)
    x=copy(Q.xbase);x[Q.keep].+=z[1:nk]
    v=Q.full_rhs(z,w).v
    A=zeros(nraw,nraw);B=zeros(nraw,78);C=zeros(78,nraw);G=copy(m.net.Y)
    for i in 1:10
        vi=(2(i+29)-1):(2(i+29));u=v[vi];sx=m.sgidx[i];p=m.sp[i];ns=length(sx)
        js=ForwardDiff.jacobian(vcat(x[sx],u)) do xu
            p.controlled ? R.SG.rhs(view(xu,1:ns),view(xu,ns+1:ns+2),p) :
                           R.SG.rhs_uncontrolled(view(xu,1:ns),view(xu,ns+1:ns+2),p)
        end
        A[sx,sx].=js[:,1:ns];B[sx,vi].=js[:,ns+1:ns+2]
        C[vi,sx].+=(1-m.rho[i]).*ForwardDiff.jacobian(xs->R.SG.output(xs,u,p),x[sx])
        D,_=R.norton(x[sx],p);G[vi,vi].+=(1-m.rho[i]).*D
        A[last(sx),Q.ig-1]-=p.omega_base
        gx=m.gfidx[i];gp=m.gp[i]
        jg=ForwardDiff.jacobian(xu->dq_local_rhs(view(xu,1:9),view(xu,10:11),gp,m.kp[i],m.ki[i],w[78+i]),vcat(x[gx],u))
        A[gx,gx].=jg[:,1:9];B[gx,vi].=jg[:,10:11]
        A[gx[3],Q.ig-1]-=m.sp[end].omega_base
        sn,cs=sincos(x[gx[3]]);iq=x[gx[6]];id=x[gx[7]];rho=m.rho[i]
        C[vi,gx[3]].+=rho.*[-sn*id-cs*iq,cs*id-sn*iq]
        C[vi,gx[6]].+=rho.*[-sn,cs]
        C[vi,gx[7]].+=rho.*[cs,sn]
    end
    vx=-(G\C);red=A+B*vx
    J=zeros(Q.nx,Q.nx);J[1:nk,1:nk].=red[Q.keep,Q.keep]
    ref=findfirst(==(Q.ig-1),Q.keep)
    for bus in 1:39
        ur,ui=v[2bus-1:2bus];phase=(-ui.*vx[2bus-1,Q.keep]+ur.*vx[2bus,Q.keep])/(ur^2+ui^2)
        ie=nk+bus;ic=nk+39+bus
        J[ie,1:nk].=phase/Q.tau_f;J[ie,ref]-=m.sp[end].omega_base;J[ie,ie]=-1/Q.tau_f
        J[ic,1:nk].=phase/(2pi*Q.tau_f*Q.tau_r);J[ic,ie]=-1/(2pi*Q.tau_f*Q.tau_r);J[ic,ic]=-1/Q.tau_r
    end
    J
end
