module AnalyticContinuation

using LinearAlgebra

export control_authority, replacement_direction, box_limited_direction

function control_authority(Jrho::AbstractMatrix,JK::AbstractMatrix,Rinv::AbstractMatrix;rtol=1e-10)
    A=JK*Rinv*transpose(JK)
    Ap=pinv(A;rtol)
    P=A*Ap
    Aperp=(I(size(A,1))-P)*Jrho
    Qpll=transpose(Jrho)*Ap*Jrho
    return (A=A,Ap=Ap,projector=P,uncompensable=Aperp,Qpll=Qpll,
        singular_values=svdvals(A),rank=rank(A;rtol))
end

"Exact projected replacement tangent and minimum-normalized-gain compensation."
function replacement_direction(p,Jrho,JK,Rinv,Qgraph;
        epsilon=1e-6,eta=1.0,gamma=1.0,free_rho=trues(length(p)))
    auth=control_authority(Jrho,JK,Rinv)
    Qeff=epsilon*I(length(p))+eta*auth.Qpll+gamma*Qgraph
    indices=findall(free_rho)
    isempty(indices) && return (drho=zeros(length(p)),dK=zeros(size(JK,2)),authority=auth,Qeff=Qeff)
    Q=Matrix(Qeff[indices,indices]);v=Q\p[indices]
    Aperp=auth.uncompensable[:,indices]
    if norm(Aperp)>1e-10
        qap=Q\transpose(Aperp)
        v-=qap*(pinv(Aperp*qap)*(Aperp*v))
    end
    drho=zeros(length(p));drho[indices].=v
    scale=maximum(abs.(drho))
    scale>0 && (drho./=scale)
    dK=-Rinv*transpose(JK)*auth.Ap*(Jrho*drho)
    return (drho=drho,dK=dK,authority=auth,Qeff=Qeff)
end

"Set outward tangent components to zero at frozen box faces and recompute."
function box_limited_direction(p,Jrho,JK,Rinv,Qgraph,rho,Kp,Ki,domain;tol=1e-8)
    n=length(rho);free=trues(n);Rwork=Matrix(Rinv)
    for i in 1:n
        if Kp[i]<=domain["Kp_nom"]*domain["Kp_factor_min"]+tol ||
           Kp[i]>=domain["Kp_nom"]*domain["Kp_factor_max"]-tol
            Rwork[i,i]=0
        end
        if Ki[i]<=domain["Ki_nom"]*domain["Ki_factor_min"]+tol ||
           Ki[i]>=domain["Ki_nom"]*domain["Ki_factor_max"]-tol
            Rwork[n+i,n+i]=0
        end
    end
    result=nothing
    for _ in 1:3n+1
        result=replacement_direction(p,Jrho,JK,Rwork,Qgraph;free_rho=free,
            epsilon=domain["epsilon_metric"],eta=domain["eta_PLL_metric"],
            gamma=domain["gamma_graph_metric"])
        outward=[i for i in 1:n if free[i] && ((rho[i]>=1-tol && result.drho[i]>1e-12) ||
            (rho[i]<=tol && result.drho[i]<-1e-12))]
        gain_out=Int[]
        for i in 1:n
            if Rwork[i,i]>0 && ((Kp[i]<=domain["Kp_nom"]*domain["Kp_factor_min"]+tol && result.dK[i]<-1e-12) ||
                (Kp[i]>=domain["Kp_nom"]*domain["Kp_factor_max"]-tol && result.dK[i]>1e-12))
                push!(gain_out,i)
            end
            j=n+i
            if Rwork[j,j]>0 && ((Ki[i]<=domain["Ki_nom"]*domain["Ki_factor_min"]+tol && result.dK[j]<-1e-12) ||
                (Ki[i]>=domain["Ki_nom"]*domain["Ki_factor_max"]-tol && result.dK[j]>1e-12))
                push!(gain_out,j)
            end
        end
        isempty(outward) && isempty(gain_out) && return merge(result,(free_gains=[Rwork[i,i]>0 for i in 1:2n],))
        free[outward].=false
        for i in gain_out;Rwork[i,i]=0;end
    end
    error("box tangent active set did not converge")
end

end
