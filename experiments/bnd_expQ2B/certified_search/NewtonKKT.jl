function newton_kkt(a,x0;maxiter=20,label="NEWTON")
    x=copy(x0);v=L.multi_evaluate(a,x;derivatives=true);history=NamedTuple[];status="MAX_ITERATIONS"
    for it in 1:maxiter
        au=audit(a,x,v;tol=1e-6)
        active=[k for (k,m) in zip(au.ids,au.mu) if k>0 && m>1e-10]
        fixed=findall((x.<1e-7).|(x.>1-1e-7))
        free=setdiff(1:length(x),fixed)
        J=v.Jac[active,free];sc=max.(sqrt.(sum(abs2,J;dims=2))[:],1e-12)
        Js=J./sc;mu=pinv(J';rtol=1e-11)*(-a.c[free])
        # Keep multiplier signs as an explicit gate; SQP must release a bound
        # or constraint if the equality KKT branch requires a negative multiplier.
        grad=a.c[free]+J'*mu
        println("NEWTON $it J=$(v.J) primal=$(au.primal) stationarity=$(au.stationarity) free=$(length(free)) active=$active mu=$mu");flush(stdout)
        push!(history,(;iteration=it,J_MW=v.J,primal=au.primal,stationarity=au.stationarity,
            complementarity=au.complementarity,active=join(v.names[active],";"),free=join(free,";"),minimum_multiplier=minimum(mu;init=Inf)))
        CSV.write(joinpath(L.OUT,"$(label)_HISTORY.csv"),DataFrame(history))
        save_point("$(label)_CURRENT.toml",a,x,v,au)
        au.primal<2e-8 && au.stationarity<1e-6 && (status="FIRST_ORDER_CANDIDATE";break)
        H=zeros(length(free),length(free));mudense=zeros(length(v.g));mudense[active]=mu
        for (jj,j) in enumerate(free)
            h=min(j<=length(a.support) ? 1e-5 : .001,.025/max(abs(v.Jac[2,j]),1e-9),x[j]/3,(1-x[j])/3)
            h<1e-9 && error("Free coordinate $j needs an explicit bound event")
            xp=copy(x);xm=copy(x);xp[j]+=h;xm[j]-=h
            vp=L.multi_evaluate(a,xp;derivatives=true);vm=L.multi_evaluate(a,xm;derivatives=true)
            H[:,jj]=((L.matched_jacobian(vp,v)-L.matched_jacobian(vm,v))'*mudense)[free]/(2h)
        end
        H=(H+H')/2;nf=length(free);na=length(active)
        K=[H Js';Js zeros(na,na)];rhs=-vcat(grad,v.g[active]./sc)
        step=pinv(K;rtol=1e-12)*rhs;dx=zeros(length(x));dx[free]=step[1:nf]
        dmu=step[nf+1:end]./sc
        fraction=min(1.,.12/max(norm(dx,Inf),1e-15))
        for j in free
            dx[j]>0 && (fraction=min(fraction,(1-x[j])/dx[j]))
            dx[j]<0 && (fraction=min(fraction,-x[j]/dx[j]))
        end
        merit0=norm(vcat(grad,v.g[active]./sc))
        accepted=false
        for power in 0:12
            t=fraction*2. ^(-power);xt=clamp.(x+t*dx,0.,1.)
            vt=try L.multi_evaluate(a,xt;derivatives=true) catch;continue end
            mt=mu+t*dmu;Jt=L.matched_jacobian(vt,v)[active,free]
            gt=matched_values(vt,v)[active];merit=norm(vcat(a.c[free]+Jt'*mt,gt./sc))
            if merit<(1-1e-3*t)*merit0 && maximum(vt.g)<max(.01,2au.primal)
                x=xt;v=vt;accepted=true;break
            end
        end
        if !accepted;status="KKT_LINE_SEARCH_BLOCKED";break;end
    end
    au=audit(a,x,v);save_point("$(label)_FINAL.toml",a,x,v,au;status)
    println("NEWTON_FINAL $status J=$(v.J) KKT=$(au.stationarity) primal=$(au.primal)");flush(stdout)
    (;a,x,v,audit=au,status)
end
