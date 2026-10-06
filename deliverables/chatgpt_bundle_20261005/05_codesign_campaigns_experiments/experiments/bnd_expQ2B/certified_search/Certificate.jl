function local_certificate(a,x;label="LOCAL_CERTIFICATE")
    v=L.multi_evaluate(a,x;derivatives=true);au=audit(a,x,v;tol=2e-7)
    data=Dict{String,Any}("retained_MW"=>v.J,"primal"=>au.primal,"stationarity"=>au.stationarity,
        "complementarity"=>au.complementarity,"active_constraints"=>au.names,"multipliers"=>au.mu,
        "LICQ_rank"=>au.rank,"active_count"=>au.active_count,"active_Jacobian_condition"=>cond(au.As),
        "alpha"=>v.alpha,"beta_observed"=>v.beta,"Fpeak"=>v.Fpeak,"Rpeak"=>v.Rpeak,
        "Fsteady"=>maximum(abs.(v.dc)),"formal_certificate"=>false,"global_certified"=>false,
        "SOSC"=>"NOT_TESTED_FIRST_ORDER_FAILED","local_KKT_certificate"=>false)
    if au.primal<2e-8 && au.stationarity<1e-6 && au.complementarity<1e-8 && au.rank==au.active_count
        strong=findall(au.mu.>1e-10);Z=nullspace(au.As[strong,:];rtol=1e-10)
        mudense=zeros(length(v.g))
        for (id,mu) in zip(au.ids,au.mu);id>0 && (mudense[id]=mu);end
        Hs=Matrix{Float64}[];valid=true
        for h in (2e-5,1e-5,5e-6)
            H=zeros(size(Z,2),size(Z,2))
            for j in axes(Z,2)
                dir=Z[:,j];dir[abs.(dir).<1e-12].=0.
                xp=x+h*dir;xm=x-h*dir
                if minimum(xp)<-1e-12||maximum(xp)>1+1e-12||minimum(xm)<-1e-12||maximum(xm)>1+1e-12
                    valid=false;break
                end
                vp=L.multi_evaluate(a,clamp.(xp,0.,1.);derivatives=true)
                vm=L.multi_evaluate(a,clamp.(xm,0.,1.);derivatives=true)
                H[:,j]=Z'*((L.matched_jacobian(vp,v)-L.matched_jacobian(vm,v))'*mudense)/(2h)
            end
            push!(Hs,(H+H')/2)
        end
        if valid
            delta=max(opnorm(Hs[1]-Hs[2]),opnorm(Hs[2]-Hs[3]));mineig=minimum(eigvals(Symmetric(Hs[3]));init=Inf)
            data["projected_Hessian_eigenvalues"]=eigvals(Symmetric(Hs[3]))
            data["Hessian_step_variation"]=delta;data["critical_subspace_dimension"]=size(Z,2)
            data["weak_multiplier_count"]=count(au.mu.<=1e-10)
            data["SOSC"]=(mineig>max(1e-8,10delta)) ? "NUMERICALLY_POSITIVE_ON_SUFFICIENT_CRITICAL_SUBSPACE" : "INCONCLUSIVE"
            data["local_KKT_certificate"]=(mineig>max(1e-8,10delta))
        else
            data["SOSC"]="WEAK_BOUND_CRITICAL_CONE_REQUIRES_ONE_SIDED_AUDIT"
        end
    end
    open(joinpath(L.OUT,"$label.toml"),"w") do io;TOML.print(io,data);end
    CSV.write(joinpath(L.OUT,"$(label)_MULTIPLIERS.csv"),DataFrame(constraint=au.names,multiplier=au.mu,value=au.values))
    println(data);flush(stdout)
    data
end
