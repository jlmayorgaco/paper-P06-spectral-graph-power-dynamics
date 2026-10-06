function feasible_roundoff_representative(a,x0)
    x=copy(x0);ref=L.multi_evaluate(a,x;derivatives=true);au=audit(a,x,ref)
    active=[k for (k,mu) in zip(au.ids,au.mu) if k>0 && mu>1e-10]
    targets=[k==2 ? -1e-7 : -1e-8 for k in active]
    for _ in 1:6
        v=L.multi_evaluate(a,x;derivatives=true);g=matched_values(v,ref)[active]-targets
        maximum(abs.(g))<2e-10 && break
        J=L.matched_jacobian(v,ref)[active,1:length(a.support)]
        scales=sqrt.(sum(abs2,J;dims=2))[:]
        x[1:length(a.support)].-=pinv(J./scales;rtol=1e-12)*(g./scales)
    end
    v=L.multi_evaluate(a,x;derivatives=true);au=audit(a,x,v)
    save_point("FEASIBLE_NUMERICAL_REPRESENTATIVE.toml",a,x,v,au;status="ROUNDING_INTERIOR_TO_ORIGINAL_CONSTRAINTS")
    println((retained_MW=v.J,cost_increment_MW=v.J-ref.J,primal=au.primal,stationarity=au.stationarity,beta=v.beta,F=v.Fpeak));flush(stdout)
    (;a,x,v,audit=au)
end
function robust_full_band(A;beta=L.BETA,maxnodes=50000)
    As=Matrix(A+.05I);n=size(A,1);started=time();tail=norm(As)+beta+1
    # High-precision quadratic Hermitian pencil near the almost-active peak.
    # H(w)=As'As + iw(As-As') +(w^2-beta^2)I. Its three degree-two
    # Bernstein coefficients being positive definite covers a whole interval.
    # One high precision congruence avoids catastrophic cancellation in A'A.
    # R=(-As)^-1; K=beta R; R'H(w)R=I-K'K+i(w/beta)(K-K')
    # +(w/beta)^2 K'K. Congruence preserves positive definiteness.
    Ab=BigFloat.(As);bb=BigFloat(beta);Rb=inv(-Ab);Kb=bb*Rb;Qb=Kb'*Kb
    inverse_residual_big=Float64(norm(I+Ab*Rb))
    K=Float64.(Kb);Q=Float64.(Qb);nK=opnorm(K);nQ=opnorm(Q)
    minimum_psd_margin=Ref(Inf)
    function big_psd(a,b)
        for (w,product) in ((a,a^2),((a+b)/2,a*b),(b,b^2))
            H=Matrix(I-Q+im*(w/beta)*(K-K')+(product/beta^2)*Q)
            allowance=128n*eps()*(1+nQ+2abs(w/beta)*nK+abs(product/beta^2)*nQ)
            low=eigmin(Hermitian(H))-allowance
            low>0 || return false
            minimum_psd_margin[]=min(minimum_psd_margin[],low)
        end
        true
    end
    near=.01;stack=[(0.,near,0)];near_nodes=0;depthmax=0
    while !isempty(stack)
        aa,zz,depth=pop!(stack);near_nodes+=1;depthmax=max(depthmax,depth)
        if big_psd(aa,zz);continue;end
        if near_nodes>maxnodes||depth>=60
            return Dict("status"=>"INCOMPLETE_NEAR_FREQUENCY_COVERAGE","near_nodes"=>near_nodes,"formal"=>false)
        end
        mid=(aa+zz)/2;push!(stack,(aa,mid,depth+1),(mid,zz,depth+1))
        near_nodes%100==0 && (println("ROBUST_NEAR nodes=$near_nodes pending=$(length(stack))");flush(stdout))
    end
    # Away from the active peak, validated residual bounds on inverse solves
    # plus the 1-Lipschitz singular-value bound cover each entire interval.
    function lower(w)
        M=Matrix(im*w*I-As);R=inv(M)
        residual=norm(I-M*R)+8n*eps()*norm(abs.(M)*abs.(R))
        max(0.,(1-residual)/(opnorm(R)*(1+16n*eps())))
    end
    function local_psd(aa,zz)
        c=(aa+zz)/2;M=Matrix(im*c*I-As);R=inv(M)
        # Account for the nonexact inverse in the congruence coefficients.
        residual=norm(I-M*R)+8n*eps()*norm(abs.(M)*abs.(R))
        Q=R'*R;nR=opnorm(R);nQ=opnorm(Q)
        for (delta,product) in ((aa-c,(aa-c)^2),(0.,(aa-c)*(zz-c)),(zz-c,(zz-c)^2))
            H=Matrix(I+im*delta*(R-R')+(product-beta^2)*Q)
            allowance=128n*eps()*(1+2abs(delta)*nR+(abs(product)+beta^2)*nQ)+2residual+residual^2+2abs(delta)*residual*nR
            eigmin(Hermitian(H))>allowance || return false
        end
        true
    end
    grid=sort(unique(vcat(near,10. .^range(log10(near),log10(tail),length=100),tail)))
    vals=lower.(grid);stack2=[(grid[k],grid[k+1],vals[k],vals[k+1],0) for k in 1:length(grid)-1]
    far_nodes=length(grid);unresolved=0
    while !isempty(stack2)
        aa,zz,fa,fb,depth=pop!(stack2)
        local_psd(aa,zz) && continue
        max(fa,fb)-(zz-aa)>=beta && continue
        (fa+fb-(zz-aa))/2>=beta && continue
        if far_nodes>=maxnodes||depth>=60;unresolved+=1;continue;end
        mid=(aa+zz)/2;fm=lower(mid);far_nodes+=1
        push!(stack2,(aa,mid,fa,fm,depth+1),(mid,zz,fm,fb,depth+1))
        far_nodes%1000==0 && (println("ROBUST_FAR nodes=$far_nodes pending=$(length(stack2))");flush(stdout))
    end
    Dict("status"=>unresolved==0 ? "NUMERICALLY_CERTIFIED_FULL_BAND" : "INCOMPLETE_FAR_FREQUENCY_COVERAGE",
        "beta_lower_tested"=>beta,"near_nodes"=>near_nodes,"far_nodes"=>far_nodes,
        "unresolved_intervals"=>unresolved,"tail_start_rad_s"=>tail,"near_cutoff_rad_s"=>near,
        "precision_bits_near"=>precision(BigFloat),"rounding"=>"nearest; NOT outward-rounded",
        "inverse_residual_big"=>inverse_residual_big,"minimum_congruence_PSD_margin"=>minimum_psd_margin[],
        "formal"=>false,"elapsed_s"=>time()-started,"maximum_depth_near"=>depthmax)
end
function time_full_audit(v;limitF=.5,limitR=.5,maxnodes=200000)
    started=time();rows=NamedTuple[];lam=v.lam
    # Real slowest mode must dominate the signed tail beyond a computable T.
    i=argmax(real.(lam));abs(imag(lam[i]))<1e-10 || error("Tail proof needs a real leading pole")
    tailtime=60.
    function tail_ok(ch,kind,t)
        co=v.R[ch,:].*L.tailcoeff(lam,kind==:F ? .5 : 1.,kind)
        shift=kind==:F ? .5 : 1.
        if kind==:F
            lead=real(co[i]);rest=sum(abs(co[j])*exp((real(lam[j])-real(lam[i]))*(t-shift)) for j in eachindex(lam) if j!=i)
            size=100sum(abs.(co).*exp.(real.(lam)*(t-shift)))
            return abs(v.dc[ch])+size<limitF || (v.dc[ch]<0 && lead>rest && size<abs(v.dc[ch]) && abs(v.dc[ch])<limitF)
        else
            return 100sum(abs.(co).*exp.(real.(lam)*(t-shift)))<limitR
        end
    end
    while !all(tail_ok(ch,kind,tailtime) for ch in 1:10 for kind in (:F,:R))
        tailtime*=1.5;tailtime>10000 && error("No verified monotone tail")
    end
    for kind in (:F,:R),ch in 1:10
        lim=kind==:F ? limitF : limitR;stack=[(0.,.5),(.5,1.),(1.,tailtime)];nodes=0;uppermax=0.;failed=false
        while !isempty(stack)
            aa,bb=pop!(stack);cc=(aa+bb)/2;h=(bb-aa)/2;nodes+=1
            f=L.output(v,cc,kind,ch);fp=L.output_dt(v,cc,kind,ch)
            # Analytic bound on the second time derivative on this segment.
            if aa>=(kind==:F ? .5 : 1.)
                second=100sum(abs.(v.R[ch,:].*lam.^2 .*L.tailcoeff(lam,aa,kind)))
            else
                ts=kind==:F ? [(1.,aa),(-1.,aa-.5)] : [(1.,aa),(-2.,aa-.5),(1.,aa-1.)]
                second=100sum(abs(c)*sum(abs.(v.R[ch,:]).*exp.(real.(lam)*t)) for (c,t) in ts if t>=0)/(kind==:F ? .5 : .25)
            end
            ub=abs(f)+h*abs(fp)+h^2*second/2+2e-11
            if ub<=lim;uppermax=max(uppermax,ub);continue;end
            if nodes>=maxnodes||h<1e-10;failed=true;break;end
            push!(stack,(aa,cc),(cc,bb))
        end
        push!(rows,(;kind=String(kind),bus=ch+29,nodes,upper_bound=uppermax,tail_start_s=tailtime,pass=!failed))
    end
    CSV.write(joinpath(L.OUT,"ALL_TIME_BOUNDS.csv"),DataFrame(rows))
    (;pass=all(r.pass for r in rows),nodes=sum(r.nodes for r in rows),tailtime,elapsed_s=time()-started,
        formal=false,scope="Analytic derivative and signed modal-tail bounds, Float64 numerical enclosure allowance")
end
