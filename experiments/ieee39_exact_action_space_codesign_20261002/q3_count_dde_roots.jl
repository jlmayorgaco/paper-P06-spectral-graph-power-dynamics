using CSV, DataFrames, LinearAlgebra, TOML

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const DC=include(joinpath(ROOT,"experiments","delay_dressed_replacement_frontier_20261002","DelayCharacteristic.jl"))
const R=DC.R

function contour(gamma,radius,step)
    corners=ComplexF64[gamma-im*radius,radius-im*radius,radius+im*radius,gamma+im*radius,gamma-im*radius]
    z=ComplexF64[]
    for k in 1:4
        a=corners[k];b=corners[k+1];n=max(2,ceil(Int,abs(b-a)/step))
        append!(z,[a+(b-a)*(j/n) for j in 0:n-1])
    end
    z
end

"Count roots in Re(s)>gamma with a finite contour from a matrix-norm root bound.
The arithmetic/phase count is Float64 numerical evidence, not a directed-rounding certificate.
"
function balance_diagonal(A;max_sweeps=100)
    d=ones(Float64,size(A,1)); n=length(d)
    for _ in 1:max_sweeps
        changed=false
        for i in 1:n
            r=0.0;c=0.0
            for j in 1:n
                i==j && continue
                r+=abs(A[i,j])*d[j]/d[i]
                c+=abs(A[j,i])*d[i]/d[j]
            end
            (r==0 || c==0) && continue
            f=2.0^round(Int,0.5*log2(r/c))
            if f!=1.0 && isfinite(d[i]*f)
                d[i]*=f;changed=true
            end
        end
        !changed && break
    end
    d
end

function balanced_action_model(L)
    Ac=L.A0+L.B*transpose(L.C); C=transpose(L.C)
    d=balance_diagonal(Ac)
    Ab=(Ac.*transpose(d))./reshape(d,:,1)
    Bb=L.B./reshape(d,:,1)
    Cb=C.*transpose(d)
    (;Ac=Ab,B=Bb,C=Cb,scale=d,original_norm=opnorm(Ac,2),balanced_norm=opnorm(Ab,2))
end

function transfer_tail_bound(Ac,B,C,tau,gamma)
    nA=opnorm(Ac,2); cB=opnorm(C*B,2); cA=opnorm(C*Ac,2); bnorm=opnorm(B,2)
    emax=maximum(1+exp(-gamma*tau[j]) for j in eachindex(tau))
    f(r)=emax*(cB/r+cA*bnorm/(r*(r-nA)))
    lo=nA+max(1.0,1e-10*nA); hi=max(2lo,lo+1.0)
    while f(hi)>=1 && hi<1e12; hi*=2; end
    f(hi)<1 || return (;Ac,nA,cB,cA,bnorm,emax,tail_value=f(hi),radius=Inf)
    for _ in 1:80
        mid=(lo+hi)/2
        f(mid)<1 ? (hi=mid) : (lo=mid)
    end
    (;Ac,B,C,nA,cB,cA,bnorm,emax,tail_value=f(hi),radius=hi+max(1e-6,1e-12*hi))
end

function count_roots(L,tau;gamma=-0.05,max_points=200000,max_refine=6)
    bal=balanced_action_model(L)
    tb=transfer_tail_bound(bal.Ac,bal.B,bal.C,tau,gamma)
    Ac=tb.Ac;radius=tb.radius
    info=(;original_Ac_norm=bal.original_norm,balanced_Ac_norm=bal.balanced_norm,
        diagonal_scale_ratio=maximum(bal.scale)/minimum(bal.scale),CB_norm=tb.cB,CA_norm=tb.cA,B_norm=tb.bnorm)
    isfinite(radius) || return merge((;status="INDETERMINATE_NO_FINITE_TAIL_BOUND",count=missing,gamma,
        bound=tb.tail_value,radius,points=0,max_phase_increment=Inf,min_logabs_small_factor=Inf,refinements=0),info)
    F=schur(ComplexF64.(Ac));T=F.T;Z=F.Z;mu=diag(T)
    left=bal.C*Z;right=adjoint(Z)*bal.B
    maxdelay=maximum(tau;init=0.0)
    step=maxdelay==0 ? min(20.0,radius/100) : min(20.0,pi/(4maxdelay))
    z=contour(gamma,radius,step)
    length(z)<=max_points || return merge((;status="INDETERMINATE_CONTOUR_TOO_LARGE",count=missing,gamma,
        bound=tb.tail_value,radius,points=length(z),max_phase_increment=Inf,
        min_logabs_small_factor=Inf,refinements=0),info)
    cache=Dict{ComplexF64,Tuple{ComplexF64,Float64}}();minlog=Inf;maxjump=0.0;maxdepth=0
    function phase_at(s)
        haskey(cache,s) && return cache[s]
        length(cache)>=max_points && error("adaptive contour point budget exceeded")
        Q=UpperTriangular(s*I-T);X=Q\right
        E=Diagonal(exp.(-s.*tau));small=Matrix{ComplexF64}(I,length(tau),length(tau))-(E-I)*left*X
        ld,sgn=logabsdet(small);minlog=min(minlog,ld)
        (iszero(sgn)||!isfinite(ld)) && error("reduced determinant singular/nonfinite on contour")
        ph=exp(im*(sum(angle.(s .- mu))+angle(sgn)))
        cache[s]=(ph,ld);(ph,ld)
    end
    function integrate_segment(a,pa,b,pb,depth)
        dphi=angle(pb/pa);maxjump=max(maxjump,abs(dphi));maxdepth=max(maxdepth,depth)
        if abs(dphi)<=pi/3;return dphi;end
        depth>=max_refine && error("phase increment remains unresolved after adaptive bisection")
        mid=(a+b)/2;pm=phase_at(mid)[1]
        integrate_segment(a,pa,mid,pm,depth+1)+integrate_segment(mid,pm,b,pb,depth+1)
    end
    try
        ph=[phase_at(s)[1] for s in z]
        wind=0.0
        for j in eachindex(z)
            k=mod1(j+1,length(z));wind+=integrate_segment(z[j],ph[j],z[k],ph[k],0)
        end
        winding=wind/(2pi);cnt=round(Int,winding)
        return merge((;status="NUMERICAL_CONTOUR_COUNT",count=cnt,gamma,bound=tb.tail_value,radius,
            points=length(cache),max_phase_increment=maxjump,min_logabs_small_factor=minlog,
            refinements=maxdepth,winding_unrounded=winding),info)
    catch err
        msg=sprint(showerror,err)
        status=occursin("budget exceeded",msg) ? "INDETERMINATE_CONTOUR_TOO_LARGE" : "INDETERMINATE_CONTOUR_RESOLUTION"
        return merge((;status,count=missing,gamma,bound=tb.tail_value,radius,points=length(cache),
            max_phase_increment=maxjump,min_logabs_small_factor=minlog,refinements=maxdepth,
            error=msg),info)
    end
end

function main()
    seed=TOML.parsefile(joinpath(@__DIR__,"seed_uniform_875.toml"))
    ctx=R.N.design_context(ROOT);rho=Float64.(seed["rho"]);kp=Float64.(seed["Kp"]);ki=Float64.(seed["Ki"])
    L=DC.linearization(ctx,rho,kp,ki)
    rows=NamedTuple[]
    for tau_ms in (0.0,20.0,40.0)
        tau=fill(tau_ms/1000,length(L.Ai))
        if tau_ms==0
            vals=eigvals(L.A0+L.B*transpose(L.C));cnt=count(x->real(x)>-0.05,vals)
            row=(;case_id="Q0_seed",tau_ms,status="NUMERICALLY_VALIDATED_ODE_COUNT",roots_right_of_margin=cnt,
                unverified_winding_index=missing,margin_s_inv=-0.05,critical_real_s_inv=maximum(real,vals),norm_bound=missing,contour_radius=missing,
                contour_points=0,max_phase_increment=0.0,min_logabs_small_factor=missing,refinements=0,
                original_Ac_norm=opnorm(L.A0+L.B*transpose(L.C),2),balanced_Ac_norm=missing,
                diagonal_scale_ratio=missing,CB_norm=missing,CA_norm=missing,B_norm=missing,
                contour_error="",certified=false,notes="tau=0 reduces to finite ODE; complete eigenvalue list counted")
        else
            println("Q3_COUNT_START tau_ms=",tau_ms);flush(stdout)
            max_points=isempty(ARGS) ? 200000 : parse(Int,ARGS[1])
            r=count_roots(L,tau;gamma=-0.05,max_points,max_refine=14)
            # A Float64 winding index is diagnostic only until independently refined
            # characteristic roots validate it. Never publish it as a root count.
            row=(;case_id="Q0_seed",tau_ms,status="INDETERMINATE",
                roots_right_of_margin=missing,
                unverified_winding_index=(r.status=="NUMERICAL_CONTOUR_COUNT" ? r.count : missing),
                margin_s_inv=-0.05,
                critical_real_s_inv=missing,norm_bound=r.bound,contour_radius=r.radius,contour_points=r.points,
                max_phase_increment=r.max_phase_increment,min_logabs_small_factor=r.min_logabs_small_factor,
                refinements=r.refinements,original_Ac_norm=r.original_Ac_norm,balanced_Ac_norm=r.balanced_Ac_norm,
                diagonal_scale_ratio=r.diagonal_scale_ratio,CB_norm=r.CB_norm,CA_norm=r.CA_norm,B_norm=r.B_norm,
                contour_error=hasproperty(r,:error) ? r.error : "",certified=false,
                notes="Unverified Float64 winding index; no independent root validation or directed-rounding proof")
        end
        push!(rows,row);CSV.write(joinpath(@__DIR__,"TABLE_Q03_DDE_ROOT_COUNTS.csv"),DataFrame(rows))
        println("Q3_COUNT_DONE ",row);flush(stdout)
    end
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    main()
end
