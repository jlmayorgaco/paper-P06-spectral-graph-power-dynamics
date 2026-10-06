using CSV, DataFrames, LinearAlgebra, TOML

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const N=include(joinpath(ROOT,"src","bnd_model_expN","PDExactDesignN.jl"))
const SIGMA=0.05

function geometry(ctx,m,kp,ki)
    n=m.n_dynamic;ny=m.n_algebraic;out=NamedTuple[]
    for bus in 30:39
        i=bus-29
        row=only(eachrow(m.state_map[(m.state_map.bus.==bus).&(m.state_map.kind.=="GFL"),:]))
        firstix=Int(row.first);lastix=Int(row.last);ix=firstix:lastix
        op=N.trim_gfl(ctx,bus);theta=op.x[3];ur=op.u[1];ui=op.u[2]
        cx=zeros(n);cx[firstix+2]=-cos(theta)*ur-sin(theta)*ui
        cy=zeros(ny);cy[2bus-1]=-sin(theta);cy[2bus]=cos(theta)
        bp=zeros(n);bi=zeros(n);bp[firstix+3]=1/op.pars.pll_tau;bi[firstix+4]=1.0
        b=kp[i].*bp.+ki[i].*bi
        push!(out,(;bus,i,ix,bp,bi,b,cx,cy))
    end
    out
end

function descriptor_dde_pencil(ctx,rho,kp,ki,tau,s,m=N.descriptor(ctx,rho,kp,ki))
    n=m.n_dynamic;ny=m.n_algebraic;geom=geometry(ctx,m,kp,ki)
    A0=copy(m.A);B0=copy(m.B)
    for g in geom
        A0.-=g.b*transpose(g.cx)
        B0.-=g.b*transpose(g.cy)
    end
    Ac=ComplexF64.(A0);Bc=ComplexF64.(B0)
    for g in geom
        fac=exp(-s*tau[g.i])
        Ac.+=fac.*(g.b*transpose(g.cx))
        Bc.+=fac.*(g.b*transpose(g.cy))
    end
    topL=Matrix{ComplexF64}(s*I(n)-Ac)
    topR=-Bc
    bottomL=ComplexF64.(m.C);bottomR=ComplexF64.(m.Gy)
    T=[topL topR;bottomL bottomR]
    schur=topL-topR*(bottomR\bottomL)
    Cport=m.Gy\m.C
    cbar=[g.cx-vec(transpose(g.cy)*Cport) for g in geom]
    Ared0=copy(m.Ared)
    for (g,cr) in zip(geom,cbar);Ared0.-=g.b*transpose(cr);end
    direct=Matrix{ComplexF64}(s*I(n)-Ared0)
    for (g,cr) in zip(geom,cbar);direct.-=exp(-s*tau[g.i]).*(g.b*transpose(cr));end
    schurerr=norm(schur-direct)/max(norm(schur),eps())
    (;T,m,geom,n,ny,schurerr)
end

function roots_of_quadratic(trM,detM)
    if abs(detM)<1e-13
        return abs(trM)<1e-13 ? ComplexF64[] : ComplexF64[-1/trM]
    end
    disc=sqrt(complex(trM^2-4detM))
    ComplexF64[(-trM+disc)/(2detM),(-trM-disc)/(2detM)]
end

function validate_q5(ctx,rho,kp,ki)
    buses=collect(30:39);freqs=[0.017,0.5,3.0];delays=[0.0,0.020];h=0.01
    rows=NamedTuple[];nrow=0
    for tau0 in delays, bus in buses, f in freqs
        i=bus-29;s=-SIGMA+2pi*im*f;tau=fill(tau0,10);rho0=copy(rho);eps0=1-rho0[i]
        m0=N.descriptor(ctx,rho0,kp,ki);T0=descriptor_dde_pencil(ctx,rho0,kp,ki,tau,s,m0).T
        rp=copy(rho0);rm=copy(rho0);rp[i]+=h;rm[i]-=h
        mp=N.descriptor(ctx,rp,kp,ki);mm=N.descriptor(ctx,rm,kp,ki)
        Tp=descriptor_dde_pencil(ctx,rp,kp,ki,tau,s,mp).T
        Tm=descriptor_dde_pencil(ctx,rm,kp,ki,tau,s,mm).T
        dTdE=-(Tp-Tm)/(2h)
        n=m0.n_dynamic;ny=m0.n_algebraic;rr=n .+ (2bus-1:2bus)
        U=zeros(ComplexF64,n+ny,2);U[rr[1],1]=1;U[rr[2],2]=1
        W=adjoint(U)*dTdE
        row_support_res=norm(dTdE-U*W)/max(norm(dTdE),eps())
        affine_res=norm(Tp+Tm-2T0)/max(norm(T0),eps())
        singrow=svdvals(dTdE[rr,:]); rank2=count(x->x>1e-10*maximum(singrow),singrow)
        Tref=T0-eps0.*(U*W)
        H=W*(Tref\U);trH=tr(H);detH=det(H);candidates=roots_of_quadratic(trH,detH)
        for (q,z) in enumerate(candidates)
            physical=abs(imag(z))<1e-8 && 0<=real(z)<=1
            full_res=NaN
            if physical
                rs=copy(rho0);rs[i]=1-real(z)
                ms=N.descriptor(ctx,rs,kp,ki)
                Ts=descriptor_dde_pencil(ctx,rs,kp,ki,tau,s,ms).T
                full_res=minimum(svdvals(Ts))/max(opnorm(Ts),eps())
            end
            pol_res=abs(1+z*trH+z^2*detH)
            nrow+=1
            push!(rows,(;bus,tau_ms=1000*tau0,target_frequency_Hz=f,target_real_s_inv=-SIGMA,
                root_index=q,epsilon_root=real(z),epsilon_imag=imag(z),rho_root=physical ? 1-real(z) : NaN,
                physical_root_0_to_1=physical,rank_update_numerical=rank2,
                update_row_support_relative_residual=row_support_res,affine_second_difference_relative_residual=affine_res,
                schur_reduction_relative_residual=descriptor_dde_pencil(ctx,rho0,kp,ki,tau,s,m0).schurerr,
                quadratic_root_residual=pol_res,full_descriptor_relative_singular_residual=full_res,
                trace_N=trH,det_N=detH))
        end
    end
    CSV.write(joinpath(@__DIR__,"TABLE_Q05_RHO_CLOSED_FORM_VALIDATION.csv"),DataFrame(rows))
    (;rows,buses,freqs,delays)
end

function validate_q6(ctx,rho,kp,ki)
    rho2=copy(rho);kp2=copy(kp);ki2=copy(ki)
    for i in 1:10
        rho2[i]+=0.005*(-1)^i
        kp2[i]*=1+0.02*cos(i)
        ki2[i]*=1+0.02*sin(i)
    end
    delays=[0.0,0.020];freqs=[0.017,0.5,3.0];rows=NamedTuple[]
    for tau0 in delays,f in freqs
        s=-SIGMA+2pi*im*f;tau=fill(tau0,10)
        m0=N.descriptor(ctx,rho,kp,ki);m1=N.descriptor(ctx,rho2,kp2,ki2)
        p0=descriptor_dde_pencil(ctx,rho,kp,ki,tau,s,m0);p1=descriptor_dde_pencil(ctx,rho2,kp2,ki2,tau,s,m1)
        dT=p1.T-p0.T;n=m0.n_dynamic;ny=m0.n_algebraic;geom=p0.geom
        U=zeros(ComplexF64,n+ny,30);H=zeros(ComplexF64,30,n+ny)
        for g in geom
            db=(kp2[g.i]-kp[g.i]).*g.bp.+(ki2[g.i]-ki[g.i]).*g.bi
            U[1:n,g.i]=db
            H[g.i,:].=-exp(-s*tau[g.i]).*vcat(g.cx,g.cy)
        end
        for i in 1:10
            bus=29+i
            for q in 1:2
                row=n+2bus-2+q;col=10+2i-2+q
                U[row,col]=1.0;H[col,:]=dT[row,:]
            end
        end
        reconstruction=U*H
        err=norm(dT-reconstruction)/max(norm(dT),eps())
        sv=svdvals(dT);rk=count(x->x>1e-9*maximum(sv),sv)
        A0m1=copy(m1.A);A0m0=copy(m0.A);B0m1=copy(m1.B);B0m0=copy(m0.B)
        for g in p1.geom;A0m1.-=g.b*transpose(g.cx);B0m1.-=g.b*transpose(g.cy);end
        for g in p0.geom;A0m0.-=g.b*transpose(g.cx);B0m0.-=g.b*transpose(g.cy);end
        coreerr=max(norm(A0m1-A0m0)/max(norm(A0m0),eps()),norm(B0m1-B0m0)/max(norm(B0m0),eps()))
        push!(rows,(;tau_ms=1000*tau0,target_frequency_Hz=f,full_descriptor_dimension=n+ny,
            analytical_action_dimension_upper_bound=30,numerical_difference_rank=rk,
            explicit_action_factorization_relative_residual=err,PLL_removed_core_invariance_error=coreerr,
            Schur_reduction_relative_residual=p0.schurerr,
            similarity_or_gauge_removed=false,claim_class="exact augmented descriptor update; numerical factorization check"))
    end
    CSV.write(joinpath(@__DIR__,"TABLE_Q06_MASTER_ACTION_SPACE.csv"),DataFrame(rows))
    rows
end

function main()
    seed=TOML.parsefile(joinpath(@__DIR__,"seed_uniform_875.toml"))
    ctx=N.design_context(ROOT);rho=Float64.(seed["rho"]);kp=Float64.(seed["Kp"]);ki=Float64.(seed["Ki"])
    q5=validate_q5(ctx,rho,kp,ki);q6=validate_q6(ctx,rho,kp,ki)
    q5status=Dict("status"=>"COMPLETED_CONDITIONAL_RANK_TWO_TEST",
        "tests"=>length(q5.rows),"candidate_buses"=>q5.buses,"tau_ms"=>1000 .*q5.delays,
        "target_frequencies_Hz"=>q5.freqs,
        "max_update_row_support_relative_residual"=>maximum(getproperty.(q5.rows,:update_row_support_relative_residual)),
        "max_affine_second_difference_relative_residual"=>maximum(getproperty.(q5.rows,:affine_second_difference_relative_residual)),
        "physical_conditional_roots"=>count(getproperty.(q5.rows,:physical_root_0_to_1)),
        "max_full_descriptor_boundary_residual"=>(isempty(filter(isfinite,getproperty.(q5.rows,:full_descriptor_relative_singular_residual))) ? "not_applicable_no_physical_real_roots" : maximum(filter(isfinite,getproperty.(q5.rows,:full_descriptor_relative_singular_residual)))),
        "max_schur_reduction_relative_residual"=>maximum(getproperty.(q5.rows,:schur_reduction_relative_residual)))
    q6status=Dict("status"=>"COMPLETED_NUMERICAL_MASTER_ACTION_FACTOR_TEST",
        "tests"=>length(q6),"action_rank_upper_bound"=>30,
        "max_factorization_relative_residual"=>maximum(getproperty.(q6,:explicit_action_factorization_relative_residual)),
        "max_core_invariance_error"=>maximum(getproperty.(q6,:PLL_removed_core_invariance_error)),
        "max_schur_reduction_relative_residual"=>maximum(getproperty.(q6,:Schur_reduction_relative_residual)),
        "max_numerical_difference_rank"=>maximum(getproperty.(q6,:numerical_difference_rank)))
    open(joinpath(@__DIR__,"Q5_STATUS.toml"),"w") do io;TOML.print(io,q5status);end
    open(joinpath(@__DIR__,"Q6_STATUS.toml"),"w") do io;TOML.print(io,q6status);end
    println("Q5 ",q5status);println("Q6 ",q6status)
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
