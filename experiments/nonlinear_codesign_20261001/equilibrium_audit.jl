include("ReducedDAE.jl")
using .ReducedDAE, LinearAlgebra, CSV, DataFrames, TOML, Test, SHA
const R=ReducedDAE
BLAS.set_num_threads(1);const CTX=R.N.design_context(R.ROOT)
function phase_symmetry(x,m)
    q=zeros(length(x));Dq=zeros(length(x),length(x))
    for i in 1:10
        !isempty(m.sgidx[i]) && (q[last(m.sgidx[i])]=1.)
        if !isempty(m.gfidx[i])
            ix=m.gfidx[i];q[ix[3]]=1.;q[ix[6]]=x[ix[7]];q[ix[7]]=-x[ix[6]]
            Dq[ix[6],ix[7]]=1.;Dq[ix[7],ix[6]]=-1.
        end
    end
    q,Dq
end
function relative_equilibrium(m;initial=nothing)
    x=initial===nothing ? copy(m.x0) : copy(initial.x)
    omega=initial===nothing ? 0. : initial.omega;nx=length(x)
    gauge=last(last(filter(!isempty,m.sgidx)));resid=Inf;iterations=0
    function residual(x,omega)
        f=zeros(nx);R.rhs!(f,x,m,0.);q,Dq=phase_symmetry(x,m)
        vcat(f-omega*q,x[gauge]-m.x0[gauge])
    end
    for k in 1:30
        iterations=k;rr=residual(x,omega);resid=norm(rr,Inf)
        resid<1e-7 && break
        de=R.derivatives(x,m);q,Dq=phase_symmetry(x,m)
        J=zeros(nx+1,nx+1);J[1:nx,1:nx]=de.Fx-omega*Dq;J[1:nx,end]=-q;J[end,gauge]=1.
        step=J\rr;ok=false
        for ls in 0:14
            a=2.0^-ls;trial=x-a*step[1:nx];w=omega-a*step[end]
            if norm(residual(trial,w),Inf)<resid
                x=trial;omega=w;ok=true;break
            end
        end
        ok || break
    end
    rr=residual(x,omega);resid=norm(rr,Inf)
    q,Dq=phase_symmetry(x,m);Q=nullspace(reshape(q/norm(q),1,:));A=R.derivatives(x,m).Fx-omega*Dq
    alpha=maximum(real.(eigvals(Q'*A*Q)));v=R.voltage(x,m;allbus=true)
    # Explicit check that changing PLL gains preserves the same relative equilibrium.
    alternative=merge(m,(;kp=clamp.(1.3m.kp,CTX.kpmin,CTX.kpmax),ki=clamp.(.7m.ki,CTX.kimin,CTX.kimax)))
    ff=zeros(nx);R.rhs!(ff,x,alternative,0.)
    gain_residual=norm(ff-omega*q,Inf)
    (;x,omega,resid,iterations,alpha,gain_residual,Vmin=minimum(hypot.(v[1:2:end],v[2:2:end])),Vmax=maximum(hypot.(v[1:2:end],v[2:2:end])))
end
function main()
    path=isempty(ARGS) ? nothing : abspath(ARGS[1])
    d=path===nothing ? Dict("rho"=>fill(.8,10),"Kp"=>fill(R.N.K0P,10),"Ki"=>fill(R.N.K0I,10)) : TOML.parsefile(path)
    rows=NamedTuple[]
    for bus in (8,16,29),delta in (100.,-100.)
        result=nothing
        # Follow the post-event equilibrium branch. A direct large Newton step
        # can cross an actuator clamp and make its Jacobian singular even when
        # the desired nonsaturated equilibrium is regular.
        for fraction in range(0.,1.;length=11)
            m=R.model(CTX,d["rho"],d["Kp"],d["Ki"];bus,delta=delta*fraction,dc_convention=Symbol(get(d,"dc_convention","legacy")))
            result=relative_equilibrium(m;initial=result)
            result.resid<1e-6 || error("Equilibrium continuation failed at fraction $fraction")
        end
        row=(;bus,delta,frequency_offset_Hz=result.omega/(2pi),residual=result.resid,
            changed_PLL_residual=result.gain_residual,physical_alpha=result.alpha,
            result.Vmin,result.Vmax,iterations=result.iterations)
        push!(rows,row);println("RELATIVE_EQUILIBRIUM ",row);flush(stdout)
        @test result.resid<1e-6
        @test result.gain_residual<1e-6
    end
    label=length(ARGS)>1 ? ARGS[2] : (path===nothing ? "equilibrium_anchor" : "equilibrium_final")
    dest=joinpath(R.OUT,label*".csv")
    CSV.write(dest,DataFrame(rows))
    if path!==nothing
        open(joinpath(R.OUT,label*"_provenance.toml"),"w") do io
            TOML.print(io,Dict("candidate_sha256"=>bytes2hex(sha256(read(path))),
                "scope"=>"Numerical relative-equilibrium solves and local quotient spectra. No region-of-attraction certificate."))
        end
    end
end
if abspath(PROGRAM_FILE)==(@__FILE__);main();end
