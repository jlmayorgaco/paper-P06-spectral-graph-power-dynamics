using LinearAlgebra, CSV, DataFrames, TOML, SciMLBase, OrdinaryDiffEqRosenbrock
include(joinpath(@__DIR__,"..","nonlinear_codesign_20261001","ReducedDAE.jl"))
const R=ReducedDAE
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
BLAS.set_num_threads(1)

"True constant-delay method of steps; no rational delay approximation."
function integrate_case(ctx,case_path,amplitude;reltol=1e-9,abstol=1e-11,dtmax=.0025,tag="standard")
    c=TOML.parsefile(joinpath(case_path,"case.toml"));key=c["design"];tau=c["tau"]
    d=TOML.parsefile(joinpath(@__DIR__,"models",key,"design.toml"))
    m=R.model(ctx,Float64.(d["rho"]),Float64.(d["Kp"]),Float64.(d["Ki"]);dc_convention=:physical_supply)
    tab=CSV.read(joinpath(case_path,"mode.csv"),DataFrame)
    v=ComplexF64.(tab.v_real.+im.*tab.v_imag);w=ComplexF64.(tab.w_real.+im.*tab.w_imag)
    s=c["real"]+im*c["imag"]
    xstar=copy(m.x0);n=length(xstar)
    B=Matrix(CSV.read(joinpath(@__DIR__,"models",key,"B.csv"),DataFrame))
    sols=Any[];ends=Float64[]
    function history(t)
        t<=0 && return xstar+real.(amplitude.*v.*exp(s*t))
        isempty(sols) && error("future history requested")
        j=clamp(searchsortedfirst(ends,t-1e-13),1,length(sols))
        t<=ends[end]+1e-10 || error("method-of-steps requested uncomputed history")
        sols[j](clamp(t,sols[j].t[1],sols[j].t[end]))
    end
    function fun!(dx,x,p,t)
        R.rhs!(dx,x,m,t)
        if tau>0
            xp=history(t-tau);up=R.voltage(xp,m);u=R.voltage(x,m)
            for i in 1:10
                ix=m.gfidx[i];th=x[ix[3]];tp=xp[ix[3]]
                e=-sin(th)*u[2i-1]+cos(th)*u[2i]
                ep=-sin(tp)*up[2i-1]+cos(tp)*up[2i]
                dx[ix[4]]+=m.kp[i]/m.gp[i].pll_tau*(ep-e)
                dx[ix[5]]+=m.ki[i]*(ep-e)
            end
        end
        nothing
    end
    function jac!(J,x,p,t)
        de=R.derivatives(x,m);J.=de.Fx
        if tau>0
            C=zeros(10,n)
            for i in 1:10
                ix=m.gfidx[i][3];th=x[ix]
                C[i,:].=-sin(th).*de.vx[2i-1,:]+cos(th).*de.vx[2i,:]
                C[i,ix]+=-cos(th)*de.v[2i-1]-sin(th)*de.v[2i]
            end
            J.-=B*C
        end
        nothing
    end
    t0=0.;x=history(0.);horizon=3.
    fn=ODEFunction(fun!;jac=jac!)
    while t0<horizon-1e-12
        t1=min(horizon,tau>0 ? t0+tau : horizon)
        prob=ODEProblem(fn,x,(t0,t1))
        sol=solve(prob,Rodas5P();reltol,abstol,dtmax,dense=true,save_everystep=true,maxiters=1000000)
        SciMLBase.successful_retcode(sol) || error("integration failed $(sol.retcode)")
        push!(sols,sol);push!(ends,t1);x=copy(sol.u[end]);t0=t1
    end
    rows=NamedTuple[];maxerr=0.;maxref=0.;maxpll=0.
    pidx=[ix[3] for ix in m.gfidx];oidx=[ix[end-1] for ix in m.sgidx]
    for t in 0:.002:horizon
        y=history(t)-xstar;ref=real.(amplitude.*v.*exp(s*t));z=dot(w,y);zr=dot(w,ref)
        maxerr=max(maxerr,norm(y[pidx]-ref[pidx]));maxref=max(maxref,norm(ref[pidx]));maxpll=max(maxpll,norm(y[pidx],Inf))
        push!(rows,(;time=t,modal_real=real(z),modal_imag=imag(z),reference_real=real(zr),reference_imag=imag(zr),
            pll_error_norm=norm(y[pidx]-ref[pidx]),pll_reference_norm=norm(ref[pidx]),
            max_pll_angle=norm(y[pidx],Inf),sg_speed_error_norm=norm(y[oidx]-ref[oidx])))
    end
    out=joinpath(@__DIR__,"nonlinear");mkpath(out)
    label=basename(case_path)*"_a"*string(amplitude)*"_"*tag
    CSV.write(joinpath(out,label*".csv"),DataFrame(rows))
    (;case=basename(case_path),design=key,tau_ms=1000tau,amplitude,tag,
        predicted_alpha=real(s),predicted_frequency_hz=imag(s)/(2pi),
        max_pll_absolute_error=maxerr,max_pll_reference_norm=maxref,relative_trajectory_error=maxerr/maxref,
        max_pll_angle=maxpll,intervals=length(sols),reltol,abstol,dtmax,solver="METHOD_OF_STEPS_RODAS5P",
        status="LOCAL_NONLINEAR_DDE_TEST_NOT_LARGE_EVENT_SECURITY")
end

function main()
    ctx=R.N.design_context(ROOT);rows=NamedTuple[]
    for name in sort(readdir(joinpath(@__DIR__,"nonlinear_inputs")))
        path=joinpath(@__DIR__,"nonlinear_inputs",name)
        for amplitude in (1e-5,1e-4)
            row=integrate_case(ctx,path,amplitude);push!(rows,row)
            CSV.write(joinpath(@__DIR__,"TABLE_07_NONLINEAR_DDE.csv"),DataFrame(rows))
            println("NONLINEAR_DONE ",row);flush(stdout)
            if row.relative_trajectory_error>.02
                fine=integrate_case(ctx,path,amplitude;reltol=1e-10,abstol=1e-12,dtmax=.00125,tag="refined")
                push!(rows,fine);CSV.write(joinpath(@__DIR__,"TABLE_07_NONLINEAR_DDE.csv"),DataFrame(rows))
                println("NONLINEAR_REFINED ",fine);flush(stdout)
            end
        end
    end
end
main()
