using LinearAlgebra, CSV, DataFrames, TOML, ForwardDiff, SHA, Pkg
include(joinpath(@__DIR__,"..","graph_gsp_codesign_20261003","DelayedEvents.jl"))
module ContourOracle
include(joinpath(@__DIR__,"..","analytical_delay_codesign_mega_20261002","m3_a_trace_integral.jl"))
end
const R=DelayedEvents.R
const ROOT=R.ROOT
const OUT=@__DIR__
BLAS.set_num_threads(1)
write_matrix(path,A)=CSV.write(path,DataFrame(A,:auto))
read_matrix(name)=Matrix(CSV.read(joinpath(OUT,"..","graph_gsp_codesign_20261003","model",name*".csv"),DataFrame))

function main()
    ctx=R.N.design_context(ROOT)
    version_rows=[(;package="Julia",version=string(VERSION))]
    for (_,info) in Pkg.dependencies()
        info.name in ("PowerDynamics","NetworkDynamics","ForwardDiff","SciMLBase","OrdinaryDiffEqRosenbrock") || continue
        push!(version_rows,(;package=info.name,version=string(info.version)))
    end
    CSV.write(joinpath(OUT,"PACKAGE_VERSIONS.csv"),DataFrame(version_rows))
    names=("Adev","Bv","Cs","Cf","Ds","Y","Etheta","Hv","Bp","Bi")
    param=Dict(name=>read_matrix(name) for name in names)
    parity=NamedTuple[];spectra=NamedTuple[]
    for id in ("baseline","analytic","fixed")
        started=time();path=joinpath(OUT,"designs",id*".toml");d=TOML.parsefile(path)
        rho,kp,ki=Float64.(d["rho"]),Float64.(d["Kp"]),Float64.(d["Ki"])
        m=R.model(ctx,rho,kp,ki;dc_convention=:physical_supply)
        f=zeros(length(m.x0));R.rhs!(f,m.x0,m,0.)
        A=ForwardDiff.jacobian(m.x0) do x
            dx=similar(x);R.rhs!(dx,x,m,0.);dx
        end
        C=ForwardDiff.jacobian(x->DelayedEvents.detector(x,m),m.x0)
        B=DelayedEvents.injection(m);n=length(m.x0)
        weight=repeat(1 .-rho,inner=2)
        Vk=-(param["Y"]+Diagonal(weight)*param["Ds"])\
            (Diagonal(weight)*param["Cs"]+Diagonal(1 .-weight)*param["Cf"])
        Bpar=param["Bp"]*Diagonal(kp)+param["Bi"]*Diagonal(ki)
        Cpar=param["Etheta"]+param["Hv"]*Vk
        Apar=param["Adev"]+param["Bv"]*Vk+Bpar*Cpar
        jacerr=norm(A-Apar)/norm(A);deterr=norm(C-Cpar)/norm(C)
        g,_=R.rotation_generator(m.x0,m;jacobian=false);g./=norm(g)
        gaugeerr=norm(A*g)/norm(A);cgauge=norm(C*g)
        pass=norm(f,Inf)<=1e-7 && jacerr<=1e-8 && deterr<=1e-8 &&
             gaugeerr<1e-10 && cgauge<1e-10
        push!(parity,(;design=id,equilibrium_inf=norm(f,Inf),jacobian_relative=jacerr,
            detector_relative=deterr,gauge_relative=gaugeerr,detector_gauge=cgauge,
            GFL_MW=dot(ctx.power,rho),SG_MW=dot(ctx.power,1 .-rho),pass,
            design_sha256=bytes2hex(sha256(read(path)))))
        CSV.write(joinpath(OUT,"TABLE_06_INDEPENDENT_MODEL.csv"),DataFrame(parity))
        dest=joinpath(OUT,"independent",id);mkpath(dest)
        for (name,mat) in (("A",A),("A0",A-B*C),("B",B),("C",C))
            write_matrix(joinpath(dest,name*".csv"),mat)
        end
        println("INDEPENDENT_PARITY ",last(parity));flush(stdout)
        pass || error("BLOCKED_INDEPENDENT_PARITY")
        Q=nullspace(reshape(g,1,:))
        L=(;A0=Q'*(A-B*C)*Q,B=Q'*B,C=Matrix(transpose(C*Q)))
        println("CONTOUR_START ",id);flush(stdout)
        row=try
            c=ContourOracle.trace_count(L,fill(.04,10);step=10.,max_depth=22)
            resolved=c.integer_distance<.01 && c.quadrature_error_estimate<.01 &&
                abs(c.phase_estimate-c.nearest)<.01 && c.largest_phase_step<pi/2
            (;design=id,status=resolved ? (c.nearest==0 ? "PASS_NUMERICAL" : "FAIL_ROOTS_RIGHT_OF_MARGIN") : "INDETERMINATE",
              count=c.nearest,count_real=real(c.estimate),count_imag=imag(c.estimate),
              phase_count=c.phase_estimate,integer_distance=c.integer_distance,
              quadrature_error=c.quadrature_error_estimate,largest_phase_step=c.largest_phase_step,
              min_sigma=c.min_sigma_small,radius=c.radius,evaluations=c.evaluations,
              seconds=time()-started,error="")
        catch err
            (;design=id,status="INDETERMINATE",count=-1,count_real=NaN,count_imag=NaN,
              phase_count=NaN,integer_distance=NaN,quadrature_error=NaN,
              largest_phase_step=NaN,min_sigma=NaN,radius=NaN,evaluations=0,
              seconds=time()-started,error=sprint(showerror,err))
        end
        push!(spectra,row);CSV.write(joinpath(OUT,"TABLE_07_COMPLETE_SPECTRAL_REGION.csv"),DataFrame(spectra))
        println("CONTOUR_DONE ",row);flush(stdout)
    end
    only(filter(x->x.design=="baseline",spectra)).status=="PASS_NUMERICAL" ||
        error("BLOCKED_BASELINE_SPECTRAL_PARITY")
    only(filter(x->x.design=="analytic",spectra)).status=="PASS_NUMERICAL" || return
    d=TOML.parsefile(joinpath(OUT,"designs","analytic.toml"))
    rho,kp,ki=Float64.(d["rho"]),Float64.(d["Kp"]),Float64.(d["Ki"])
    all(ctx.kpmin .<=kp.<=ctx.kpmax) && all(ctx.kimin .<=ki.<=ctx.kimax) || error("GAIN_LIMIT")
    events=NamedTuple[];mkpath(joinpath(OUT,"nonlinear"))
    for (bus,delta) in DelayedEvents.CASES
        println("EVENT_START analytic ",bus," ",delta);flush(stdout)
        m=R.model(ctx,rho,kp,ki;bus,delta,dc_convention=:physical_supply)
        row=try
            met,r=DelayedEvents.adaptive(m;tau=.04)
            pass=met.Fpeak_Hz<=.5 && met.Rpeak_Hz_s<=.5 && met.Vmin>=.9 &&
                 met.Vmax<=1.1 && met.limiter_fraction>=.002
            R.metrics(m,r;dt=.01,window=.5,monitor_buses=collect(1:39),
                savepath=joinpath(OUT,"nonlinear","bus$(bus)_$(Int(delta))_trajectory.csv"))
            (;design="analytic",bus,delta,F=met.Fpeak_Hz,R=met.Rpeak_Hz_s,Vmin=met.Vmin,
              Vmax=met.Vmax,slack=met.limiter_fraction,runtime=met.runtime_s,pass,error="")
        catch err
            (;design="analytic",bus,delta,F=NaN,R=NaN,Vmin=NaN,Vmax=NaN,slack=NaN,
              runtime=NaN,pass=false,error=sprint(showerror,err))
        end
        push!(events,row);CSV.write(joinpath(OUT,"TABLE_08_NONLINEAR_EVENTS.csv"),DataFrame(events))
        println("EVENT_DONE ",row);flush(stdout)
    end
    println("VALIDATION_COMPLETE fully_feasible=",all(x.pass for x in events));flush(stdout)
end
main()
