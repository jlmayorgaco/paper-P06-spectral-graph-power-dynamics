using LinearAlgebra, CSV, DataFrames, TOML, ForwardDiff, SHA, Pkg
include(joinpath(@__DIR__,"..","graph_gsp_codesign_20261003","DelayedEvents.jl"))
module ContourOracle
include(joinpath(@__DIR__,"..","analytical_delay_codesign_mega_20261002","m3_a_trace_integral.jl"))
end
const R=DelayedEvents.R
const OUT=@__DIR__
BLAS.set_num_threads(1)
read_matrix(name)=Matrix(CSV.read(joinpath(OUT,"..","graph_gsp_codesign_20261003","model",name*".csv"),DataFrame))
function main()
    ctx=R.N.design_context(R.ROOT)
    versions=[(;package="Julia",version=string(VERSION))]
    for (_,info) in Pkg.dependencies()
        info.name in ("PowerDynamics","NetworkDynamics","ForwardDiff","SciMLBase","OrdinaryDiffEqRosenbrock") || continue
        push!(versions,(;package=info.name,version=string(info.version)))
    end
    CSV.write(joinpath(OUT,"JULIA_VERSIONS.csv"),DataFrame(versions))
    names=("Adev","Bv","Cs","Cf","Ds","Y","Etheta","Hv","Bp","Bi")
    param=Dict(name=>read_matrix(name) for name in names)
    parity=NamedTuple[];spectra=NamedTuple[]
    for id in ("anchor","single30","single37","joint","complex_pair","corrected")
        path=joinpath(OUT,"designs",id*".toml");d=TOML.parsefile(path)
        rho,kp,ki=Float64.(d["rho"]),Float64.(d["Kp"]),Float64.(d["Ki"])
        all(ctx.kpmin .<=kp.<=ctx.kpmax) && all(ctx.kimin .<=ki.<=ctx.kimax) || error("GAIN_LIMIT")
        m=R.model(ctx,rho,kp,ki;dc_convention=:physical_supply)
        f=zeros(length(m.x0));R.rhs!(f,m.x0,m,0.)
        A=ForwardDiff.jacobian(m.x0) do x
            dx=similar(x);R.rhs!(dx,x,m,0.);dx
        end
        C=ForwardDiff.jacobian(x->DelayedEvents.detector(x,m),m.x0)
        B=DelayedEvents.injection(m);w=repeat(1 .-rho,inner=2)
        Vk=-(param["Y"]+Diagonal(w)*param["Ds"])\(Diagonal(w)*param["Cs"]+Diagonal(1 .-w)*param["Cf"])
        Bp=param["Bp"]*Diagonal(kp)+param["Bi"]*Diagonal(ki)
        Cp=param["Etheta"]+param["Hv"]*Vk
        Ap=param["Adev"]+param["Bv"]*Vk+Bp*Cp
        g,_=R.rotation_generator(m.x0,m;jacobian=false);g./=norm(g)
        jacerr=norm(A-Ap)/norm(A);deterr=norm(C-Cp)/norm(C)
        pass=norm(f,Inf)<1e-7 && jacerr<1e-8 && deterr<1e-8 && norm(C*g)<1e-10
        push!(parity,(;design=id,equilibrium=norm(f,Inf),jacobian=jacerr,detector=deterr,pass,
            GFL_MW=dot(ctx.power,rho),SG_MW=dot(ctx.power,1 .-rho),sha256=bytes2hex(sha256(read(path)))))
        CSV.write(joinpath(OUT,"TABLE_10_JULIA_PARITY.csv"),DataFrame(parity));pass || error("PARITY")
        Q=nullspace(reshape(g,1,:));L=(;A0=Q'*(A-B*C)*Q,B=Q'*B,C=Matrix(transpose(C*Q)))
        println("SPECTRAL_START ",id);flush(stdout)
        row=try
            c=ContourOracle.trace_count(L,fill(.04,10);step=10.,max_depth=22)
            resolved=c.integer_distance<.01 && c.quadrature_error_estimate<.01 && abs(c.phase_estimate-c.nearest)<.01 && c.largest_phase_step<pi/2
            (;design=id,status=resolved ? (c.nearest==0 ? "PASS_NUMERICAL" : "FAIL_MARGIN") : "INDETERMINATE",
                count=c.nearest,count_real=real(c.estimate),count_imag=imag(c.estimate),phase_count=c.phase_estimate,
                quadrature_error=c.quadrature_error_estimate,radius=c.radius,evaluations=c.evaluations,error="")
        catch err
            (;design=id,status="INDETERMINATE",count=-1,count_real=NaN,count_imag=NaN,phase_count=NaN,
                quadrature_error=NaN,radius=NaN,evaluations=0,error=sprint(showerror,err))
        end
        push!(spectra,row);CSV.write(joinpath(OUT,"TABLE_11_FULL_SPECTRUM.csv"),DataFrame(spectra))
        println("SPECTRAL_DONE ",row);flush(stdout)
    end
    only(filter(x->x.design=="corrected",spectra)).status=="PASS_NUMERICAL" || return
    d=TOML.parsefile(joinpath(OUT,"designs","corrected.toml"))
    rho,kp,ki=Float64.(d["rho"]),Float64.(d["Kp"]),Float64.(d["Ki"])
    rows=NamedTuple[];mkpath(joinpath(OUT,"nonlinear"))
    for (bus,delta) in DelayedEvents.CASES
        println("EVENT_START ",bus," ",delta);flush(stdout)
        m=R.model(ctx,rho,kp,ki;bus,delta,dc_convention=:physical_supply)
        row=try
            met,r=DelayedEvents.adaptive(m;tau=.04,horizon=60.,dtmax=.01,tol=1e-9)
            pass=met.Fpeak_Hz<=.5 && met.Rpeak_Hz_s<=.5 && met.Vmin>=.9 && met.Vmax<=1.1 && met.limiter_fraction>=.002
            R.metrics(m,r;dt=.01,window=.5,monitor_buses=collect(1:39),savepath=joinpath(OUT,"nonlinear","bus$(bus)_$(Int(delta))_trajectory.csv"))
            (;design="corrected",bus,delta,F=met.Fpeak_Hz,R=met.Rpeak_Hz_s,Vmin=met.Vmin,Vmax=met.Vmax,slack=met.limiter_fraction,runtime=met.runtime_s,pass,error="")
        catch err
            (;design="corrected",bus,delta,F=NaN,R=NaN,Vmin=NaN,Vmax=NaN,slack=NaN,runtime=NaN,pass=false,error=sprint(showerror,err))
        end
        push!(rows,row);CSV.write(joinpath(OUT,"TABLE_12_NONLINEAR_EVENTS.csv"),DataFrame(rows))
        println("EVENT_DONE ",row);flush(stdout)
    end
end
main()
