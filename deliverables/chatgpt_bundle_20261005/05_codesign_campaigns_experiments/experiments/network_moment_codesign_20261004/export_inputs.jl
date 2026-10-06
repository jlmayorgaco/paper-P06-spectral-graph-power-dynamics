using LinearAlgebra, CSV, DataFrames, TOML, ForwardDiff
include(joinpath(@__DIR__,"..","graph_gsp_codesign_20261003","DelayedEvents.jl"))
const R=DelayedEvents.R
BLAS.set_num_threads(1)
const OUT=@__DIR__
function main()
    ctx=R.N.design_context(R.ROOT)
    d=TOML.parsefile(joinpath(OUT,"..","interaction_decision_20261004","designs","corrected.toml"))
    rho,kp,ki=Float64.(d["rho"]),Float64.(d["Kp"]),Float64.(d["Ki"])
    m=R.model(ctx,rho,kp,ki;dc_convention=:physical_supply)
    x=m.x0;de=R.derivatives(x,m);B=DelayedEvents.injection(m)
    e,C,_=DelayedEvents.detector(x,m,de)
    A0=de.Fx-B*C
    nx=length(x);Fu=zeros(nx,20)
    for i=1:10
        vi=(2i-1):(2i);u=de.v[vi]
        si=m.sgidx[i];sp=m.sp[i]
        Fu[si,vi]=ForwardDiff.jacobian(v->sp.controlled ? R.SG.rhs(x[si],v,sp) : R.SG.rhs_uncontrolled(x[si],v,sp),u)
        gi=m.gfidx[i];gp=m.gp[i]
        Fu[gi,vi]=ForwardDiff.jacobian(v->R.gfl_rhs(x[gi],v,gp,kp[i],ki[i],:physical_supply),u)
    end
    lift=vcat(m.net.lift,Matrix{Float64}(I,20,20));vv=lift*de.v
    phase=zeros(39,78)
    for i=1:39
        ur,ui=vv[2i-1:2i];phase[i,2i-1]=-ui/(ur^2+ui^2);phase[i,2i]=ur/(ur^2+ui^2)
    end
    O=phase*lift*de.vx
    G=copy(m.net.Y)
    for i=1:10
        D,_=R.norton(x[m.sgidx[i]],m.sp[i]);ix=(2*(29+i)-1):(2*(29+i))
        G[ix,ix]+=(1-rho[i])*D
    end
    inputs=zeros(nx,3);feed=zeros(39,3);edirect=zeros(10,3);checks=NamedTuple[]
    for (j,bus) in enumerate([8,16,29])
        net=R.network(ctx,bus,0.);Yp=zeros(78,78)
        Yp[2bus-1,2bus-1]=Yp[2bus,2bus]=-1/(100net.vset2)
        dv=-(G\(Yp*vv));inputs[:,j]=Fu*dv[59:78];feed[:,j]=phase*dv
        for i=1:10
            theta=x[m.gfidx[i][3]];edirect[i,j]=-sin(theta)*dv[59+2(i-1)]+cos(theta)*dv[60+2(i-1)]
        end
        h=1e-3;mp=R.model(ctx,rho,kp,ki;bus,delta=h,dc_convention=:physical_supply)
        mm=R.model(ctx,rho,kp,ki;bus,delta=-h,dc_convention=:physical_supply)
        fp=zeros(nx);fm=zeros(nx);R.rhs!(fp,x,mp,0.);R.rhs!(fm,x,mm,0.)
        fd=(fp-fm)/(2h)
        push!(checks,(;bus,input_relative_error=norm(fd-inputs[:,j])/norm(inputs[:,j])))
    end
    q,_=R.rotation_generator(x,m;jacobian=false)
    mkpath(joinpath(OUT,"model"))
    for (name,mat) in [("A0",A0),("C",C),("Binput",inputs-B*edirect),("detector_input",edirect),
                       ("phase_state",O),("phase_input",feed),("gauge",reshape(q,:,1))]
        CSV.write(joinpath(OUT,"model",name*".csv"),DataFrame(mat,:auto))
    end
    CSV.write(joinpath(OUT,"TABLE_01_INPUT_PARITY.csv"),DataFrame(checks))
    println(checks)
end
main()
