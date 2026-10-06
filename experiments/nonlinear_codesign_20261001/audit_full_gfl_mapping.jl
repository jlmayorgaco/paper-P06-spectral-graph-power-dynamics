include("ReducedDAE.jl")
using .ReducedDAE, Random, SHA, TOML, LinearAlgebra
const R=ReducedDAE
const OUT=joinpath(R.OUT,"full_gfl_sector_contract")
ctx=R.N.design_context(R.ROOT)
rng=MersenneTwister(20261001)
maxrhs=0.;maxenergy=0.;cases=0
for bus in 30:39
    op=R.N.trim_gfl(ctx,bus);p=op.pars;V0=norm(op.u)
    for sample in 1:100
        reference=5randn(rng);delta=.1randn(rng);omega=randn(rng);nu=randn(rng)
        theta=reference+delta;x=copy(op.x)
        x[1:2].+=.001randn(rng,2);x[3]=theta;x[4]=omega;x[5]=randn(rng)
        id=op.x[8]+.2randn(rng);iq=p.iset_q+.2randn(rng)
        x[7]=cos(theta)*id-sin(theta)*iq;x[6]=sin(theta)*id+cos(theta)*iq
        x[8]+=.1randn(rng);x[9]+=.1randn(rng)
        du=.1randn(rng,2)
        # External voltage perturbation expressed in the PLL frame.
        upll=[V0*cos(delta),-V0*sin(delta)]+du
        u=[cos(theta) -sin(theta);sin(theta) cos(theta)]*upll
        pdc=.1randn(rng);pp=merge(p,(;Pdc=p.Pdc+pdc))
        kp=10+70rand(rng);ki=60+500rand(rng)
        dx=R.gfl_rhs(x,u,pp,kp,ki,:physical_supply)
        ed=p.dc_kp*(x[9]-p.Vdc)+x[8]-id;eq=p.iset_q-iq
        vid=p.cc_kp*ed+p.cc_ki*x[2];viq=p.cc_kp*eq+p.cc_ki*x[1]
        a=p.omega_base/p.Xf;freq=p.omega_base*p.omega_frame+omega
        iddot=a*(vid-upll[1]-p.Rf*id)+freq*iq
        iqdot=a*(viq-upll[2]-p.Rf*iq)-freq*id
        predicted=[omega-nu,(x[5]+kp*upll[2]-omega)/p.pll_tau,ki*upll[2],
                   eq,ed,iddot,iqdot,p.dc_ki*(x[9]-p.Vdc),
                   (pp.Pdc-vid*id-viq*iq)/(p.Cdc*x[9])]
        observed=[dx[3]-nu,dx[4],dx[5],dx[1],dx[2],
                  cos(theta)*dx[7]+sin(theta)*dx[6]+omega*iq,
                  -sin(theta)*dx[7]+cos(theta)*dx[6]-omega*id,dx[8],dx[9]]
        global maxrhs=max(maxrhs,norm(predicted-observed,Inf))
        edot=p.Cdc*x[9]*dx[9]+p.Xf/p.omega_base*(id*iddot+iq*iqdot)
        supply=pp.Pdc-dot(upll,[id,iq])-p.Rf*(id^2+iq^2)
        global maxenergy=max(maxenergy,abs(edot-supply))
        global cases+=1
    end
end
@assert maxrhs<1e-8
@assert maxenergy<1e-10
out=Dict("status"=>"FULL_INSTALLED_GFL_COORDINATES_AND_ENERGY_IDENTITY_VERIFIED",
         "cases"=>cases,"max_absolute_rhs_error"=>maxrhs,"max_energy_identity_error"=>maxenergy,
         "ReducedDAE_sha256"=>bytes2hex(sha256(read(joinpath(@__DIR__,"ReducedDAE.jl")))),
         "script_sha256"=>bytes2hex(sha256(read(@__FILE__))))
open(joinpath(OUT,"julia_mapping_audit.toml"),"w") do io;TOML.print(io,out);end
TOML.print(stdout,out)
