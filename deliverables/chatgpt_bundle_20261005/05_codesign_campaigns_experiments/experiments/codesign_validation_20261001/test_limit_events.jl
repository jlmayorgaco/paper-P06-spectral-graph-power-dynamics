using SciMLBase,OrdinaryDiffEqRosenbrock,Test
module NetworkDynamics
module SII
variable_symbols(nw)=nw
end
end
VIndex(bus,name)="$(bus):$(name)"
source=read(joinpath(@__DIR__,"validate_pd.jl"),String)
a=first(findfirst("function limiter_callbacks(",source))
b=first(findnext("function main()",source,a))-1
include_string(Main,source[a:b])
@testset "Existing limiter event realization and release" begin
  nw=[VIndex(30,:ctrld_gen₊gov₊xg1),VIndex(30,:ctrld_gen₊avr₊vr)]
  state=Dict(VIndex(30,:ctrld_gen₊gov₊V_min)=>0.,VIndex(30,:ctrld_gen₊gov₊V_max)=>1.,
    VIndex(30,:ctrld_gen₊avr₊vr_min)=>-1.,VIndex(30,:ctrld_gen₊avr₊vr_max)=>1.)
  rho=ones(10);rho[1]=.5
  for mode in (:upper,:lower)
    events=NamedTuple[];cbs=limiter_callbacks(nw,state,rho,events)
    function f!(du,u,p,t)
      target=mode==:upper ? (t<2 ? 2. : 0.) : (t<2 ? -1. : .5)
      du[1]=((u[1]>1 && target>u[1]) || (u[1]<0 && target<u[1])) ? 0. : target-u[1]
      du[2]=-u[2]
    end
    initial=mode==:upper ? [0.,0.] : [.5,0.]
    sol=solve(ODEProblem(f!,initial,(0.,5.)),Rodas5P();callback=CallbackSet(cbs...),
      tstops=[2.],abstol=1e-10,reltol=1e-10,saveat=.01)
    @test SciMLBase.successful_retcode(sol.retcode)
    @test length(events)==1
    @test abs(events[1].time_s-(mode==:upper ? log(2) : log(1.5)))<1e-7
    expected=mode==:upper ? exp(-3.) : .5*(1-exp(-3.))
    @test abs(sol.u[end][1]-expected)<1e-7
    # Dense output near a nonsmooth crossing has integration/interpolation
    # error in addition to the 1e-10 side-selection nudge.
    @test maximum(x[1] for x in sol.u)<=1+2e-8
    @test minimum(x[1] for x in sol.u)>=-2e-8
  end
end
