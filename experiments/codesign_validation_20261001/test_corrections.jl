using LinearAlgebra,Test
module LinearSecurity
using LinearAlgebra
reduced_step_model(ctx,m,rho,bus;kwargs...)=(;A=reshape([-1.0],1,1),B=[0.0],C_bus=zeros(10,1),D_bus=zeros(10))
function load_input_vector(ctx,m,bus;kwargs...)
    y=zeros(78);y[60]=-ctx.jump;y
end
end
include("FiniteWindow.jl")
@testset "Physical phase units and causal window" begin
  for jump in (-.02,.02)
    ctx=(net=(voltage=ones(ComplexF64,39),),jump=jump)
    m=(Gy=Matrix{Float64}(I,78,78),)
    t=collect(-100:400).*.01
    sig=FiniteWindow.phase_step(ctx,m,ones(10),t;disturbance_MW=1.,gauge_vector=x->zeros(1))
    metric=only(FiniteWindow.window_metrics(sig;windows=(.5,),dt_s=.01,horizon_s=4.))
    @test sig.phase_jump_rad[1]≈jump
    @test metric.F_peak_Hz≈abs(jump)/(2pi*.5)
    @test metric.R_peak_Hz_s≈abs(jump)/(2pi*.5^2)
  end
end
module JacobianHarness
const FREQUENCY_LIMIT=.5
const ROCOF_LIMIT=.5
decode(ctx,s,y)=(y[1:length(s)],zeros(10),zeros(10))
_bounds(ctx,s)=(zeros(length(s)+20),ones(length(s)+20))
_dc_frequency(ctx,s,y,d)=3y[1]
end
source=read(joinpath(@__DIR__,"CoreDesign.jl"),String)
a=first(findfirst("function constraint_jacobian(",source))
b=first(findnext("function _solve_qp(",source,a))-1
include_string(JacobianHarness,source[a:b])
@testset "Steady plateau envelope derivative" begin
  ev=(g=[-10.,-10.,0.,0.],rho=ones(10),Fpeak=.6,Finf=.6)
  jac=JacobianHarness.constraint_jacobian(nothing,[30],fill(.2,21),ev)
  @test jac.J[3,1]≈6.
  @test jac.J[4,1]≈6.
end
