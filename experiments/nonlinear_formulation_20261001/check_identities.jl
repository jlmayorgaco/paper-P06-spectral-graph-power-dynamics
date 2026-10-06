using LinearAlgebra, Random, Test, TOML, SHA, ForwardDiff
BLAS.set_num_threads(1)
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const DEST=joinpath(ROOT,"reports","nonlinear_formulation_20261001")
mkpath(DEST)
const SRC=joinpath(ROOT,"src","bnd_model_expN","PDExactDesignN.jl")
source=read(SRC,String)
a=first(findfirst("function gfl_rhs(",source))
b=first(findnext("function gfl_jacobians(",source,a))-1
include_string(Main,source[a:b],SRC*":gfl_rhs")

# Exact coordinate lift on c^2+s^2=1. The DC division is retained here only
# for comparing vector fields; the polynomial residual multiplies by Cdc*vdc.
lift(x)=[x[1],x[2],cos(x[3]),sin(x[3]),x[4],x[5],x[6],x[7],x[8],x[9]]
function lifted_rhs(y,u,p,kp,ki)
  gq,gd,c,s,w,wi,ifi,ifr,vdi,vdc=y
  ur,ui=u;id=c*ifr+s*ifi;iq=-s*ifr+c*ifi
  ed=(p.Vdc-vdc)*p.dc_kp+vdi-id;eq=p.iset_q-iq
  vd=p.cc_kp*ed+p.cc_ki*gd;vq=p.cc_kp*eq+p.cc_ki*gq
  vr=c*vd-s*vq;vi=s*vd+c*vq;eqpll=-s*ur+c*ui
  [eq,ed,-s*w,c*w,(wi+kp*eqpll-w)/p.pll_tau,ki*eqpll,
   (p.omega_base/p.Xf)*(vi-ui-p.Rf*ifi-p.omega_frame*p.Xf*ifr),
   (p.omega_base/p.Xf)*(vr-ur-p.Rf*ifr+p.omega_frame*p.Xf*ifi),
   (p.Vdc-vdc)*p.dc_ki,(vd*id+vq*iq-p.Pdc)/(p.Cdc*vdc)]
end

function main()
  rng=MersenneTwister(20261001);report=Dict{String,Any}()
  p=(;Rf=.01,Xf=.03,omega_base=2pi*60,omega_frame=1.,pll_tau=1/(300*2pi),
    cc_kp=.3,cc_ki=120.,Cdc=1.25,Vdc=2.5,dc_kp=90.,dc_ki=750.,iset_q=-.2,Pdc=1.)
  maxlift=0.;maxcircle=0.;maxdc=0.
  @testset "Exact GFL lift and nonlinear DC energy identity" begin
    for j in 1:200
      x=[.01randn(rng),.01randn(rng),2pi*rand(rng)-pi,10randn(rng),5randn(rng),
        randn(rng),randn(rng),1+.1randn(rng),2.5+.1randn(rng)]
      ph=2pi*rand(rng);u=(.9+.2rand(rng)).*[cos(ph),sin(ph)]
      kp=8+117rand(rng);ki=60+920rand(rng)
      y=lift(x);f=gfl_rhs(x,u,p,kp,ki);dy=lifted_rhs(y,u,p,kp,ki)
      pushed=ForwardDiff.jacobian(lift,x)*f
      err=norm(dy-pushed,Inf)/max(1.,norm(pushed,Inf));maxlift=max(maxlift,err)
      circ=abs(2y[3]*dy[3]+2y[4]*dy[4]);maxcircle=max(maxcircle,circ)
      c,s=y[3:4];id=c*x[7]+s*x[6];iq=-s*x[7]+c*x[6]
      ed=(p.Vdc-x[9])*p.dc_kp+x[8]-id;eq=p.iset_q-iq
      pac=(p.cc_kp*ed+p.cc_ki*x[2])*id+(p.cc_kp*eq+p.cc_ki*x[1])*iq
      dc=abs(p.Cdc*x[9]*f[9]-(pac-p.Pdc));maxdc=max(maxdc,dc)
      @test err<2e-13
      @test circ<1e-12
      @test dc<1e-10
    end
  end
  report["gfl_lift"]=Dict("samples"=>200,"max_relative_vector_field_error"=>maxlift,
    "max_circle_derivative"=>maxcircle,"max_dc_energy_residual"=>maxdc,
    "source_sha256"=>bytes2hex(sha256(read(SRC))))

  # Illustrative lossless graph. This is not the IEEE-39 or the full SG model.
  B=zeros(4,5);edges=[(1,2),(2,3),(3,4),(4,1),(1,3)]
  for (e,(i,j)) in enumerate(edges);B[i,e]=1;B[j,e]=-1;end
  w=[1.,1.3,.8,1.1,.7];W=Diagonal(w);L=B*W*B'
  spec=eigen(Symmetric(L));U=spec.vectors;lam=spec.values
  delta0=[.2,-.1,.05,-.15];eta0=B'*delta0;BU=B'*U
  network(d)=B*(w.*sin.(B'*d))
  maxbasis=0.;maxratio=0.
  @testset "Exact graph coordinates and cubic Taylor remainder" begin
    for j in 1:200
      dd=(.01+.5rand(rng)).*randn(rng,4);h=B'*dd
      z=U'*(delta0+dd)
      basiserr=norm(U'*network(U*z)-U'*network(delta0+dd))
      maxbasis=max(maxbasis,basiserr)
      cubic=sin.(eta0)+cos.(eta0).*h-.5sin.(eta0).*h.^2-(cos.(eta0).*h.^3)./6
      residual=norm(network(delta0+dd)-B*(w.*cubic))
      bound=opnorm(B*W)*norm(abs.(h).^4)/24
      ratio=residual/max(bound,eps());maxratio=max(maxratio,ratio)
      @test basiserr<1e-12
      @test residual<=bound+1e-13
    end
    # Nonlinear cross-mode coefficient, even though U diagonalizes the graph L.
    q=zeros(4,4,4)
    for i in 1:4,j in 1:4,k in 1:4
      q[i,j,k]=-sum(w.*sin.(eta0).*BU[:,i].*BU[:,j].*BU[:,k])
    end
    cross=maximum(abs(q[i,j,k]) for i in 2:4,j in 2:4,k in 2:4 if !(i==j==k))
    @test cross>1e-3
    D=Diagonal([.1,.6,.7,1.3]);DH=U'*D*U
    offratio=norm(DH-Diagonal(diag(DH)))/norm(DH)
    @test offratio>.1
    report["graph"]=Dict("model"=>"illustrative four-node lossless graph, not IEEE-39",
      "max_coordinate_error"=>maxbasis,"max_error_over_Taylor_bound"=>maxratio,
      "max_cross_mode_quadratic_coefficient"=>cross,"heterogeneous_gain_offdiagonal_fraction"=>offratio)
  end

  @testset "Polynomial graph locality, common mode, frequency units" begin
    signal=randn(rng,4);coef=[.8,.2,.03]
    centralized=(coef[1]*I+coef[2]*L+coef[3]*L^2)*signal
    distributed=coef[1]*signal+coef[2]*(L*signal)+coef[3]*(L*(L*signal))
    spectral=U*((coef[1].+coef[2].*lam.+coef[3].*lam.^2).*(U'*signal))
    @test norm(centralized-distributed)<1e-12
    @test norm(centralized-spectral)<1e-12
    @test norm(L*ones(4))<1e-12
    @test norm((.8I+.2L)*ones(4)-.8ones(4))<1e-12
    T=.5;t=4.;freq=.35;phase(t)=2pi*freq*t
    @test abs((phase(t)-phase(t-T))/(2pi*T)-freq)<1e-14
    @test abs((phase(t)-2phase(t-T)+phase(t-2T))/(2pi*T^2))<1e-14
  end
  report["scope"]="Algebraic identities and a toy graph verified. No nonlinear IEEE-39 controller, ROA, SOS certificate, or optimum computed."
  report["status"]="PASS"
  open(joinpath(DEST,"identity_audit.toml"),"w") do io;TOML.print(io,report);end
  println("IDENTITY_AUDIT ",report)
end
main()
