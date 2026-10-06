using LinearAlgebra, Random, Test, TOML, SHA
BLAS.set_num_threads(1)
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
include(joinpath(ROOT,"src","bnd_model_expN","PDExactDesignN.jl"))
const N=PDExactDesignN
const SG=N.CollectiveModel.AnalyticSG

"""Exact affine stator map at arbitrary machine states, in repository signs."""
function norton(x,p)
  xm=length(x)==12 ? x[7:12] : x
  pq,pd,ed,eq,w,delta=xm
  ad=(p.xdpp-p.xls)/(p.xdp-p.xls);aq=(p.xqpp-p.xls)/(p.xqp-p.xls)
  cd=ad*eq+(1-ad)*pd;cq=-aq*ed+(1-aq)*pq
  sn,cs=sincos(delta);T=[sn -cs;cs sn]
  R=[p.rs -w*p.xqpp;w*p.xdpp p.rs]
  scale=p.rating_mva/p.system_base_mva
  D=-scale*T'*(R\T);c=scale*T'*(R\(w.*[-cq,cd]))
  D,c
end

function main()
 ctx=N.design_context(ROOT);rng=MersenneTwister(73510)
 baseline=collect(Iterators.flatten(([real(v),imag(v)] for v in ctx.net.voltage)))
 sg=[copy(ctx.net.sg[b].op.x) for b in 30:39]
 gfl=[copy(N.trim_gfl(ctx,b).x) for b in 30:39]
 function assemble(xs,xg,rho)
   G=copy(ctx.net.y_static);h=zeros(78)
   for k in 1:10
     bus=k+29;ix=(2bus-1):(2bus);p=ctx.net.sg[bus].op.parameters
     D,c=norton(xs[k],p)
     G[ix,ix].+=(1-rho[k]).*D
     h[ix].+=(1-rho[k]).*c+rho[k].*xg[k][[7,6]]
   end
   G,h
 end
 # Positive symmetric graph weights from off-diagonal line-admittance blocks.
 W=zeros(39,39)
 for i in 1:39,j in i+1:39
   a=hypot(ctx.net.y_static[2i-1,2j-1],ctx.net.y_static[2i,2j-1])
   b=hypot(ctx.net.y_static[2j-1,2i-1],ctx.net.y_static[2j,2i-1])
   W[i,j]=W[j,i]=(a+b)/2
 end
 L=Diagonal(vec(sum(W,dims=2)))-W;U=eigen(Symmetric(L)).vectors
 T=kron(U,Matrix{Float64}(I,2,2))
 maxtrim=0.;maxaffine=0.;maxkcl=0.;maxsens=0.;maxgraph=0.;minsigma=Inf
 maxfinite=0.;maxdet=0.
 @testset "Exact nonlinear IEEE39 algebraic elimination" begin
   for rho in (fill(.4,10),TOML.parsefile(joinpath(ROOT,"reports","codesign_validation_20261001","candidate_improved.toml"))["rho"])
     G,h=assemble(sg,gfl,rho);v=-(G\h)
     err=norm(v-baseline,Inf);maxtrim=max(maxtrim,err)
     @test err<1e-9
   end
   for trial in 1:40
     xs=deepcopy(sg);xg=deepcopy(gfl);rho=.2 .+.7rand(rng,10)
     for k in 1:10
       xs[k].+=.002randn(rng,length(xs[k]));xg[k][6:7].+=.01randn(rng,2)
       bus=k+29;p=ctx.net.sg[bus].op.parameters;D,c=norton(xs[k],p)
       u=randn(rng,2);exact=SG.output(xs[k],u,p)
       err=norm(D*u+c-exact,Inf)/max(1.,norm(exact,Inf));maxaffine=max(maxaffine,err)
       @test err<1e-12
     end
     G,h=assemble(xs,xg,rho);v=-(G\h);r=ctx.net.y_static*v
     minsigma=min(minsigma,minimum(svdvals(G)))
     for k in 1:10
       bus=k+29;ix=(2bus-1):(2bus)
       r[ix].+=(1-rho[k]).*SG.output(xs[k],v[ix],ctx.net.sg[bus].op.parameters)+rho[k].*xg[k][[7,6]]
     end
     resid=norm(r,Inf);maxkcl=max(maxkcl,resid)
     @test resid<1e-10
     k=1+mod(trial-1,10);bus=k+29;ix=(2bus-1):(2bus)
     injection=zeros(78);injection[ix]=xg[k][[7,6]]-SG.output(xs[k],v[ix],ctx.net.sg[bus].op.parameters)
     vp=-(G\injection)
     hfd=1e-5;rp=copy(rho);rm=copy(rho);rp[k]+=hfd;rm[k]-=hfd
     Gp,hp=assemble(xs,xg,rp);Gm,hm=assemble(xs,xg,rm)
     fd=(-(Gp\hp)+(Gm\hm))/(2hfd)
     err=norm(vp-fd,Inf)/max(norm(vp,Inf),1e-9);maxsens=max(maxsens,err)
     @test err<1e-5
     # Exact finite rho update, at fixed internal device states.
     delta=.03;D,c=norton(xs[k],ctx.net.sg[bus].op.parameters)
     S=zeros(78,2);S[ix,:]=Matrix{Float64}(I,2,2)
     GS=G\S;Z=S'*GS;M2=I-delta*D*Z
     mismatch=xg[k][[7,6]]-SG.output(xs[k],v[ix],ctx.net.sg[bus].op.parameters)
     updated=v-delta*GS*(M2\mismatch)
     rf=copy(rho);rf[k]+=delta;Gf,hf=assemble(xs,xg,rf)
     actual=-(Gf\hf);err=norm(updated-actual,Inf);maxfinite=max(maxfinite,err)
     @test err<1e-10
     logG,signG=logabsdet(G);logF,signF=logabsdet(Gf);logM,signM=logabsdet(M2)
     derr=abs(logF-logG-logM);maxdet=max(maxdet,derr)
     @test derr<1e-10 && signF==signG*signM
     vh=-((T'*G*T)\(T'*h));err=norm(T*vh-v,Inf);maxgraph=max(maxgraph,err)
     @test err<1e-10
   end
 end
 d=Dict("status"=>"PASS","random_states"=>40,"stator_affinity_checks"=>400,
   "max_frozen_voltage_error_pu"=>maxtrim,"max_relative_stator_affinity_error"=>maxaffine,
   "max_nonlinear_KCL_residual_pu"=>maxkcl,"max_relative_rho_derivative_error"=>maxsens,
   "max_graph_coordinate_voltage_error_pu"=>maxgraph,"minimum_sampled_sigma_G"=>minsigma,
   "max_finite_rho_update_voltage_error_pu"=>maxfinite,"max_determinant_identity_log_error"=>maxdet,
   "source_sha256"=>Dict(f=>bytes2hex(sha256(read(joinpath(ROOT,f)))) for f in
      ["src/bnd_model_expN/PDExactDesignN.jl","src/bnd_design/AnalyticSG.jl","src/bnd_design_e/CollectiveModel.jl"]),
   "scope"=>"Exact stator algebra and constant-impedance network, random nearby states. Sampled sigma is not a region certificate. No controller or nonlinear optimum computed.")
 out=joinpath(ROOT,"reports","nonlinear_formulation_20261001");mkpath(out)
 open(joinpath(out,"nonlinear_network_audit.toml"),"w") do io;TOML.print(io,d);end
 println("NONLINEAR_NETWORK_AUDIT ",d)
end
main()
