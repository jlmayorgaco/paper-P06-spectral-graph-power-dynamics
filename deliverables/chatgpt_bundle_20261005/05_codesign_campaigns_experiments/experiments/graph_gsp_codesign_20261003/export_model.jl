using LinearAlgebra, ForwardDiff, CSV, DataFrames, TOML, SHA
include(joinpath(@__DIR__,"..","nonlinear_codesign_20261001","ReducedDAE.jl"))
const R=ReducedDAE
const ROOT=R.ROOT
BLAS.set_num_threads(1)
write_matrix(name,A)=CSV.write(joinpath(@__DIR__,"model",name*".csv"),DataFrame(A,:auto))
ctx=R.N.design_context(ROOT);rho=fill(.875,10);kp=fill(.9R.N.K0P,10);ki=fill(R.N.K0I,10)
m=R.model(ctx,rho,kp,ki;dc_convention=:physical_supply);n=length(m.x0);mkpath(joinpath(@__DIR__,"model"))
open(joinpath(@__DIR__,"baseline.toml"),"w") do io
    TOML.print(io,Dict("rho"=>rho,"Kp"=>kp,"Ki"=>ki,"tau"=>.04))
end
Adev=zeros(n,n);Bv=zeros(n,20);Cs=zeros(20,n);Cf=zeros(20,n);Ds=zeros(20,20)
Etheta=zeros(10,n);Hv=zeros(10,20);Bp=zeros(n,10);Bi=zeros(n,10)
v=R.voltage(m.x0,m)
for i in 1:10
    vi=(2i-1):(2i);u=v[vi];ix=m.sgidx[i];xx=m.x0[ix];p=m.sp[i];nn=length(ix)
    jf=ForwardDiff.jacobian(vcat(xx,u)) do z
        p.controlled ? R.SG.rhs(view(z,1:nn),view(z,nn+1:nn+2),p) : R.SG.rhs_uncontrolled(view(z,1:nn),view(z,nn+1:nn+2),p)
    end
    Adev[ix,ix].=jf[:,1:nn];Bv[ix,vi].=jf[:,nn+1:nn+2]
    Cs[vi,ix].=ForwardDiff.jacobian(x->R.SG.output(x,u,p),xx)
    Ds[vi,vi].=first(R.norton(xx,p))
    ix=m.gfidx[i];xx=m.x0[ix];pg=m.gp[i]
    jg=ForwardDiff.jacobian(z->R.gfl_rhs(view(z,1:9),view(z,10:11),pg,0.,0.,:physical_supply),vcat(xx,u))
    Adev[ix,ix].=jg[:,1:9];Bv[ix,vi].=jg[:,10:11]
    Cf[2i-1,ix[7]]=1;Cf[2i,ix[6]]=1
    th=xx[3];Hv[i,vi].=[-sin(th),cos(th)];Etheta[i,ix[3]]=-cos(th)*u[1]-sin(th)*u[2]
    Bp[ix[4],i]=1/pg.pll_tau;Bi[ix[5],i]=1
end
for (name,A) in (("Adev",Adev),("Bv",Bv),("Cs",Cs),("Cf",Cf),("Ds",Ds),("Y",m.net.reduced),
    ("Etheta",Etheta),("Hv",Hv),("Bp",Bp),("Bi",Bi),("x0",reshape(m.x0,:,1)))
    write_matrix(name,A)
end
ports=DataFrame(bus=30:39,state_index=[ix[end-1]-1 for ix in m.sgidx],M=[2*m.sp[i].inertia*(1-rho[i])*m.sp[i].rating_mva for i=1:10],
    rho=rho,Kp=kp,Ki=ki,P0=ctx.power,pll_angle_index=[ix[3]-1 for ix in m.gfidx],pll_frequency_index=[ix[4]-1 for ix in m.gfidx])
CSV.write(joinpath(@__DIR__,"model","ports.csv"),ports)
g,_=R.rotation_generator(m.x0,m;jacobian=false);write_matrix("gauge",reshape(g,:,1))
for (j,rr) in enumerate((rho,[.84+.005i for i in 1:10]))
    mm=R.model(ctx,rr,kp,ki;dc_convention=:physical_supply)
    write_matrix("A_check$j",R.derivatives(mm.x0,mm).Fx)
    write_matrix("rho_check$j",reshape(rr,:,1))
end
cp(joinpath(ROOT,"experiments","delay_dressed_replacement_frontier_20261002","GRAPH_LP_LOSSLESS_PORT.csv"),joinpath(@__DIR__,"model","L.csv");force=true)
println("PARAMETRIC_MODEL_EXPORTED states=",n);flush(stdout)
