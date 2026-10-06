include(joinpath(pwd(),"src","bnd_design_g","BNDDesignG.jl"))
using .BNDDesignG, LinearAlgebra
net=BNDDesignG.CollectiveModel.frozen_network(pwd())
buses=sort(collect(keys(net.sg))); n=length(buses)
nom=BNDDesignG.PhysicalData.nominal_pll_gains()
kp=fill(nom.Kp,n); ki=fill(nom.Ki,n)
println("bus,eps,dim,alpha,positive_poles,nearest_to_zero")
for epsv in (1e-4,0.01,0.1,0.5,1.0)
 for i in 1:n
  rho=ones(n); rho[i]=1-epsv
  A=BNDDesignG.CollectiveModel.mixed_jacobian(net,rho,kp,ki).Ared
  λ=eigvals(A); phys=[z for z in λ if abs(z)>1e-6]
  println(buses[i],",",epsv,",",size(A,1),",",maximum(real.(phys)),",",count(z->real(z)>1e-8,phys),",",phys[argmin(abs.(phys))])
 end
end
