include(joinpath(pwd(),"src","bnd_design_g","BNDDesignG.jl"))
using .BNDDesignG, LinearAlgebra, CSV, DataFrames
net=BNDDesignG.CollectiveModel.frozen_network(pwd())
buses=sort(collect(keys(net.sg))); n=length(buses)
nom=BNDDesignG.PhysicalData.nominal_pll_gains()
kp=fill(nom.Kp,n); ki=fill(nom.Ki,n); rho=ones(n)
s0=BNDDesignG.zero_structure(net,rho,kp,ki)
println("ZERO_CLASS=",s0.classification," jordan_chain_residual=",s0.generalized_vector_residual)
println("SINGULAR_VALUES_NEAR_ZERO=",s0.singular_values[end-8:end])
println("QUOTIENT_NEAR_POLES=",sort(s0.physical,by=abs)[1:8])
rows=NamedTuple[]
for e in [1e-8,3e-8,1e-7,3e-7,1e-6,3e-6,1e-5,3e-5,1e-4,3e-4,1e-3]
 r=copy(rho); r[1]=1-e
 sp=BNDDesignG.spectrum_with_gauge(net,r,kp,ki)
 nearest=sort(sp.physical,by=abs)[1:min(4,length(sp.physical))]
 push!(rows,(epsilon=e,lambda1_real=real(nearest[1]),lambda1_imag=imag(nearest[1]),
  lambda2_real=real(nearest[2]),lambda2_imag=imag(nearest[2]),alpha=sp.spectral_abscissa))
 println("eps=",e," near=",nearest[1:min(2,length(nearest))]," alpha=",sp.spectral_abscissa)
end
CSV.write("reports/experiment_G/tables/TABLE_G01_jordan_scaling_probe.csv",DataFrame(rows))
