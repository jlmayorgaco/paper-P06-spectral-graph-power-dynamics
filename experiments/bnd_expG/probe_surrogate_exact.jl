include(joinpath(pwd(),"src","bnd_design_g","BNDDesignG.jl"))
using .BNDDesignG, LinearAlgebra
root=pwd(); net=BNDDesignG.CollectiveModel.frozen_network(root)
buses=sort(collect(keys(net.sg))); n=length(buses); nom=BNDDesignG.PhysicalData.nominal_pll_gains()
kpf=fill(nom.Kp,n); kif=fill(nom.Ki,n)
eps=zeros(n); eps[findfirst(==(38),buses)]=0.001; eps[findfirst(==(39),buses)]=0.12
s=BNDDesignG.spectrum_with_gauge(net,1 .- eps,kpf,kif)
println("alpha=",s.spectral_abscissa," critical=",s.critical," support=",[(b,eps[i]) for (i,b) in enumerate(buses) if eps[i]>0])
r=BNDDesignG.robustness_certificate(s.quotient_matrix,3.4e-6,0.05;scale=1.0)
println("mrob=",r.m_rob," beta*",r.beta_star," SGpass=",r.small_gain_pass," robustalpha=",s.spectral_abscissa+r.m_rob)

