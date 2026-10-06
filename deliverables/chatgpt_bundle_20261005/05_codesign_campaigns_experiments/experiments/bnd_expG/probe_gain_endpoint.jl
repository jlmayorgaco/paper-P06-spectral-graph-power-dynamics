include(joinpath(pwd(),"src","bnd_design_g","BNDDesignG.jl"))
using .BNDDesignG, LinearAlgebra
root=pwd(); net=BNDDesignG.CollectiveModel.frozen_network(root)
buses=sort(collect(keys(net.sg))); n=length(buses)
nom=BNDDesignG.PhysicalData.nominal_pll_gains()
for (label,kpf,kif) in (("nominal",1.0,1.0),("max",4.0,4.0),("min",0.25,0.25))
 kp=fill(kpf*nom.Kp,n); ki=fill(kif*nom.Ki,n)
 s=BNDDesignG.zero_structure(net,ones(n),kp,ki)
 println(label," class=",s.classification," physical_alpha=",s.spectral_abscissa,
  " near=",s.physical_nearzero," gauge_res=",s.gauge_residual,
  " max_gfl_gain=",maximum(real.(s.physical)))
end
