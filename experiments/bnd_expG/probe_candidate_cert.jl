include(joinpath(pwd(),"src","bnd_design_g","BNDDesignG.jl"))
using .BNDDesignG, LinearAlgebra
net=BNDDesignG.CollectiveModel.frozen_network(pwd())
n=length(net.sg); nom=BNDDesignG.PhysicalData.nominal_pll_gains()
e=[0.0,0.4,1.0,1.0,1.0,1.0,1.0,0.0,1.0,1.0]
kp=fill(nom.Kp,n); ki=fill(nom.Ki,n)
s=BNDDesignG.spectrum_with_gauge(net,1 .- e,kp,ki)
r=BNDDesignG.robustness_certificate(s.quotient_matrix,1.6991206999182038e-6,0.05;scale=1.0)
println("alpha=",s.spectral_abscissa," critical=",s.critical," condition=",s.critical_nonnormality)
println("modal_bound=",r.modal_residue_bound," mrob=",r.m_rob," robust_alpha=",s.spectral_abscissa+r.m_rob)
println("beta_star=",r.beta_star," exact_small_gain_pass=",r.small_gain_pass," exact_margin=",r.small_gain_margin," omega=",r.omega_peak)
