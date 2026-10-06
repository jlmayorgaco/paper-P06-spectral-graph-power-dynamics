include(joinpath(pwd(),"src","bnd_design_g","BNDDesignG.jl"))
using .BNDDesignG, LinearAlgebra
net=BNDDesignG.CollectiveModel.frozen_network(pwd())
n=length(net.sg); g=BNDDesignG.PhysicalData.nominal_pll_gains()
sp=BNDDesignG.spectrum_with_gauge(net,zeros(n),fill(g.Kp,n),fill(g.Ki,n)); m=sp.model
rp=BNDDesignG.robustness_certificate(sp.quotient_matrix,0.0,0.05;scale=1.0)
println("ALPHA=",maximum(real.(eigvals(sp.quotient_matrix))))
println("BETA_STAR=",rp.beta_star," RESOLVENT_PEAK=",rp.resolvent_peak,
 " MODAL_R=",rp.modal_residue_bound," OMEGA=",rp.omega_peak)
println("H=",[(b,net.sg[b].op.parameters.inertia*net.sg[b].op.parameters.rating_mva) for b in sort(collect(keys(net.sg)))])

