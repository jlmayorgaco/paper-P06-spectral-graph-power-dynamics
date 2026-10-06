include(joinpath(pwd(),"src","bnd_design_g","BNDDesignG.jl"))
using .BNDDesignG, LinearAlgebra
root=pwd(); net=BNDDesignG.CollectiveModel.frozen_network(root)
buses=sort(collect(keys(net.sg))); n=length(buses); nom=BNDDesignG.PhysicalData.nominal_pll_gains()
kp=fill(nom.Kp,n); ki=fill(nom.Ki,n); P=BNDDesignG.PhysicalData.discover_generators(root).SG_dispatch_initial_MW
h=[net.sg[b].op.parameters.inertia*net.sg[b].op.parameters.rating_mva for b in buses]
# Reference robust radius is frozen from the all-SG analytical model.
base=BNDDesignG.spectrum_with_gauge(net,zeros(n),kp,ki)
rb=BNDDesignG.robustness_certificate(base.quotient_matrix,0.0,0.05;scale=1.0)
spare=-base.spectral_abscissa-0.05
beta=0.5*min(rb.beta_star,spare/(2rb.modal_residue_bound))
println("reference beta=",beta," base modal R=",rb.modal_residue_bound)
bR=60*100/(2*0.5); eps=ones(n)
function metrics(e)
 s=BNDDesignG.spectrum_with_gauge(net,1 .- e,kp,ki)
 mr=BNDDesignG.modal_residue_bound(s.quotient_matrix).value*beta
 (;s, mr, feasible=s.spectral_abscissa+mr<=-0.05+1e-8 && dot(h,e)>=bR-1e-8,
   retained=dot(P,e))
end
println("allSG=",metrics(eps).s.spectral_abscissa," mrob=",metrics(eps).mr)
for sweep in 1:2
 changed=false
 for i in 1:n
   best=eps[i]; bestcost=metrics(eps).retained
   for x in 0.0:0.1:1.0
     trial=copy(eps); trial[i]=x; mt=metrics(trial)
     if mt.feasible && dot(P,trial)<bestcost-1e-7
       best=x; bestcost=dot(P,trial)
     end
   end
   if best < eps[i]-1e-8
     eps[i]=best; changed=true
     println("sweep=",sweep," bus=",buses[i]," eps=",best," retained=",bestcost," alpha=",metrics(eps).s.spectral_abscissa," mrob=",metrics(eps).mr)
   end
 end
 changed || break
end
println("FINAL_EPS=",eps," retained=",dot(P,eps)," metrics=",metrics(eps).s.spectral_abscissa," mrob=",metrics(eps).mr)
