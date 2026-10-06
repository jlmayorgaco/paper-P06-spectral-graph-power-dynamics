include(joinpath(pwd(),"src","bnd_design_g","BNDDesignG.jl"))
using .BNDDesignG, LinearAlgebra, CSV, DataFrames
root=pwd()
net=BNDDesignG.CollectiveModel.frozen_network(root)
buses=sort(collect(keys(net.sg)))
ng=length(buses)
nominal=BNDDesignG.PhysicalData.nominal_pll_gains()
kp=fill(nominal.Kp,ng); ki=fill(nominal.Ki,ng)
patterns=(
("allSG",zeros(ng),kp,ki),
("allGFL_nominal",ones(ng),kp,ki),
("allGFL_pattern_A",ones(ng),kp.*[0.25,0.5,1,2,4,0.75,1.5,3,0.4,2.5],ki.*[4,3,2,1,0.25,2.5,0.4,1.5,0.75,0.5]),
("allGFL_pattern_B",ones(ng),kp.*[1.7,0.3,3.2,0.6,2.8,0.45,1.3,3.8,0.9,2.1],ki.*[0.4,2.3,0.7,3.5,1.1,4,0.3,1.8,2.7,0.55]))
for (name,rho,kp0,ki0) in patterns
 s=BNDDesignG.zero_structure(net,rho,kp0,ki0)
 println(name," n=",size(s.model.Ared,1)," alpha=",s.spectral_abscissa," gauge=",s.gauge_pole,
 " gauge_res=",s.gauge_residual," cls=",s.classification," near=",s.physical_nearzero,
 " nullity=",s.geometric_multiplicity," svmin=",last(s.singular_values))
end
d=BNDDesignG.PhysicalData.discover_generators(root)
println("GEN_BUSES=",buses," total_actual_MW=",sum(d.SG_dispatch_initial_MW)," KCL=",net.no_load_kcl_inf)
