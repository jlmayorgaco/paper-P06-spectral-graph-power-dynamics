include("oracle.jl")
d=TOML.parsefile(joinpath(OUT,"candidate_improved.toml"))
ep=Float64.(d["epsilon"]);kp=Float64.(d["Kp"]);ki=Float64.(d["Ki"])
rows=NamedTuple[];options=Any[]
for factor in (1.0,1.03,1.06,1.10,1.15), floor in (.35,.5,.7,1.0)
  # A tuning heuristic in model units, not a grid-coupled damping certificate.
  kk=max.(kp,2floor.*sqrt.(ki));ee=min.(ep.*factor,1-1e-5)
  v=oracle(CoreDesign.encode(CTX,ee,kk,ki);dt=.01,horizon=90.)
  feasible=maximum(v.g)<=0
  row=(;factor,pll_floor_heuristic=floor,J=v.J,alpha=v.alpha,beta=v.beta,
    F16=v.peaks[1],F8=v.peaks[2],F29=v.peaks[3],Rmax=maximum(v.rocofs),feasible)
  push!(rows,row);println("REPAIR_SCREEN ",row);flush(stdout)
  feasible && floor>=.7 && push!(options,(;v,factor,floor))
end
CSV.write(joinpath(OUT,"repair_damping_screen.csv"),DataFrame(rows))
isempty(options) && error("No repaired linear-feasible candidate")
sort!(options;by=x->(x.v.J,-x.floor))
chosen=first(options);v=chosen.v
save_candidate("candidate_repaired.toml",v.ep,v.kp,v.ki;extra=Dict{String,Any}(
  "alpha"=>v.alpha,"beta_sampled"=>v.beta,"Fpeak_Hz"=>v.peaks,"Rpeak_Hz_s"=>v.rocofs,
  "retention_factor_from_improved"=>chosen.factor,"pll_floor_heuristic"=>chosen.floor,
  "claim"=>"retuned after nonlinear bus8 failure; independent validation pending",
  "local_KKT_certified"=>false,"global_optimum_certified"=>false))
CSV.write(joinpath(OUT,"parameters_repaired.csv"),DataFrame(bus=30:39,rho=1 .-v.ep,
  epsilon=v.ep,retained_SG_MW=CTX.power.*v.ep,Kp=v.kp,Ki=v.ki))
println("REPAIRED_SELECTED ",chosen.factor," ",chosen.floor," ",v.J);flush(stdout)
