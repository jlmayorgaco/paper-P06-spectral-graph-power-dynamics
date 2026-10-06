include("oracle.jl")
function full_band(A,beta_req;maxnodes=30000)
    As=A+0.05I;n=size(A,1);anorm=opnorm(As)
    # A conservative working precision allowance, not directed interval arithmetic.
    allowance=64eps(Float64)*n*max(anorm,1.0)
    f(w)=minimum(svdvals(im*w*I-As))
    wmax=anorm+beta_req+1.0
    grid=sort!(unique(vcat(0.0,10.0.^range(-5,log10(wmax),length=100),
      abs.(imag.(eigvals(A))))))
    values=f.(grid);upper=minimum(values);omega=grid[argmin(values)]
    stack=[(grid[j],grid[j+1],values[j],values[j+1],0) for j in 1:length(grid)-1]
    nodes=length(grid);accepted=Float64[];unresolved=0;maxdepth=0
    while !isempty(stack)
      a,b,fa,fb,depth=pop!(stack)
      lb=max(0.0,(fa+fb-(b-a))/2-allowance)
      if lb>=beta_req
        push!(accepted,lb);continue
      end
      if min(fa,fb)<beta_req-allowance
        return Dict("status"=>"REJECTED_POINTWISE","beta_observed_upper"=>upper,"omega"=>omega,
          "nodes"=>nodes,"formal_certificate"=>false)
      end
      if nodes>=maxnodes || depth>=60
        unresolved+=1;push!(accepted,lb);continue
      end
      m=(a+b)/2;fm=f(m);nodes+=1;maxdepth=max(maxdepth,depth+1)
      if fm<upper;upper=fm;omega=m;end
      push!(stack,(a,m,fa,fm,depth+1),(m,b,fm,fb,depth+1))
    end
    Dict("status"=>unresolved==0 ? "FULL_BAND_NUMERICAL_LOWER_BOUND" : "INCOMPLETE_INTERVAL_BOUND",
      "beta_lower_bound"=>minimum(accepted),"beta_observed_upper"=>upper,"omega"=>omega,
      "beta_req"=>beta_req,"unresolved_intervals"=>unresolved,"nodes"=>nodes,"max_depth"=>maxdepth,
      "roundoff_allowance"=>allowance,"A_shift_norm"=>anorm,"tail_lower"=>wmax-anorm,
      "outward_rounded"=>false,"formal_certificate"=>false)
end

function main()
  seed=abspath(ARGS[1]);label=length(ARGS)>1 ? ARGS[2] : "improved"
  vals=parse.(Float64,split(strip(read(seed,String)),','))
  # Preserve a feasible operating margin after optimization, separately labelled.
  factor=length(ARGS)>2 ? parse(Float64,ARGS[3]) : 1.01
  vals[1:10]=min.(vals[1:10].*factor,1-1e-5)
  v=oracle(vals;dt=.005,horizon=120.)
  extra=Dict{String,Any}("design_label"=>label,"post_search_retention_factor"=>factor,
    "alpha"=>v.alpha,"beta_sampled"=>v.beta,"Fpeak_Hz"=>v.peaks,"Rpeak_Hz_s"=>v.rocofs,
    "Fsteady_Hz"=>v.steady,"peak_times_s"=>v.pt,"peak_buses"=>v.pb,
    "frequency_window_s"=>.5,"design_event_buses"=>[16,8,29],"design_Pset_step_MW"=>100.,
    "local_KKT_certified"=>false,"global_optimum_certified"=>false,
    "claim"=>"candidate with explicit post-search margin; independent validation pending")
  save_candidate("candidate_"*label*".toml",v.ep,v.kp,v.ki;extra)
  CSV.write(joinpath(OUT,"parameters_"*label*".csv"),DataFrame(bus=30:39,rho=1 .-v.ep,
    epsilon=v.ep,retained_SG_MW=CTX.power.*v.ep,Kp=v.kp,Ki=v.ki))
  println("FROZEN ",label," J=",v.J," alpha=",v.alpha," F=",v.peaks," R=",v.rocofs);flush(stdout)
  rows=NamedTuple[]
  for (name,kp,ki) in (("tuned",v.kp,v.ki),("nominal",fill(N.K0P,10),fill(N.K0I,10)),
      ("common_tuned_mean",fill(sum(v.kp)/10,10),fill(sum(v.ki)/10,10)))
    vv=oracle(CoreDesign.encode(CTX,v.ep,kp,ki);dt=.01,horizon=90.)
    push!(rows,(;gains=name,J=vv.J,alpha=vv.alpha,beta_sampled=vv.beta,
      F16=vv.peaks[1],F8=vv.peaks[2],F29=vv.peaks[3],max_R=maximum(vv.rocofs)))
  end
  CSV.write(joinpath(OUT,"gain_comparison_"*label*".csv"),DataFrame(rows))
  CSV.write(joinpath(OUT,"Aq_"*label*".csv"),DataFrame(v.A,:auto))
  t=@elapsed cert=full_band(v.A,CoreDesign.BETA_REQ)
  cert["elapsed_s"]=t;cert["candidate_sha256"]=bytes2hex(sha256(read(joinpath(OUT,"candidate_"*label*".toml"))))
  save_toml("robustness_"*label*".toml",cert)
  println("ROBUSTNESS ",cert);flush(stdout)
end
if abspath(PROGRAM_FILE)==(@__FILE__)
  main()
end
