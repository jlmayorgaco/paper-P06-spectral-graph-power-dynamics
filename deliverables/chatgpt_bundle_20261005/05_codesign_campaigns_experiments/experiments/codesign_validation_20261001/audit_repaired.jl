include("freeze_and_audit.jl")
d=TOML.parsefile(joinpath(OUT,"candidate_repaired.toml"))
ep=Float64.(d["epsilon"]);kp=Float64.(d["Kp"]);ki=Float64.(d["Ki"])
rows=NamedTuple[]
for (name,p,i) in (("tuned",kp,ki),("nominal",fill(N.K0P,10),fill(N.K0I,10)),
    ("common_tuned_mean",fill(sum(kp)/10,10),fill(sum(ki)/10,10)))
  v=oracle(CoreDesign.encode(CTX,ep,p,i);dt=.0025,horizon=120.)
  push!(rows,(;gains=name,J=v.J,alpha=v.alpha,beta_sampled=v.beta,
    F16=v.peaks[1],F8=v.peaks[2],F29=v.peaks[3],max_R=maximum(v.rocofs)))
  if name=="tuned"
    t=@elapsed cert=full_band(v.A,CoreDesign.BETA_REQ)
    cert["elapsed_s"]=t;cert["candidate_sha256"]=bytes2hex(sha256(read(joinpath(OUT,"candidate_repaired.toml"))))
    save_toml("robustness_repaired.toml",cert)
    println("REPAIRED_ROBUSTNESS ",cert);flush(stdout)
  end
end
CSV.write(joinpath(OUT,"gain_comparison_repaired.csv"),DataFrame(rows))
