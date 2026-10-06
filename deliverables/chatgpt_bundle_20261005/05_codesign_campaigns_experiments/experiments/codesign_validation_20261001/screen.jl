include("setup.jl")
rows=NamedTuple[]
for factor in (1.0,1.05,1.15), bus in (16,8,29)
    e=min.(factor.*ORIGINAL["epsilon"],0.999)
    kp=ORIGINAL["Kp"];ki=ORIGINAL["Ki"];rho=1 .-e
    sp=N.spectrum(CTX,rho,kp,ki)
    Aq=sp.quotient'*sp.model.Ared*sp.quotient
    beta=CoreDesign.robust_beta_sampled(Aq,sp.lambda)
    t=@elapsed tm=FiniteWindow.design_metrics(CTX,sp.model,rho;
      load_bus=bus,disturbance_MW=100.0,windows=(0.2,0.5,1.0,2.0),
      dt_s=0.01,horizon_s=60.0,gauge_vector=N.gauge_vector)
    for m in tm.metrics
      push!(rows,merge((;factor,event_bus=bus,J=dot(CTX.power,e),alpha=sp.alpha,
        beta_sampled=beta.beta,runtime_s=t),m))
    end
    row=only(filter(x->x.window_s==0.5,tm.metrics))
    println("SCREEN ",factor," bus=",bus," J=",dot(CTX.power,e)," alpha=",sp.alpha,
      " beta=",beta.beta," F=",row.F_peak_Hz," R=",row.R_peak_Hz_s," seconds=",t);flush(stdout)
    CSV.write(joinpath(OUT,"screen_corrected.csv"),DataFrame(rows))
end
save_toml("provenance.toml",Dict("julia"=>string(VERSION),"model_sha"=>ORIGINAL["model_sha"],
  "corrected_phase_sha"=>bytes2hex(sha256(read(joinpath(@__DIR__,"FiniteWindow.jl")))),
  "corrected_core_sha"=>bytes2hex(sha256(read(joinpath(@__DIR__,"CoreDesign.jl")))),
  "candidate_input_sha"=>bytes2hex(sha256(read(joinpath(@__DIR__,"frozen","original_candidate.toml")))),
  "Project_sha"=>bytes2hex(sha256(read(joinpath(ROOT,"Project.toml")))),
  "Manifest_sha"=>bytes2hex(sha256(read(joinpath(ROOT,"Manifest.toml"))))))
