if !isdefined(Main,:CTX);include("setup.jl");end
label=isempty(ARGS) ? "improved" : ARGS[1]
d=TOML.parsefile(joinpath(OUT,"candidate_"*label*".toml"));rho=d["rho"]
m=N.descriptor(CTX,rho,d["Kp"],d["Ki"]);rows=NamedTuple[]
for bus in (16,8,29)
  r=LinearSecurity.reduced_step_model(CTX,m,rho,bus;gauge_vector=N.gauge_vector)
  E=eigen(r.A);lam=E.values;V=E.vectors;z=V\r.B;res=(r.C_bus*V).*reshape(z,1,:)
  all(real.(lam).<0) || error("unstable tail")
  T=.5;t0=120.
  fss=100 .* real.(r.D_bus-(r.C_bus*V)*(z./lam))
  fc=abs.(res.*reshape(expm1.(lam*T)./(T.*lam.^2),1,:))
  rc=abs.(res.*reshape(expm1.(lam*T).^2 ./(T^2 .*lam.^2),1,:))
  ftail=100 .* (fc*exp.(real.(lam).*(t0-T)))
  rtail=100 .* (rc*exp.(real.(lam).*(t0-2T)))
  fine=FiniteWindow.design_metrics(CTX,m,rho;load_bus=bus,disturbance_MW=100.,windows=(.5,),
    dt_s=.0025,horizon_s=120.,gauge_vector=N.gauge_vector)
  ff=only(fine.metrics)
  row=(;event_bus=bus,horizon_s=t0,dt_s=.0025,Fpeak_refined_Hz=ff.F_peak_Hz,
    Rpeak_refined_Hz_s=ff.R_peak_Hz_s,tail_F_upper_Hz=maximum(abs.(fss)+ftail),
    tail_R_upper_Hz_s=maximum(rtail),
    pass=maximum(abs.(fss)+ftail)<.5 && maximum(rtail)<.5 && ff.F_peak_Hz<.5 && ff.R_peak_Hz_s<.5)
  push!(rows,row);println("TAIL_AUDIT ",row);flush(stdout)
end
CSV.write(joinpath(OUT,"linear_peak_and_tail_audit"*(label=="improved" ? "" : "_"*label)*".csv"),DataFrame(rows))
