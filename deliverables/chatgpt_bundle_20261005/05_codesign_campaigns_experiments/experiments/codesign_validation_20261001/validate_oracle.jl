if !isdefined(Main,:CTX);include("setup.jl");end
rows=NamedTuple[]
targets=isempty(ARGS) ? [("original",joinpath(@__DIR__,"frozen","original_candidate.toml")),
                    ("improved",joinpath(OUT,"candidate_improved.toml"))] : [(ARGS[1],joinpath(OUT,"candidate_"*ARGS[1]*".toml"))]
for (label,path) in targets
  d=TOML.parsefile(path);rho=d["rho"];m=N.descriptor(CTX,rho,d["Kp"],d["Ki"])
  for bus in (16,8,29)
    r=LinearSecurity.reduced_step_model(CTX,m,rho,bus;gauge_vector=N.gauge_vector)
    dv=-(m.Gy\LinearSecurity.load_input_vector(CTX,m,bus))
    jump=[imag(complex(dv[2k-1],dv[2k])/CTX.net.voltage[k]) for k in 30:39]
    n=size(r.A,1);aug=zeros(n+11,n+11)
    aug[1:n,1:n]=r.A;aug[1:n,end]=r.B
    aug[n+1:n+10,1:n]=2pi.*r.C_bus;aug[n+1:n+10,end]=2pi.*r.D_bus
    initial=vcat(zeros(n),jump,1.)
    times=[0.,.0001,.2,.5,.50001,1.,2.,10.,60.]
    modal=FiniteWindow.phase_step(CTX,m,rho,times;load_bus=bus,disturbance_MW=1.,gauge_vector=N.gauge_vector)
    exact=hcat([(exp(aug*t)*initial)[n+1:n+10] for t in times]...)
    err=maximum(abs.(modal.phase_rad-exact));rel=err/max(maximum(abs,exact),1e-12)
    row=(;candidate=label,event_bus=bus,max_abs_phase_error_rad=err,relative_phase_error=rel,
      modal_eigenvector_condition=modal.eigenvector_condition,pass=rel<1e-6)
    push!(rows,row);println("ORACLE_VALIDATION ",row);flush(stdout)
  end
end
CSV.write(joinpath(OUT,"independent_matrix_exponential"*(isempty(ARGS) ? "" : "_"*ARGS[1])*".csv"),DataFrame(rows))
all(r.pass for r in rows) || error("modal versus matrix exponential mismatch")
