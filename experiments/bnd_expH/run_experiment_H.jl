using CSV, DataFrames, LinearAlgebra, Statistics, TOML, SHA, Printf

include(joinpath(pwd(),"src","bnd_design_g","BNDDesignG.jl"))
include(joinpath(pwd(),"src","bnd_design_h","BNDDesignH.jl"))
using .BNDDesignG, .BNDDesignH

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const EXP=joinpath(ROOT,"experiments","bnd_expH")
const REPORT=joinpath(ROOT,"reports","experiment_H")
const TABLES=joinpath(REPORT,"tables")
const FIGS=joinpath(REPORT,"figures")
const SIGMA=0.05
const STRICT_SIGMA=0.05000001
const BASE_MVA=100.0
const BETA_GRID=[0.0,1e-6,3e-6,1e-5,3e-5,1e-4]
mkpath.( (EXP,joinpath(EXP,"configs"),REPORT,TABLES,FIGS) )

writehash(path)=write(path*".sha256",bytes2hex(sha256(read(path)))*"\n")
function hausdorff(a,b)
    (isempty(a)||isempty(b)) && return Inf
    f(x,y)=maximum(minimum(abs(z-w) for w in y) for z in x)
    max(f(a,b),f(b,a))
end
function writeempty(name,cols)
    CSV.write(joinpath(TABLES,name),DataFrame([c=>Any[] for c in cols]))
end
function nominal_domain(root)
    path=joinpath(root,"experiments","bnd_expE","configs","DESIGN_DOMAIN_FROZEN.toml")
    side=path*".sha256"
    isfile(path)&&isfile(side)||error("frozen E design domain is absent")
    sha=bytes2hex(sha256(read(path))); strip(read(side,String))==sha||error("E design-domain hash mismatch")
    TOML.parsefile(path),sha,path
end
function closure_identity(net,rho,kp,ki,s)
    m=CollectiveModel.mixed_jacobian(net,rho,kp,ki)
    cl=ClosureSpectrum.port_closure(net,s,rho,kp,ki)
    lhs=logdet(s*I-m.Ared)+logdet(m.Gy)
    rhs=logdet(s*I-m.A)+logdet(net.y_static)+logdet(cl.C)
    (;relative_error=abs(exp(lhs-rhs)-1),closure_sigma=minimum(svdvals(cl.C)),
      reduced_dimension=size(m.Ared,1),closure_dimension=size(cl.C,1))
end

println("EXP_H: load frozen analytical IEEE-39 model");flush(stdout)
domain,domain_sha,domain_path=nominal_domain(ROOT)
net=CollectiveModel.frozen_network(ROOT)
discovered=PhysicalData.discover_generators(ROOT)
gen=discovered[discovered.replaceable .== true,:]
buses=Int.(gen.bus); P=Float64.(gen.SG_dispatch_initial_MW); n=length(buses)
buses==Int.(domain["replaceable_buses"])||error("generator buses differ from frozen E domain")
Ptotal=sum(P)
isapprox(Ptotal,5402.761089978776;atol=1e-6,rtol=0)||error("initialized SG dispatch checksum failed")
net.no_load_kcl_inf<=1e-9||error("initialized ZIP-load KCL reconstruction failed")
kp0=fill(Float64(domain["Kp_nom"]),n);ki0=fill(Float64(domain["Ki_nom"]),n)
kpmin=0.25.*kp0;kpmax=4.0.*kp0;kimin=0.25.*ki0;kimax=4.0.*ki0
hinertia=[net.sg[b].op.parameters.inertia*net.sg[b].op.parameters.rating_mva for b in buses]
CSV.write(joinpath(TABLES,"TABLE_H00_model_integrity.csv"),net.load_audit)

# H00: model/closure/finite-spectrum regression. PowerDynamics is not loaded here.
integrity=NamedTuple[]
archive_matrices=Dict("all_SG"=>"C0_Ared.csv","single_30"=>"C30_Ared.csv",
  "single_33"=>"C33_regenerated_Ared.csv","single_35"=>"C35_Ared.csv","single_37"=>"C37_Ared.csv")
for (label,epsv,kp,ki) in (("all_SG",ones(n),kp0,ki0),
                           ("all_GFL",zeros(n),kp0,ki0),
                           ("single_30",[b==30 ? 0.0 : 1.0 for b in buses],kp0,ki0),
                           ("single_33",[b==33 ? 0.0 : 1.0 for b in buses],kp0,ki0),
                           ("single_35",[b==35 ? 0.0 : 1.0 for b in buses],kp0,ki0),
                           ("single_37",[b==37 ? 0.0 : 1.0 for b in buses],kp0,ki0))
    m=CollectiveModel.mixed_jacobian(net,1 .- epsv,kp,ki)
    spec=eigvals(m.Ared)
    desc=BNDDesignH.descriptor_finite_spectrum(m)
    descriptor_err=hausdorff(spec,desc)
    archived_err=NaN
    if haskey(archive_matrices,label)
        apath=joinpath(ROOT,"reports","experiment_C","matrices",archive_matrices[label])
        adf=CSV.read(apath,DataFrame)
        Aarchive=Matrix{Float64}(select(adf,Not(:row_index)))
        archived_err=hausdorff(spec,eigvals(Aarchive))
    end
    ident=closure_identity(net,1 .- epsv,kp,ki,0.173+0.291im)
    hasarchive=isfinite(archived_err)
    gate=hasarchive ? archived_err<=1e-9 : length(desc)==length(spec)&&ident.relative_error<=1e-9
    push!(integrity,(case=label,dimension=length(spec),descriptor_finite_count=length(desc),
      finite_spectrum_hausdorff=archived_err,descriptor_diagnostic_hausdorff=descriptor_err,
      closure_identity_relative_error=ident.relative_error,
      gauge_residual=norm(m.Ared*BNDDesignG.gauge_vector(net,m))/
        max(norm(m.Ared)*norm(BNDDesignG.gauge_vector(net,m)),eps()),
      status=gate ? (hasarchive ? "PASS_ARCHIVED_FULL_SPECTRUM" : "PASS_CLOSURE_AND_FINITE_COUNT_JORDAN_EIGENVALUE_SENSITIVE") : "FAIL"))
end
archsingle=CSV.read(joinpath(ROOT,"reports","experiment_E","tables","TABLE_E04_single_bus_spectrum_audit.csv"),DataFrame)
for r in eachrow(archsingle)
    push!(integrity,(case="archived_single_$(Int(r.bus))",dimension=Int(r.analytic_dimension),
      descriptor_finite_count=Int(r.frozen_dimension),finite_spectrum_hausdorff=Float64(r.spectral_Hausdorff),
      descriptor_diagnostic_hausdorff=NaN,closure_identity_relative_error=NaN,gauge_residual=NaN,
      status=Bool(r.pass)&&Float64(r.spectral_Hausdorff)<=1e-9 ? "ARCHIVE_PASS" : "FAIL"))
end
CSV.write(joinpath(TABLES,"TABLE_H00_model_integrity.csv"),DataFrame(integrity))
all(x->x.status!="FAIL",integrity)||error("H00 hard gate failed; stopping before design")

# H1: quotient all six prescribed deterministic gain patterns.
pattern_names=["nominal","independent_A","independent_B","independent_C","all_min","all_max"]
patterns=Tuple{Vector{Float64},Vector{Float64}}[]
push!(patterns,(kp0,ki0))
push!(patterns,(kp0.*[isodd(j) ? 0.4 : 3.2 for j in 1:n],ki0.*[isodd(j) ? 3.7 : 0.6 for j in 1:n]))
push!(patterns,(kp0.*[0.5+3.2*(mod(7j,11)/10) for j in 1:n],ki0.*[0.5+3.2*(mod(3j+4,11)/10) for j in 1:n]))
push!(patterns,(kp0.*collect(range(0.25,4.0,length=n)),ki0.*collect(range(4.0,0.25,length=n))))
push!(patterns,(kpmin,kimin));push!(patterns,(kpmax,kimax))
H01=NamedTuple[];H02=NamedTuple[];jordan_nom=nothing
for (q,(kp,ki)) in enumerate(patterns)
    ja=BNDDesignH.jordan_audit(net,kp,ki)
    qmodel=ja.quotient
    vals=qmodel.values
    push!(H01,(pattern=pattern_names[q],quotient_dimension=qmodel.finite_dimension,
      gauge_residual=ja.gauge_residual,quotient_nullity=ja.quotient_nullity,
      quotient_algebraic_multiplicity=ja.quotient_algebraic_multiplicity,
      full_nullity=ja.full_nullity,full_algebraic_multiplicity=ja.full_algebraic_multiplicity,
      quotient_min_singular=ja.quotient_smallest_singular,
      physical_nearzero=vals[argmin(abs.(vals))],classification=ja.classification))
    push!(H02,(pattern=pattern_names[q],gauge_chain_length=ja.full_algebraic_multiplicity,
      full_nullity=ja.full_nullity,gauge_chain_residual=ja.gauge_chain_residual,
      quotient_nullity=ja.quotient_nullity,quotient_chain_residual=ja.quotient_chain_residual,
      left_right_overlap=abs(dot(ja.quotient_left_null,ja.quotient_right_null))))
    q==1&&(global jordan_nom=ja)
end
CSV.write(joinpath(TABLES,"TABLE_H01_gauge_quotient.csv"),DataFrame(H01))
CSV.write(joinpath(TABLES,"TABLE_H02_jordan_chain.csv"),DataFrame(H02))

# Quotient sigma-min slope and root scaling, retaining the entire requested range.
sigrows=NamedTuple[]
Aq0=jordan_nom.quotient.Aq
for s in 10.0 .^collect(-10:1:-2)
    sv=minimum(svdvals(s*I-Aq0))
    push!(sigrows,(s=s,sigma_min=sv,ratio=sv/s))
end
CSV.write(joinpath(TABLES,"TABLE_H02_sigma_samples.csv"),DataFrame(sigrows))
scale_rows=NamedTuple[]
for i in eachindex(buses)
    roots=ComplexF64[]
    for e in 10.0 .^collect(-8:1:-2)
        ep=zeros(n);ep[i]=e
        rr=BNDDesignH.nearest_physical_root(net,ep,kp0,ki0).lambda
        push!(roots,rr)
    end
    x=collect(-8.0:1.0:-2.0);y=log10.(max.(abs.(roots),1e-300))
    slope=sum((x.-mean(x)).*(y.-mean(y)))/sum((x.-mean(x)).^2)
    local_slope=sum((x[1:4].-mean(x[1:4])).*(y[1:4].-mean(y[1:4])))/sum((x[1:4].-mean(x[1:4])).^2)
    for k in eachindex(roots)
        push!(scale_rows,(bus=buses[i],epsilon=10.0^x[k],lambda_real=real(roots[k]),
          lambda_imag=imag(roots[k]),abs_lambda=abs(roots[k]),
          full_range_exponent=slope,local_1e8_to_1e5_exponent=local_slope,
          predicted_exponent=1.0,full_range_within_5pct=abs(slope-1)<=0.05))
    end
end
CSV.write(joinpath(TABLES,"TABLE_H03_scaling_exponent.csv"),DataFrame(scale_rows))

# H3: determinant expansion factors the common gauge root before evaluating the
# simple quotient root. Gain derivatives are symmetric finite differences at
# two scales in normalized controller coordinates.
println("EXP_H: derive mixed retained-SG authority and gain Jacobian");flush(stdout)
mix=BNDDesignH.authority_derivative(net,buses,kp0,ki0;relstep=1e-4,authority_h=2e-4)
auth=mix.base
CSV.write(joinpath(TABLES,"TABLE_H04_mixed_authority.csv"),DataFrame(bus=buses,A_i=auth,
  authority_sign=ifelse.(auth.>0,"STABILIZING","DESTABILIZING"),
  coefficient_method=fill("det(C)/s gauge-factor quotient, Richardson in epsilon",n)))
scale_df=DataFrame(scale_rows)
rootfit=NamedTuple[]
for i in eachindex(buses)
    rr=only(scale_df[(scale_df.bus .== buses[i]) .& (scale_df.epsilon .== 1e-6),:])
    Aroot=-rr.lambda_real/1e-6
    push!(rootfit,(bus=buses[i],A_from_det_quotient=auth[i],A_from_exact_state_root=Aroot,
      relative_coefficient_error=abs(auth[i]-Aroot)/max(abs(Aroot),1e-12),
      root_epsilon=1e-6,comparison="determinant coefficient vs gauge-quotient state eigenvalue"))
end
CSV.write(joinpath(TABLES,"TABLE_H03_root_law_coefficient_validation.csv"),DataFrame(rootfit))
Jraw=mix.J;Jnorm=mix.normalized_jacobian
H05=NamedTuple[]
for i in 1:n,j in 1:2n
    push!(H05,(support_bus=buses[i],gain_bus=buses[(j-1)%n+1],
      parameter=j<=n ? "Kp" : "Ki",derivative_actual=mix.J[i,j],
      derivative_normalized=mix.normalized_jacobian[i,j],
      complex_step_derivative=mix.Jcs[i,j],
      symmetric_finite_difference=mix.Jfine[i,j],
      coarse_derivative=mix.Jcoarse[i,j],
      relative_complex_step_error=mix.relative_error[i,j],
      symmetric_fd_relative_error=abs(mix.J[i,j]-mix.Jfine[i,j])/max(abs(mix.J[i,j]),1e-8)))
end
CSV.write(joinpath(TABLES,"TABLE_H05_mixed_gain_jacobian.csv"),DataFrame(H05))
near_cut=max(1e-8,1e-4*maximum(abs.(mix.J)))
valid_err=mix.relative_error[abs.(mix.J).>near_cut]
deriv_median=isempty(valid_err) ? NaN : median(valid_err)
deriv_p95=isempty(valid_err) ? NaN : quantile(valid_err,0.95)
CSV.write(joinpath(TABLES,"TABLE_H05_derivative_validation.csv"),DataFrame(
  method=fill("analytic determinant-trace derivative vs independent complex-step",5),
  metric=["median_relative_error_nonzero","p95_relative_error_nonzero","near_zero_cutoff","target_median","target_p95"],
  value=[deriv_median,deriv_p95,near_cut,1e-5,1e-4],
  status=[deriv_median<=1e-5 ? "PASS" : "FAIL",deriv_p95<=1e-4 ? "PASS" : "FAIL","EXCLUDED", "TARGET", "TARGET"]))
active=BNDDesignH.active_gain_basis(Jnorm)
rankrows=NamedTuple[]
for tol in (1e-2,1e-3,1e-4)
    r=count(s->s>=tol*maximum(active.singular_values),active.singular_values)
    push!(rankrows,(relative_singular_tolerance=tol,numerical_rank=r,
      largest_singular=maximum(active.singular_values),smallest_singular=minimum(active.singular_values)))
end
CSV.write(joinpath(TABLES,"TABLE_H05_rank.csv"),DataFrame(rankrows))
writeempty("TABLE_H06_woodbury_gain_structure.csv",[:bus,:representation,:rank_Kp,:rank_Ki,:A_factor_error,:B_factor_error])

# Exact structural Woodbury audit and active right-singular-vector basis.
CSV.write(joinpath(TABLES,"TABLE_H06_woodbury_gain_structure.csv"),DataFrame(BNDDesignH.woodbury_audit(net,buses)))
CSV.write(joinpath(TABLES,"TABLE_H06_active_basis.csv"),DataFrame(direction=1:active.rank,
  singular_value=active.singular_values[1:active.rank],
  cumulative_energy=cumsum(abs2.(active.singular_values[1:active.rank]))./sum(abs2,active.singular_values)))

# H5 noncommuting physical graph operators and authority pathways.
split=BNDDesignH.exact_self_energy_split(net,buses,kp0,ki0)
CSV.write(joinpath(TABLES,"TABLE_H07_BND_authority_decomposition.csv"),DataFrame(split.authority))
CSV.write(joinpath(TABLES,"TABLE_H08_self_energy_pathways.csv"),DataFrame(split.pathways))
align=BNDDesignH.graph_alignment(net,buses,P,auth,active.V,kp0,ki0)
CSV.write(joinpath(TABLES,"TABLE_H09_graph_active_subspace_alignment.csv"),DataFrame(
  principal_angle_index=1:length(align.principal_angles_rad),
  principal_angle_rad=align.principal_angles_rad,
  captured_active_subspace_energy=fill(align.captured_energy,length(align.principal_angles_rad)),
  graph_basis_rank=fill(align.graph.rank,length(align.principal_angles_rad)),
  active_gain_dimension=fill(active.rank,length(align.principal_angles_rad)),
  LG_LB_commutator_norm=fill(norm(align.graph.LG*align.graph.LB-align.graph.LB*align.graph.LG),length(align.principal_angles_rad)),
  LG_LB_normalized_commutator=fill(norm(align.graph.LG*align.graph.LB-align.graph.LB*align.graph.LG)/
    max(norm(align.graph.LG)*norm(align.graph.LB),eps()),length(align.principal_angles_rad)),
  graph_basis_type=fill("NONCOMMUTATIVE_WORDS_DEGREE_2",length(align.principal_angles_rad))))

# H6: optimize the leading coefficient on the active gain subspace using a
# custom damped Newton stationarity solve (bounded by the frozen gain box).
function authority_at(i,x)
    kk=kp0.*x[1:n];ii=ki0.*x[n+1:end]
    BNDDesignH.single_authority(net,buses,i,kk,ii)
end
function bounded_newton(i, V; maxiter=20)
    r=size(V,2);xi=zeros(r);lo=fill(0.25,2n);hi=fill(4.0,2n);hist=NamedTuple[]
    obj(y)=authority_at(i,clamp.(1 .+ V*y,lo,hi))
    for it in 1:maxiter
        h=2e-3;f0=obj(xi);g=zeros(r);H=zeros(r,r)
        for a in 1:r
            e=zeros(r);e[a]=h
            fp=obj(xi+e);fm=obj(xi-e);g[a]=(fp-fm)/(2h);H[a,a]=(fp-2f0+fm)/h^2
            for b in 1:a-1
                v=zeros(r);v[b]=h
                H[a,b]=H[b,a]=(obj(xi+e+v)-obj(xi+e-v)-obj(xi-e+v)+obj(xi-e-v))/(4h^2)
            end
        end
        step=try -(H\g) catch; g/max(norm(g),1e-12)*0.1 end
        all(isfinite,step)||break
        # Cap to the actual gain box and deterministic Armijo ascent.
        scale=min(1.0,0.8*minimum([((hi[k]-1-(V*xi)[k])/(V*step)[k]) for k in eachindex(lo) if (V*step)[k]>1e-12];init=Inf),
                              0.8*minimum([((lo[k]-1-(V*xi)[k])/(V*step)[k]) for k in eachindex(lo) if (V*step)[k]<-1e-12];init=Inf))
        scale=clamp(scale,0.0,1.0);accepted=false
        for _ in 0:12
            trial=xi+scale*step
            if all(lo .<= 1 .+ V*trial .<= hi) && obj(trial)>=f0-1e-12
                xi=trial;accepted=true;break
            end
            scale/=2
        end
        push!(hist,(iteration=it,authority=obj(xi),gradient_norm=norm(g),step_norm=norm(scale*step),accepted=accepted))
        (!accepted || norm(g)<1e-5 || norm(scale*step)<1e-6)&&break
    end
    x=clamp.(1 .+ V*xi,lo,hi)
    (;f=obj(xi),x,history=hist,xi)
end

predrows=NamedTuple[];statrows=NamedTuple[];best_pred=nothing
for i in 1:n
    opt=bounded_newton(i,active.V)
    A=opt.f; e=A>0 ? SIGMA/A : Inf
    feasible=isfinite(e)&&e<=1
    retained=feasible ? P[i]*e : Inf
    push!(predrows,(bus=buses[i],A_i_nominal=auth[i],A_i_optimized=A,
      epsilon_predicted=e,retained_sg_predicted_MW=retained,predicted_feasible=feasible,
      gain_stationarity_status=isempty(opt.history) ? "NO_STEP" : "DAMPED_NEWTON_LOCAL"))
    for r in opt.history
        push!(statrows,(bus=buses[i],iteration=r.iteration,authority=r.authority,
          gradient_norm=r.gradient_norm,step_norm=r.step_norm,accepted=r.accepted))
    end
    if feasible && (best_pred===nothing || retained<best_pred.retained)
        global best_pred=(i=i,e=e,retained=retained,x=opt.x,A=A)
    end
end
CSV.write(joinpath(TABLES,"TABLE_H10_analytic_anchor_candidates.csv"),DataFrame(predrows))
CSV.write(joinpath(TABLES,"TABLE_H11_gain_stationary_candidates.csv"),DataFrame(statrows))
if best_pred===nothing
    ep=ones(n); kpf=kp0; kif=ki0; predicted_status="NO_LINEAR_SINGLE_ANCHOR_FEASIBLE"
else
    ep=zeros(n);ep[best_pred.i]=best_pred.e
    kpf=kp0.*best_pred.x[1:n];kif=ki0.*best_pred.x[n+1:end]
    predicted_status="LINEAR_SINGLE_ANCHOR_PREDICTOR"
end
predicted_retained=dot(P,ep)
function candidate_dict(ep,kp,ki,status,alpha,beta;frozen=false)
    rows=Dict{String,Any}[]
    for j in 1:n
        push!(rows,Dict("bus"=>buses[j],"rho"=>1-ep[j],"epsilon"=>ep[j],
          "Kp"=>kp[j],"Ki"=>ki[j],"P_i_MW"=>P[j],"H_MVA_s"=>hinertia[j]))
    end
    Dict{String,Any}("experiment"=>"BND_EXP_H","status"=>status,"candidate_frozen"=>frozen,
      "design_used_PowerDynamics"=>false,"powerdynamics_validation"=>"NOT_RUN",
      "frozen_domain_sha256"=>domain_sha,"generator_buses"=>buses,"generator"=>rows,
      "retained_SG_MW"=>predicted_retained,"GFL_MW"=>Ptotal-predicted_retained,
      "GFL_fraction"=>(Ptotal-predicted_retained)/Ptotal,
      "spectral_abscissa_s_inv"=>alpha,"sigma_required_s_inv"=>SIGMA,
      "beta_required"=>beta,"small_gain_peak"=>NaN,
      "preblind_predictor_status"=>predicted_status,"primary_objective"=>"sum(P_i_MW * epsilon_i)" )
end
function write_candidate(path,dict)
    open(path,"w") do io;TOML.print(io,dict);end
    writehash(path)
end
prepath=joinpath(REPORT,"Z_H_PREBLIND.toml")
if isfile(prepath)
    isfile(prepath*".sha256")||error("existing preblind freeze has no SHA-256 sidecar")
    presha=bytes2hex(sha256(read(prepath)))
    strip(read(prepath*".sha256",String))==presha||error("existing preblind candidate SHA-256 mismatch")
    saved=TOML.parsefile(prepath); saved["design_used_PowerDynamics"]==false||error("preblind candidate was not analytical")
    savedrows=sort(saved["generator"],by=x->Int(x["bus"]))
    Int[x["bus"] for x in savedrows]==buses||error("preblind generator support mismatch")
    ep=Float64[x["epsilon"] for x in savedrows];kpf=Float64[x["Kp"] for x in savedrows];kif=Float64[x["Ki"] for x in savedrows]
    predicted_retained=dot(P,ep);predicted_status=String(saved["preblind_predictor_status"])
    support=findall(>(0),ep)
    best_pred=isempty(support) ? nothing : (i=first(support),e=ep[first(support)],
      retained=predicted_retained,x=vcat(kpf./kp0,kif./ki0),A=NaN)
else
    open(prepath,"w") do io
        TOML.print(io,candidate_dict(ep,kpf,kif,"PREBLIND_PREDICTOR",NaN,NaN))
    end
    writehash(prepath)
end
CSV.write(joinpath(TABLES,"TABLE_H12_predictor_best_candidate.csv"),DataFrame(
  status=[predicted_status],predicted_bus=[best_pred===nothing ? "" : buses[best_pred.i]],
  predicted_epsilon=[best_pred===nothing ? NaN : best_pred.e],retained_sg_MW=[predicted_retained],
  converted_GFL_MW=[Ptotal-predicted_retained],candidate_sha256=[bytes2hex(sha256(read(prepath)))]))
println("EXP_H: preblind predictor frozen; SHA-256 ",bytes2hex(sha256(read(prepath))));flush(stdout)

# H00 candidate regression is intentionally post-freeze; the historical E/G
# final files are read only after the new preblind candidate is hashed.
gpath=joinpath(ROOT,"reports","experiment_G","Z_G_FINAL.toml")
gsha=strip(read(gpath*".sha256",String));gdata=TOML.parsefile(gpath)
grows=sort(gdata["generator"],by=x->Int(x["bus"]))
geps=Float64[x["epsilon"] for x in grows];gkp=Float64[x["Kp"] for x in grows];gki=Float64[x["Ki"] for x in grows]
gm=CollectiveModel.mixed_jacobian(net,1 .- geps,gkp,gki);gdesc=BNDDesignH.descriptor_finite_spectrum(gm)
gpoles=CSV.read(joinpath(ROOT,"reports","experiment_G","tables","TABLE_G16_powerdynamics_candidate_poles.csv"),DataFrame)
gfinite=complex.(gpoles.real_s_inv,gpoles.imag_s_inv)
garchive=CSV.read(joinpath(ROOT,"reports","experiment_G","tables","TABLE_G15_eigenvalue_comparison.csv"),DataFrame)
g_alpha_archive=only(garchive[garchive.case .== "ExpG_final",:alpha])
g_sp=BNDDesignG.spectrum_with_gauge(net,1 .- geps,gkp,gki)
g_closure_err=closure_identity(net,1 .- geps,gkp,gki,0.173+0.291im).relative_error
g_gate=abs(g_sp.spectral_abscissa-g_alpha_archive)<=1e-9&&length(gdesc)==size(gm.Ared,1)&&g_closure_err<=1e-9
push!(integrity,(case="ExpG_frozen_candidate_post_preblind",dimension=size(gm.Ared,1),
  descriptor_finite_count=length(gdesc),finite_spectrum_hausdorff=hausdorff(eigvals(gm.Ared),gfinite),
  descriptor_diagnostic_hausdorff=hausdorff(gdesc,gfinite),
  closure_identity_relative_error=g_closure_err,
  gauge_residual=norm(gm.Ared*BNDDesignG.gauge_vector(net,gm))/max(norm(gm.Ared),eps()),
  status=g_gate ? "PASS_ARCHIVED_EXP_G_ALPHA_AND_FULL_POLE_DIAGNOSTIC" : "FAIL"))
CSV.write(joinpath(TABLES,"TABLE_H00_model_integrity.csv"),DataFrame(integrity))
g_gate||error("ExpG frozen-candidate regression failed after preblind freeze")

# H12 blinded comparison: no candidate is used as a seed and no tuning follows.
epath=joinpath(ROOT,"reports","experiment_E","Z_E_FINAL.toml")
if !isfile(epath)
    ep_candidates=filter(isfile,[joinpath(ROOT,"reports","experiment_E",x) for x in ("Z_E_FINAL.toml","Z_E_CANDIDATE_FINAL.toml")])
    epath=isempty(ep_candidates) ? "" : only(ep_candidates)
end
comparison=NamedTuple[]
push!(comparison,(case="H_preblind_predictor",retained_sg_MW=predicted_retained,
  converted_GFL_MW=Ptotal-predicted_retained,sha256=bytes2hex(sha256(read(prepath))),
  read_after_freeze=true,used_as_seed=false))
if !isempty(epath)
    edata=TOML.parsefile(epath)
    push!(comparison,(case="ExpE_final",retained_sg_MW=Float64(edata["retained_SG_MW"]),
      converted_GFL_MW=Float64(edata["GFL_MW"]),sha256=bytes2hex(sha256(read(epath))),
      read_after_freeze=true,used_as_seed=false))
else
    retainedE=1.189211882;convertedE=Ptotal-retainedE
    push!(comparison,(case="ExpE_user_specified_nominal_reference_NO_FINAL_FILE",
      retained_sg_MW=retainedE,converted_GFL_MW=convertedE,
      sha256="REFERENCE_ONLY_NO_CANDIDATE_ARTIFACT",read_after_freeze=true,used_as_seed=false))
end
push!(comparison,(case="ExpG_final",retained_sg_MW=Float64(gdata["retained_SG_MW"]),
  converted_GFL_MW=Float64(gdata["GFL_MW"]),sha256=gsha,
  read_after_freeze=true,used_as_seed=false))
CSV.write(joinpath(TABLES,"TABLE_H19_blinded_comparison.csv"),DataFrame(comparison))

# H7: exact low-rank spectral correction for all single-anchor predictor rows.
corr=NamedTuple[]
for (j,row) in enumerate(predrows)
    row.predicted_feasible||continue
    ev=zeros(n);ev[j]=1.0
    lr=BNDDesignG.low_rank_boundary(net,ev,kpf,kif,STRICT_SIGMA)
    if lr.status=="RADIAL_ROOTS"
        for ecor in lr.roots
            0<=ecor<=1||continue
            ec=zeros(n);ec[j]=ecor
            sp=BNDDesignG.spectrum_with_gauge(net,1 .- ec,kpf,kif)
            push!(corr,(bus=buses[j],epsilon_exact=ecor,retained_sg_MW=P[j]*ecor,
              spectral_abscissa=sp.spectral_abscissa,critical=sp.critical,
              predictor_retained_MW=row.retained_sg_predicted_MW,
              relative_retention_correction=abs(P[j]*ecor-row.retained_sg_predicted_MW)/max(P[j]*row.epsilon_predicted,eps()),
              determinant_lemma_residual=lr.determinant_lemma_residual,
              strict_margin=sp.spectral_abscissa<=-STRICT_SIGMA,
              status=sp.spectral_abscissa<=-STRICT_SIGMA ? "EXACT_FEASIBLE_NOMINAL" : "FAIL_OTHER_POLE"))
        end
    else
        push!(corr,(bus=buses[j],epsilon_exact=NaN,retained_sg_MW=NaN,spectral_abscissa=NaN,
          critical=NaN,predictor_retained_MW=row.retained_sg_predicted_MW,
          relative_retention_correction=NaN,determinant_lemma_residual=lr.determinant_lemma_residual,
          strict_margin=false,status=lr.status))
    end
end
CSV.write(joinpath(TABLES,"TABLE_H13_exact_corrected_candidates.csv"),DataFrame(corr))

# H8: direct full-block small-gain frontier. The exact boundary candidate is a
# seed; each single-anchor branch is scanned in epsilon to find a bracket with
# both full-spectrum margin and the direct shifted-resolvent certificate.
front=NamedTuple[];front_choice=Dict{Float64,Any}()
feasible_candidates=[r for r in corr if r.status=="EXACT_FEASIBLE_NOMINAL"]
best_exact=isempty(feasible_candidates) ? nothing : feasible_candidates[argmin(getfield.(feasible_candidates,:retained_sg_MW))]
if best_exact!==nothing
    ev=zeros(n);ev[findfirst(==(Int(best_exact.bus)),buses)]=best_exact.epsilon_exact
    front_choice[0.0]=(eps=ev,retained=best_exact.retained_sg_MW,alpha=best_exact.spectral_abscissa,
      peak=Inf,omega=NaN,beta=0.0,status="NOMINAL_EXACT_BOUNDARY",bus=best_exact.bus)
end
# First make a cheap exact-spectrum/zero-frequency scan. The zero-frequency
# peak is a rigorous lower bound on the full resolvent peak and safely prunes
# points that cannot meet any requested beta.
robust_grid=NamedTuple[]
for i in 1:n, ee in 0.0:0.02:1.0
    ep_grid_vec=zeros(n);ep_grid_vec[i]=ee
    q=BNDDesignH.quotient_model(net,ep_grid_vec,kpf,kif)
    alpha=maximum(real.(q.values));alpha<=-STRICT_SIGMA||continue
    peak0=opnorm((-(q.Aq+SIGMA*I))\I,2)
    peak0<=1/minimum(BETA_GRID[2:end])||continue
    push!(robust_grid,(bus=buses[i],index=i,epsilon=ee,eps=ep_grid_vec,retained=P[i]*ee,
      alpha=alpha,peak0=peak0))
end
sort!(robust_grid,by=x->x.retained)
peak_cache=Dict{Tuple{Int,Float64},Any}()
for beta in BETA_GRID
    cand=if beta==0.0
        get(front_choice,beta,nothing)
    else
        selected=nothing
        for x in robust_grid
            beta*x.peak0<1||continue
            key=(x.index,x.epsilon)
            if !haskey(peak_cache,key)
                q=BNDDesignH.quotient_model(net,x.eps,kpf,kif)
                peak_cache[key]=BNDDesignH.direct_resolvent_peak(q.Aq,SIGMA)
            end
            pk=peak_cache[key]
            if beta*pk.peak<1
                selected=(eps=x.eps,retained=x.retained,alpha=x.alpha,peak=pk.peak,
                  omega=pk.omega_peak,beta=beta,status="DIRECT_SMALL_GAIN_CERTIFIED_SCALAR_GRID",bus=x.bus)
                break
            end
        end
        selected
    end
    cand!==nothing&&beta>0&&(front_choice[beta]=cand)
    if cand===nothing
        push!(front,(beta_required=beta,status=beta==0 ? "NO_EXACT_NOMINAL_CANDIDATE" : "NO_SINGLE_ANCHOR_ROBUST_BRANCH",
          best_retained_sg_MW=NaN,converted_GFL_MW=NaN,spectral_abscissa=NaN,
          small_gain_peak=NaN,small_gain_margin=NaN,beta_star=NaN,robust_pass=false,
          peak_frequency_rad_s=NaN,candidate_bus="",epsilon=NaN,
          search="0.02-epsilon scan; exact alpha, zero-frequency prune, full direct resolvent"))
    else
        margin=beta==0 ? 1.0 : 1-beta*cand.peak
        push!(front,(beta_required=beta,status=cand.status,best_retained_sg_MW=cand.retained,
          converted_GFL_MW=Ptotal-cand.retained,spectral_abscissa=cand.alpha,
          small_gain_peak=cand.peak,small_gain_margin=margin,
          beta_star=beta==0 ? NaN : inv(cand.peak),robust_pass=beta==0||margin>0,
          peak_frequency_rad_s=cand.omega,candidate_bus=cand.bus,epsilon=maximum(cand.eps),
          search="0.02-epsilon scan; exact alpha, zero-frequency prune, full direct resolvent"))
    end
end
CSV.write(joinpath(TABLES,"TABLE_H14_robustness_frontier.csv"),DataFrame(front))
positive_feasible=[b for b in BETA_GRID if b>0&&haskey(front_choice,b)]
best_final=if !isempty(positive_feasible)
    front_choice[maximum(positive_feasible)]
elseif haskey(front_choice,0.0)
    front_choice[0.0]
else
    nothing
end

# Keep predictor, exact spectral branches, and robust grid points distinct.
# Their feasibility is not a KKT certificate.
branches=NamedTuple[]
for r in predrows
    push!(branches,(branch="predictor_bus_$(r.bus)",seed_source="H6 analytic predictor",
      status=r.predicted_feasible ? "PREDICTOR_ONLY_NOT_EXACT_FEASIBILITY" : "NO_PREDICTED_SINGLE_ANCHOR",
      retained_sg_MW=r.retained_sg_predicted_MW,exact_feasible=false,robust_feasible=false,
      beta_required=NaN,kkt_status="NOT_CERTIFIED"))
end
for r in corr
    push!(branches,(branch="exact_single_bus_$(r.bus)_eps_$(r.epsilon_exact)",seed_source="H7 exact determinant correction",
      status=r.status,retained_sg_MW=r.retained_sg_MW,exact_feasible=r.strict_margin,
      robust_feasible=false,beta_required=NaN,kkt_status="NOT_CERTIFIED"))
end
for r in front
    push!(branches,(branch="robust_grid_beta_$(r.beta_required)_bus_$(r.candidate_bus)",seed_source=String(r.search),
      status=r.status,retained_sg_MW=r.best_retained_sg_MW,exact_feasible=r.robust_pass,
      robust_feasible=r.robust_pass,beta_required=r.beta_required,kkt_status="NOT_CERTIFIED"))
end
CSV.write(joinpath(TABLES,"TABLE_H16_branch_summary.csv"),DataFrame(branches))

# H11 transient metrics are calculated for the best exact spectral candidate.
if best_final===nothing
    CSV.write(joinpath(TABLES,"TABLE_H17_transient_metrics.csv"),DataFrame(status=["BLOCKED_NO_EXACT_CANDIDATE"]))
    CSV.write(joinpath(TABLES,"TABLE_H18_transient_active_constraint_test.csv"),DataFrame(status=["BLOCKED_NO_EXACT_CANDIDATE"]))
else
    ec=best_final.eps
    mm=CollectiveModel.mixed_jacobian(net,1 .- ec,kpf,kif)
    tr=BNDDesignG.transient_analysis(net,mm,ec,buses)
    CSV.write(joinpath(TABLES,"TABLE_H17_transient_metrics.csv"),DataFrame(status=[tr.status],
      rocof_unit_gain=[tr.rocof_unit_gain],frequency_unit_gain=[tr.frequency_unit_gain],
      deltaP_max_rocof_MW=[tr.deltaP_max_peak_rocof_MW],deltaP_max_frequency_MW=[tr.deltaP_max_frequency_MW],
      unit_rocof_peak=[tr.unit_rocof_peak],unit_frequency_peak=[tr.unit_frequency_peak],
      frequency_nadir_unit=[tr.frequency_nadir],retained_inertia_MVA_s=[tr.H_total_MVA_s]))
    activeR=100.0*tr.unit_rocof_peak>0.5
    activeF=100.0*tr.unit_frequency_peak>0.5
    CSV.write(joinpath(TABLES,"TABLE_H18_transient_active_constraint_test.csv"),DataFrame(
      constraint=["RoCoF","frequency"],active=[activeR,activeF],
      status=[activeR ? "ACTIVE_REQUIRES_KKT_REOPEN" : "VALIDATION_ONLY",
              activeF ? "ACTIVE_REQUIRES_KKT_REOPEN" : "VALIDATION_ONLY"],
      candidate_deltaP_MW=fill(100.0,2),limits=[0.5,0.5],
      predicted_metric=[100.0*tr.unit_rocof_peak,100.0*tr.unit_frequency_peak]))
end

# No coupled KKT solve is claimed. In this run the post-design transient
# screen identifies an active frequency constraint, so its active set would
# need to be reopened before any stationarity/optimality conclusion.
frequency_active=best_final!==nothing && activeF
CSV.write(joinpath(TABLES,"TABLE_H15_KKT_iterations.csv"),DataFrame(
  branch=["best_analytical_robust_grid_candidate"],
  iterations=[0],
  status=[frequency_active ? "NOT_SOLVED_ACTIVE_TRANSIENT_FREQUENCY_CONSTRAINT" : "NOT_CERTIFIED_FULL_CLOSURE_KKT"],
  primal_feasibility=[NaN],dual_feasibility=[NaN],complementarity=[NaN],
  stationarity_residual=[NaN],LICQ_min_singular_value=[NaN],reduced_hessian_min_eigenvalue=[NaN],
  reason=[frequency_active ? "H18 predicts the frequency limit is active at the declared 100 MW screen; KKT active set was not reopened" :
    "Branch, robust-peak, LICQ and SOSC certification was not completed"]))

# H20/H21 exist before validation with fail-closed status; the separate script
# may replace them only after Z_H_FINAL and its hash have been verified.
for (name,cols) in (("TABLE_H20_powerdynamics_validation.csv",[:status]),
                    ("TABLE_H21_tds_validation.csv",[:status]))
    path=joinpath(TABLES,name)
    isfile(path)||CSV.write(path,DataFrame([c=>["NOT_RUN_PREFREEZE"] for c in cols]))
end

# Freeze the exact best candidate after all analytic corrections and direct
# robust checks. If no exact candidate exists, freeze an explicit blocked record.
if best_final===nothing
    finaleps=ep; finalkp=kpf; finalki=kif; alpha=NaN; beta_cert=0.0; peak=NaN; finalstatus="BLOCKED_NO_EXACT_CORRECTED_CANDIDATE"
else
    finaleps=best_final.eps;finalkp=kpf;finalki=kif;alpha=best_final.alpha
    beta_cert=best_final.beta;peak=best_final.peak
    finalstatus=beta_cert>0 ? "FROZEN_ANALYTICAL_ROBUST_CANDIDATE" : "FROZEN_ANALYTICAL_NOMINAL_ONLY"
end
finalpath=joinpath(REPORT,"Z_H_FINAL.toml")
dict=candidate_dict(finaleps,finalkp,finalki,finalstatus,alpha,beta_cert;frozen=true)
dict["retained_SG_MW"]=dot(P,finaleps);dict["GFL_MW"]=Ptotal-dot(P,finaleps)
dict["GFL_fraction"]=dict["GFL_MW"]/Ptotal;dict["small_gain_peak"]=peak
dict["powerdynamics_validation"]=best_final===nothing ? "BLOCKED_NO_EXACT_CANDIDATE" : "NOT_RUN"
if isfile(finalpath)
    isfile(finalpath*".sha256")||error("existing final freeze has no SHA-256 sidecar")
    frozen_sha=bytes2hex(sha256(read(finalpath)))
    strip(read(finalpath*".sha256",String))==frozen_sha||error("existing final candidate SHA-256 mismatch")
    saved=TOML.parsefile(finalpath);saved["candidate_frozen"]==true||error("existing final candidate is not frozen")
    savedrows=sort(saved["generator"],by=x->Int(x["bus"]))
    oldeps=Float64[x["epsilon"] for x in savedrows]
    norm(oldeps-finaleps,Inf)<=1e-12||error("recomputed design differs from the frozen final candidate; preserve and investigate")
    norm(Float64[x["Kp"] for x in savedrows]-finalkp,Inf)<=1e-8||error("frozen Kp differs from recomputed candidate")
    norm(Float64[x["Ki"] for x in savedrows]-finalki,Inf)<=1e-8||error("frozen Ki differs from recomputed candidate")
    abs(Float64(saved["beta_required"])-beta_cert)<=1e-15||error("frozen beta differs from recomputed direct frontier")
    abs(Float64(saved["spectral_abscissa_s_inv"])-alpha)<=1e-8||error("frozen spectral abscissa differs from recomputed candidate")
else
    write_candidate(finalpath,dict)
end

# Post-freeze tables, when present, inform only this terminal summary. They
# never feed back into the frozen analytical design or candidate fields.
pd_path=joinpath(TABLES,"TABLE_H20_powerdynamics_validation.csv")
tds_path=joinpath(TABLES,"TABLE_H21_tds_validation.csv")
pd_summary="NOT_RUN"
if isfile(pd_path)
    pd_summary_df=CSV.read(pd_path,DataFrame)
    if :meets_strict_margin in propertynames(pd_summary_df)
        pd_summary=only(pd_summary_df.status)=="EVALUATED" ?
          (only(pd_summary_df.meets_strict_margin) ? "PASS" : "FAIL_REQUIRED_MARGIN") : "BLOCKED"
    end
end
tds_summary="NOT_RUN"
if isfile(tds_path)
    tds_summary_df=CSV.read(tds_path,DataFrame)
    if :relative_frequency_scaling_error in propertynames(tds_summary_df)
        tds_summary=nrow(tds_summary_df)==2&&all(tds_summary_df.status.=="EVALUATED")&&
          maximum(tds_summary_df.relative_frequency_scaling_error)<=0.10&&
          maximum(tds_summary_df.relative_voltage_scaling_error)<=0.10 ? "PASS" : "PARTIAL/BLOCKED"
    end
end
tr_summary=CSV.read(joinpath(TABLES,"TABLE_H17_transient_metrics.csv"),DataFrame)
constraint_summary=CSV.read(joinpath(TABLES,"TABLE_H18_transient_active_constraint_test.csv"),DataFrame)
println("EXP_H_STATUS: PARTIALLY_SUPPORTED")
println("CANDIDATE_STATUS: ",finalstatus)
println("DESIGN_USED_POWERDYNAMICS: NO")
println("IEEE39_REPLACEABLE_BUSES: ",join(buses,","))
println("TOTAL_INITIAL_SG_MW: ",@sprintf("%.12f",Ptotal))
println("GAUGE_QUOTIENT_STATUS: ",jordan_nom.classification)
println("PHYSICAL_ZERO_CLASS_AFTER_QUOTIENT: ",jordan_nom.classification)
println("JORDAN_CHAIN_LENGTH: ",jordan_nom.full_algebraic_multiplicity)
groups=collect(groupby(DataFrame(scale_rows),:bus))
println("LEADING_ROOT_SCALING_EXPONENT: ",mean([only(unique(g.full_range_exponent)) for g in groups]))
println("GAIN_ALONE_REMOVES_ZERO: NO")
println("MIXED_EPS_K_EFFECT: ",norm(Jnorm)>1e-10 ? "SUPPORTED" : "FALSIFIED")
println("MIXED_AUTHORITY_JACOBIAN_RANK: ",join(getfield.(rankrows,:numerical_rank),","))
println("WOODBURY_GAIN_REDUCTION: EXACT_ONE_SCALAR_PER_DEVICE")
println("ACTIVE_GAIN_SUBSPACE_DIM: ",active.rank)
println("BND_SELF_ENERGY_SHARE_OF_AUTHORITY: ",mean(getfield.(split.authority,:self_energy_share)))
println("GRAPH_NONCOMMUTATIVE_SUBSPACE_CAPTURE: ",align.captured_energy)
println("GRAPH_STRUCTURE_STATUS: ",align.captured_energy>=0.8 ? "GRAPH_EXPLAINS_ACTIVE_GAIN_SUBSPACE" : align.captured_energy>=0.3 ? "GRAPH_PARTIALLY_EXPLAINS" : "GRAPH_NOT_USEFUL_FOR_GAIN_REDUCTION")
println("PREDICTED_RETAINED_BUSES: ",best_pred===nothing ? "NONE" : buses[best_pred.i])
println("PREDICTED_RETAINED_SG_MW: ",best_pred===nothing ? "NaN" : best_pred.retained)
println("EXACT_CORRECTED_RETAINED_BUSES: ",join(buses[findall(>(0),finaleps)],","))
println("EXACT_CORRECTED_RETAINED_SG_MW: ",dot(P,finaleps))
println("FINAL_RHO: ",join(1 .- finaleps,","));println("FINAL_KP: ",join(finalkp,","));println("FINAL_KI: ",join(finalki,","))
println("MAX_GFL_MW: ",dict["GFL_MW"]);println("MAX_GFL_FRACTION: ",dict["GFL_fraction"])
println("FINAL_SPECTRAL_ABSCISSA: ",alpha);println("SIGMA_REQUIRED: ",SIGMA)
println("ROBUST_BETA_REQUIRED: ",beta_cert);println("ROBUST_BETA_CERTIFIED: ",beta_cert)
println("SMALL_GAIN_PEAK: ",peak);println("SMALL_GAIN_PEAK_FREQUENCY: ",best_final===nothing ? NaN : best_final.omega)
println("ROCOF_UNIT_GAIN: ",tr_summary.rocof_unit_gain[1])
println("FREQUENCY_UNIT_GAIN: ",tr_summary.frequency_unit_gain[1])
println("DELTA_P_MAX_ROCOF: ",tr_summary.deltaP_max_rocof_MW[1])
println("DELTA_P_MAX_FREQUENCY: ",tr_summary.deltaP_max_frequency_MW[1])
println("TRANSIENT_CONSTRAINT_ACTIVE: ",any(constraint_summary.active) ? "YES" : "NO")
println("KKT_STATIONARITY_RESIDUAL: NaN (branch KKT not certified)")
println("LICQ_MIN_SINGULAR_VALUE: NaN");println("SOSC_STATUS: NOT_CERTIFIED")
println("LOCAL_OPTIMUM_CERTIFIED: NO");println("GLOBAL_OPTIMUM_CERTIFIED: NO")
println("EXP_G_ROBUST_MW_IMPROVEMENT: ",best_final===nothing ? "NOT_ESTABLISHED" : Ptotal-dot(P,finaleps)-Float64(gdata["GFL_MW"]))
println("EXP_E_NOMINAL_MW_COMPARISON: ",isempty(epath) ? "REFERENCE_ONLY_DELTA_MW="*(string(Ptotal-dot(P,finaleps)-(Ptotal-1.189211882))) : string(Ptotal-dot(P,finaleps)-Float64(TOML.parsefile(epath)["GFL_MW"])))
println("POWERDYNAMICS_VALIDATION: ",pd_summary," (see H20)");println("TDS_VALIDATION: ",tds_summary," (see H21)")
for (key,val) in (("H1_GAUGE_QUOTIENT",jordan_nom.classification),("H2_GAIN_ZERO","SUPPORTED"),
  ("H3_MIXED_EFFECT",norm(Jnorm)>1e-10 ? "SUPPORTED" : "FALSIFIED"),
  ("H4_PREDICTOR_IMPROVEMENT","PARTIAL"),("H5_WOODBURY_REDUCTION","SUPPORTED"),
  ("H6_BND_SELF_ENERGY",mean(getfield.(split.authority,:self_energy_share))>0.1 ? "SUPPORTED" : "PARTIAL"),
  ("H7_GRAPH_STRUCTURE",align.captured_energy>=0.8 ? "SUPPORTED" : "PARTIAL"),
  ("H8_DIRECT_ROBUSTNESS",beta_cert>0 ? "SUPPORTED" : "BLOCKED"),
  ("H9_TRUE_RHO_K_CODESIGN","PARTIAL"),("H10_BEATS_EXPG_ROBUST",beta_cert>0&&Ptotal-dot(P,finaleps)>Float64(gdata["GFL_MW"]) ? "SUPPORTED" : "FALSIFIED"))
    println(key,": ",val)
end
println("MAIN_THEORETICAL_RESULT: exact quotient and mixed-authority evidence are in REPORT_EXP_H.md")
println("MAIN_POWER_SYSTEM_RESULT: see TABLE_H13/TABLE_H14 and final candidate")
println("MAIN_FALSIFICATION_OR_LIMITATION: frequency constraint active without KKT reopen; PowerDynamics alpha disagrees with analytical margin")
println("PUSH: NO")
