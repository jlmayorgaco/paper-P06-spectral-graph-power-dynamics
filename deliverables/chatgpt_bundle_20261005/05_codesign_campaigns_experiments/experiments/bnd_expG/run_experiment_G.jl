using CSV, DataFrames, LinearAlgebra, TOML, SHA, Statistics, Printf
include(joinpath(pwd(),"src","bnd_design_g","BNDDesignG.jl"))
using .BNDDesignG

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const EXP=joinpath(ROOT,"experiments","bnd_expG")
const REPORT=joinpath(ROOT,"reports","experiment_G")
const TABLES=joinpath(REPORT,"tables")
const FIGS=joinpath(REPORT,"figures")
const SIGMA=0.05
const BASE_MVA=100.0
const F0=60.0
const DESIGN_DELTA_P_MW=100.0
const DESIGN_ROCOF_LIMIT=0.5
const DESIGN_FREQUENCY_LIMIT=0.5
const EVENT_BUS=16
const ROBUST_SCALE=1.0 # additive state-matrix uncertainty scale, s^-1
const RESUME_FROZEN="--resume-frozen" in ARGS

mkpath(TABLES); mkpath(FIGS)
const EXISTING_FINAL_PATH=joinpath(REPORT,"Z_G_FINAL.toml")
const EXISTING_FINAL_SHA=isfile(EXISTING_FINAL_PATH) ? BNDDesignG.file_sha256(EXISTING_FINAL_PATH) : ""
isfile(EXISTING_FINAL_PATH) && !RESUME_FROZEN &&
    error("Z_G_FINAL.toml already exists; pass --resume-frozen to verify it without overwriting")
RESUME_FROZEN && isempty(EXISTING_FINAL_SHA) &&
    error("--resume-frozen requires an existing frozen Z_G_FINAL.toml")
RESUME_FROZEN && strip(read(EXISTING_FINAL_PATH*".sha256",String))!=EXISTING_FINAL_SHA &&
    error("existing frozen candidate SHA-256 sidecar does not verify")

function max_nearest_error(a,b)
    if isempty(a)||isempty(b)
        return Inf
    end
    f(x,y)=maximum(minimum(abs(z-w) for w in y) for z in x)
    max(f(a,b),f(b,a))
end

function current_git_commit(root)
    ref=joinpath(root,".git","refs","heads","research","expG-robust-spectral-codesign")
    isfile(ref) && return strip(read(ref,String))
    head=strip(read(joinpath(root,".git","HEAD"),String))
    startswith(head,"ref: ") || return head
    path=joinpath(root,".git",split(replace(head,"ref: "=>""),'/')...)
    isfile(path) ? strip(read(path,String)) : "UNRESOLVED"
end

function source_model_hash(root)
    paths=[
      "src/bnd_design_g/BNDDesignG.jl","src/bnd_design_g/Endpoint.jl",
      "src/bnd_design_g/Robustness.jl","src/bnd_design_g/Transients.jl",
      "src/bnd_design_g/GraphAblation.jl","src/bnd_design_e/CollectiveModel.jl",
      "src/bnd_design_e/ClosureSpectrum.jl","src/bnd_design_e/PhysicalData.jl",
      "src/bnd_design/AnalyticSG.jl","src/bnd_design/AnalyticGFLPLL.jl",
      "experiments/bnd_expG/run_experiment_G.jl"]
    bytes=UInt8[]
    for p in paths
        append!(bytes,read(joinpath(root,p))); append!(bytes,codeunits("\n"))
    end
    return bytes2hex(sha256(bytes)),paths
end

function write_hash(path)
    write(path*".sha256",BNDDesignG.file_sha256(path)*"\n")
end

function closure_audit(net,rho,kp,ki,s)
    m=CollectiveModel.mixed_jacobian(net,rho,kp,ki)
    pc=ClosureSpectrum.port_closure(net,s,rho,kp,ki)
    lhs=logdet(s*I-m.Ared)+logdet(m.Gy)
    rhs=logdet(s*I-m.A)+logdet(net.y_static)+logdet(pc.C)
    return (;point=s,relative_logdet_residual=abs(exp(lhs-rhs)-1),
        closure_sigma_min=minimum(svdvals(pc.C)),
        reduced_schur_residual=norm((s*I-m.Ared)-(s*I-m.A+m.B*(m.Gy\m.C)))/
            max(norm(s*I-m.Ared),eps()),
        condition_static=cond(net.y_static),condition_Gy=cond(m.Gy))
end

function finite_difference_gamma(net,buses,kp,ki;h=1e-4)
    n=length(buses); rho0=ones(n)
    base=spectrum_with_gauge(net,rho0,kp,ki)
    λ0=base.physical[argmin(abs.(base.physical))]
    rows=NamedTuple[]
    for i in 1:n
        estimates=Float64[]
        for step in (h,h/2)
            r=copy(rho0); r[i]=1-step
            sp=spectrum_with_gauge(net,r,kp,ki)
            λ=sp.physical[argmin(abs.(sp.physical .- λ0))]
            push!(estimates,-real((λ-λ0)/step))
        end
        gamma=estimates[end]
        rel=abs(estimates[2]-estimates[1])/max(abs(estimates[2]),1e-8)
        push!(rows,(gamma=gamma,coarse_gamma=estimates[1],fine_gamma=estimates[2],
            step_convergence_relative_error=rel,selected_endpoint_pole=λ0,
            finite_difference_method="full finite-state spectrum after exact gauge quotient"))
    end
    return base,DataFrame(rows)
end

function gain_endpoint_derivatives(net,buses,kp,ki;step_rel=1e-4)
    n=length(buses); rho=ones(n)
    base=spectrum_with_gauge(net,rho,kp,ki)
    λ0=base.physical[argmin(abs.(base.physical))]
    rows=NamedTuple[]
    for group in (:kp,:ki), i in 1:n
        estimates=Float64[]
        x=group===:kp ? kp : ki
        for factor in (step_rel,step_rel/2)
            step=factor*max(abs(x[i]),1.0)
            xp=copy(x); xm=copy(x); xp[i]+=step; xm[i]-=step
            ap=group===:kp ? spectrum_with_gauge(net,rho,xp,ki) : spectrum_with_gauge(net,rho,kp,xp)
            am=group===:kp ? spectrum_with_gauge(net,rho,xm,ki) : spectrum_with_gauge(net,rho,kp,xm)
            lp=ap.physical[argmin(abs.(ap.physical .- λ0))]
            lm=am.physical[argmin(abs.(am.physical .- λ0))]
            push!(estimates,real((lp-lm)/(2step)))
        end
        push!(rows,(parameter=String(group),bus=buses[i],dRe_lambda_dK=estimates[end],
            kappa=-estimates[end],coarse=estimates[1],fine=estimates[end],
            relative_step_error=abs(estimates[2]-estimates[1])/max(abs(estimates[2]),1e-8),
            endpoint_mode_status="DEFECTIVE_GAUGE_QUOTIENT_FD"))
    end
    return base,DataFrame(rows)
end

function self_energy_decomposition(net,buses,kp,ki,gamma)
    n=length(buses); rho=ones(n)
    cl=ClosureSpectrum.port_closure(net,0.0,rho,kp,ki;derivatives=true)
    F=svd(cl.C); u=F.U[:,end]; v=F.V[:,end]
    den=dot(u,cl.Cs*v)
    rows=NamedTuple[]
    if abs(den)<1e-12
        for i in 1:n
            push!(rows,(bus=buses[i],gamma_total=gamma[i],gamma_direct=NaN,
                gamma_collective=NaN,closure_slope=abs(den),closure_sigma_min=F.S[end],
                split_status="DEGENERATE_CLOSURE_JORDAN"))
        end
        paths=NamedTuple[]
    else
        for i in 1:n
            pi=2i-1:2i
            dCeps=-cl.d_rho[i]
            dCdirect=zeros(ComplexF64,size(cl.C))
            dCdirect[pi,pi].=dCeps[pi,pi]
            gd=real(dot(u,dCdirect*v)/den)
            gt=real(dot(u,dCeps*v)/den)
            push!(rows,(bus=buses[i],gamma_total=gamma[i],gamma_direct=gd,
                gamma_collective=gt-gd,closure_slope=abs(den),closure_sigma_min=F.S[end],
                split_status="NEP_FIRST_ORDER_LOCAL_VS_OFFDIAGONAL"))
        end
        paths=NamedTuple[]
        for i in 1:n, j in 1:n
            i==j && continue
            push!(paths,(from_bus=buses[i],to_bus=buses[j],
                port_coupling_norm=norm(cl.G[2i-1:2i,2j-1:2j]),
                normalized_by_row=norm(cl.G[2i-1:2i,2j-1:2j])/
                    max(norm(cl.G[2i-1:2i,:]),eps())))
        end
    end
    return DataFrame(rows),DataFrame(paths),F.S[end],abs(den)
end

function all_candidates_tables(P,gamma,h,bS,bR,buses)
    singles=NamedTuple[]
    for i in eachindex(P)
        x=NaN; cost=Inf; feasible=false
        if gamma[i]>0 && h[i]>0
            x=max(bS/gamma[i],bR/h[i])
            feasible=0<=x<=1
            feasible && (cost=P[i]*x)
        end
        push!(singles,(bus=buses[i],gamma=gamma[i],h_MVA_s=h[i],P_MW=P[i],
            epsilon=x,retained_sg_MW=cost,feasible=feasible))
    end
    pairs=NamedTuple[]
    for i in 1:length(P)-1, j in i+1:length(P)
        D=gamma[i]*h[j]-gamma[j]*h[i]
        ei=NaN; ej=NaN; cost=Inf; feasible=false
        if abs(D)>1e-12
            ei=(bS*h[j]-bR*gamma[j])/D
            ej=(bR*gamma[i]-bS*h[i])/D
            feasible=0<=ei<=1 && 0<=ej<=1
            feasible && (cost=P[i]*ei+P[j]*ej)
        end
        push!(pairs,(bus_i=buses[i],bus_j=buses[j],Dij=D,epsilon_i=ei,epsilon_j=ej,
            retained_sg_MW=cost,feasible=feasible,
            spectral_slack=feasible ? gamma[i]*ei+gamma[j]*ej-bS : NaN,
            rocof_slack=feasible ? h[i]*ei+h[j]*ej-bR : NaN))
    end
    DataFrame(singles),DataFrame(pairs)
end

function coordinate_search(net,buses,P,h,kp,ki,beta,sigma,bR;order=collect(eachindex(P)),
                           gridstep=0.1,max_sweeps=3)
    n=length(P); e=ones(n); evaluations=0
    function metrics(x)
        sp=spectrum_with_gauge(net,1 .- x,kp,ki)
        mr=modal_residue_bound(sp.quotient_matrix;scale=ROBUST_SCALE).value*beta
        evaluations+=1
        inertia=dot(h,x)
        return (;sp,mrob=mr,inertia,retained=dot(P,x),
            feasible=isfinite(mr)&&sp.spectral_abscissa+mr<=-sigma+1e-9&&inertia>=bR-1e-8)
    end
    history=NamedTuple[]
    for sweep in 1:max_sweeps
        changed=false
        for i in order
            old=e[i]; best=old; bestcost=dot(P,e); found=false
            for x in 0.0:gridstep:1.0
                x>=old-1e-10 && continue
                trial=copy(e); trial[i]=x; mt=metrics(trial)
                if mt.feasible && mt.retained<bestcost-1e-8
                    best=x; bestcost=mt.retained; found=true
                end
            end
            if found
                e[i]=best; changed=true
                mt=metrics(e)
                push!(history,(sweep=sweep,bus=buses[i],epsilon=best,
                    retained_sg_MW=mt.retained,spectral_abscissa=mt.sp.spectral_abscissa,
                    m_rob=mt.mrob,inertia_MVA_s=mt.inertia,source="coarse coordinate grid"))
            end
        end
        changed || break
    end
    # Refine each remaining fractional anchor along its coordinate using a
    # deterministic bracketed boundary correction; no general-purpose optimizer.
    for i in eachindex(e)
        (1e-8<e[i]<1-1e-8) || continue
        xs=collect(range(0,e[i],length=11)); feasible_x=Float64[]
        for x in xs
            trial=copy(e); trial[i]=x; metrics(trial).feasible && push!(feasible_x,x)
        end
        isempty(feasible_x) && continue
        hi=minimum(feasible_x); hi>=e[i]-1e-12 && continue
        lo=maximum(x for x in xs if x<hi)
        for _ in 1:10
            mid=(lo+hi)/2; trial=copy(e); trial[i]=mid
            if metrics(trial).feasible; hi=mid else lo=mid end
        end
        e[i]=hi; mt=metrics(e)
        push!(history,(sweep=max_sweeps+1,bus=buses[i],epsilon=hi,
            retained_sg_MW=mt.retained,spectral_abscissa=mt.sp.spectral_abscissa,
            m_rob=mt.mrob,inertia_MVA_s=mt.inertia,source="deterministic coordinate bisection"))
    end
    final=metrics(e)
    return (;eps=e,rho=1 .- e,metrics=final,history=DataFrame(history),
        evaluations,status=final.feasible ? "FEASIBLE_COORDINATE_CANDIDATE" : "NO_FEASIBLE_COORDINATE_CANDIDATE")
end

function robust_feasibility(net,buses,P,h,eps,kp,ki,beta,sigma,bR)
    sp=spectrum_with_gauge(net,1 .- eps,kp,ki)
    mr=modal_residue_bound(sp.quotient_matrix;scale=ROBUST_SCALE).value*beta
    (;spectrum=sp,mrob=mr,retained=dot(P,eps),inertia=dot(h,eps),
      feasible=sp.spectral_abscissa+mr<=-sigma+1e-9 && dot(h,eps)>=bR-1e-8)
end


# ---------- G0: freeze analytical inputs and the uncertainty convention ----------
println("EXP_G: load frozen analytical domain"); flush(stdout)
domain_path=joinpath(ROOT,"experiments","bnd_expE","configs","DESIGN_DOMAIN_FROZEN.toml")
domain_hash_path=domain_path*".sha256"
isfile(domain_path)&&isfile(domain_hash_path)||error("frozen ExpE design-domain inputs are missing")
domain_sha=BNDDesignG.file_sha256(domain_path)
strip(read(domain_hash_path,String))==domain_sha||error("frozen design-domain SHA-256 mismatch")
domain=TOML.parsefile(domain_path)
net=CollectiveModel.frozen_network(ROOT)
discovered=PhysicalData.discover_generators(ROOT)
replaceable=discovered[discovered.replaceable .== true,:]
buses=Int.(replaceable.bus); n=length(buses)
buses==Int.(domain["replaceable_buses"])||error("discovered replaceable buses disagree with frozen domain")
buses==collect(30:39)||error("expected IEEE-39 replacement set changed; refusing silent hardcoding")
P=Float64.(replaceable.SG_dispatch_initial_MW)
totalP=sum(P)
isapprox(totalP,5402.761089978776;atol=1e-6,rtol=0)||error("reconstructed actual SG dispatch integrity check failed")
net.no_load_kcl_inf<=1e-9||error("frozen analytical operating point fails no-load KCL")
h=[net.sg[b].op.parameters.inertia*net.sg[b].op.parameters.rating_mva for b in buses]
kpnom=fill(Float64(domain["Kp_nom"]),n); kinom=fill(Float64(domain["Ki_nom"]),n)
kpmin=Float64[Float64(domain["Kp_factor_min"])*x for x in kpnom]
kpmax=Float64[Float64(domain["Kp_factor_max"])*x for x in kpnom]
kimin=Float64[Float64(domain["Ki_factor_min"])*x for x in kinom]
kimax=Float64[Float64(domain["Ki_factor_max"])*x for x in kinom]
CSV.write(joinpath(TABLES,"TABLE_G00_domain_reconstruction.csv"),replaceable)
CSV.write(joinpath(TABLES,"TABLE_G00_initialized_load_reconstruction.csv"),net.load_audit)

# The all-SG model provides a frozen reference point for a dimensionless
# full-block state-matrix uncertainty radius. No physical percentage is claimed.
sg_ref=spectrum_with_gauge(net,zeros(n),kpnom,kinom)
sg_rob0=robustness_certificate(sg_ref.quotient_matrix,0.0,SIGMA;scale=ROBUST_SCALE)
spare_ref=-sg_ref.spectral_abscissa-SIGMA
spare_ref>0||error("all-SG frozen reference does not meet the requested nominal margin")
beta_modal_cap=spare_ref/max(2sg_rob0.modal_residue_bound,eps())
beta_target=0.5*min(sg_rob0.beta_star,beta_modal_cap)
beta_target>0&&isfinite(beta_target)||error("could not establish a positive normalized robust radius")
mrob_ref=beta_target*sg_rob0.modal_residue_bound
uncertainty=Dict(
 "experiment"=>"BND_EXP_G",
 "uncertainty_model"=>"additive reduced-state matrix full block: A_delta=A+scale*Delta",
 "delta_norm"=>"spectral norm <= beta",
 "state_matrix_scale_s_inv"=>ROBUST_SCALE,
 "beta_units"=>"dimensionless after the declared 1 s^-1 matrix scale",
 "beta_design"=>beta_target,
 "beta_policy"=>"half the smaller of all-SG exact resolvent radius and half-margin modal-residue radius",
 "beta_star_all_sg_exact_resolvent"=>sg_rob0.beta_star,
 "modal_residue_all_sg"=>sg_rob0.modal_residue_bound,
 "modal_margin_all_sg_s_inv"=>spare_ref,
 "sigma_req_s_inv"=>SIGMA,
 "structured_physical_uncertainty"=>"local normalized state-matrix perturbation envelope only; not a physical percentage standard",
 "load_setpoint_variations"=>"excluded; no unequilibrated load-setpoint uncertainty is asserted",
 "rocof_scenario_deltaP_MW"=>DESIGN_DELTA_P_MW,
 "rocof_limit_Hz_s"=>DESIGN_ROCOF_LIMIT,
 "frequency_excursion_limit_Hz"=>DESIGN_FREQUENCY_LIMIT,
 "event_load_bus"=>EVENT_BUS,
 "base_MVA"=>BASE_MVA,
 "base_frequency_Hz"=>F0,
 "frozen_domain_sha256"=>domain_sha)
uncertainty_path=joinpath(EXP,"configs","UNCERTAINTY_SET_FROZEN.toml")
open(uncertainty_path,"w") do io; TOML.print(io,uncertainty); end
write_hash(uncertainty_path)
uncertainty_sha=BNDDesignG.file_sha256(uncertainty_path)
config_static=Dict("sigma_req_s_inv"=>SIGMA,"beta_target"=>beta_target,
 "beta_policy"=>uncertainty["beta_policy"],"deltaP_MW"=>DESIGN_DELTA_P_MW,
 "rocof_limit_Hz_s"=>DESIGN_ROCOF_LIMIT,"frequency_limit_Hz"=>DESIGN_FREQUENCY_LIMIT,
 "event_bus"=>EVENT_BUS,"base_MVA"=>BASE_MVA,"base_frequency_Hz"=>F0,
 "design_uses_PowerDynamics"=>false,"gain_bounds"=>"frozen ExpE independent rectangles")
config_path=joinpath(EXP,"configs","DESIGN_SCENARIO_FROZEN.toml")
open(config_path,"w") do io; TOML.print(io,config_static); end
write_hash(config_path)

# Frozen all-SG spectrum comparison to ExpC is a read-only integrity check.
C0path=joinpath(ROOT,"reports","experiment_C","matrices","C0_Ared.csv")
if isfile(C0path)
    c0df=CSV.read(C0path,DataFrame)
    A_C0=Matrix{Float64}(c0df[:,2:end])
    haus=max_nearest_error(eigvals(sg_ref.model.Ared),eigvals(A_C0))
    baseline_status=haus<=1e-5 ? "PASS" : "FAIL"
else
    haus=NaN; baseline_status="BLOCKED_MISSING_ARCHIVED_C0"
end
baseline_table=DataFrame(check=["all-SG analytical spectrum vs archived ExpC C0"],
    analytic_dimension=[size(sg_ref.model.Ared,1)],archived_dimension=[isfile(C0path) ? size(A_C0,1) : 0],
    spectral_Hausdorff_error=[haus],no_load_KCL_inf=[net.no_load_kcl_inf],
    actual_SG_total_MW=[totalP],status=[baseline_status])
CSV.write(joinpath(TABLES,"TABLE_G00_analytical_model_integrity.csv"),baseline_table)

# ---------- G1: all-GFL zero/gauge structure ----------
println("EXP_G: classify all-GFL endpoint"); flush(stdout)
patA_kp=kpnom.*[0.25,0.5,1.0,2.0,4.0,0.75,1.5,3.0,0.4,2.5]
patA_ki=kinom.*[4.0,3.0,2.0,1.0,0.25,2.5,0.4,1.5,0.75,0.5]
patB_kp=kpnom.*[1.7,0.3,3.2,0.6,2.8,0.45,1.3,3.8,0.9,2.1]
patB_ki=kinom.*[0.4,2.3,0.7,3.5,1.1,4.0,0.3,1.8,2.7,0.55]
patterns=[("nominal",kpnom,kinom),("deterministic_A",patA_kp,patA_ki),
          ("deterministic_B",patB_kp,patB_ki)]
g1rows=NamedTuple[]; svrows=NamedTuple[]; nzrows=NamedTuple[]
for (label,kp,ki) in patterns
    z=zero_structure(net,ones(n),kp,ki)
    leftres=norm(adjoint(z.model.Ared)*z.left_null)/
        max(norm(z.model.Ared)*norm(z.left_null),eps())
    push!(g1rows,(gain_pattern=label,class=z.classification,
        allGFL_dynamic_states=z.model.n_dynamic,raw_finite_poles=length(z.lambda),
        gauge_detected=z.gauge_detected,gauge_pole_real=real(z.gauge_pole),
        gauge_pole_imag=imag(z.gauge_pole),uniform_angle_gauge_residual=z.gauge_residual,
        physical_nearzero_real=real(z.physical_nearzero),
        physical_nearzero_imag=imag(z.physical_nearzero),
        physical_spectral_abscissa=z.spectral_abscissa,
        algebraic_multiplicity_estimate=z.algebraic_multiplicity_estimate,
        geometric_multiplicity_estimate=z.geometric_multiplicity,
        generalized_vector_residual=z.generalized_vector_residual,
        left_null_residual=leftres,smallest_singular_value=last(z.singular_values),
        near_zero_tolerance=z.near_zero_tolerance))
    for (k,v) in enumerate(z.singular_values)
        push!(svrows,(gain_pattern=label,index=k,singular_value=v))
    end
    for (k,v) in enumerate(sort(z.physical,by=abs)[1:min(12,length(z.physical))])
        push!(nzrows,(gain_pattern=label,index=k,lambda_real=real(v),lambda_imag=imag(v),abs_lambda=abs(v)))
    end
end
g1=DataFrame(g1rows); CSV.write(joinpath(TABLES,"TABLE_G01_all_gfl_zero_structure.csv"),g1)
CSV.write(joinpath(TABLES,"TABLE_G01_zero_singular_values.csv"),DataFrame(svrows))
CSV.write(joinpath(TABLES,"TABLE_G01_near_zero_poles.csv"),DataFrame(nzrows))

# ---------- G2: isolated analytical PLL seed, explicitly clipped for initialization ----------
println("EXP_G: derive isolated PLL reference"); flush(stdout)
pll=pll_reference(net,buses,kpmin,kpmax,kimin,kimax)
CSV.write(joinpath(TABLES,"TABLE_G02_isolated_pll_reference.csv"),pll)
kpseed=Float64.(pll.kp_used_seed); kiseed=Float64.(pll.ki_used_seed)

# The isolated triple-root seed lies beyond ExpE gain boxes in this case.
# Check it in a deterministic inertia-efficient fractional-port probe; if that
# probe exposes an unstable finite pole, use frozen nominal gains as allowed.
bR=F0*DESIGN_DELTA_P_MW/(2DESIGN_ROCOF_LIMIT)
seed_anchor=argmax(h./P)
seed_eps_probe=clamp(bR/h[seed_anchor],1e-4,1.0)
epsprobe=zeros(n); epsprobe[seed_anchor]=seed_eps_probe
seed_probe=spectrum_with_gauge(net,1 .- epsprobe,kpseed,kiseed)
seed_probe_pass=seed_probe.spectral_abscissa<=1e-8
working_kp=seed_probe_pass ? kpseed : kpnom
working_ki=seed_probe_pass ? kiseed : kinom
seed_status=seed_probe_pass ? "PROJECTED_ISOLATED_SEED_ACCEPTED" :
    "PROJECTED_SEED_REJECTED_BY_FRACTIONAL_NETWORK_SPECTRUM; NOMINAL_FALLBACK"
println("EXP_G seed gate: ",seed_status," (probe bus ",buses[seed_anchor],
    ", epsilon=",seed_eps_probe,", alpha=",seed_probe.spectral_abscissa,")"); flush(stdout)

# ---------- G3: collective SG authority after exact gauge quotient ----------
println("EXP_G: measure collective spectral authority"); flush(stdout)
allgfl,authority_fd=finite_difference_gamma(net,buses,working_kp,working_ki)
gamma=Float64.(authority_fd.fine_gamma)
direct,paths,closure_svmin,closure_slope=self_energy_decomposition(
    net,buses,working_kp,working_ki,gamma)
authority=DataFrame(bus=buses,P_i_MW=P,gamma=gamma,gamma_per_MW=gamma./P,
    gamma_coarse=authority_fd.coarse_gamma,gamma_fine=authority_fd.fine_gamma,
    derivative_step_convergence_error=authority_fd.step_convergence_relative_error)
authority=leftjoin(authority,direct[:,[:bus,:gamma_direct,:gamma_collective,:split_status]],
    on=:bus)
CSV.write(joinpath(TABLES,"TABLE_G03_sg_spectral_authority.csv"),authority)
CSV.write(joinpath(TABLES,"TABLE_G04_self_energy_pathways.csv"),paths)
println("G3 finite-difference gamma: ",gamma); flush(stdout)

# ---------- G4: gain authority on the gauge-free all-GFL physical branch ----------
println("EXP_G: measure PLL gain authority"); flush(stdout)
gain_base,gains=gain_endpoint_derivatives(net,buses,working_kp,working_ki)
CSV.write(joinpath(TABLES,"TABLE_G05_gain_authority.csv"),gains)
gain_history=DataFrame(iteration=[0],stage=["initialization"],
    seed_status=[seed_status],active_gfl_buses=[join(buses[findall(>(0),1 .- zeros(n))],":")],
    normalized_gain_correction=[0.0],active_bounds=[""],
    note=["gain authority is reported on the all-GFL physical zero branch; final correction is re-evaluated on the exact mixed spectrum"])
# ---------- G5: local normalized uncertainty and small-gain certificate ----------
rob_allgfl=robustness_certificate(allgfl.quotient_matrix,beta_target,SIGMA;scale=ROBUST_SCALE)
CSV.write(joinpath(TABLES,"TABLE_G07_uncertainty_parameters.csv"),
    DataFrame(parameter=["unstructured state matrix"],normalization=["1 s^-1 spectral-norm scale"],
      beta=[beta_target],physical_range=["NOT_CLAIMED"],derivative_basis=["full reduced-state resolvent"]))
# ---------- G6: exact finite single/pair enumeration of the linearized surrogate ----------
mrob_gfl=beta_target*modal_residue_bound(allgfl.quotient_matrix;scale=ROBUST_SCALE).value
a0=allgfl.spectral_abscissa
bS=a0+SIGMA+mrob_ref
surrogate=enumerate_surrogate(P,gamma,h,bS,bR)
single_table,pair_table=all_candidates_tables(P,gamma,h,bS,bR,buses)
CSV.write(joinpath(TABLES,"TABLE_G08_all_single_anchor_candidates.csv"),single_table)
CSV.write(joinpath(TABLES,"TABLE_G09_all_pair_anchor_candidates.csv"),pair_table)
if surrogate.feasible
    survrow=DataFrame(status=[surrogate.status],retained_sg_MW=[surrogate.cost],
      support_buses=[join(buses[surrogate.support],":")],
      epsilon=[join(surrogate.eps,",")],rho=[join(1 .- surrogate.eps,",")],
      bS_s_inv=[bS],mrob_reference_s_inv=[mrob_ref],bR_MVA_s=[bR],
      enumerated_feasible_vertices=[nrow(surrogate.candidates)],
      theorem="two unsaturated anchors at a vertex with two active linear inequalities")
else
    survrow=DataFrame(status=[surrogate.status],retained_sg_MW=[Inf],support_buses=[""],
      epsilon=[""],rho=[""],bS_s_inv=[bS],mrob_reference_s_inv=[mrob_ref],
      bR_MVA_s=[bR],enumerated_feasible_vertices=[0],
      theorem="two unsaturated anchors at a vertex with two active linear inequalities")
end
CSV.write(joinpath(TABLES,"TABLE_G10_surrogate_global_solution.csv"),survrow)
surrogate_exact=if surrogate.feasible
    robust_feasibility(net,buses,P,h,surrogate.eps,working_kp,working_ki,
        beta_target,SIGMA,bR)
else
    nothing
end
surrogate_exact_alpha=surrogate_exact===nothing ? NaN : surrogate_exact.spectrum.spectral_abscissa
surrogate_exact_mrob=surrogate_exact===nothing ? NaN : surrogate_exact.mrob

# ---------- G7: deterministic explicit full-closure coordinate iteration ----------
println("EXP_G: run deterministic exact-spectrum coordinate design"); flush(stdout)
orders=[collect(1:n),reverse(collect(1:n)),sortperm(gamma./P;rev=true)]
order_names=["bus_order","reverse_bus_order","authority_per_MW_order"]
designs=NamedTuple[]; order_rows=NamedTuple[]
for q in eachindex(orders)
    println("  coordinate order ",order_names[q]); flush(stdout)
    des=coordinate_search(net,buses,P,h,working_kp,working_ki,beta_target,SIGMA,bR;
        order=orders[q],gridstep=0.1,max_sweeps=3)
    push!(designs,des)
    CSV.write(joinpath(TABLES,"TABLE_G06_gain_correction_history_order$(q).csv"),des.history)
    push!(order_rows,(order=order_names[q],status=des.status,retained_sg_MW=des.metrics.retained,
        spectral_abscissa=des.metrics.sp.spectral_abscissa,m_rob=des.metrics.mrob,
        robustified_abscissa=des.metrics.sp.spectral_abscissa+des.metrics.mrob,
        inertia_MVA_s=des.metrics.inertia,coordinate_evaluations=des.evaluations,
        epsilon=join(des.eps,","),retained_buses=join(buses[findall(>(1e-10),des.eps)],":")))
end
coord_table=DataFrame(order_rows)
CSV.write(joinpath(TABLES,"TABLE_G06_coordinate_searches.csv"),coord_table)
valid=findall(d->d.metrics.feasible,designs)
isempty(valid)&&error("all-SG fallback was feasible but every deterministic coordinate path failed")
bestindex=valid[argmin([designs[i].metrics.retained for i in valid])]
best=designs[bestindex]
epsfinal=copy(best.eps); rhofinal=1 .- epsfinal
kpfinal=copy(working_kp); kifinal=copy(working_ki)
final_exact=robust_feasibility(net,buses,P,h,epsfinal,kpfinal,kifinal,beta_target,SIGMA,bR)
final_exact.feasible||error("coordinate candidate lost its exact finite-spectrum robust feasibility")
println("EXP_G coordinate candidate: retained MW=",final_exact.retained,
    " alpha=",final_exact.spectrum.spectral_abscissa," modal mrob=",final_exact.mrob,
    " support=",buses[findall(>(1e-10),epsfinal)]); flush(stdout)

# Declared design-parameter sweeps. These are independent deterministic runs;
# they never overwrite the selected main candidate.
println("EXP_G: sweep beta and RoCoF scenario parameters"); flush(stdout)
sensrows=NamedTuple[]
for factor in (0.1,0.25,0.5,1.0,2.0)
    des=coordinate_search(net,buses,P,h,working_kp,working_ki,beta_target*factor,SIGMA,bR;
        order=orders[1],gridstep=0.2,max_sweeps=2)
    push!(sensrows,(sweep="beta",parameter=factor,parameter_value=beta_target*factor,
      status=des.status,retained_SG_MW=des.metrics.retained,max_GFL_MW=totalP-des.metrics.retained,
      spectral_abscissa=des.metrics.sp.spectral_abscissa,m_rob=des.metrics.mrob,
      robustified_abscissa=des.metrics.sp.spectral_abscissa+des.metrics.mrob,
      retained_buses=join(buses[findall(>(1e-9),des.eps)],":"),epsilon=join(des.eps,",")))
end
for Rmax in (0.25,0.5,1.0,2.0)
    req=F0*DESIGN_DELTA_P_MW/(2Rmax)
    des=coordinate_search(net,buses,P,h,working_kp,working_ki,beta_target,SIGMA,req;
        order=orders[1],gridstep=0.2,max_sweeps=2)
    push!(sensrows,(sweep="rocof",parameter=Rmax,parameter_value=Rmax,
      status=des.status,retained_SG_MW=des.metrics.retained,max_GFL_MW=totalP-des.metrics.retained,
      spectral_abscissa=des.metrics.sp.spectral_abscissa,m_rob=des.metrics.mrob,
      robustified_abscissa=des.metrics.sp.spectral_abscissa+des.metrics.mrob,
      retained_buses=join(buses[findall(>(1e-9),des.eps)],":"),epsilon=join(des.eps,",")))
end
design_sensitivity=DataFrame(sensrows)
CSV.write(joinpath(TABLES,"TABLE_G12_design_sensitivity_beta_rocof.csv"),design_sensitivity)

# The gain active-set correction is explicitly evaluated. If the candidate is
# already feasible, the minimum-energy correction is exactly zero.
spectral_residual=max(0.0,final_exact.spectrum.spectral_abscissa+
    SIGMA+final_exact.mrob)
push!(gain_history,(iteration=1,stage="exact mixed closure",
    seed_status=seed_status,active_gfl_buses=join(buses[findall(>(1e-10),rhofinal)],":"),
    normalized_gain_correction=0.0,active_bounds="",
    note=spectral_residual<=1e-8 ? "no gain correction required at the feasible coordinate candidate" :
        "gain correction requires a finite-difference active-pole Jacobian; design marked partial"))
CSV.write(joinpath(TABLES,"TABLE_G06_gain_correction_history.csv"),gain_history)

# ---------- G8: exact collective-closure corrector ----------
println("EXP_G: exact low-rank closure and algebraic audits"); flush(stdout)
closure_test=0.7+1.3im
cl_audit=closure_audit(net,rhofinal,kpfinal,kifinal,closure_test)
final_closure=ClosureSpectrum.port_closure(net,final_exact.spectrum.critical,
    rhofinal,kpfinal,kifinal)
state_eig_res=norm(final_exact.spectrum.model.Ared*final_exact.spectrum.critical_right-
    final_exact.spectrum.critical*final_exact.spectrum.critical_right)/
    max(norm(final_exact.spectrum.model.Ared)*norm(final_exact.spectrum.critical_right),eps())
lowrank=if surrogate.feasible && length(surrogate.support)<=2
    low_rank_boundary(net,surrogate.eps,kpfinal,kifinal,SIGMA+final_exact.mrob;
        omega=imag(final_exact.spectrum.critical))
else
    (;status="SKIPPED_SURROGATE_SUPPORT_NOT_ONE_OR_TWO",roots=Float64[],
      determinant_residual=NaN,determinant_lemma_residual=NaN)
end
closuretab=DataFrame(check=["full-state/collective determinant identity at a regular complex point",
                              "critical exact closure residual",
                              "critical full-state eigenvector residual",
                              "low-rank determinant lemma"],
    residual=[cl_audit.relative_logdet_residual,minimum(svdvals(final_closure.C)),
              state_eig_res,hasproperty(lowrank,:determinant_lemma_residual) ?
              lowrank.determinant_lemma_residual : NaN],
    tolerance=[1e-9,1e-7,1e-8,1e-8],
    status=[cl_audit.relative_logdet_residual<=1e-9 ? "PASS" : "FAIL",
      minimum(svdvals(final_closure.C))<=1e-7 ? "PASS" : "FAIL",
      state_eig_res<=1e-8 ? "PASS" : "FAIL",
      hasproperty(lowrank,:determinant_lemma_residual)&&lowrank.determinant_lemma_residual<=1e-8 ?
        "PASS" : lowrank.status])
CSV.write(joinpath(TABLES,"TABLE_G00_closure_identity_audit.csv"),closuretab)
lowranktab=DataFrame(status=[lowrank.status],
    surrogate_anchors=[hasproperty(lowrank,:anchors) ? join(buses[lowrank.anchors],":") : ""],
    predicted_boundary_s_inv=[-(SIGMA+final_exact.mrob)],
    radial_roots=[hasproperty(lowrank,:roots) ? join(lowrank.roots,",") : ""],
    determinant_lemma_residual=[hasproperty(lowrank,:determinant_lemma_residual) ?
        lowrank.determinant_lemma_residual : NaN],
    exact_corrected_candidate_spectral_abscissa=[final_exact.spectrum.spectral_abscissa],
    correction_status=["surrogate low-rank root is a predictor only; final candidate uses complete state spectrum"])
CSV.write(joinpath(TABLES,"TABLE_G08_low_rank_boundary_correction.csv"),lowranktab)

# ---------- G9: full finite-spectrum and deterministic local audit ----------
println("EXP_G: deterministic local perturbation audit"); flush(stdout)
radii=[1e-4,3e-4,1e-3,3e-3,1e-2,3e-2]
auditrows=NamedTuple[]; global tested=0; global improving=0
for r in radii
    for i in eachindex(epsfinal)
        epsfinal[i]>0 || continue
        trial=copy(epsfinal); trial[i]=max(0.0,trial[i]-min(r,epsfinal[i]))
        mt=robust_feasibility(net,buses,P,h,trial,kpfinal,kifinal,beta_target,SIGMA,bR)
        global tested+=1
        improve=mt.feasible && mt.retained<final_exact.retained-1e-8
        improve && (global improving+=1)
        push!(auditrows,(radius=r,direction="reduce_retained_SG",bus_from=buses[i],
          bus_to=missing,feasible=mt.feasible,retained_sg_MW=mt.retained,
          robustified_abscissa=mt.spectrum.spectral_abscissa+mt.mrob,
          objective_improves=improve,status=improve ? "COUNTEREXAMPLE" : "NO_IMPROVEMENT"))
    end
    active=findall(>(1e-8),rhofinal)
    for i in active, group in (:kp,:ki), sign in (-1.0,1.0)
        kp2=copy(kpfinal); ki2=copy(kifinal)
        if group===:kp
            kp2[i]=clamp(kp2[i]+sign*r*kpnom[i],kpmin[i],kpmax[i])
        else
            ki2[i]=clamp(ki2[i]+sign*r*kinom[i],kimin[i],kimax[i])
        end
        (kp2==kpfinal && ki2==kifinal) && continue
        mt=robust_feasibility(net,buses,P,h,epsfinal,kp2,ki2,beta_target,SIGMA,bR)
        global tested+=1
        push!(auditrows,(radius=r,direction="gain_$(group)",bus_from=buses[i],
          bus_to=missing,feasible=mt.feasible,retained_sg_MW=mt.retained,
          robustified_abscissa=mt.spectrum.spectral_abscissa+mt.mrob,
          objective_improves=false,status=mt.feasible ? "FEASIBLE_GAIN_PERTURBATION" : "INFEASIBLE"))
    end
end
# SG relocation tests at each preregistered radius, restricted to deterministic
# source/destination pairs with potential primary-objective reduction.
for r in radii
    pairs=Tuple{Int,Int}[]
    sources=sort(findall(>(1e-8),epsfinal),by=i->-P[i])[1:min(3,count(>(1e-8),epsfinal))]
    targets=sort(findall(<(1-1e-8),epsfinal),by=i->P[i])[1:min(3,count(<(1-1e-8),epsfinal))]
    for i in sources,j in targets
        i==j && continue
        P[i]>P[j] && epsfinal[i]>=r && epsfinal[j]<=1-r && push!(pairs,(i,j))
    end
    for (i,j) in unique(pairs)
        trial=copy(epsfinal); trial[i]-=r; trial[j]+=r
        mt=robust_feasibility(net,buses,P,h,trial,kpfinal,kifinal,beta_target,SIGMA,bR)
        global tested+=1
        improve=mt.feasible && mt.retained<final_exact.retained-1e-8
        improve&&(global improving+=1)
        push!(auditrows,(radius=r,direction="SG_relocation",bus_from=buses[i],
          bus_to=buses[j],feasible=mt.feasible,retained_sg_MW=mt.retained,
          robustified_abscissa=mt.spectrum.spectral_abscissa+mt.mrob,
          objective_improves=improve,status=improve ? "COUNTEREXAMPLE" :
            (mt.feasible ? "FEASIBLE_NO_OBJECTIVE_IMPROVEMENT" : "INFEASIBLE")))
    end
end
auditdf=DataFrame(auditrows)
CSV.write(joinpath(TABLES,"TABLE_G09_deterministic_local_audit.csv"),auditdf)
local_audit_status=improving==0 ? "PASS_ON_TESTED_RETENTION_AND_RELOCATION_DIRECTIONS" : "FAIL"
# This audit is deliberately not upgraded to strict KKT/SOSC certification.
kkt_status="NOT_CERTIFIED_FULL_CLOSURE_KKT_SOSC"

# ---------- G10: transients from the exact same reduced state model ----------
println("EXP_G: derive modal and time-domain disturbance response"); flush(stdout)
transient=try
    transient_analysis(net,final_exact.spectrum.model,epsfinal,buses;event_bus=EVENT_BUS,
      system_base_mva=BASE_MVA,f0=F0,rocof_limit=DESIGN_ROCOF_LIMIT,
      frequency_limit=DESIGN_FREQUENCY_LIMIT)
catch err
    (;status="BLOCKED_TRANSIENT_EXCEPTION",error=sprint(showerror,err))
end
if hasproperty(transient,:modal_residues)
    CSV.write(joinpath(TABLES,"TABLE_G11_modal_transient_residues.csv"),transient.modal_residues)
    CSV.write(joinpath(TABLES,"TABLE_G12_unit_disturbance_response.csv"),
      DataFrame(time_s=transient.times,frequency_Hz_per_MW=transient.frequency,
        rocof_Hz_s_per_MW=transient.rocof))
    caprow=(status=transient.status,rocof_unit_gain_Hz_s_per_MW=transient.rocof_unit_gain,
      frequency_unit_gain_Hz_per_MW=transient.frequency_unit_gain,
      unit_peak_rocof_Hz_s_per_MW=transient.unit_rocof_peak,
      unit_peak_frequency_Hz_per_MW=transient.unit_frequency_peak,
      frequency_nadir_per_MW=transient.frequency_nadir,
      frequency_nadir_time_s=transient.frequency_nadir_time,
      H_total_MVA_s=transient.H_total_MVA_s,
      deltaP_max_initial_rocof_MW=transient.deltaP_max_initial_rocof_MW,
      deltaP_max_peak_rocof_MW=transient.deltaP_max_peak_rocof_MW,
      deltaP_max_frequency_MW=transient.deltaP_max_frequency_MW,
      deltaP_max_total_MW=transient.deltaP_max_total_MW,
      direct_modal_rocof_relative_error=transient.direct_modal_rocof_relative_error,
      direct_modal_frequency_relative_error=transient.direct_modal_frequency_relative_error,
      dominant_mode=transient.dominant_mode,
      dominant_pair_tail_bound=transient.dominant_pair_tail_bound,
      dominant_pair_justified=transient.dominant_pair_justified,
      rocof_limit_Hz_s=DESIGN_ROCOF_LIMIT,frequency_limit_Hz=DESIGN_FREQUENCY_LIMIT)
else
    CSV.write(joinpath(TABLES,"TABLE_G12_unit_disturbance_response.csv"),
      DataFrame(status=[transient.status]))
    CSV.write(joinpath(TABLES,"TABLE_G11_modal_transient_residues.csv"),
      DataFrame(status=[transient.status],error=[transient.error]))
    caprow=(status=transient.status,rocof_unit_gain_Hz_s_per_MW=NaN,
      frequency_unit_gain_Hz_per_MW=NaN,unit_peak_rocof_Hz_s_per_MW=NaN,
      unit_peak_frequency_Hz_per_MW=NaN,frequency_nadir_per_MW=NaN,
      frequency_nadir_time_s=NaN,H_total_MVA_s=dot(h,epsfinal),
      deltaP_max_initial_rocof_MW=NaN,deltaP_max_peak_rocof_MW=NaN,
      deltaP_max_frequency_MW=NaN,deltaP_max_total_MW=NaN,
      direct_modal_rocof_relative_error=NaN,direct_modal_frequency_relative_error=NaN,
      dominant_mode=missing,dominant_pair_tail_bound=NaN,dominant_pair_justified=false,
      rocof_limit_Hz_s=DESIGN_ROCOF_LIMIT,frequency_limit_Hz=DESIGN_FREQUENCY_LIMIT)
end
capacity=DataFrame([caprow])
CSV.write(joinpath(TABLES,"TABLE_G12_transient_capacity.csv"),capacity)
CSV.write(joinpath(TABLES,"TABLE_G11_direct_modal_validation.csv"),
    hasproperty(transient,:checked_times) ?
      DataFrame(time_s=transient.checked_times,rocof_direct=transient.checked_direct_rocof,
        rocof_modal=transient.checked_modal_rocof,frequency_direct=transient.checked_direct_frequency,
        frequency_modal=transient.checked_modal_frequency) :
      DataFrame(status=[transient.status]))

# ---------- G11 optional graph-compressed ablation; it cannot affect the design ----------
gf=graph_features(net,buses,P,gamma)
graphdf=DataFrame(bus=buses)
for j in 1:length(gf.names); graphdf[!,Symbol(gf.names[j])]=gf.raw_features[:,j]; end
CSV.write(joinpath(TABLES,"TABLE_G11_graph_features.csv"),graphdf)
ablation=DataFrame(controller_class=["uniform gains","single-L graph feature class",
    "noncommutative LG/LB feature class","fully free local gains"],
  parameter_count=[2,2,2gf.rank,2n],retained_sg_MW=[missing,missing,missing,final_exact.retained],
  robust_margin=[missing,missing,missing,-(final_exact.spectrum.spectral_abscissa+final_exact.mrob)],
  deltaP_max_MW=[missing,missing,missing,caprow.deltaP_max_total_MW],
  gain_effort=[missing,missing,missing,norm(vcat((kpfinal.-kpnom)./kpnom,(kifinal.-kinom)./kinom))],
  status=["OPTIONAL_ABLATION_NOT_RUN","OPTIONAL_ABLATION_NOT_RUN",
          "OPTIONAL_ABLATION_NOT_RUN","MAIN_DESIGN"],
  graph_feature_rank=fill(gf.rank,4),simultaneous_diagonalization_claim=fill(false,4))
CSV.write(joinpath(TABLES,"TABLE_G11_graph_controller_ablation.csv"),ablation)

# ---------- G13 pre-blind freeze: no ExpE result has influenced the design ----------
retained=dot(P,epsfinal); gflmw=totalP-retained
rob_sg=robustness_certificate(sg_ref.quotient_matrix,beta_target,SIGMA;scale=ROBUST_SCALE)
rob_candidate=robustness_certificate(final_exact.spectrum.quotient_matrix,beta_target,SIGMA;scale=ROBUST_SCALE)
robrows=NamedTuple[]
for (label,r) in (("all_SG_reference",rob_sg),("all_GFL_nominal",rob_allgfl),
                  ("ExpG_exact_coordinate_candidate",rob_candidate))
    push!(robrows,(case=label,beta=r.beta,beta_star_exact_small_gain=r.beta_star,
      exact_small_gain_peak=r.resolvent_peak,omega_peak=r.omega_peak,
      exact_small_gain_pass=r.small_gain_pass,small_gain_margin=r.small_gain_margin,
      modal_residue_bound=r.modal_residue_bound,m_rob=r.m_rob,
      nominal_alpha=r.nominal_abscissa,robustified_alpha=r.nominal_abscissa+r.m_rob,
      sigma_required=SIGMA,status=r.status))
end
CSV.write(joinpath(TABLES,"TABLE_G07_robustness_certificate.csv"),DataFrame(robrows))
CSV.write(joinpath(TABLES,"TABLE_G07_resolvent_sigma_min.csv"),
  DataFrame(omega=rob_candidate.frequency,sigma_min=rob_candidate.sigma_min))
betasweep=sort(unique([0.0,beta_target/4,beta_target/2,beta_target,
    min(2beta_target,rob_candidate.beta_star/2),rob_candidate.beta_star]))
CSV.write(joinpath(TABLES,"TABLE_G07_beta_sweep.csv"),DataFrame(beta=betasweep,
  small_gain_margin=1 .- betasweep.*rob_candidate.resolvent_peak,
  certificate_pass=(betasweep.*rob_candidate.resolvent_peak .< 1)))
modelhash,sourcepaths=source_model_hash(ROOT)
inputdir=joinpath(ROOT,"reports","experiment_D","inputs")
input_hashes=Dict(f=>BNDDesignG.file_sha256(joinpath(inputdir,f)) for f in
    ("bus.csv","machine.csv","load.csv","branch.csv","avr.csv","gov.csv"))
generators=[Dict("bus"=>b,"epsilon"=>epsfinal[i],"rho"=>rhofinal[i],
    "Kp"=>kpfinal[i],"Ki"=>kifinal[i],"P_SG_initial_MW"=>P[i],
    "H_MVA_s"=>h[i],"retained_SG_MW"=>epsfinal[i]*P[i],
    "converted_GFL_MW"=>rhofinal[i]*P[i],
    "Kp_min"=>kpmin[i],"Kp_max"=>kpmax[i],"Ki_min"=>kimin[i],"Ki_max"=>kimax[i])
    for (i,b) in enumerate(buses)]
preblind=Dict("experiment"=>"EXP_G","status"=>"FROZEN_PREBLIND",
    "candidate_frozen"=>true,"design_used_PowerDynamics"=>false,
    "design_method"=>"finite single/pair authority surrogate plus deterministic exact-spectrum coordinate correction",
    "surrogate_status"=>surrogate.status,"full_closure_local_optimality"=>kkt_status,
    "surrogate_global_support_buses"=>surrogate.feasible ? buses[surrogate.support] : Int[],
    "exact_coordinate_order"=>order_names[bestindex],
    "exact_coordinate_candidate_status"=>best.status,
    "retained_SG_MW"=>retained,"GFL_MW"=>gflmw,"total_initial_SG_MW"=>totalP,
    "spectral_abscissa_s_inv"=>final_exact.spectrum.spectral_abscissa,
    "sigma_required_s_inv"=>SIGMA,"beta_design"=>beta_target,
    "beta_star_exact_small_gain"=>rob_candidate.beta_star,
    "exact_small_gain_pass"=>rob_candidate.small_gain_pass,
    "exact_small_gain_margin"=>rob_candidate.small_gain_margin,
    "modal_robust_margin_s_inv"=>final_exact.mrob,
    "robustified_spectral_abscissa_s_inv"=>final_exact.spectrum.spectral_abscissa+final_exact.mrob,
    "critical_pole_real"=>real(final_exact.spectrum.critical),
    "critical_pole_imag"=>imag(final_exact.spectrum.critical),
    "critical_pole_condition_number"=>final_exact.spectrum.critical_nonnormality,
    "rocof_inertia_requirement_MVA_s"=>bR,"retained_inertia_MVA_s"=>final_exact.inertia,
    "uncertainty_config_sha256"=>uncertainty_sha,"analytical_model_sha256"=>modelhash,
    "frozen_input_sha256"=>input_hashes,"git_commit"=>current_git_commit(ROOT),
    "generator"=>generators)
preblind_path=joinpath(REPORT,"Z_G_CANDIDATE_PREBLIND.toml")
if RESUME_FROZEN
    isfile(preblind_path) && isfile(preblind_path*".sha256") || error("resume requires the frozen pre-blind candidate and SHA sidecar")
    frozen_preblind=TOML.parsefile(preblind_path)
    frozen_modelhash=String(frozen_preblind["analytical_model_sha256"])
    candidate_check=copy(preblind)
    candidate_check["analytical_model_sha256"]=frozen_modelhash
    isequal(frozen_preblind,candidate_check) || error("recomputed pre-blind candidate differs from the frozen candidate")
    preblind=frozen_preblind
    modelhash=frozen_modelhash
    preblind_hash=BNDDesignG.file_sha256(preblind_path)
    strip(read(preblind_path*".sha256",String))==preblind_hash || error("pre-blind candidate SHA-256 sidecar does not verify")
else
    open(preblind_path,"w") do io; TOML.print(io,preblind); end
    write_hash(preblind_path)
    preblind_hash=BNDDesignG.file_sha256(preblind_path)
end
println("EXP_G_PREBLIND_FREEZE_SHA256: ",preblind_hash); flush(stdout)

# ---------- G12: blinded comparison, strictly after candidate and tables are frozen ----------
println("EXP_G: read provisional ExpE benchmark after pre-blind freeze"); flush(stdout)
eprov=joinpath(ROOT,"reports","experiment_E","PROVISIONAL_ANALYTIC_RESULT.json")
function json_numeric(raw,keys)
    for key in keys
        m=match(Regex("\\\""*key*"\\\"\\s*:\\s*([-+]?(?:[0-9]+\\.?[0-9]*|\\.[0-9]+)(?:[eE][-+]?[0-9]+)?)","i"),raw)
        m===nothing || return parse(Float64,m.captures[1])
    end
    NaN
end
function json_text(raw,keys)
    for key in keys
        m=match(Regex("\\\""*key*"\\\"\\s*:\\s*\\\"([^\\\"]+)\\\"","i"),raw)
        m===nothing || return m.captures[1]
    end
    ""
end
if isfile(eprov)
    eraw=read(eprov,String)
    e_mw=json_numeric(eraw,["retained_SG_MW","retained_sg_MW","retained_MW","retained_mw"])
    e_gfl=json_numeric(eraw,["GFL_MW","gfl_MW","converted_GFL_MW","gfl_mw"])
    e_alpha=json_numeric(eraw,["spectral_abscissa","alpha","alpha_s_inv"])
    ebus=json_numeric(eraw,["anchor_bus","bus"])
    ehash=BNDDesignG.file_sha256(eprov)
    exp_e_read="YES"
else
    e_mw=NaN; e_gfl=NaN; e_alpha=NaN; ebus=NaN; ehash=""; exp_e_read="NO"
end
bus38match=isfinite(ebus) ? (38 in buses[findall(>(1e-8),epsfinal)])==(round(Int,ebus)==38) : missing
comparison=DataFrame(metric=["ExpG retained SG MW","ExpG converted GFL MW","ExpG spectral abscissa",
    "ExpE provisional retained SG MW","ExpE provisional GFL MW","ExpE provisional alpha",
    "ExpG retained-MW difference vs ExpE","ExpE anchor bus","ExpG support includes bus 38",
    "ExpE provisional source hash","candidate preblind hash"],
 value=[string(retained),string(gflmw),string(final_exact.spectrum.spectral_abscissa),
    string(e_mw),string(e_gfl),string(e_alpha),string(isfinite(e_mw) ? retained-e_mw : NaN),
    isfinite(ebus) ? string(round(Int,ebus)) : "NOT_PARSED",
    ismissing(bus38match) ? "NOT_PARSED" : string(bus38match),ehash,preblind_hash])
CSV.write(joinpath(TABLES,"TABLE_G14_blinded_ExpE_comparison.csv"),comparison)


# ---------- G13: immutable final analytical candidate ----------
finalcand=copy(preblind)
finalcand["status"]="FROZEN"
finalcand["expE_blinded_comparison_status"]=exp_e_read=="YES" ? "READ_AFTER_PREBLIND_FREEZE" : "NOT_CHECKED"
finalcand["expE_provisional_retained_SG_MW"]=e_mw
finalcand["expE_provisional_GFL_MW"]=e_gfl
finalcand["expE_provisional_alpha"]=e_alpha
finalcand["expE_provisional_anchor_bus"]=isfinite(ebus) ? round(Int,ebus) : -1
finalcand["expE_source_sha256"]=ehash
finalcand["preblind_candidate_sha256"]=preblind_hash
finalcand["powerdynamics_validation"]="NOT_RUN"
finalcand["tds_validation"]="NOT_RUN"
final_path=joinpath(REPORT,"Z_G_FINAL.toml")
if RESUME_FROZEN
    isequal(TOML.parsefile(final_path),finalcand) || error("recomputed candidate differs from frozen Z_G_FINAL.toml")
    BNDDesignG.file_sha256(final_path)==EXISTING_FINAL_SHA || error("frozen candidate SHA changed during resume")
else
    open(final_path,"w") do io; TOML.print(io,finalcand); end
    write_hash(final_path)
end
final_hash=BNDDesignG.file_sha256(final_path)

gen_table=DataFrame(bus=buses,P_i_MW=P,H_i_MVA_s=h,
  epsilon_retained=epsfinal,rho_converted=rhofinal,Kp=kpfinal,Ki=kifinal,
  Kp_min=kpmin,Kp_max=kpmax,Ki_min=kimin,Ki_max=kimax,
  retained_SG_MW=P.*epsfinal,converted_GFL_MW=P.*rhofinal,
  H_times_S_MVA_s=h.*epsfinal,gain_active=rhofinal.>0)
CSV.write(joinpath(TABLES,"TABLE_G13_final_per_generator_design.csv"),gen_table)
CSV.write(joinpath(TABLES,"TABLE_G13_rho_Kp_Ki_by_bus.csv"),gen_table)
eigenmap=DataFrame(case=["all_SG","all_GFL","surrogate_candidate","ExpG_final"],
  alpha=[sg_ref.spectral_abscissa,allgfl.spectral_abscissa,
    surrogate_exact===nothing ? NaN : surrogate_exact.spectrum.spectral_abscissa,
    final_exact.spectrum.spectral_abscissa],
  robustified_alpha=[sg_ref.spectral_abscissa+mrob_ref,
    allgfl.spectral_abscissa+mrob_gfl,
    surrogate_exact===nothing ? NaN : surrogate_exact.spectrum.spectral_abscissa+surrogate_exact.mrob,
    final_exact.spectrum.spectral_abscissa+final_exact.mrob],
  sigma_required=fill(SIGMA,4),critical_real=[NaN,real(allgfl.critical),
    surrogate_exact===nothing ? NaN : real(surrogate_exact.spectrum.critical),
    real(final_exact.spectrum.critical)],
  critical_imag=[NaN,imag(allgfl.critical),
    surrogate_exact===nothing ? NaN : imag(surrogate_exact.spectrum.critical),
    imag(final_exact.spectrum.critical)])
CSV.write(joinpath(TABLES,"TABLE_G15_eigenvalue_comparison.csv"),eigenmap)

# Acceptance/falsification ledger. A failed or unavailable gate is retained.
g1status=all(g1.gauge_detected) ? "PASS" : "FAIL"
derivstatus=maximum(authority.derivative_step_convergence_error)<=1e-4 ? "PASS" : "PARTIAL"
robuststatus=rob_candidate.small_gain_pass &&
    final_exact.spectrum.spectral_abscissa+final_exact.mrob<=-SIGMA+1e-9 ? "PASS" : "FAIL"
transientstatus=hasproperty(transient,:direct_modal_rocof_relative_error) &&
    transient.direct_modal_rocof_relative_error<=1e-6 &&
    transient.direct_modal_frequency_relative_error<=1e-6 ? "PASS" : "PARTIAL"
gate_rows=DataFrame(gate=["A analytical endpoint reproduction","B all-GFL zero/gauge classification",
  "C retained-SG derivative validation","D finite single/pair surrogate enumeration",
  "E exact full finite-spectrum strict feasibility","F robust certificate",
  "G transient/direct matrix exponential agreement","H pre-blind freeze before ExpE read",
  "I independent PowerDynamics validation","J nonlinear TDS validation"],
  status=[baseline_status,g1status,derivstatus,surrogate.status,
    final_exact.feasible ? "PASS" : "FAIL",robuststatus,transientstatus,
    "PASS","NOT_RUN","NOT_RUN"],
  evidence=["TABLE_G00_analytical_model_integrity.csv",
    "TABLE_G01_all_gfl_zero_structure.csv and Jordan scaling probe",
    "TABLE_G03_sg_spectral_authority.csv and TABLE_G05_gain_authority.csv",
    "TABLE_G08_all_single_anchor_candidates.csv / TABLE_G09_all_pair_anchor_candidates.csv",
    "TABLE_G13_final_per_generator_design.csv","TABLE_G07_robustness_certificate.csv",
    "TABLE_G11_direct_modal_validation.csv","pre-blind candidate and SHA-256",
    "independent validation pending","independent validation pending"])
CSV.write(joinpath(TABLES,"TABLE_G00_acceptance_gate_ledger.csv"),gate_rows)
gstatus="PARTIAL" # PD/TDS and a strict KKT/SOSC certificate remain pending.

sha_rows=NamedTuple[]
for p in (uncertainty_path,config_path,preblind_path,final_path)
    push!(sha_rows,(path=relpath(p,ROOT),sha256=BNDDesignG.file_sha256(p)))
end
for (f,hsh) in input_hashes
    push!(sha_rows,(path=joinpath("reports","experiment_D","inputs",f),sha256=hsh))
end
for (p,hsh) in (("experiments/bnd_expE/configs/DESIGN_DOMAIN_FROZEN.toml",domain_sha),
                ("reports/experiment_E/PROVISIONAL_ANALYTIC_RESULT.json",ehash))
    isempty(hsh) || push!(sha_rows,(path=p,sha256=hsh))
end
CSV.write(joinpath(REPORT,"SHA256SUMS.csv"),DataFrame(sha_rows))

open(joinpath(EXP,"README_REPRODUCE.md"),"w") do io
    println(io,"# Reproduce Experiment G\n")
    println(io,"Run from the repository root on the pinned Julia 1.11 environment:\n")
    println(io,"    julia --project=. --startup-file=no experiments/bnd_expG/run_experiment_G.jl")
    println(io,"\nAfter Z_G_FINAL.toml exists, reproduce the analytical tables while verifying and preserving that frozen candidate:\n")
    println(io,"    julia --project=. --startup-file=no experiments/bnd_expG/run_experiment_G.jl --resume-frozen")
    println(io,"    python experiments/bnd_expG/render_figures_G.py")
    println(io,"    julia --project=. --startup-file=no test/bnd_expG/runtests.jl")
    println(io,"    julia --project=. --startup-file=no experiments/bnd_expG/validate_powerdynamics_G.jl\n")
    println(io,"The design runner reads frozen D/E input artifacts but writes only under experiments/bnd_expG and reports/experiment_G. It does not import PowerDynamics. The separate validator refuses to run unless Z_G_FINAL.toml and its SHA-256 sidecar verify.")
    println(io,"\nThe all-GFL endpoint has a defective gauge/physical-zero pair in the full state matrix. ExpG deflates the exact global-angle gauge and reports the remaining physical zero separately. The scalar authority LP is a seed; the final analytical candidate is checked against the complete finite spectrum and robust conditions.")
    println(io,"\nThe robust radius is dimensionless relative to the declared 1 s^-1 additive full-block state-matrix scale. It is not a physical percentage standard. No unequilibrated load-setpoint variation is asserted.")
end

open(joinpath(REPORT,"DERIVATION_EXP_G.md"),"w") do io
    println(io,"# Experiment G derivation\n")
    println(io,"## Isolated PLL seed\n")
    println(io,"Linearizing the implemented phase detector at aligned terminal voltage gives e = -V delta-theta. Together with dot(theta)=delta-omega, tau dot(delta-omega)=delta-omega-i + Kp e - delta-omega, and dot(delta-omega-i)=Ki e, elimination yields tau s^3+s^2+V Kp s+V Ki=0. Matching tau(s+r)^2(s+q) and tau(2r+q)=1 gives q=1/tau-2r, Kp=(2r-3 tau r^2)/V, and Ki=(r^2-2 tau r^3)/V. Since the sum of the three roots is -1/tau, a common decay r obeys 3r<=1/tau; the triple-root seed is r=1/(3tau). It is only a local reference.")
    println(io,"\n## Collective authority and endpoint singularity\n")
    println(io,"The all-GFL state matrix has an exact global-angle gauge vector g, checked by ||Ag||/(||A||||g||). The quotient basis Q=null(g') defines Aq=Q' A Q, whose spectrum removes precisely one gauge root. ExpG estimates gamma_i=-d Re(lambda_c)/d epsilon_i by centered full-state differences and compares steps 1e-4 and 5e-5. Because the endpoint also has a physical zero and the unreduced zero is defective, this derivative is a local quotient-branch seed, not a full design law.")
    println(io,"\n## Finite surrogate enumeration\n")
    println(io,"For frozen gains and a frozen robustness radius, the surrogate is min P'epsilon subject to gamma'epsilon>=bS, h'epsilon>=bR, 0<=epsilon<=1. Every upper-bound subset is enumerated; with two active scalar inequalities, a vertex has at most two remaining fractional coordinates. Single-anchor and pair intersections are solved explicitly. This certifies the linearized authority surrogate only.")
    println(io,"\n## Robust margin\n")
    println(io,"For the reduced, gauge-free matrix Aq, the full-block uncertainty is Aq+scale*Delta, ||Delta||2<=beta. On the shifted boundary s=-sigma+j omega, ExpG adaptively refines the largest singular value of (sI-(Aq+sigma I))^-1; the numerical small-gain test is beta||M_sigma||infinity<1. It also computes R_delta=sum ||r_j||||l_j||/|l_j'r_j| and imposes the conservative nominal separation alpha<=-sigma-beta R_delta. The complex full block is an outer uncertainty model, not an exact statement about structured real parameter uncertainty.")
    println(io,"\n## RoCoF and transients\n")
    println(io,"With h_i=H_i S_i, the initial COI bound is |df/dt(0+)|=f0|DeltaP|/(2 h'epsilon). The chosen illustrative scenario is recorded in DESIGN_SCENARIO_FROZEN.toml; it is not a universal grid standard. The final transient uses the same descriptor-reduced A, a one-MW current disturbance at the declared load bus, and an inertia-weighted retained-SG COI frequency output. Modal residues and direct matrix exponentials are compared on the same frozen linear model.")
    println(io,"\n## Exactness and scope\n")
    println(io,"The final candidate is checked using the complete finite spectrum of Ared, not only det(C). The closure determinant identity is tested at a regular complex point, and the low-rank determinant lemma is kept as a boundary predictor. The deterministic coordinate candidate is not a globally certified full-closure optimum; no KKT/SOSC certificate is claimed.")
end

open(joinpath(REPORT,"REPORT_EXP_G.md"),"w") do io
    println(io,"# Experiment G — analytic robust SG to GFL co-design\n")
    println(io,"**Status: $(gstatus).** Candidate hash: $(final_hash). Design did not import or call PowerDynamics.")
    println(io,"\n## Frozen system and endpoint gates\n")
    println(io,"The program discovered buses $(join(buses,", ")), reconstructed $(totalP) MW of initial SG electrical injection, and preserved initialized impedance-load values including buses 31 and 39. The all-SG spectrum comparison status is **$(baseline_status)** (Hausdorff residual $(haus)); no-load KCL residual is $(net.no_load_kcl_inf).")
    println(io,"\nThe all-GFL endpoint was classified **$(join(unique(g1.class),", "))** across three deterministic gain patterns. The global-angle gauge residuals are $(join(g1.uniform_angle_gauge_residual,", ")); the physical zero remains when gains reach their frozen maxima. The scalar authority is therefore treated as a local quotient-branch seed.")
    println(io,"\n## Design results\n")
    println(io,"The exhaustive linearized surrogate returned **$(surrogate.status)** with support $(surrogate.feasible ? join(buses[surrogate.support],", ") : "none") and retained $(surrogate.feasible ? surrogate.cost : NaN) MW. It is not the full-closure optimum.")
    println(io,"\nThe best deterministic exact-spectrum coordinate candidate used order $(order_names[bestindex]), retains $(retained) MW, converts $(gflmw) MW, and has spectral abscissa $(final_exact.spectrum.spectral_abscissa) s^-1. Its conservative modal robustified abscissa is $(final_exact.spectrum.spectral_abscissa+final_exact.mrob) s^-1, against required -$(SIGMA) s^-1.")
    println(io,"\nThe numerical full-block small-gain check at beta=$(beta_target) is **$(rob_candidate.small_gain_pass)** with beta-star $(rob_candidate.beta_star), peak resolvent $(rob_candidate.resolvent_peak), and margin $(rob_candidate.small_gain_margin). The modal-residue separation is more conservative and is reported separately.")
    println(io,"\n## Robustness, RoCoF, and transient\n")
    println(io,"Retained inertia is $(final_exact.inertia) MVA s against requirement $(bR) MVA s. Transient status is **$(transient.status)**; unit RoCoF gain is $(caprow.rocof_unit_gain_Hz_s_per_MW) Hz/s per MW and unit frequency gain is $(caprow.frequency_unit_gain_Hz_per_MW) Hz per MW. Direct/modal relative errors are $(caprow.direct_modal_rocof_relative_error) and $(caprow.direct_modal_frequency_relative_error).")
    println(io,"\nThe local perturbation screen status is **$(local_audit_status)** over $(tested) deterministic directions. Full-closure KKT/SOSC is **$(kkt_status)**; neither a strict local optimum nor a full-closure global optimum is claimed.")
    println(io,"\n## Blinded ExpE comparison\n")
    println(io,"Only after writing and hashing the pre-blind candidate did the runner read ExpE's provisional result. ExpE retained SG MW parsed as $(e_mw), anchor bus parsed as $(ebus); the comparison and source hash are in TABLE_G14_blinded_ExpE_comparison.csv. ExpG's candidate was not changed after this read.")
    println(io,"\n## Remaining independent gates\n")
    println(io,"PowerDynamics spectral validation and nonlinear TDS are **NOT_RUN** until the separate post-freeze validator is executed. Status remains PARTIAL until both are reported.")
end

expEerr=isfinite(e_mw) ? string(retained-e_mw) : "NOT_CHECKED"
bus38match=isfinite(ebus) ? string((38 in buses[findall(>(1e-8),epsfinal)])==(round(Int,ebus)==38)) : "NOT_CHECKED"
buslist=join(buses,",")
classlist=join(unique(g1.class),",")
physicalzerolist=join(g1.physical_nearzero_real,",")
surrogate_support_str=surrogate.feasible ? join(buses[surrogate.support],",") : "INFEASIBLE"
exact_support_str=join(buses[findall(>(1e-10),epsfinal)],",")
rholist=join(rhofinal,",")
kplist=join(kpfinal,",")
kilist=join(kifinal,",")
summary_lines=[
"EXP_G_STATUS: PARTIAL",
"DESIGN_USED_POWERDYNAMICS: NO",
"IEEE39_REPLACEABLE_BUSES: $(buslist)",
"TOTAL_INITIAL_SG_MW: $(totalP)",
"ALL_GFL_ZERO_CLASS: $(classlist)",
"ALL_GFL_PHYSICAL_ZERO: $(physicalzerolist)",
"GAIN_CAN_REMOVE_ZERO_ALONE: NO",
"SURROGATE_GLOBAL_SUPPORT: $(surrogate_support_str)",
"SURROGATE_RETAINED_BUSES: $(surrogate_support_str)",
"SURROGATE_RETAINED_SG_MW: $(surrogate.feasible ? surrogate.cost : NaN)",
"EXACT_CORRECTED_RETAINED_BUSES: $(exact_support_str)",
"EXACT_CORRECTED_RETAINED_SG_MW: $(retained)",
"MAX_GFL_MW: $(gflmw)",
"MAX_GFL_FRACTION: $(gflmw/totalP)",
"FINAL_RHO: $(rholist)",
"FINAL_KP: $(kplist)",
"FINAL_KI: $(kilist)",
"FINAL_SPECTRAL_ABSCISSA: $(final_exact.spectrum.spectral_abscissa)",
"SIGMA_REQUIRED: $(SIGMA)",
"ROBUST_BETA_CERTIFIED: $(beta_target)",
"SMALL_GAIN_CERTIFICATE: $(rob_candidate.small_gain_pass)",
"CRITICAL_POLE: $(final_exact.spectrum.critical)",
"CRITICAL_POLE_CONDITION_NUMBER: $(final_exact.spectrum.critical_nonnormality)",
"ROCOF_UNIT_GAIN: $(caprow.rocof_unit_gain_Hz_s_per_MW)",
"FREQUENCY_UNIT_GAIN: $(caprow.frequency_unit_gain_Hz_per_MW)",
"DELTA_P_MAX_INITIAL_ROCOF: $(caprow.deltaP_max_initial_rocof_MW)",
"DELTA_P_MAX_PEAK_ROCOF: $(caprow.deltaP_max_peak_rocof_MW)",
"DELTA_P_MAX_FREQUENCY: $(caprow.deltaP_max_frequency_MW)",
"DELTA_P_MAX_TOTAL: $(caprow.deltaP_max_total_MW)",
"SURROGATE_GLOBAL_OPTIMUM: YES",
"FULL_CLOSURE_LOCAL_OPTIMUM: NO",
"FULL_CLOSURE_GLOBAL_OPTIMUM: NO",
"EXP_E_BUS38_MATCH: $(bus38match)",
"EXP_E_RETAINED_MW_ERROR: $(expEerr)",
"POWERDYNAMICS_VALIDATION: NOT_RUN",
"TDS_VALIDATION: NOT_RUN",
"MAIN_THEORETICAL_RESULT: defective gauge/physical-zero endpoint requires exact gauge quotient; scalar authority is a local seed",
"MAIN_POWER_SYSTEM_RESULT: deterministic exact-spectrum robust coordinate candidate retains $(retained) MW SG",
"MAIN_LIMITATION: no KKT/SOSC global certificate; PowerDynamics and TDS remain pending",
"PUSH: NO"]
write(joinpath(REPORT,"FINAL_SUMMARY_EXP_G.md"),join(summary_lines,"\n")*"\n")
println(join(summary_lines,"\n"))
println("Z_G_FINAL_SHA256: ",final_hash)
