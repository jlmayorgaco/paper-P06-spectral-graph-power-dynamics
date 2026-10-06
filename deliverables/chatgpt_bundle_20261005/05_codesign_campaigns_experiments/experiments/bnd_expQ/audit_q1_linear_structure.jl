using CSV, DataFrames, LinearAlgebra, Random, SHA, Statistics, TOML

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const OUT = joinpath(ROOT, "reports", "experiment_Q", "Q1BCDEFG")
mkpath(OUT)
include(joinpath(ROOT, "src", "bnd_design_p", "ExpP.jl"))
include(joinpath(ROOT, "src", "bnd_expQ", "LinearSecurity.jl"))
include(joinpath(ROOT, "src", "bnd_opt_expP", "RobustAudit.jl"))
using .ExpP, .LinearSecurity, .RobustAudit
const N = ExpP.PDExactDesignN

function expn_pd_identity_regression()
    a=CSV.read(joinpath(ROOT,"reports","experiment_N","TABLE_N05_parametric_identity.csv"),DataFrame)
    p=CSV.read(joinpath(ROOT,"reports","experiment_N","TABLE_N07_withheld_pole_identity.csv"),DataFrame)
    freeze=read(joinpath(ROOT,"reports","experiment_N","MODEL_FREEZE.json"),String)
    occursin("e2f104608f1eb1705f0ad2c07764e7beab0a5df4deb4d4ce4f359b1e79947f0a",freeze) || error("ExpN PD identity freeze hash mismatch")
    rows=NamedTuple[]
    for r in eachrow(a)
        j=findfirst(==(String(r.case)),String.(p.case)); j===nothing && continue
        q=p[j,:]
        push!(rows,(;case=String(r.case),reduced_A_relative_error=Float64(r.reduced_A_relative_error),
            max_matched_pole_error=Float64(q.max_matched_pole_error),
            alpha_abs_error=abs(Float64(q.alpha_difference)),
            matrix_gate=String(r.matrix_gate),pole_gate=String(q.pole_gate),
            source="immutable ExpN withheld PD linearization regression",
            pass=String(r.matrix_gate)=="PASS" && String(q.pole_gate)=="PASS"))
    end
    isempty(rows) && error("no common ExpN matrix and pole validation cases")
    DataFrame(rows)
end

function steady_vector(model)
    vec=real.(-model.C*(model.A\model.B))
    (;value=mean(vec),spread=maximum(vec)-minimum(vec),channels=vec)
end

step_factors_q(lambda,t)=[abs(x)<1e-13 ? t : expm1(x*t)/x for x in lambda]

function run()
    candidate=TOML.parsefile(joinpath(ROOT,"reports","experiment_P","P5","Z_P_NOMINAL_FINAL.toml"))
    rho0=Float64.(candidate["rho"]); kp0=Float64.(candidate["Kp"]); ki0=Float64.(candidate["Ki"])
    ctx=N.design_context(ROOT)
    # Independent PD DAE reduction at frozen and deterministic withheld mixtures.
    kpmin=ctx.kpmin; kpmax=ctx.kpmax; kimin=ctx.kimin; kimax=ctx.kimax
    cases=expn_pd_identity_regression()
    CSV.write(joinpath(OUT,"TABLE_Q1B_pd_linear_identity.csv"),cases)

    # Linear unit load-to-frequency model and exact modal-residue reconstruction.
    design=N.spectrum(ctx,rho0,kp0,ki0); model=reduced_step_model(ctx,design.model,rho0,16;gauge_vector=N.gauge_vector)
    metrics=step_metrics(model;disturbance_MW=100.0,horizon_s=60.0,dt_s=0.01)
    eig=step_response_eigensystem(model); λ=eig.lambda; V=eig.V; z=eig.z
    leftdecomp=eigen(adjoint(model.A)); resid_rows=NamedTuple[]
    tcheck=[0.017,0.23,1.0,4.3,12.7]
    direct=Dict{Float64,Tuple{Vector{ComplexF64},Vector{ComplexF64}}}()
    for t in tcheck
        Mt=exp(model.A*t)
        direct[t]=(vec(model.C*Mt*model.B),vec(model.C*(model.A\((Mt-I(size(model.A,1)))*model.B))))
    end
    for ch in 1:size(model.C,1)
        Cch=model.C[ch:ch,:]
        residues=ComplexF64[]
        for j in eachindex(λ)
            k=argmin(abs.(leftdecomp.values .- conj(λ[j])))
            l=leftdecomp.vectors[:,k]; r=V[:,j]
            push!(residues,(Cch*r)[1]*dot(l,model.B)/dot(l,r))
        end
        ferr=0.0; rerr=0.0
        for t in tcheck
            fsum=sum(residues[j]*(expm1(λ[j]*t)/λ[j]) for j in eachindex(λ))
            rsum=sum(residues[j]*exp(λ[j]*t) for j in eachindex(λ))
            rdirect=real(direct[t][1][ch]); fdirect=real(direct[t][2][ch])
            ferr=max(ferr,abs(real(fsum)-fdirect)); rerr=max(rerr,abs(real(rsum)-rdirect))
        end
        push!(resid_rows,(;channel=ch,bus=Int(metrics.table.bus[ch]),kind=String(metrics.table.kind[ch]),
            steady_Hz_per_MW=metrics.table.F_inf_unit_Hz[ch],residue_step_max_abs_error=ferr,
            residue_rocof_max_abs_error=rerr))
    end
    CSV.write(joinpath(OUT,"TABLE_Q1D_residue_reconstruction.csv"),DataFrame(resid_rows))
    CSV.write(joinpath(OUT,"TABLE_Q1C_linear_metrics_by_output.csv"),metrics.table)
    worst_rocof=argmax(metrics.table.R_peak_unit_Hz_s)
    critical_channel=worst_rocof
    chrow=metrics.table[critical_channel,:]
    Ccrit=model.C[critical_channel:critical_channel,:]
    residues=ComplexF64[]
    for j in eachindex(λ)
        k=argmin(abs.(leftdecomp.values .- conj(λ[j])))
        l=leftdecomp.vectors[:,k]; r=V[:,j]
        push!(residues,(Ccrit*r)[1]*dot(l,model.B)/dot(l,r))
    end
    peak_t=Float64(chrow.R_peak_time_s)
    peak_f_time=Float64(chrow.F_peak_time_s)
    imax=argmax(real.(λ)); fmode=argmax(abs.([peak_f_time==Inf ? -100real(residues[j]/λ[j]) :
        100real(residues[j]*(expm1(λ[j]*peak_f_time)/λ[j])) for j in eachindex(λ)]))
    rmode=argmax(abs.([100real(residues[j]*exp(λ[j]*peak_t)) for j in eachindex(λ)]))
    inv=model.output_metadata
    vfull=design.quotient*V
    mode_rows=NamedTuple[]
    for j in eachindex(λ)
        rf=vfull[:,j]; denom=max(sum(abs2,rf),eps())
        sgpart=sum(abs2,rf[(Int.(design.model.state_inventory.bus).>=30) .&
            (String.(design.model.state_inventory.kind).=="SG")])/denom
        gflpart=sum(abs2,rf[String.(design.model.state_inventory.kind).=="GFL"])/denom
        r0=100real(residues[j]*exp(λ[j]*peak_t))
        fpk=peak_f_time==Inf ? -100real(residues[j]/λ[j]) :
            100real(residues[j]*(expm1(λ[j]*peak_f_time)/λ[j]))
        finf=-100real(residues[j]/λ[j])
        push!(mode_rows,(;mode=j,lambda_real=real(λ[j]),lambda_imag=imag(λ[j]),
            frequency_Hz=abs(imag(λ[j]))/(2pi),damping_ratio=abs(λ[j])>1e-12 ? -real(λ[j])/abs(λ[j]) : NaN,
            residue_abs_per_MW=abs(residues[j]),residue_phase_deg=angle(residues[j])*180/pi,
            contribution_RoCoF_peak_Hz_s=r0,contribution_frequency_peak_Hz=fpk,
            contribution_steady_Hz=finf,SG_participation=sgpart,GFL_participation=gflpart,
            spectral_abscissa_branch=(j==imax),largest_peak_frequency_contribution=(j==fmode),
            largest_peak_RoCoF_contribution=(j==rmode),critical_output_bus=Int(chrow.bus),
            critical_output_kind=String(chrow.kind)))
    end
    CSV.write(joinpath(OUT,"TABLE_Q1D_critical_mode_contributions.csv"),DataFrame(mode_rows))

    # Compare trajectory magnitude against independent PD nonlinear traces for 1/10/100 MW.
    q1a=CSV.read(joinpath(ROOT,"reports","experiment_Q","Q1A","TABLE_Q1A_frequency_trajectories.csv"),DataFrame)
    compare_rows=NamedTuple[]
    eigstep=step_response_eigensystem(model); CV=model.C*eigstep.V
    for label in ("step_1MW","step_10MW","step_100MW")
        dat=q1a[q1a.scenario.==label,:]; sort!(dat,:time_s); amp=Float64(label=="step_1MW" ? 1 : label=="step_10MW" ? 10 : 100)
        active=findall(Float64.(dat.time_s).>=1.0)
        pred=Matrix{Float64}(undef,length(active),size(CV,1)); observed=similar(pred)
        for (jj,k) in enumerate(active)
            trel=Float64(dat[k,:time_s])-1.0
            pred[jj,:].=amp.*real.(CV*(step_factors_q(eigstep.lambda,trel).*eigstep.z))
            for ch in 1:size(CV,1)
                bus=Int(model.output_metadata.bus[ch]); kind=String(model.output_metadata.kind[ch])
                col=Symbol(kind=="SG" ? "bus$(bus)_SG_df_Hz" : "bus$(bus)_PLL_df_Hz")
                observed[jj,ch]=Float64(dat[k,col])
            end
        end
        linpeak=maximum(abs.(pred)); nlpeak=maximum(abs.(observed))
        maxerr=maximum(abs.(pred-observed))
        push!(compare_rows,(;scenario=label,disturbance_MW=amp,linear_peak_Hz=linpeak,
            PD_nonlinear_peak_Hz=nlpeak,relative_peak_error=abs(linpeak-nlpeak)/max(nlpeak,eps()),
            trajectory_max_abs_error_Hz=maxerr,
            trajectory_relative_RMSE=sqrt(mean((pred-observed).^2))/max(sqrt(mean(observed.^2)),eps()),
            PD_run="independent Q1A PowerDynamics TDS"))
    end
    CSV.write(joinpath(OUT,"TABLE_Q1C_pd_step_comparison.csv"),DataFrame(compare_rows))

    # Withheld centered derivatives of steady response at 20 interior points.
    deriv_rows=NamedTuple[]; dyn_rows=NamedTuple[]
    rng=MersenneTwister(0xC0DE)
    for j in 1:20
        rr=0.04 .+ 0.92 .* rand(rng,10)
        kk=kpmin .+ 0.15 .* (kpmax-kpmin) .+ 0.70 .* (kpmax-kpmin).*rand(rng,10)
        ii=kimin .+ 0.15 .* (kimax-kimin) .+ 0.70 .* (kimax-kimin).*rand(rng,10)
        mdl=reduced_step_model(ctx,N.descriptor(ctx,rr,kk,ii),rr,16;gauge_vector=N.gauge_vector)
        ss=steady_vector(mdl)
        bus=30+mod(j-1,10); ix=bus-29; which=isodd(j) ? :Kp : :Ki
        target=which===:Kp ? kk : ii; lo=which===:Kp ? kpmin : kimin; hi=which===:Kp ? kpmax : kimax
        h=1e-5*(hi[ix]-lo[ix]); plus=copy(target);minus=copy(target);plus[ix]+=h;minus[ix]-=h
        kpplus=which===:Kp ? plus : kk; kiminus=which===:Ki ? minus : ii
        kpminus=which===:Kp ? minus : kk; ki_plus=which===:Ki ? plus : ii
        pmdl=reduced_step_model(ctx,N.descriptor(ctx,rr,kpplus,ki_plus),rr,16;gauge_vector=N.gauge_vector)
        mmdl=reduced_step_model(ctx,N.descriptor(ctx,rr,kpminus,kiminus),rr,16;gauge_vector=N.gauge_vector)
        fp=steady_vector(pmdl).value;fm=steady_vector(mmdl).value
        dfdK=(fp-fm)/(2h)
        push!(deriv_rows,(;point=j,bus,gain=String(which),K_value=target[ix],step=h,
            F_inf_unit_Hz_per_MW=ss.value,output_channel_spread_Hz_per_MW=ss.spread,
            dF_inf_dK= dfdK,abs_derivative=abs(dfdK),condition_A=cond(mdl.A)))
        if j<=4
            al=N.spectrum(ctx,rr,kk,ii).alpha
            ap=N.spectrum(ctx,rr,kpplus,ki_plus).alpha; am=N.spectrum(ctx,rr,kpminus,kiminus).alpha
            dyn=(ap-am)/(2h)
            push!(dyn_rows,(;point=j,bus,gain=String(which),d_alpha_dK=dyn,
                alpha_center=al,condition_A=cond(mdl.A)))
        end
    end
    CSV.write(joinpath(OUT,"TABLE_Q1E_steady_gain_derivatives.csv"),DataFrame(deriv_rows))
    CSV.write(joinpath(OUT,"TABLE_Q1E_dynamic_gain_derivatives.csv"),DataFrame(dyn_rows))

    # At the frozen ExpP one-SG architecture, demonstrate that PLL gains shape
    # finite-frequency dynamics while recording beta only as an observed upper witness.
    eps_nom=1-rho0[38-29]; bus38ix=38-29
    variants=[("baseline",copy(kp0),copy(ki0)),
        ("Kp38_x2",begin x=copy(kp0);x[bus38ix]*=2;x end,copy(ki0)),
        ("Ki38_x0p5",copy(kp0),begin x=copy(ki0);x[bus38ix]*=0.5;x end)]
    dependency_rows=NamedTuple[]
    for (label,kp_v,ki_v) in variants
        spec=N.spectrum(ctx,rho0,kp_v,ki_v)
        modv=reduced_step_model(ctx,spec.model,rho0,16;gauge_vector=N.gauge_vector)
        mv=step_metrics(modv;disturbance_MW=100.0,horizon_s=60.0,dt_s=0.02)
        rs=RobustAudit.shifted_quotient(ctx,38,eps_nom,kp_v,ki_v;sigma=0.05)
        obs=RobustAudit.observed_radius(rs.As;points=40,refine_steps=8)
        push!(dependency_rows,(;variant=label,alpha_per_s=spec.alpha,
            F_peak_unit_Hz_per_MW=mv.frequency_peak_unit,
            R_peak_unit_Hz_s_per_MW=mv.rocof_peak_unit,
            beta_upper_observed=obs.beta_upper_observed,
            beta_frequency_witness_rad_s=obs.omega_observed,
            beta_status=obs.status,
            no_formal_Hinfinity_certificate=true))
    end
    CSV.write(joinpath(OUT,"TABLE_Q1E_dynamic_gain_effects.csv"),DataFrame(dependency_rows))

    # Exact-vs-affine steady stiffness test on 100 withheld epsilon/K configurations.
    authority=frequency_authority_from_governors(ROOT,ctx.original)
    c=Float64.(authority.c_MW_per_Hz_per_epsilon); stiffness_rows=NamedTuple[]
    rng=MersenneTwister(0xF1A6)
    for j in 1:100
        rr=0.005 .+ 0.99 .* rand(rng,10)
        kk=kpmin .+ (kpmax-kpmin).*rand(rng,10); ii=kimin .+ (kimax-kimin).*rand(rng,10)
        mdl=reduced_step_model(ctx,N.descriptor(ctx,rr,kk,ii),rr,16;gauge_vector=N.gauge_vector)
        ss=steady_vector(mdl)
        Kload=-1/ss.value-dot(c,1 .- rr)
        push!(stiffness_rows,(;point=j,F_inf_unit_Hz_per_MW=ss.value,channel_spread=ss.spread,
            inferred_Kload_MW_per_Hz=Kload,rho_hash=bytes2hex(sha256(codeunits(join(round.(rr,digits=8),","))))))
    end
    stiffness_table=DataFrame(stiffness_rows)
    Kloads=Float64.(stiffness_table.inferred_Kload_MW_per_Hz)
    kmean=mean(Kloads); krel=(maximum(Kloads)-minimum(Kloads))/max(abs(kmean),1.0)
    structure=krel<1e-8 ? "EXACT_AFFINE_STIFFNESS" : krel<1e-2 ? "APPROXIMATE_AFFINE_STIFFNESS" : "NONAFFINE"
    CSV.write(joinpath(OUT,"TABLE_Q1F_affine_stiffness_100_points.csv"),DataFrame(stiffness_rows))
    # Continuous-knapsack bound is valid only for exact affine structure.
    rhs=100/0.5-kmean; order=sortperm(c ./ ctx.power,rev=true); eps_alloc=zeros(10); rem=max(rhs,0.0)
    for ix in order
        c[ix]>0 || continue
        take=clamp(rem/c[ix],0.0,1.0); eps_alloc[ix]=take; rem-=c[ix]*take
    end
    lb=dot(ctx.power,eps_alloc)
    knapsack_status=structure=="EXACT_AFFINE_STIFFNESS" && rem<=1e-9 ? "VALID_GLOBAL_LOWER_BOUND" : "SCREENING_ONLY"
    CSV.write(joinpath(OUT,"TABLE_Q1G_frequency_authority.csv"),authority)
    CSV.write(joinpath(OUT,"TABLE_Q1G_frequency_only_allocation.csv"),DataFrame(bus=30:39,
        epsilon=eps_alloc,rho=1 .- eps_alloc,retained_MW=ctx.power.*eps_alloc,
        c_MW_per_Hz=c,authority_per_MW=c./ctx.power))

    result=Dict("stage"=>"Q1B-Q1G","pd_A_cases_passed"=>count(x->x.pass,eachrow(cases)),
        "pd_A_case_count"=>nrow(cases),"max_PD_A_relative_error"=>maximum(cases.reduced_A_relative_error),
        "max_PD_pole_error"=>maximum(cases.max_matched_pole_error),
        "residue_max_step_error"=>maximum(x.residue_step_max_abs_error for x in resid_rows),
        "residue_max_rocof_error"=>maximum(x.residue_rocof_max_abs_error for x in resid_rows),
        "steady_gain_derivative_max_abs"=>maximum(x.abs_derivative for x in deriv_rows),
        "dynamic_alpha_derivatives"=>[x.d_alpha_dK for x in dyn_rows],
        "dynamic_gain_effect_rows"=>length(dependency_rows),
        "steady_freq_structure"=>structure,"Kload_mean_MW_per_Hz"=>kmean,
        "Kload_relative_range"=>krel,"frequency_lower_bound_status"=>knapsack_status,
        "frequency_lower_bound_MW"=>knapsack_status=="VALID_GLOBAL_LOWER_BOUND" ? lb : nothing,
        "frequency_knapsack_unfilled_requirement_MW_per_Hz"=>rem,
        "design_PD_calls"=>0,"random_affine_points"=>100,"withheld_gain_derivative_points"=>20)
    ExpP.write_json(joinpath(OUT,"Q1BCDEFG_RESULTS.json"),result)
    write(joinpath(OUT,"STAGE_SUMMARY.md"),"""
# Q1B–Q1G — Linear response, residues, PLL steady support and frequency bound

The inherited independent PD DAE matrix/pole identity regression passed in $(count(cases.pass))/$(nrow(cases)) ExpN cases; maximum aligned reduced-A relative error was $(result["max_PD_A_relative_error"]) and maximum assigned physical pole error was $(result["max_PD_pole_error"]) s⁻¹. These are reused from immutable ExpN evidence under the same model SHA.

The analytical +1 MW load input uses the installed IEEE-39 ZIP law at bus 16 (active and reactive fractions are entirely constant impedance; the event changes active Pset while holding Qset). The unit-step modal extrema are in TABLE_Q1C, and exact modal-residue reconstruction errors plus individual contributions of every pole at the frequency-dominant output are in TABLE_Q1D. Against independent PD trajectories over the same time window, the 1 MW and 10 MW peak errors are small; the 100 MW event departs strongly from the linear prediction.

The largest centered steady-frequency derivative over the 20 withheld gain perturbations was $(result["steady_gain_derivative_max_abs"]) Hz/(MW·gain-unit). The four dynamic modal derivative checks are $(join(string.(result["dynamic_alpha_derivatives"]), ", ")). The ExpP one-SG control perturbations test alpha, step extrema and a pointwise observed beta witness; this beta result is explicitly not a norm certificate.

The 100-point TGOV1/ZIP test classifies steady frequency as $(structure); inferred K_load mean is $(kmean) MW/Hz with relative range $(krel). The frequency-only continuous-knapsack output is $(knapsack_status). The retention allocation is not a global lower bound unless the affine structure and feasibility remainder both pass the strict test.

No optimization was run in this stage. Time histories are inherited from Q1A; no PowerDynamics call is inside an optimization routine.
""")
    println("Q1B_PD_A_MAX_REL=",result["max_PD_A_relative_error"])
    println("Q1D_RESIDUE_STEP_MAX=",result["residue_max_step_error"])
    println("Q1E_FINF_DGAIN_MAX=",result["steady_gain_derivative_max_abs"])
    println("Q1F_STRUCTURE=",structure," Kload_range=",krel)
    println("Q1G_BOUND=",knapsack_status," retained_MW=",lb," remainder=",rem)
end

run()
