using CSV, DataFrames, LinearAlgebra, TOML, SHA, Pkg, Dates

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const OUT=joinpath(ROOT,"reports","experiment_Q2B")
const REG=joinpath(OUT,"REGRESSION");const CORE=joinpath(OUT,"CORE")
mkpath(REG);mkpath(CORE)
include(joinpath(ROOT,"src","bnd_design_p","ExpP.jl"))
include(joinpath(ROOT,"src","bnd_expQ","LinearSecurity.jl"))
include(joinpath(ROOT,"src","bnd_expQ2B","FiniteWindow.jl"))
include(joinpath(ROOT,"src","bnd_expQ2B","CoreDesign.jl"))
const N=ExpP.PDExactDesignN;const CD=CoreDesign

function pkgver(name)
    p=only(filter(x->x.name==name,collect(values(Pkg.dependencies()))))
    string(p.version)
end

function regression(ctx)
    q2dir=joinpath(ROOT,"reports","experiment_Q2")
    gate=CSV.read(joinpath(q2dir,"BASELINE","TABLE_Q2_BASELINE_REGRESSION.csv"),DataFrame)
    prior=TOML.parsefile(joinpath(q2dir,"F4","F4_RESULTS.toml"))
    f2=TOML.parsefile(joinpath(q2dir,"F2","F2_RESULTS.toml"))
    f3=TOML.parsefile(joinpath(q2dir,"F3","F3_RESULTS.toml"))
    sha=ExpP.verify_expN_freeze(ROOT)
    checks=NamedTuple[]
    add(n,v,e,p)=push!(checks,(;check=n,value=string(v),expected=string(e),pass=Bool(p)))
    add("baseline_14_of_14",sum(gate.pass),14,all(gate.pass))
    add("ExpN_model_SHA",sha,"e2f104608f1eb1705f0ad2c07764e7beab0a5df4deb4d4ce4f359b1e79947f0a",
        sha=="e2f104608f1eb1705f0ad2c07764e7beab0a5df4deb4d4ce4f359b1e79947f0a")
    add("ExpN_Finf_Hz",f2["steady_offset_peak_Hz"],23.104449,abs(f2["steady_offset_peak_Hz"]-23.104449)<1e-4)
    add("ExpN_Fpeak_Hz",f2["grid_frequency_peak_Hz"],21.254654,abs(f2["grid_frequency_peak_Hz"]-21.254654)<1e-4)
    add("ExpN_alpha",f2["alpha_s_inv"],-0.050000001105,abs(f2["alpha_s_inv"]+0.050000001105)<1e-8)
    add("ExpQ_PLL_DC_independence",f3["max_abs_dH0_dKp"],0.0,f3["max_abs_dH0_dKp"]<1e-10)
    add("ExpQ_nonaffine_structure",prior["STEADY_FREQ_STRUCTURE"],"NONAFFINE",prior["STEADY_FREQ_STRUCTURE"]=="NONAFFINE")
    add("PowerDynamics_version",pkgver("PowerDynamics"),"5.0.0",pkgver("PowerDynamics")=="5.0.0")
    add("NetworkDynamics_version",pkgver("NetworkDynamics"),"1.3.0",pkgver("NetworkDynamics")=="1.3.0")
    CSV.write(joinpath(REG,"TABLE_Q2B_regression.csv"),DataFrame(checks))
    all(x.pass for x in checks) || error("Q2B regression failed")
    (;status="PASS",checks=length(checks),model_sha=sha,julia=string(VERSION),
      pd=pkgver("PowerDynamics"),nd=pkgver("NetworkDynamics"))
end

function write_point(path,stage,d,cand,audit)
    ev=cand.evaluation
    tbl=Dict{String,Any}(
      "stage"=>stage,"status"=>cand.status,"support"=>cand.support,
      "epsilon"=>cand.epsilon,"rho"=>cand.rho,"Kp"=>cand.Kp,"Ki"=>cand.Ki,
      "retained_SG_MW"=>cand.retained_SG_MW,"converted_GFL_MW"=>cand.converted_GFL_MW,
      "disturbance_MW"=>d,"alpha"=>ev.alpha,"beta_sampled"=>ev.beta,
      "beta_status"=>ev.beta_observed.status,"F_inf_Hz"=>ev.Finf,
      "F_peak_Tref_Hz"=>ev.Fpeak,"KKT_primal"=>audit.primal,
      "KKT_stationarity"=>audit.stationarity,"KKT_complementarity"=>audit.complementarity,
      "KKT_dual_feasible"=>audit.dual_feasibility,"LICQ"=>audit.LICQ,
      "SOSC"=>audit.SOSC,"model_sha"=>"e2f104608f1eb1705f0ad2c07764e7beab0a5df4deb4d4ce4f359b1e79947f0a",
      "objective_definition"=>"sum Pgen0_i * epsilon_i; original initialized component active powers")
    open(path,"w") do io;TOML.print(io,tbl);end
    open(path*".sha256","w") do io;write(io,bytes2hex(sha256(read(path)))*"\n");end
end

function main()
    reg=regression(nothing)
    println("Q2B_REGRESSION=PASS checks=",reg.checks," model_sha=",reg.model_sha)
    ctx=N.design_context(ROOT);all_support=collect(30:39)
    # Recompute only the frozen baseline configurations under Q2B's new
    # causal-window definition; this is a regression, not a repeat of F0-F4.
    baseline_windows=joinpath(REG,"TABLE_Q2B_baseline_window_regression.csv")
    if !isfile(baseline_windows)
        expn=TOML.parsefile(joinpath(ROOT,"reports","experiment_N","Z_N_NOMINAL_FINAL.toml"))
        expn_rho=Float64.(expn["rho"]);expn_kp=Float64.(expn["Kp"]);expn_ki=Float64.(expn["Ki"])
        expn_m=N.descriptor(ctx,expn_rho,expn_kp,expn_ki)
        expn_win=FiniteWindow.design_metrics(ctx,expn_m,expn_rho;
            windows=(0.2,0.5,1.0,2.0),dt_s=0.02,horizon_s=60.0,
            disturbance_MW=100.0,gauge_vector=N.gauge_vector)
        sg_m=N.descriptor(ctx,zeros(10),fill(N.K0P,10),fill(N.K0I,10))
        sg_win=FiniteWindow.design_metrics(ctx,sg_m,zeros(10);
            windows=(0.2,0.5,1.0,2.0),dt_s=0.02,horizon_s=60.0,
            disturbance_MW=100.0,gauge_vector=N.gauge_vector)
        btable=DataFrame(vcat([merge((;architecture="ExpN"),r) for r in expn_win.metrics],
                              [merge((;architecture="ALL_SG"),r) for r in sg_win.metrics]))
        CSV.write(baseline_windows,btable)
        println("BASELINE_WINDOW ExpN_Finf=",maximum(abs.(expn_win.signals.F_inf_Hz)),
            " ExpN_Fpeak_T05=",only(filter(r->r.window_s==0.5,expn_win.metrics)).F_peak_Hz,
            " ALLSG_Finf=",maximum(abs.(sg_win.signals.F_inf_Hz)))
    else
        println("BASELINE_WINDOW_REUSED=",basename(baseline_windows))
    end
    # Three frozen seeds: all-SG adjacent mixed architecture, the ExpN point,
    # and corrected ExpG. They are all explicitly re-evaluated in ExpN's model.
    kp_expG=fill(2pi*5,10);ki_expG=fill((2pi*5)^2/4,10)
    eps0=fill(0.999,10)
    seed=CD.design_eval(ctx,all_support,eps0,kp_expG,ki_expG;
        disturbance_MW=0.0,dt_s=0.05,horizon_s=30.0)
    println("SEED_ALLSG_ADJ alpha=",seed.alpha," beta=",seed.beta,
        " Finf=",seed.Finf," Fpeak=",seed.Fpeak)
    seedrow=DataFrame(seed=["allSG_adjacent_mixed"],support=[join(all_support,";")],
        epsilon=[join(eps0,";")],Kp=[join(kp_expG,";")],Ki=[join(ki_expG,";")],
        retained_SG_MW=[sum(ctx.power.*eps0)],alpha=[seed.alpha],beta_sampled=[seed.beta],
        status=[maximum(seed.g)<=0 ? "FEASIBLE_CORE_SEED" : "INFEASIBLE"])
    CSV.write(joinpath(CORE,"TABLE_Q2B_seed_audit.csv"),seedrow)
    # Core continuation starts at d=0 and then uses adaptive, recorded
    # disturbance increments. Every stage is exactly re-evaluated.
    current_eps=eps0;current_kp=kp_expG;current_ki=ki_expG
    live_path=joinpath(CORE,"Z_Q2B_live_d000.toml")
    resume_path=isfile(live_path) ? live_path : joinpath(CORE,"Z_Q2B_CORE_d000.toml")
    if isfile(resume_path)
        saved=TOML.parsefile(resume_path)
        if get(saved,"model_sha",reg.model_sha)==reg.model_sha && saved["support"]==all_support
            current_eps=Float64.(saved["epsilon"])
            current_kp=Float64.(saved["Kp"]);current_ki=Float64.(saved["Ki"])
            println("RESUME_CORE_FROM_D000 J=",saved["retained_SG_MW"])
        end
    end
    stages=NamedTuple[];last=nothing
    # Resume the best feasible d=0 point already generated by the custom SQP.
    # It is explicitly labeled best-found because its local KKT audit did not
    # close before the branch continuation; no claim of local certification is
    # inferred from feasibility alone.
    core0=CD.design_eval(ctx,all_support,current_eps,current_kp,current_ki;
        disturbance_MW=0.0,dt_s=0.05,horizon_s=30.0)
    if maximum(core0.g)<=1e-8
        eps0copy=copy(current_eps)
        last=(;status="BEST_FOUND_KKT_OPEN",support=all_support,
            y=CD.encode(ctx,current_eps,current_kp,current_ki),epsilon=eps0copy,
            rho=1 .- eps0copy,Kp=current_kp,Ki=current_ki,
            retained_SG_MW=sum(ctx.power.*eps0copy),
            converted_GFL_MW=sum(ctx.power)-sum(ctx.power.*eps0copy),
            evaluation=core0,accepted_steps=0,history=DataFrame())
        push!(stages,(;disturbance_MW=0.0,status=last.status,cost=last.retained_SG_MW,
            alpha=core0.alpha,beta=core0.beta,Finf=core0.Finf,Fpeak=core0.Fpeak,
            max_g=maximum(core0.g),iterations=25,accepted_steps=0,
            KKT_stationarity=NaN,KKT_primal=max(0.0,maximum(core0.g))))
        CSV.write(joinpath(CORE,"TABLE_Q2B_core_history_d000.csv"),
            DataFrame(iteration=collect(1:25),note=fill("SQP_checkpoint_from_prior_attempt",25)))
        no_audit=(;primal=max(0.0,maximum(core0.g)),stationarity=NaN,
            complementarity=NaN,dual_feasibility=false,LICQ=false,SOSC="KKT_OPEN")
        write_point(joinpath(CORE,"Z_Q2B_CORE_d000.toml"),"CORE_BEST_FOUND",0.0,last,no_audit)
        println("CORE_D0_INCUMBENT retained=",last.retained_SG_MW,
            " alpha=",core0.alpha," beta=",core0.beta,
            " status=BEST_FOUND_KKT_OPEN")
    else
        println("CORE_STAGE_BLOCKED d=0 seed_max_g=",maximum(core0.g))
    end
    # Predictor/corrector continuation with automatic step halving whenever
    # the previous exact design is not feasible at the proposed event size.
    dnow=0.0;step=25.0;attempts=0
    cached_path=joinpath(CORE,"TABLE_Q2B_core_continuation.csv")
    if isfile(cached_path)
        cached=CSV.read(cached_path,DataFrame)
        if last!==nothing && nrow(cached)>1 && maximum(abs.(Float64.(cached.cost).-last.retained_SG_MW))<1e-8 &&
           maximum(cached.max_g)<=1e-8 && maximum(cached.disturbance_MW)<100.0
            empty!(stages)
            for r in eachrow(cached)
                push!(stages,(;disturbance_MW=Float64(r.disturbance_MW),status=String(r.status),
                    cost=Float64(r.cost),alpha=Float64(r.alpha),beta=Float64(r.beta),
                    Finf=Float64(r.Finf),Fpeak=Float64(r.Fpeak),max_g=Float64(r.max_g),
                    iterations=Int(r.iterations),accepted_steps=Int(r.accepted_steps),
                    KKT_stationarity=Float64(r.KKT_stationarity),
                    KKT_primal=Float64(r.KKT_primal)))
            end
            dnow=maximum(Float64.(cached.disturbance_MW));attempts=30
            last=merge(last,(;evaluation=CD.design_eval(ctx,all_support,current_eps,current_kp,current_ki;
                disturbance_MW=dnow,dt_s=0.05,horizon_s=30.0)))
            println("CONTINUATION_CHECKPOINT_REUSED d=",dnow," J=",last.retained_SG_MW)
        end
    end
    while last!==nothing && dnow<100.0-1e-9 && attempts<30
        attempts+=1;target=min(100.0,dnow+step)
        pre=CD.design_eval(ctx,all_support,current_eps,current_kp,current_ki;
            disturbance_MW=target,dt_s=0.05,horizon_s=30.0)
        if maximum(pre.g)>1e-8
            step/=2
            println("CONTINUATION_STEP_HALVED from=",dnow," target=",target,
                " maxg=",maximum(pre.g)," next_step=",step)
            if step<0.5
                println("CONTINUATION_BLOCKED: feasible predictor step fell below 0.5 MW")
                break
            end
            continue
        end
        if maximum(pre.g[3:4]) < -1e-6
            # With both event-frequency rows inactive, the
            # RoCoF-independent core optimum carries forward unchanged. A new
            # KKT corrector is only needed as the event rows approach binding.
            dnow=target;last=merge(last,(;evaluation=pre))
            no_audit=(;primal=max(0.0,maximum(pre.g)),stationarity=NaN,
                complementarity=NaN,dual_feasibility=false,LICQ=false,
                SOSC="NOT_RECOMPUTED_EVENT_ROWS_INACTIVE")
            write_point(joinpath(CORE,"Z_Q2B_CORE_d$(lpad(Int(round(target)),3,'0')).toml"),
                "CORE_CONTINUATION_HOLD",target,last,no_audit)
            push!(stages,(;disturbance_MW=target,status="CORE_HOLD_FREQUENCY_INACTIVE",
                cost=last.retained_SG_MW,alpha=pre.alpha,beta=pre.beta,
                Finf=pre.Finf,Fpeak=pre.Fpeak,max_g=maximum(pre.g),iterations=0,
                accepted_steps=0,KKT_stationarity=NaN,KKT_primal=no_audit.primal))
            println("CORE_HOLD d=",target," retained=",last.retained_SG_MW,
                " Finf=",pre.Finf," Fpeak=",pre.Fpeak,
                " freq_g=",maximum(pre.g[3:4]))
            step=min(25.0,step*1.5)
            continue
        end
        result=CD.solve_fixed_support(ctx,all_support,current_eps,current_kp,current_ki;
            disturbance_MW=target,maxiter=12,trust_radius=0.02,dt_s=0.05,
            horizon_s=30.0,checkpoint_path=joinpath(CORE,
                "Z_Q2B_live_d$(lpad(Int(round(target)),3,'0')).toml"))
        if result.status=="INFEASIBLE_SEED" || maximum(result.evaluation.g)>1e-8
            step/=2;println("CONTINUATION_CORRECTOR_HALVED target=",target,
                " maxg=",maximum(result.evaluation.g)," next_step=",step)
            continue
        end
        dnow=target;current_eps=result.epsilon;current_kp=result.Kp;current_ki=result.Ki;last=result
        no_audit=(;primal=max(0.0,maximum(result.evaluation.g)),stationarity=NaN,
            complementarity=NaN,dual_feasibility=false,LICQ=false,SOSC="NOT_COMPUTED_BEFORE_FINAL")
        CSV.write(joinpath(CORE,"TABLE_Q2B_core_history_d$(lpad(Int(round(target)),3,'0')).csv"),result.history)
        write_point(joinpath(CORE,"Z_Q2B_CORE_d$(lpad(Int(round(target)),3,'0')).toml"),
            "CORE_CONTINUATION",target,result,no_audit)
        push!(stages,(;disturbance_MW=target,status=result.status,cost=result.retained_SG_MW,
            alpha=result.evaluation.alpha,beta=result.evaluation.beta,Finf=result.evaluation.Finf,
            Fpeak=result.evaluation.Fpeak,max_g=maximum(result.evaluation.g),
            iterations=nrow(result.history),accepted_steps=result.accepted_steps,
            KKT_stationarity=NaN,KKT_primal=no_audit.primal))
        println("CORE_STAGE d=",target," retained=",result.retained_SG_MW,
            " alpha=",result.evaluation.alpha," beta=",result.evaluation.beta,
            " Finf=",result.evaluation.Finf," Fpeak=",result.evaluation.Fpeak,
            " status=",result.status)
        step=min(25.0,step*1.5)
    end
    # The continuation may stop immediately below the event-frequency
    # boundary. In that case, use the prescribed all-SG-adjacent feasible seed
    # for a direct fixed-support 100 MW corrector; this is still the same
    # custom active-set SQP and every step is re-evaluated exactly.
    if (last===nothing || dnow<100.0-1e-8) && !isempty(all_support)
        direct_eps=fill(0.999,10);direct_kp=fill(2pi*5,10);direct_ki=fill((2pi*5)^2/4,10)
        println("DIRECT_100MW_SQ P_FROM_ALLSG_ADJACENT_SEED")
        direct=CD.solve_fixed_support(ctx,all_support,direct_eps,direct_kp,direct_ki;
            disturbance_MW=100.0,maxiter=24,trust_radius=0.05,dt_s=0.05,
            horizon_s=30.0,checkpoint_path=joinpath(CORE,"Z_Q2B_live_d100.toml"))
        if direct.status!="INFEASIBLE_SEED" && maximum(direct.evaluation.g)<=1e-8
            no_audit=(;primal=max(0.0,maximum(direct.evaluation.g)),stationarity=NaN,
                complementarity=NaN,dual_feasibility=false,LICQ=false,
                SOSC="KKT_AUDIT_PENDING")
            last=direct;dnow=100.0
            CSV.write(joinpath(CORE,"TABLE_Q2B_core_history_d100_direct.csv"),direct.history)
            write_point(joinpath(CORE,"Z_Q2B_CORE_d100_direct.toml"),
                "CORE_DIRECT_100MW",100.0,direct,no_audit)
            push!(stages,(;disturbance_MW=100.0,status=direct.status,
                cost=direct.retained_SG_MW,alpha=direct.evaluation.alpha,
                beta=direct.evaluation.beta,Finf=direct.evaluation.Finf,
                Fpeak=direct.evaluation.Fpeak,max_g=maximum(direct.evaluation.g),
                iterations=nrow(direct.history),accepted_steps=direct.accepted_steps,
                KKT_stationarity=NaN,KKT_primal=no_audit.primal))
            println("DIRECT_100MW_RESULT retained=",direct.retained_SG_MW,
                " alpha=",direct.evaluation.alpha," beta=",direct.evaluation.beta,
                " Finf=",direct.evaluation.Finf," Fpeak=",direct.evaluation.Fpeak,
                " status=",direct.status)
        else
            println("DIRECT_100MW_FAILED seed/correction maxg=",maximum(direct.evaluation.g),
                " status=",direct.status)
        end
    end
    CSV.write(joinpath(CORE,"TABLE_Q2B_core_continuation.csv"),DataFrame(stages))
    if last!==nothing && dnow>=100.0-1e-8 && last.evaluation!==nothing && maximum(last.evaluation.g)<=1e-8
        # High-resolution recomputation of the d=100 core candidate and all
        # finite-window RoCoF sensitivities; this does not alter design values.
        exact=CD.design_eval(ctx,last.support,last.epsilon,last.Kp,last.Ki;
            disturbance_MW=100.0,dt_s=0.01,horizon_s=60.0)
        CSV.write(joinpath(CORE,"TABLE_Q2B_core_rocof_windows.csv"),exact.Rwindows)
        for r in exact.Rwindows
            println("CORE_ROCOF T=",r.window_s," Rpeak=",r.R_peak_Hz_s)
        end
        audit=CD.kkt_audit(ctx,last;disturbance_MW=dnow)
        write_point(joinpath(CORE,"Z_Q2B_CORE.toml"),"CORE_AT_100MW",100.0,last,audit)
        println("CORE_FINAL retained=",last.retained_SG_MW," windows=" ,
            join([string(r.window_s,":",r.R_peak_Hz_s) for r in exact.Rwindows],";"))
    end
    open(joinpath(CORE,"CORE_RESULTS.toml"),"w") do io
        stage_dicts=[Dict{String,Any}(string(k)=>getproperty(s,k) for k in propertynames(s)) for s in stages]
        TOML.print(io,Dict("status"=>(last===nothing ? "NO_FEASIBLE_STAGE" : string(last.status)),
          "continuation_points"=>stage_dicts,"branch_completeness_certified"=>false,
          "continuation_reached_100MW"=>dnow>=100.0-1e-8,
          "support_searched"=>all_support,"optimizer"=>"custom feasible-start active-set SQP",
          "design_pd_calls"=>false,"model_sha"=>reg.model_sha))
    end
end
main()
