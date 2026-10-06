using LinearAlgebra, CSV, DataFrames, SHA, TOML, Statistics
BLAS.set_num_threads(1)
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
include(joinpath(ROOT,"src","bnd_design_k","BNDDesignK.jl"))
using .BNDDesignK
const OUT=joinpath(ROOT,"reports","experiment_K")
const TABLES=joinpath(OUT,"tables")
mkpath(TABLES)
const TARGET=-0.050001
const BETAS=[0.0,1e-6,1.6991206999182038e-6,3e-6,1e-5,3e-5,1e-4]
ctx=design_context(ROOT)
length(ctx.buses)==10 || error("wrong generator count")
abs(sum(ctx.power)-5402.761089978776)<1e-6 || error("SG dispatch integrity failure")
catalog=CSV.read(joinpath(TABLES,"TABLE_K01_support_catalog.csv"),DataFrame)
nrow(catalog)==1024 || error("support catalog incomplete")

function freeze(path,record)
    isfile(path) && error("candidate already frozen: $path")
    open(path,"w") do io; TOML.print(io,record;sorted=true) end
    digest=bytes2hex(sha256(read(path)))
    write(path*".sha256",digest*"\n")
    println("FROZEN: ",basename(path)," SHA-256=",digest)
    digest
end

function candidate_record(ev,kind)
    Dict("experiment"=>"EXP_K","status"=>"FROZEN_ANALYTICAL_CANDIDATE",
        "kind"=>kind,"blind_integrity"=>"CONTAMINATED_PRIOR_EXP_E_G_SUMMARY_READ",
        "global_optimality_certified"=>false,"design_used_PowerDynamics"=>false,
        "support_mask"=>ev.support,"support_buses"=>support_buses(ctx,ev.support),
        "epsilon"=>ev.epsilon,"rho"=>1 .- ev.epsilon,
        "Kp"=>ev.kp,"Ki"=>ev.ki,"P_initial_MW"=>ctx.power,
        "retained_SG_MW"=>ev.retained_mw,"GFL_MW"=>ev.gfl_mw,
        "spectral_abscissa_s_inv"=>ev.alpha,
        "required_spectral_abscissa_s_inv"=>-BNDDesignK.SIGMA_REQUIRED,
        "active_poles_real"=>real.(ev.lambda[ev.active]),
        "active_poles_imag"=>imag.(ev.lambda[ev.active]),
        "gauge_residual"=>ev.gauge_residual,
        "source_model_sha256"=>bytes2hex(sha256(read(joinpath(ROOT,"src","bnd_design_k","BNDDesignK.jl")))))
end

history=NamedTuple[]
bestrows=NamedTuple[]
seeds=Dict{Int,Any}()
anchor_levels=[0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.35,0.5,0.7,0.9,0.99]
for mask in 0:1023
    best=nothing
    for level in anchor_levels
        ee=support_epsilon(ctx,mask,level)
        ev=evaluate(ctx,ee)
        feasible=ev.alpha<=TARGET
        push!(history,(support_mask=mask,phase="uniform_seed",iteration=0,
            epsilon_level=level,retained_SG_MW=ev.retained_mw,
            alpha_s_inv=ev.alpha,feasible=feasible,
            note="complete finite quotient spectrum"))
        if feasible && (best===nothing || ev.retained_mw<best.retained_mw)
            best=ev
        end
    end
    best!==nothing && (seeds[mask]=best)
    mask%128==0 && (println("K design seed mask $mask / 1023");flush(stdout))
end

# Explicit one-coordinate active-boundary corrector. The global spectral
# abscissa is recomputed at every trial, so pole switches are honored. The
# sampling grid detects nonmonotone feasible components before bisection.
function reduce_coordinate(ev,i;max_sweeps=1)
    x=ev.epsilon[i]
    x>1e-8 || return ev
    samples=sort(unique(vcat(1e-8,x.*[0.05,0.1,0.2,0.35,0.5,0.7,0.85,1.0])))
    samples=[s for s in samples if s<=x]
    vals=Float64[]
    for s in samples
        ee=copy(ev.epsilon);ee[i]=s
        push!(vals,evaluate(ctx,ee,ev.kp,ev.ki).alpha)
    end
    feasible=findall(a->a<=TARGET,vals)
    isempty(feasible) && return ev
    k=first(feasible)
    lo=samples[k]
    if k>1
        a=samples[k-1];b=lo
        for _ in 1:32
            mid=(a+b)/2
            ee=copy(ev.epsilon);ee[i]=mid
            if evaluate(ctx,ee,ev.kp,ev.ki).alpha<=TARGET; b=mid else a=mid end
        end
        lo=b
    end
    if lo<x-1e-9
        ee=copy(ev.epsilon);ee[i]=lo
        trial=evaluate(ctx,ee,ev.kp,ev.ki)
        if trial.alpha<=TARGET+1e-8 && trial.retained_mw<ev.retained_mw
            return trial
        end
    end
    ev
end

function gain_authority_trial(ev)
    v=evaluate(ctx,ev.epsilon,ev.kp,ev.ki;vectors=true)
    # Equal nonnegative weights are a deterministic multi-active-mode trial,
    # not KKT multipliers or a proof of optimal controller elimination.
    gp=zeros(10);gi=zeros(10)
    for j in v.active
        g=fixed_support_gradient(ctx,v,j)
        gp .-= real.(g.kp)./length(v.active)
        gi .-= real.(g.ki)./length(v.active)
    end
    x=vcat(ev.kp,ev.ki)
    lo=vcat(ctx.kpmin,ctx.kimin);hi=vcat(ctx.kpmax,ctx.kimax)
    box=gain_box_support(vcat(gp,gi),x,lo,hi)
    winner=ev
    for step in (1/32,1/16,1/8,1/4,1/2,1.0)
        xx=x.+step.*box.delta
        trial=evaluate(ctx,ev.epsilon,xx[1:10],xx[11:20])
        if trial.alpha<winner.alpha-1e-8
            winner=trial
        end
    end
    winner
end

bestseed=sort(collect(values(seeds));by=x->x.retained_mw)[1]
println("K best uniform seed: mask=",bestseed.support," retained=",bestseed.retained_mw)
# Prioritize every seeded architecture capable of competing with the current
# incumbent under the sampled seed family; the others remain unexhausted.
ordered=sort(collect(keys(seeds));by=m->seeds[m].retained_mw)
incumbent=bestseed
optimized=0
for mask in ordered
    global incumbent, optimized
    seed=seeds[mask]
    if seed.retained_mw>max(3incumbent.retained_mw,50.0)
        push!(bestrows,(support_mask=mask,support_buses=join(support_buses(ctx,mask),":"),
            status="FEASIBLE_SEED_UNEXHAUSTED",retained_SG_MW=seed.retained_mw,
            GFL_MW=seed.gfl_mw,alpha_s_inv=seed.alpha,iterations=0))
        continue
    end
    current=seed; iterations=0
    for sweep in 1:3
        old=current.retained_mw
        for i in eachindex(ctx.buses)
            current.epsilon[i]>0 || continue
            next=reduce_coordinate(current,i)
            if next.retained_mw<current.retained_mw-1e-8
                current=next;iterations+=1
                push!(history,(support_mask=mask,phase="coordinate_corrector",iteration=iterations,
                    epsilon_level=current.epsilon[i],retained_SG_MW=current.retained_mw,
                    alpha_s_inv=current.alpha,feasible=current.alpha<=TARGET+1e-8,
                    note="bus $(ctx.buses[i]); branch switches recomputed"))
            end
        end
        old-current.retained_mw<1e-6 && break
    end
    if current.retained_mw<max(2incumbent.retained_mw,20.0)
        tuned=gain_authority_trial(current)
        if tuned.alpha<current.alpha-1e-8
            current=tuned
            for i in eachindex(ctx.buses)
                current.epsilon[i]>0 || continue
                current=reduce_coordinate(current,i)
            end
            push!(history,(support_mask=mask,phase="gain_box_authority",iteration=iterations+1,
                epsilon_level=NaN,retained_SG_MW=current.retained_mw,
                alpha_s_inv=current.alpha,feasible=current.alpha<=TARGET+1e-8,
                note="equal-weight active-mode box support trial"))
        end
    end
    optimized+=1
    push!(bestrows,(support_mask=mask,support_buses=join(support_buses(ctx,mask),":"),
        status="DETERMINISTIC_LOCAL_SEARCH_NOT_KKT_CERTIFIED",
        retained_SG_MW=current.retained_mw,GFL_MW=current.gfl_mw,
        alpha_s_inv=current.alpha,iterations=iterations))
    if current.alpha<=TARGET+1e-8 && current.retained_mw<incumbent.retained_mw
        incumbent=current
        println("K incumbent mask=",mask," retained=",incumbent.retained_mw,
            " alpha=",incumbent.alpha);flush(stdout)
    end
end
for mask in 0:1023
    haskey(seeds,mask) && continue
    push!(bestrows,(support_mask=mask,support_buses=join(support_buses(ctx,mask),":"),
        status="NO_FEASIBLE_UNIFORM_SEED_NOT_PROVEN_INFEASIBLE",
        retained_SG_MW=NaN,GFL_MW=NaN,alpha_s_inv=NaN,iterations=0))
end
sort!(bestrows;by=x->x.support_mask)
CSV.write(joinpath(TABLES,"TABLE_K07_support_branch_history.csv"),DataFrame(history))
CSV.write(joinpath(TABLES,"TABLE_K08_support_best_candidates.csv"),DataFrame(bestrows))

nominal=evaluate(ctx,incumbent.epsilon,incumbent.kp,incumbent.ki;vectors=true)
nominal.alpha<=TARGET+1e-8 || error("nominal candidate failed exact margin")
nominal_record=candidate_record(nominal,"NOMINAL_BEST_SAMPLED_ANALYTIC")
nominal_record["supports_enumerated"]=1024
nominal_record["supports_with_feasible_uniform_seed"]=length(seeds)
nominal_record["supports_locally_searched"]=optimized
nominal_record["KKT_certified"]=false
freeze(joinpath(OUT,"Z_K_NOMINAL_PREBLIND.toml"),nominal_record)
freeze(joinpath(OUT,"Z_K_NOMINAL_FINAL.toml"),nominal_record)

# The direct full-block radius is evaluated for a finite, declared candidate
# pool. No claim of robust Pareto optimality is made for unexamined interiors.
pool=Any[nominal]
for mask in ordered
    seed=seeds[mask]
    seed.retained_mw<=max(5nominal.retained_mw,100.0) || continue
    seed.alpha<=-0.06 || continue
    push!(pool,seed)
    length(pool)>=12 && break
end
unique_pool=Any[]
for c in pool
    any(x->x.support==c.support && abs(x.retained_mw-c.retained_mw)<1e-8,unique_pool) || push!(unique_pool,c)
end
radius=NamedTuple[]
for (idx,c) in enumerate(unique_pool)
    ev=evaluate(ctx,c.epsilon,c.kp,c.ki;vectors=false)
    rp=exact_beta_star(ev)
    push!(radius,(index=idx,support_mask=c.support,retained_SG_MW=c.retained_mw,
        GFL_MW=c.gfl_mw,alpha_s_inv=c.alpha,beta_star=rp.beta_star,
        small_gain_peak=rp.peak,peak_frequency_rad_s=rp.omega_peak,
        radius_status=rp.status))
    println("K radius $idx / $(length(unique_pool)): ",rp.beta_star);flush(stdout)
end
frontier=NamedTuple[]
for beta in BETAS
    eligible=[r for r in radius if r.beta_star>=beta && r.alpha_s_inv<=TARGET+1e-8]
    if isempty(eligible)
        push!(frontier,(beta_req=beta,status="NO_QUALIFIED_SAMPLED_CANDIDATE",
            support_mask=-1,retained_SG_MW=NaN,GFL_MW=NaN,
            alpha_s_inv=NaN,beta_star=NaN,small_gain_peak=NaN,
            peak_frequency_rad_s=NaN,global_certified=false))
    else
        r=sort(eligible;by=x->x.retained_SG_MW)[1]
        push!(frontier,(beta_req=beta,status="BEST_SAMPLED_CANDIDATE",
            support_mask=r.support_mask,retained_SG_MW=r.retained_SG_MW,GFL_MW=r.GFL_MW,
            alpha_s_inv=r.alpha_s_inv,beta_star=r.beta_star,
            small_gain_peak=r.small_gain_peak,peak_frequency_rad_s=r.peak_frequency_rad_s,
            global_certified=false))
    end
end
CSV.write(joinpath(TABLES,"TABLE_K09_robust_frontier.csv"),DataFrame(frontier))
for row in frontier
    row.status=="BEST_SAMPLED_CANDIDATE" || continue
    c=unique_pool[row.support_mask==nominal.support && abs(row.retained_SG_MW-nominal.retained_mw)<1e-8 ? 1 :
        findfirst(x->x.support==row.support_mask && abs(x.retained_mw-row.retained_SG_MW)<1e-8,unique_pool)]
    rec=candidate_record(evaluate(ctx,c.epsilon,c.kp,c.ki;vectors=true),"ROBUST_BEST_SAMPLED_ANALYTIC")
    rec["beta_req"]=row.beta_req
    rec["beta_star"]=row.beta_star
    rec["small_gain_peak"]=row.small_gain_peak
    rec["peak_frequency_rad_s"]=row.peak_frequency_rad_s
    label=replace(string(row.beta_req),"."=>"p","-"=>"m","+"=>"p")
    freeze(joinpath(OUT,"Z_K_ROBUST_beta_$(label).toml"),rec)
end
println("EXP_K_DESIGN_STATUS: PARTIAL_UNCERTIFIED")
println("SUPPORTS_EVALUATED: 1024")
println("SUPPORTS_WITH_FEASIBLE_SEEDS: ",length(seeds))
println("NOMINAL_BEST_SUPPORT: ",join(support_buses(ctx,nominal.support),","))
println("NOMINAL_RETAINED_SG_MW: ",nominal.retained_mw)
println("NOMINAL_GFL_MW: ",nominal.gfl_mw)
println("NOMINAL_ALPHA: ",nominal.alpha)
println("NOMINAL_GLOBAL_CERTIFIED: NO")
