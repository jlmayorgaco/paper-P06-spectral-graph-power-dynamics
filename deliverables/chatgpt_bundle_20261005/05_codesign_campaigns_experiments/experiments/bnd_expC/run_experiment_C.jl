using CSV
using DataFrames
using Dates
using LinearAlgebra
using SHA
using Statistics
using TOML
using PowerDynamics
using NetworkDynamics

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const OUT = joinpath(ROOT, "reports", "experiment_C")
const TABLES = joinpath(OUT, "tables")
const FIGURES = joinpath(OUT, "figures")
const MATRICES = joinpath(OUT, "matrices")
const CONFIG = TOML.parsefile(joinpath(@__DIR__, "configs", "expC.toml"))

include(joinpath(ROOT, "src", "bnd_graph", "BNDGraph.jl"))
include(joinpath(ROOT, "src", "pd39", "PD39.jl"))
using .BNDGraph
using .PD39

const GDS = BNDGraph.GeneralizedDynamicSelfEnergy
const GB = BNDGraph.GraphBackbone
const IP = BNDGraph.IntermodalPathways
const RP = BNDGraph.RealPoleAnalysis
const RSA = BNDGraph.RealSystemAdapter

for d in (OUT,TABLES,FIGURES,MATRICES); mkpath(d); end

const FREQS = sort!(unique!(vcat(10.0 .^ range(-3,2;length=81), collect(range(0.01,5.0;length=121)))))
const GIT_HEAD = strip(read(`git -C $ROOT rev-parse HEAD`,String))

json_escape(s) = replace(string(s), "\\"=>"\\\\", "\""=>"\\\"", "\n"=>"\\n", "\r"=>"\\r", "\t"=>"\\t")
function json_value(x)
    x === nothing && return "null"
    x === missing && return "null"
    x isa Bool && return x ? "true" : "false"
    x isa AbstractString && return "\"$(json_escape(x))\""
    x isa Symbol && return json_value(string(x))
    x isa Real && return isfinite(x) ? string(x) : "null"
    x isa Complex && return json_value("$(real(x)) + $(imag(x))im")
    x isa NamedTuple && return json_value(Dict(string(k)=>v for (k,v) in pairs(x)))
    x isa AbstractDict && return "{"*join(("\"$(json_escape(k))\":"*json_value(v) for (k,v) in sort!(collect(x);by=p->string(first(p)))),",")*"}"
    (x isa AbstractArray || x isa Tuple) && return "["*join(json_value.(collect(x)),",")*"]"
    return json_value(string(x))
end

function write_json(path,payload)
    open(path,"w") do io
        println(io,json_value(payload))
    end
end

function sha(path)
    return bytes2hex(SHA.sha256(read(path)))
end

function cache_fingerprint()
    paths=(joinpath(ROOT,"Project.toml"),joinpath(ROOT,"Manifest.toml"),
        joinpath(ROOT,"src","pd39","model.jl"),
        joinpath(ROOT,"src","bnd_graph","RealSystemAdapter.jl"),
        joinpath(ROOT,"reports","experiment_A","matrices","bus33_A_reduced.csv"),
        joinpath(ROOT,"reports","experiment_A","tables","TABLE_A03_state_partition.csv"))
    return bytes2hex(SHA.sha256(join([sha(p) for p in paths],"|")))
end

function save_case_cache(id,case; regenerated=false)
    stem=id=="C33"&&regenerated ? "C33_regenerated" : id
    matrix_csv(joinpath(MATRICES,"$(stem)_Ared.csv"),case.Ared)
    CSV.write(joinpath(MATRICES,"$(stem)_states.csv"),DataFrame(
        state_index=collect(eachindex(case.state_names)),state_name=case.state_names))
    open(joinpath(MATRICES,"$(stem)_cache.toml"),"w") do io
        TOML.print(io,Dict("schema_version"=>"1.0","fingerprint"=>cache_fingerprint(),
             "source"=>case.source,"case"=>id,"bus"=>case.bus))
    end
end

function load_case_cache(id; regenerated=false)
    stem=id=="C33"&&regenerated ? "C33_regenerated" : id
    meta_path=joinpath(MATRICES,"$(stem)_cache.toml")
    isfile(meta_path)||return nothing
    meta=TOML.parsefile(meta_path)
    get(meta,"fingerprint","")==cache_fingerprint()||return nothing
    A=RSA.matrix_from_csv(joinpath(MATRICES,"$(stem)_Ared.csv"))
    st=CSV.read(joinpath(MATRICES,"$(stem)_states.csv"),DataFrame)
    names=String.(st.state_name)
    bus=id=="C0" ? 0 : parse(Int,id[2:end])
    label=id=="C0" ? "C0_all_SG" : "$(id)_bus$(bus)_SimpleGFLDC"
    return RSA.case_from_reduced(label,bus,A,names)
end

function load_all_case_caches()
    out=Dict{String,Any}()
    for id in ("C0","C30","C33","C35","C37")
        c=load_case_cache(id;regenerated=(id=="C33"))
        c===nothing&&return nothing
        out[id]=c
    end
    return out
end

function matrix_csv(path,A)
    return RSA.write_matrix_csv(path,A)
end

function as_dataframe(rows)
    return isempty(rows) ? DataFrame() : DataFrame(rows)
end

function csv_write(name,rows)
    df=rows isa DataFrame ? rows : as_dataframe(rows)
    if isempty(names(df))
        schemas=Dict(
            "TABLE_C07_physical_pole_graph_mapping.csv"=>["case","pole_id","lambda_real","lambda_imag","frequency_hz","damping_ratio","retained_participation","dominant_graph_mode","graph_concentration","second_graph_mode","top1_top2_ratio","gauge_score","pll_state_share","dominant_states","role","q_norm"],
            "TABLE_C08_Gamma_exact.csv"=>["case","pole_id","graph_mode_k","evaluation_s_real","evaluation_s_imag","gamma_real","gamma_imag","gamma_abs","Trr_sigma_min","Trr_cond","schur_residual","determinant_factorization_error"],
            "TABLE_C09_uncoupled_modal_roots.csv"=>["case","root_id","graph_mode","root_real","root_imag","root_residual","seed_pole_id","convergence","psi_derivative_check","nearest_controller_pole","degenerate_cluster"],
            "TABLE_C10_pole_predictor.csv"=>["case","graph_mode","s0_real","s0_imag","gamma_real","gamma_imag","D_real","D_imag","pred_real","pred_imag","exact_real","exact_imag","abs_error","rel_error","real_error","freq_error_hz","eta_shift","Trr_sigma_min","Trr_cond","offdiag_ratio","graph_concentration","nearest_complementary_root_distance","matched_pole_id"],
            "TABLE_C11_exact_pathways.csv"=>["case","pole_id","root_id","evaluation_s_real","evaluation_s_imag","evaluation_kind","mode_k","mode_l","mode_m","contribution_real","contribution_imag","contribution_abs","diagonal_path","rank","gamma_real","gamma_imag","reconstruction_error"],
            "TABLE_C12_pairwise_susceptibility.csv"=>["case","mode_k","mode_l","coupling_product_abs","dynamic_detuning_abs","susceptibility_real","susceptibility_imag","susceptibility_abs","exact_diagonal_path_abs","coupling_rank","susceptibility_rank","exact_rank","t_l0_real","t_l0_imag"],
            "TABLE_C13_baseline_vs_replacement.csv"=>["replacement_bus","matched_mode","baseline_lambda","replacement_lambda","delta_real","delta_imag","baseline_gamma","replacement_gamma","delta_gamma","graph_mode_overlap","baseline_graph_mode","replacement_graph_mode","q_shape_overlap"],
            "TABLE_C14_backbone_sensitivity.csv"=>["case","physical_mode","primary_graph_mode","secondary_graph_mode","subspace_overlap","gamma_primary","gamma_secondary","predictor_error_primary","predictor_error_secondary"],
            "TABLE_C21_degenerate_block_pathways.csv"=>["case","pole_id","root_id","evaluation_s_real","evaluation_s_imag","evaluation_kind","target_mode","target_cluster","row_cluster","column_cluster","row_modes","column_modes","block_trace_real","block_trace_imag","block_frobenius","target_cluster_dimension","diagonal_block","gamma_block_frobenius","reconstruction_error","degenerate_target"])
        if haskey(schemas,name)
            df=DataFrame()
            for col in schemas[name]; df[!,Symbol(col)]=Any[]; end
        end
    end
    CSV.write(joinpath(TABLES,name),df)
    return df
end

function markdown_tables(filenames)
    open(joinpath(OUT,"TABLES_EXP_C.md"),"w") do io
        println(io,"# Experiment C tables\n")
        for name in filenames
            path=joinpath(TABLES,name)
            isfile(path) || continue
            df=CSV.read(path,DataFrame)
            println(io,"## `$(replace(name,".csv"=>""))`\n")
            if nrow(df)==0
                println(io,"No rows were generated.\n")
                continue
            end
            println(io,"| ",join(string.(names(df))," | ")," |")
            println(io,"| ",join(fill("---",ncol(df))," | ")," |")
            for i in 1:min(nrow(df),30)
                vals=[replace(string(df[i,j]),"|"=>"\\|") for j in 1:ncol(df)]
                println(io,"| ",join(vals," | ")," |")
            end
            nrow(df)>30 && println(io,"\nShowing first 30 of $(nrow(df)) rows. Full data: `tables/$name`.\n")
            println(io)
        end
    end
end

function isfinite_matrix(A)
    return all(isfinite,real.(A)) && all(isfinite,imag.(A))
end

function matrix_p95(v)
    isempty(v) && return NaN
    return quantile(Float64.(v),0.95)
end

function describe_state(name)
    m=match(r"VIndex\((\d+),\s*:(.*)\)",name)
    m===nothing && return (bus=missing,component="unparsed")
    var=m.captures[2]
    parts=split(var,'₊')
    return (bus=parse(Int,m.captures[1]),component=length(parts)>1 ? join(parts[1:end-1],"₊") : var)
end

function build_op(case)
    x=case.matrices
    return GDS.AngleOperator(x.M,x.D0,x.L0,x.Avc,x.Acc,x.Acq,x.Acv)
end

function case_analysis(case)
    op=build_op(case)
    n=size(op.M,1)
    mass=GB.mass_audit(op.M)
    mass.spd || return (case=case,op=op,mass=mass,blocked="M is not Hermitian positive definite")
    g=ones(n)
    primary=GB.primary_backbone(op.M,op.L0,g)
    secondary=GB.euclidean_backbone(op.L0,g)
    basis=GB.generalized_basis(op.M,primary.LG)
    basis2=GB.generalized_basis(op.M,secondary.LG)
    poles=RP.select_real_poles(case.Ared,case.state_names,case.partition.q,
        case.partition.v,basis.Phi,op.M;low_hz=CONFIG["engineering_band_hz"][1],
        high_hz=CONFIG["engineering_band_hz"][2],
        keep_band=CONFIG["pole_selection_low_damping_pairs"])
    return (case=case,op=op,mass=mass,primary=primary,secondary=secondary,
        basis=basis,basis2=basis2,poles=poles,blocked="")
end

function eval_graph(ctx,s)
    op=ctx.op; b=ctx.basis; LG=ctx.primary.LG
    pi=GDS.pi_q(op,s)
    psi=GDS.psi(op,LG,s)
    phat=b.Phi' * psi * b.Phi
    That=IP.graph_operator(b.Lambda,phat,s)
    Tdirect=b.Phi' * GDS.retained_operator(op,s) * b.Phi
    residual=norm(Tdirect-That)/max(norm(Tdirect),eps(Float64))
    Dhat=s*(b.Phi' * op.D0 * b.Phi)
    DLhat=b.Phi'*(op.L0-LG)*b.Phi
    Pihat=b.Phi' * pi * b.Phi
    return (Pi=pi,Psi=psi,Psihat=phat,T_hat=That,Tdirect=Tdirect,
        graph_error=residual,D0hat=Dhat,DeltaLhat=DLhat,Pihat=Pihat)
end

function eval_frequency_tables(ctx)
    c=ctx.case
    graphrows=NamedTuple[]
    metricrows=NamedTuple[]
    chigrows=NamedTuple[]
    for f in FREQS
        s=2pi*f*im
        v=eval_graph(ctx,s)
        DG=IP.effective_dissipative_operator(v.Psihat,2pi*f)
        dvals=real.(eigvals(Hermitian((DG+DG')/2)))
        chi,C=IP.commutator_metric(ctx.basis.Lambda,v.Psihat)
        nrm=norm(v.Psihat)
        off=v.Psihat-Diagonal(diag(v.Psihat))
        push!(graphrows,(case=c.label,frequency_hz=f,s_real=real(s),s_imag=imag(s),
            direct_norm=norm(v.Tdirect),reconstruction_error=v.graph_error))
        push!(metricrows,(case=c.label,frequency_hz=f,psi_norm=norm(v.Psihat),
            pi_norm=norm(v.Pihat),D0_term_norm=norm(v.D0hat),DeltaL_norm=norm(v.DeltaLhat),
            diagonal_norm=norm(Diagonal(diag(v.Psihat))),offdiagonal_norm=norm(off),
            offdiag_ratio=norm(off)/(nrm+eps(Float64)),chi_comm=chi,
            DG_eff_min=minimum(dvals),DG_eff_max=maximum(dvals),
            DG_eff_offdiag_norm=norm(DG-Diagonal(diag(DG)))))
        cg=IP.candidate_chi_G(ctx.basis.Lambda,v.Psihat,s)
        push!(chigrows,(case=c.label,frequency_hz=f,chi_G=cg.value,
            sigma_min_TD=cg.sigma_min_TD,condition_TD=cg.condition_TD,
            classification="candidate collective coupling diagnostic; not a stability certificate"))
    end
    for row in ctx.poles.selected
        s=row.lambda+1e-4*(1+im)
        v=eval_graph(ctx,s)
        push!(graphrows,(case=c.label,frequency_hz=abs(imag(s))/(2pi),s_real=real(s),s_imag=imag(s),
            direct_norm=norm(v.Tdirect),reconstruction_error=v.graph_error))
    end
    return (graph=graphrows,metrics=metricrows,chig=chigrows)
end

function select_gamma_rows(ctx)
    rows=NamedTuple[]
    for pole in ctx.poles.selected
        k=pole.dominant_graph_mode
        1<=k<=length(ctx.basis.Lambda) || continue
        s=pole.lambda+1e-4*(1+im)
        v=eval_graph(ctx,s)
        gm=IP.schur_gamma(v.T_hat,k)
        push!(rows,(case=ctx.case.label,pole_id="eig$(pole.index)",graph_mode_k=k,
            evaluation_s_real=real(s),evaluation_s_imag=imag(s),gamma_real=real(gm.gamma),
            gamma_imag=imag(gm.gamma),gamma_abs=abs(gm.gamma),Trr_sigma_min=gm.sigma_min,
            Trr_cond=gm.condition,schur_residual=gm.residual,
            determinant_factorization_error=gm.determinant_factorization_error,
            graph_concentration=pole.graph_concentration,role=pole.role))
    end
    return rows
end

function roots_and_predictors(ctx,gamma_rows)
    op=ctx.op; basis=ctx.basis; poles=ctx.poles
    roots=NamedTuple[]
    pred=NamedTuple[]
    pathwayrows=NamedTuple[]
    blockpathrows=NamedTuple[]
    susceptibilityrows=NamedTuple[]
    actual=poles.all_positive
    seedrows=filter(r -> CONFIG["engineering_band_hz"][1] <= r.frequency_hz <= CONFIG["engineering_band_hz"][2],actual)
    modeks=unique([p.dominant_graph_mode for p in poles.selected if p.dominant_graph_mode>0])
    for k in modeks
        f=s->begin
            try
                IP.modal_diagonal(basis.Lambda,basis.Phi' * GDS.psi(op,ctx.primary.LG,s)*basis.Phi,s,k)
            catch
                complex(Inf,Inf)
            end
        end
        df=s->begin
            try
                2s+(basis.Phi' * GDS.psi_derivative(op,s)*basis.Phi)[k,k]
            catch
                complex(Inf,Inf)
            end
        end
        candidates=filter(r -> r.projection[k] >= CONFIG["pole_matching_min_graph_projection"],seedrows)
        isempty(candidates) && (candidates=seedrows)
        roots_k=ComplexF64[]
        for seed in candidates
            admissible=s->abs(imag(s))/(2pi) <= 15.0 && abs(real(s)) <= 1e4 &&
                try cond(s*I-op.Acc) < 1e12 catch; false end
            nr=IP.safeguarded_newton(f,df,seed.lambda;tol=1e-10,maxiter=80,max_step=2.0,admissible=admissible)
            converged=nr.converged && CONFIG["engineering_band_hz"][1] <= abs(imag(nr.root))/(2pi) <= CONFIG["engineering_band_hz"][2]
            if converged && all(abs(nr.root-r)>1e-6 for r in roots_k)
                push!(roots_k,nr.root)
                root_id="mode$(k)_root$(length(roots_k))"
                push!(roots,(case=ctx.case.label,root_id=root_id,graph_mode=k,root_real=real(nr.root),root_imag=imag(nr.root),
                    root_residual=abs(f(nr.root)),seed_pole_id="eig$(seed.index)",convergence=true,
                    psi_derivative_check=GDS.derivative_check(op,nr.root),
                    nearest_controller_pole=minimum(abs.(nr.root .- eigvals(op.Acc))),
                    degenerate_cluster=length(first(filter(cl -> k in cl,basis.clusters)))>1))
            end
        end
        for (root_index,s0) in enumerate(roots_k)
            root_id="mode$(k)_root$(root_index)"
            deg=length(first(filter(cl -> k in cl,basis.clusters)))>1
            v=eval_graph(ctx,s0)
            bp=IP.block_pathway_matrix(v.Psihat,basis.Lambda,s0,k,basis.clusters)
            for br in bp.contributions
                push!(blockpathrows,(case=ctx.case.label,pole_id="",root_id=root_id,
                    evaluation_s_real=real(s0),evaluation_s_imag=imag(s0),evaluation_kind="diagonal_root",
                    target_mode=k,
                    target_cluster=join(bp.target_cluster,","),row_cluster=br.row_cluster,
                    column_cluster=br.column_cluster,row_modes=br.row_modes,
                    column_modes=br.column_modes,block_trace_real=real(br.block_trace),
                    block_trace_imag=imag(br.block_trace),block_frobenius=br.block_frobenius,
                    target_cluster_dimension=br.block_dimension,diagonal_block=br.diagonal_block,
                    gamma_block_frobenius=norm(bp.gamma),reconstruction_error=bp.reconstruction_error,
                    degenerate_target=length(bp.target_cluster)>1))
            end
            if deg
                # Individual-mode scalar predictors are basis-dependent inside
                # repeated eigenvalue clusters, so retain roots but do not claim them.
                continue
            end
            gamma=IP.schur_gamma(v.T_hat,k)
            derivative=2s0+(basis.Phi' * GDS.psi_derivative(op,s0)*basis.Phi)[k,k]
            shift=IP.pole_shift(gamma.gamma,derivative,s0)
            cand=filter(r -> r.projection[k] >= CONFIG["pole_matching_min_graph_projection"],actual)
            isempty(cand) && (cand=actual)
            exact=isempty(cand) ? nothing : cand[argmin([abs(r.lambda-shift.predicted) for r in cand])]
            err=exact===nothing ? NaN : abs(exact.lambda-shift.predicted)
            rel=exact===nothing ? NaN : err/max(abs(exact.lambda),eps(Float64))
            eta=abs(shift.delta)/max(abs(s0),eps(Float64))
            off=v.Psihat-Diagonal(diag(v.Psihat))
            otherroots=[complex(r.root_real,r.root_imag) for r in roots if r.case==ctx.case.label && r.graph_mode != k]
            dist=isempty(otherroots) ? NaN : minimum(abs.(s0 .- otherroots))
            push!(pred,(case=ctx.case.label,graph_mode=k,s0_real=real(s0),s0_imag=imag(s0),
                gamma_real=real(gamma.gamma),gamma_imag=imag(gamma.gamma),D_real=real(derivative),D_imag=imag(derivative),
                pred_real=real(shift.predicted),pred_imag=imag(shift.predicted),
                exact_real=exact===nothing ? NaN : real(exact.lambda),exact_imag=exact===nothing ? NaN : imag(exact.lambda),
                abs_error=err,rel_error=rel,real_error=exact===nothing ? NaN : abs(real(exact.lambda)-real(shift.predicted)),
                freq_error_hz=exact===nothing ? NaN : abs(abs(imag(exact.lambda))-abs(imag(shift.predicted)))/(2pi),
                eta_shift=eta,Trr_sigma_min=gamma.sigma_min,Trr_cond=gamma.condition,
                offdiag_ratio=norm(off)/(norm(v.Psihat)+eps(Float64)),
                graph_concentration=exact===nothing ? NaN : exact.graph_concentration,
                nearest_complementary_root_distance=dist,matched_pole_id=exact===nothing ? "" : "eig$(exact.index)"))
            # Exact pathways, block-aggregated degenerate diagnostics, and the
            # diagonal-complement susceptibility are evaluated at this root.
            pw=IP.pathway_matrix(v.Psihat,basis.Lambda,s0,k)
            ord=sortperm(abs.(pw.G[:]);rev=true)
            for (rank,idx) in enumerate(ord)
                li,mi=Tuple(CartesianIndices(pw.G)[idx])
                l=pw.modes[li]; m=pw.modes[mi]
                push!(pathwayrows,(case=ctx.case.label,pole_id="",root_id=root_id,
                    evaluation_s_real=real(s0),evaluation_s_imag=imag(s0),evaluation_kind="diagonal_root",
                    mode_k=k,mode_l=l,mode_m=m,
                    contribution_real=real(pw.G[li,mi]),contribution_imag=imag(pw.G[li,mi]),
                    contribution_abs=abs(pw.G[li,mi]),diagonal_path=l==m,rank=rank,
                    gamma_real=real(pw.gamma),gamma_imag=imag(pw.gamma),
                    reconstruction_error=pw.reconstruction_error))
            end
            pair=IP.pairwise_susceptibility(v.Psihat,basis.Lambda,s0,k)
            diagpath=Dict(pw.modes[i]=>abs(pw.G[i,i]) for i in eachindex(pw.modes))
            exactabs=[get(diagpath,r.mode,0.0) for r in pair]
            cRank=sortperm(sortperm([r.coupling_product for r in pair]))
            sRank=sortperm(sortperm([r.susceptibility for r in pair]))
            eRank=sortperm(sortperm(exactabs))
            for j in eachindex(pair)
                r=pair[j]
                push!(susceptibilityrows,(case=ctx.case.label,mode_k=k,mode_l=r.mode,
                    coupling_product_abs=r.coupling_product,dynamic_detuning_abs=r.dynamic_detuning,
                    susceptibility_real=real(r.complex_susceptibility),
                    susceptibility_imag=imag(r.complex_susceptibility),susceptibility_abs=r.susceptibility,
                    exact_diagonal_path_abs=exactabs[j],coupling_rank=cRank[j],
                    susceptibility_rank=sRank[j],exact_rank=eRank[j],t_l0_real=real(r.t_l0),t_l0_imag=imag(r.t_l0)))
            end
        end
    end
    # The exact pathway identity is also evaluated at real full-system poles,
    # so it remains testable when a diagonal-root search does not converge.
    for pole in poles.selected
        k=pole.dominant_graph_mode
        k>0||continue
        s=pole.lambda+1e-4*(1+im)
        v=eval_graph(ctx,s)
        pw=IP.pathway_matrix(v.Psihat,basis.Lambda,s,k)
        for (rank,idx) in enumerate(sortperm(abs.(pw.G[:]);rev=true))
            li,mi=Tuple(CartesianIndices(pw.G)[idx])
            l=pw.modes[li]; m=pw.modes[mi]
            push!(pathwayrows,(case=ctx.case.label,pole_id="eig$(pole.index)",root_id="",
                evaluation_s_real=real(s),evaluation_s_imag=imag(s),evaluation_kind="full_pole_offset",
                mode_k=k,mode_l=l,mode_m=m,
                contribution_real=real(pw.G[li,mi]),contribution_imag=imag(pw.G[li,mi]),
                contribution_abs=abs(pw.G[li,mi]),diagonal_path=l==m,rank=rank,
                gamma_real=real(pw.gamma),gamma_imag=imag(pw.gamma),
                reconstruction_error=pw.reconstruction_error))
        end
        bp=IP.block_pathway_matrix(v.Psihat,basis.Lambda,s,k,basis.clusters)
        for br in bp.contributions
            push!(blockpathrows,(case=ctx.case.label,pole_id="eig$(pole.index)",root_id="",
                evaluation_s_real=real(s),evaluation_s_imag=imag(s),evaluation_kind="full_pole_offset",
                target_mode=k,
                target_cluster=join(bp.target_cluster,","),row_cluster=br.row_cluster,
                column_cluster=br.column_cluster,row_modes=br.row_modes,
                column_modes=br.column_modes,block_trace_real=real(br.block_trace),
                block_trace_imag=imag(br.block_trace),block_frobenius=br.block_frobenius,
                target_cluster_dimension=br.block_dimension,diagonal_block=br.diagonal_block,
                gamma_block_frobenius=norm(bp.gamma),reconstruction_error=bp.reconstruction_error,
                degenerate_target=length(bp.target_cluster)>1))
        end
    end
    return (roots=roots,predictors=pred,pathways=pathwayrows,
        block_pathways=blockpathrows,susceptibilities=susceptibilityrows)
end

function zero_frequency_diagnostic(ctx)
    op=ctx.op; A=op.Acc; n=size(A,1)
    tols=CONFIG["zero_singularity_relative_tolerances"]
    sv=svdvals(A); sv2=svdvals(A*A)
    nullrows=NamedTuple[]
    nullities=NamedTuple[]
    F=svd(A;full=true)
    for τ in tols
        n0=count(x->x<=τ*maximum(sv),sv)
        n2=count(x->x<=τ*maximum(sv2),sv2)
        push!(nullities,(relative_tolerance=τ,nullity_A=n0,nullity_A2=n2,
            sigma_min_A=minimum(sv),sigma_min_A2=minimum(sv2)))
    end
    τ=Float64(tols[2])
    n0=last(filter(r->r.relative_tolerance==τ,nullities)).nullity_A
    right= n0==0 ? zeros(Float64,n,0) : F.V[:,(n-n0+1):n]
    left= n0==0 ? zeros(Float64,n,0) : F.U[:,(n-n0+1):n]
    names=ctx.case.state_names[ctx.case.partition.condensed]
    for j in 1:size(right,2)
        order=sortperm(abs.(right[:,j]);rev=true)
        for (rank,i) in enumerate(order[1:min(20,length(order))])
            desc=describe_state(names[i])
            push!(nullrows,(singular_vector="right_$j",state_index=ctx.case.partition.condensed[i],
                state_name=names[i],bus=desc.bus,component=desc.component,abs_weight=abs(right[i,j]),rank=rank))
        end
    end
    for j in 1:size(left,2)
        order=sortperm(abs.(left[:,j]);rev=true)
        for (rank,i) in enumerate(order[1:min(20,length(order))])
            desc=describe_state(names[i])
            push!(nullrows,(singular_vector="left_$j",state_index=ctx.case.partition.condensed[i],
                state_name=names[i],bus=desc.bus,component=desc.component,abs_weight=abs(left[i,j]),rank=rank))
        end
    end
    P=nothing
    semi_res=NaN
    residue=nothing
    semisimple=all(r->r.nullity_A==r.nullity_A2 && r.nullity_A==n0,nullities)
    if n0>0 && semisimple
        overlap=left' * right
        P=right*(overlap\left')
        semi_res=max(norm(A*P),norm(P*A))/max(norm(A)*norm(P),eps(Float64))
        residue=-op.Avc*P*op.Acq
    end
    scalings=NamedTuple[]
    dirs=("positive_real"=>0.0,"positive_imaginary"=>pi/2,"45_degree"=>pi/4)
    for (d,θ) in dirs, ab in CONFIG["zero_frequency_magnitudes"]
        s=ab*cis(θ)
        Tcc=s*I-A
        condv=try cond(Tcc) catch; Inf end
        pinorm=try norm(GDS.pi_q(op,s)) catch; NaN end
        fitrange=CONFIG["zero_frequency_fit_range"]
        push!(scalings,(s_abs=ab,direction=d,Pi_norm2=pinorm,Tcc_sigma_min=minimum(svdvals(Tcc)),
            Tcc_cond=condv,fitted_region=(fitrange[1]<=ab<=fitrange[2] && condv<CONFIG["zero_frequency_fit_cond_limit"])))
    end
    fitrows=filter(r->r.fitted_region && isfinite(r.Pi_norm2) && r.Pi_norm2>0,scalings)
    if length(fitrows)>=3
        x=-log.(getproperty.(fitrows,:s_abs)); y=log.(getproperty.(fitrows,:Pi_norm2))
        slope=dot(x.-mean(x),y.-mean(y))/sum(abs2,x.-mean(x))
    else
        slope=NaN
    end
    residuecheck=NamedTuple[]
    if semisimple && residue!==nothing
        for (d,θ) in dirs, ab in (1e-4,1e-5,1e-6,1e-7)
            s=ab*cis(θ)
            val=s*GDS.pi_q(op,s)
            err=norm(val-residue)/max(norm(residue),eps(Float64))
            push!(residuecheck,(s_abs=ab,direction=d,residue_relative_error=err,
                residue_norm=norm(residue),sPi_norm=norm(val)))
        end
    end
    top=isempty(nullrows) ? String[] : [r.state_name for r in nullrows
        if r.singular_vector=="right_1" && r.abs_weight>=1e-6]
    state="UNRESOLVED"
    if n0>0 && semisimple && length(top)>0 && !isempty(residuecheck) && matrix_p95(getproperty.(residuecheck,:residue_relative_error))<0.1
        state="EXPLAINED"
    elseif n0>0 && !isempty(scalings)
        state="PARTIALLY_EXPLAINED"
    end
    return (nullities=nullities,null_rows=nullrows,scalings=scalings,
        residue_checks=residuecheck,nullity=n0,semisimple=semisimple,
        semisimple_residual=semi_res,rank_tolerance=τ,residue=residue,observed_power=slope,
        status=state,dominant_states=top,projector=P)
end

function power_identity(op,basis,LG,s;z=nothing)
    ω=imag(s)
    z===nothing && (z=ComplexF64.(collect(1:length(basis.Lambda)) .+ im .* reverse(collect(1:length(basis.Lambda)))))
    Psi=GDS.psi(op,LG,s)
    q=basis.Phi*z
    v=im*ω*q
    p_direct=real(dot(v,Psi*q))/2
    Shat=basis.Phi'*Psi*basis.Phi
    DG=IP.effective_dissipative_operator(Shat,ω)
    p_graph=ω^2*real(dot(z,DG*z))/2
    return abs(p_direct-p_graph)/max(1.0,abs(p_direct),abs(p_graph))
end

function synthetic_sigma_special_case()
    cases=BNDGraph.SyntheticSystems.build_synthetic_cases()
    sys=first(filter(c->c.name=="S2",cases)).modal_sys
    basis=first(filter(c->c.name=="S2",cases)).basis
    σ=BNDGraph.DynamicSelfEnergy.sigma
    errs=Float64[]
    for f in (0.2,0.8,2.1,4.5,12.0)
        ω=2pi*f
        S=σ(sys,im*ω)
        Shat=basis.Phi'*S*basis.Phi
        Pih=im*ω*Shat
        DG=IP.effective_dissipative_operator(Pih,ω)
        push!(errs,norm(DG-(Shat+Shat')/2)/max(norm(Shat),eps(Float64)))
    end
    return maximum(errs)
end

function predictor_class(rows)
    finite=filter(r->isfinite(r.rel_error),rows)
    isempty(finite) && return "BLOCKED"
    low=filter(r->r.eta_shift<0.01,finite)
    mid=filter(r->0.01<=r.eta_shift<0.05,finite)
    if !isempty(low) && median(getproperty.(low,:rel_error))<=CONFIG["predictor_small_error_relative"]
        larger=filter(r->r.eta_shift>=0.01,finite)
        if !isempty(larger) && median(getproperty.(larger,:rel_error))>median(getproperty.(low,:rel_error))
            return "SUPPORTED"
        end
        return "LIMITED"
    end
    if !isempty(finite) && all(r->r.rel_error>0.05 || r.graph_concentration<0.5 || r.Trr_sigma_min<1e-8,finite)
        return "FALSIFIED"
    end
    return "LIMITED"
end

function predictor_bins(rows)
    bins=(("eta_lt_0.01",r->r.eta_shift<0.01),
          ("eta_0.01_to_0.05",r->0.01<=r.eta_shift<0.05),
          ("eta_0.05_to_0.10",r->0.05<=r.eta_shift<0.10),
          ("eta_ge_0.10",r->r.eta_shift>=0.10))
    out=Dict{String,Any}()
    for (name,pred) in bins
        group=filter(r->isfinite(r.rel_error)&&pred(r),rows)
        out[name]=Dict("count"=>length(group),"median_relative_error"=>isempty(group) ? nothing : median(getproperty.(group,:rel_error)),
            "max_relative_error"=>isempty(group) ? nothing : maximum(getproperty.(group,:rel_error)))
    end
    return out
end

function run()
    expApath=joinpath(ROOT,"reports","experiment_A","RESULTS_EXP_A.json")
    expBpath=joinpath(ROOT,"reports","experiment_B","RESULTS_EXP_B.json")
    expA=read(expApath,String); expB=read(expBpath,String)
    expApass=occursin("\"experiment\":\"BND_EXP_A\"",expA)&&occursin("\"status\":\"PASS\"",expA)
    expBpass=occursin("\"experiment\":\"BND_EXP_B_PRE\"",expB)&&occursin("\"status\":\"PASS\"",expB)
    expAprimary=RSA.load_expA_primary(ROOT)
    rebuilt=load_all_case_caches()
    if rebuilt===nothing || "--rebuild-models" in ARGS
        rebuilt=RSA.build_ieee39_cases(PD39,PowerDynamics,NetworkDynamics;
            on_case=(id,c)->save_case_cache(id,c;regenerated=(id=="C33")))
    end
    # ExpA's frozen bus-33 Ared and state map remain the primary bus-33 input.
    # A same-model reconstruction is used only to check input reproducibility.
    fresh=rebuilt["C33"]
    a3err=norm(fresh.Ared-expAprimary.Ared)/max(norm(expAprimary.Ared),eps(Float64))
    cases=Dict{String,Any}("C0"=>rebuilt["C0"])
    for bus in (30,35,37); cases["C$(bus)"]=rebuilt["C$(bus)"]; end
    cases["C33"]=RSA.case_from_reduced("C33_bus33_SimpleGFLDC",33,
        expAprimary.Ared,expAprimary.state_names)
    contexts=Dict{String,Any}()
    for id in ("C0","C30","C33","C35","C37")
        contexts[id]=case_analysis(cases[id])
    end
    all(ctx->isempty(ctx.blocked),values(contexts)) || error("At least one case has a blocked M audit")
    derivative_rows=NamedTuple[]
    unit_rows=NamedTuple[]
    for (id,ctx) in contexts
        for s in (0.3+0.9im,-0.2+1.7im,0.1+4.1im)
            condcc=cond(s*I-ctx.op.Acc)
            err=condcc<1e12 ? GDS.derivative_check(ctx.op,s) : NaN
            push!(derivative_rows,(case=ctx.case.label,s_real=real(s),s_imag=imag(s),
                cond_Tcc=condcc,relative_derivative_error=err,accepted=condcc<1e12))
        end
        for (j,pair) in enumerate(ctx.case.partition.pairs)
            gfl=pair.device=="gfl"
            push!(unit_rows,(case=ctx.case.label,coordinate_index=j,bus=pair.bus,device=pair.device,
                q_state=ctx.case.state_names[pair.index],v_state=ctx.case.state_names[ctx.case.partition.v[j]],
                q_unit="rad",v_unit=gfl ? "rad/s" : "pu on 60-Hz base",
                Ckin_value=ctx.case.partition.Ckin[j,j],
                Ckin_unit=gfl ? "(rad/s)/(rad/s)" : "(rad/s)/pu",
                M_value=ctx.op.M[j,j],
                M_unit=gfl ? "inverse PLL angle-rate conversion" : "inverse synchronous-speed angle-rate conversion",
                interpretation="second-order metric from qdot=Ckin*v; not identified as physical inertia"))
        end
    end
    cderivative=csv_write("TABLE_C19_derivative_validation.csv",derivative_rows)
    csv_write("TABLE_C20_metric_units.csv",unit_rows)
    dvals=filter(isfinite,Float64.(cderivative.relative_derivative_error))
    dmedian=isempty(dvals) ? NaN : median(dvals)
    dp95=isempty(dvals) ? NaN : quantile(dvals,0.95)

    # Provenance, source hashes, and exact input inventory.
    inputfiles=[expApath,expBpath,joinpath(ROOT,"reports","experiment_A","REPORT_EXP_A.md"),
        joinpath(ROOT,"reports","experiment_A","CLAIM_LEDGER_EXP_A.md"),
        joinpath(ROOT,"reports","experiment_A","TABLES_EXP_A.md"),
        joinpath(ROOT,"reports","experiment_B","REPORT_EXP_B.md"),
        joinpath(ROOT,"reports","experiment_B","CLAIM_LEDGER_EXP_B.md"),
        joinpath(ROOT,"reports","experiment_B","EXPA_INTERFACE_SPEC.md"),
        joinpath(ROOT,"reports","experiment_A","matrices","bus33_A_reduced.csv"),
        joinpath(ROOT,"reports","experiment_A","tables","TABLE_A03_state_partition.csv")]
    provenance=NamedTuple[]
    for path in inputfiles
        push!(provenance,(field=basename(path),value=isfile(path) ? sha(path) : "missing",
            source=relpath(path,ROOT),hash_or_revision=isfile(path) ? sha(path) : "missing"))
    end
    versions=RSA.model_versions(PowerDynamics,NetworkDynamics)
    for (field,value,source) in (("ExpA status",expApass ? "PASS" : "FAIL","RESULTS_EXP_A.json"),
        ("ExpB-pre status",expBpass ? "PASS" : "FAIL","RESULTS_EXP_B.json"),
        ("Julia version",string(VERSION),"Julia runtime"),
        ("PowerDynamics version",versions.powerdynamics,"active Manifest.toml"),
        ("NetworkDynamics version",versions.networkdynamics,"active Manifest.toml"),
        ("Git HEAD",GIT_HEAD,"git rev-parse HEAD"),
        ("BLAS",sprint(show,BLAS.get_config()),"LinearAlgebra.BLAS.get_config()"),
        ("C33 rebuilt Ared relative discrepancy",string(a3err),"ExpA export versus same-model reconstruction"))
        push!(provenance,(field=field,value=value,source=source,hash_or_revision=""))
    end
    csv_write("TABLE_C01_input_provenance.csv",provenance)

    # Mandatory zero-frequency diagnostic uses the exact ExpA condensed block.
    zero=zero_frequency_diagnostic(contexts["C33"])
    csv_write("TABLE_C02_zero_singularity_states.csv",zero.null_rows)
    csv_write("TABLE_C03_zero_frequency_scaling.csv",zero.scalings)
    open(joinpath(OUT,"ZERO_FREQUENCY_NULLITIES.csv"),"w") do io
        CSV.write(io,as_dataframe(zero.nullities))
    end

    backbone_rows=NamedTuple[]
    graphrows=NamedTuple[]
    metricrows=NamedTuple[]
    chigrows=NamedTuple[]
    gammarows=NamedTuple[]
    rootrows=NamedTuple[]
    predrows=NamedTuple[]
    pathwayrows=NamedTuple[]
    blockpathrows=NamedTuple[]
    susceptibilityrows=NamedTuple[]
    polemaprows=NamedTuple[]
    frequency_analyses=Dict{String,Any}()
    root_analyses=Dict{String,Any}()
    allcases=Dict{String,Any}()
    for id in ("C0","C30","C33","C35","C37")
        ctx=contexts[id]
        allcases[id]=Dict("M"=>ctx.op.M,"D0"=>ctx.op.D0,"L0"=>ctx.op.L0,
            "LG_primary"=>ctx.primary.LG,"LG_secondary"=>ctx.secondary.LG,
            "Phi"=>ctx.basis.Phi,"Lambda"=>ctx.basis.Lambda,
            "state_names"=>ctx.case.state_names)
        matrix_csv(joinpath(MATRICES,"$(id)_M.csv"),ctx.op.M)
        matrix_csv(joinpath(MATRICES,"$(id)_D0.csv"),ctx.op.D0)
        matrix_csv(joinpath(MATRICES,"$(id)_L0.csv"),ctx.op.L0)
        matrix_csv(joinpath(MATRICES,"$(id)_LG_primary.csv"),ctx.primary.LG)
        matrix_csv(joinpath(MATRICES,"$(id)_Phi_primary.csv"),ctx.basis.Phi)
        matrix_csv(joinpath(MATRICES,"$(id)_LG_secondary.csv"),ctx.secondary.LG)
        matrix_csv(joinpath(MATRICES,"$(id)_Phi_secondary.csv"),ctx.basis2.Phi)
        for (which,LG,bb) in (("primary",ctx.primary.LG,ctx.primary),("secondary",ctx.secondary.LG,ctx.secondary))
            bas=which=="primary" ? ctx.basis : ctx.basis2
            diag=GB.mass_audit(ctx.op.M)
            cls=isempty(bas.negative_modes) ? "LAPLACIAN_LIKE" : "INDEFINITE"
            for k in eachindex(bas.Lambda)
                push!(backbone_rows,(case=ctx.case.label,backbone=which,M_symmetry_error=diag.symmetry_error,
                    M_lambda_min=diag.lambda_min,gauge_residual=bb.gauge_residual,lambda_index=k,
                    backbone_eigenvalue=bas.Lambda[k],classification=cls))
            end
        end
        f=eval_frequency_tables(ctx)
        append!(graphrows,f.graph); append!(metricrows,f.metrics); append!(chigrows,f.chig)
        frequency_analyses[id]=f
        gr=select_gamma_rows(ctx); append!(gammarows,gr)
        ra=roots_and_predictors(ctx,gr)
        append!(rootrows,ra.roots); append!(predrows,ra.predictors)
        append!(pathwayrows,ra.pathways); append!(susceptibilityrows,ra.susceptibilities)
        append!(blockpathrows,ra.block_pathways)
        root_analyses[id]=ra
        for p in ctx.poles.selected
            sorted=sortperm(p.projection;rev=true)
            push!(polemaprows,(case=ctx.case.label,pole_id="eig$(p.index)",
                lambda_real=real(p.lambda),lambda_imag=imag(p.lambda),frequency_hz=p.frequency_hz,
                damping_ratio=p.damping_ratio,retained_participation=p.retained_participation,
                dominant_graph_mode=p.dominant_graph_mode,graph_concentration=p.graph_concentration,
                second_graph_mode=length(sorted)>1 ? sorted[2] : missing,top1_top2_ratio=p.top1_top2,
                gauge_score=p.gauge_score,pll_state_share=p.pll_state_share,dominant_states=p.dominant_states,
                role=p.role,q_norm=p.q_norm))
        end
    end
    c05=csv_write("TABLE_C05_graph_identity_residuals.csv",graphrows)
    c06=csv_write("TABLE_C06_Psi_frequency_metrics.csv",metricrows)
    c07=csv_write("TABLE_C07_physical_pole_graph_mapping.csv",polemaprows)
    c08=csv_write("TABLE_C08_Gamma_exact.csv",gammarows)
    c09=csv_write("TABLE_C09_uncoupled_modal_roots.csv",rootrows)
    c10=csv_write("TABLE_C10_pole_predictor.csv",predrows)
    c11=csv_write("TABLE_C11_exact_pathways.csv",pathwayrows)
    c12=csv_write("TABLE_C12_pairwise_susceptibility.csv",susceptibilityrows)
    c21=csv_write("TABLE_C21_degenerate_block_pathways.csv",blockpathrows)
    c04=csv_write("TABLE_C04_backbone_diagnostics.csv",backbone_rows)

    # Exact harmonic power identity and ExpB Sigma special-case recovery.
    powererrs=Dict{String,Float64}()
    for id in keys(contexts)
        powererrs[id]=power_identity(contexts[id].op,contexts[id].basis,
            contexts[id].primary.LG,2pi*0.8im)
    end
    sigma_special_error=synthetic_sigma_special_case()

    # Baseline-to-replacement comparison, with cross-case modes paired by
    # normalized Mref overlap in the common physical q-coordinate ordering.
    base=contexts["C0"]
    comparison=NamedTuple[]
    basep=base.poles.selected
    for id in ("C30","C33","C35","C37")
        ctx=contexts[id]
        overlap=GB.mode_overlap(base.basis.Phi,base.op.M,ctx.basis.Phi,ctx.op.M)
        for p in ctx.poles.selected
            candidates=filter(r->r.role==p.role,basep)
            isempty(candidates) && (candidates=basep)
            qnew=p.right_vector[ctx.case.partition.q]
            score(r)=begin
                qold=r.right_vector[base.case.partition.q]
                abs(dot(qold,qnew))/max(norm(qold)*norm(qnew),eps(Float64))
            end
            bmatch=isempty(candidates) ? nothing : candidates[argmax(score.(candidates))]
            bmode=bmatch===nothing ? 0 : bmatch.dominant_graph_mode
            kmatch=bmode==0 ? 0 : argmax(overlap[bmode,:])
            λb=bmatch===nothing ? complex(NaN,NaN) : bmatch.lambda
            # All quantities are evaluated at an offset from the actual full
            # pole; this avoids evaluating the singular full operator itself.
            function pole_mechanism(cctx,pp)
                k=pp.dominant_graph_mode
                v=eval_graph(cctx,pp.lambda+1e-4*(1+im))
                γ=IP.schur_gamma(v.T_hat,k).gamma
                χ,_=IP.commutator_metric(cctx.basis.Lambda,v.Psihat)
                ω=max(abs(imag(pp.lambda)),eps(Float64))
                DG=IP.effective_dissipative_operator(v.Psihat,ω)
                pw=IP.pathway_matrix(v.Psihat,cctx.basis.Lambda,pp.lambda+1e-4*(1+im),k)
                ij=argmax(abs.(pw.G))
                lm=Tuple(CartesianIndices(pw.G)[ij])
                pathway="$(pw.modes[lm[1]])←$(pw.modes[lm[1]]),$(pw.modes[lm[2]])"
                return (gamma=γ,psi_kk=v.Psihat[k,k],DG_kk=real(DG[k,k]),
                    chi_comm=χ,dominant_pathway=pathway)
            end
            db=bmatch===nothing ? nothing : pole_mechanism(base,bmatch)
            dg=pole_mechanism(ctx,p)
            γb=db===nothing ? complex(NaN,NaN) : db.gamma
            γg=dg.gamma
            push!(comparison,(replacement_bus=ctx.case.bus,matched_mode="$(p.role):eig$(p.index)",
                baseline_lambda=λb,replacement_lambda=p.lambda,delta_real=real(p.lambda-λb),
                delta_imag=imag(p.lambda-λb),baseline_gamma=γb,replacement_gamma=γg,
                delta_gamma=γg-γb,graph_mode_overlap=bmode==0 ? NaN : overlap[bmode,kmatch],
                baseline_graph_mode=bmode,replacement_graph_mode=p.dominant_graph_mode,
                q_shape_overlap=bmatch===nothing ? NaN : score(bmatch),
                baseline_psi_kk=db===nothing ? complex(NaN,NaN) : db.psi_kk,
                replacement_psi_kk=dg.psi_kk,delta_psi_kk=dg.psi_kk-(db===nothing ? complex(NaN,NaN) : db.psi_kk),
                baseline_DG_eff_kk=db===nothing ? NaN : db.DG_kk,replacement_DG_eff_kk=dg.DG_kk,
                baseline_chi_comm=db===nothing ? NaN : db.chi_comm,replacement_chi_comm=dg.chi_comm,
                baseline_dominant_pathway=db===nothing ? "" : db.dominant_pathway,
                replacement_dominant_pathway=dg.dominant_pathway))
        end
    end
    csv_write("TABLE_C13_baseline_vs_replacement.csv",comparison)

    # Secondary-backbone sensitivity for the same actual pole eigenvectors.
    sensrows=NamedTuple[]
    for id in ("C0","C30","C33","C35","C37")
        ctx=contexts[id]
        ov=GB.mode_overlap(ctx.basis.Phi,ctx.op.M,ctx.basis2.Phi,ctx.op.M)
        for p in ctx.poles.selected
            k=p.dominant_graph_mode
            l=argmax(ov[k,:])
            s=p.lambda+1e-4*(1+im)
            v1=eval_graph(ctx,s)
            psi2=ctx.basis2.Phi'*GDS.psi(ctx.op,ctx.secondary.LG,s)*ctx.basis2.Phi
            T2=IP.graph_operator(ctx.basis2.Lambda,psi2,s)
            γ1=IP.schur_gamma(v1.T_hat,k).gamma
            γ2=IP.schur_gamma(T2,l).gamma
            push!(sensrows,(case=ctx.case.label,physical_mode="eig$(p.index)",primary_graph_mode=k,
                secondary_graph_mode=l,subspace_overlap=ov[k,l],gamma_primary=abs(γ1),
                gamma_secondary=abs(γ2),predictor_error_primary=begin
                    rr=filter(r->r.case==ctx.case.label&&r.graph_mode==k,predrows)
                    isempty(rr) ? NaN : first(rr).rel_error
                end,predictor_error_secondary=NaN,
                primary_concentration=p.graph_concentration,
                secondary_concentration=RP.graph_projection(ctx.basis2.Phi,ctx.op.M,p.right_vector[ctx.case.partition.q]).concentration))
        end
    end
    csv_write("TABLE_C14_backbone_sensitivity.csv",sensrows)

    # Operator ablation is an algebraic diagnostic, not a physical plant.
    ablation=NamedTuple[]
    for id in ("C0","C30","C33","C35","C37")
        ctx=contexts[id]
        for f in (0.1,0.5,1.0,2.0,5.0)
            s=2pi*f*im
            op=ctx.op; b=ctx.basis
            parts=("Psi_D",s*op.D0),("Psi_DeltaL",s*op.D0+(op.L0-ctx.primary.LG)),("Psi_full",GDS.psi(op,ctx.primary.LG,s))
            for (name,Psi) in parts
                Phat=b.Phi'*Psi*b.Phi
                off=Phat-Diagonal(diag(Phat))
                k=argmax(abs.(diag(Phat)))
                T=IP.graph_operator(b.Lambda,Phat,s)
                γ=IP.schur_gamma(T,k).gamma
                push!(ablation,(case=ctx.case.label,frequency_hz=f,operator=name,
                    offdiag_ratio=norm(off)/(norm(Phat)+eps(Float64)),gamma_abs=abs(γ),
                    dynamic_component_included=name=="Psi_full"))
            end
        end
    end
    csv_write("TABLE_C17_component_ablation.csv",ablation)
    chirows=[(case=r.case,frequency_hz=r.frequency_hz,chi_G=r.chi_G,
        sigma_min_TD=r.sigma_min_TD,condition_TD=r.condition_TD,
        classification=r.classification) for r in chigrows]
    csv_write("TABLE_C18_candidate_chiG.csv",chirows)

    # Per-case exactness and predictor status.
    cross=NamedTuple[]
    for id in ("C0","C30","C33","C35","C37")
        ctx=contexts[id]
        gres=filter(:case => ==(ctx.case.label),c05).reconstruction_error
        grrows=filter(:case => ==(ctx.case.label),c08)
        gr=grrows.schur_residual
        pr=filter(r->r.case==ctx.case.label,predrows)
        pstatus=predictor_class(pr)
        prfinite=filter(r->isfinite(r.rel_error),pr)
        push!(cross,(case=id,input_pass=id=="C33" ? expApass : true,
            graph_basis_pass=ctx.basis.M_orthogonality_error<1e-10&&ctx.basis.L_diagonalization_error<1e-10,
            graph_identity_pass=!isempty(gres)&&matrix_p95(gres)<1e-9,
            Gamma_identity_pass=!isempty(gr)&&matrix_p95(gr)<1e-9,
            number_analyzed_modes=length(ctx.poles.selected),predictor_status=pstatus,
            median_predictor_error=isempty(prfinite) ? NaN : median(getproperty.(prfinite,:rel_error)),
            max_predictor_error=isempty(prfinite) ? NaN : maximum(getproperty.(prfinite,:rel_error)),
            dominant_intermodal_mode=isempty(grrows) ? missing : grrows[argmax(grrows.gamma_abs),:graph_mode_k]))
    end
    csv_write("TABLE_C15_cross_bus_summary.csv",cross)

    # Spearman test: compare susceptibility ranking with exact diagonal paths.
    correlations=Dict{String,Any}()
    for id in ("C0","C30","C33","C35","C37")
        ss=filter(r->r.case==contexts[id].case.label,susceptibilityrows)
        groups=unique(getproperty.(ss,:mode_k))
        cc=Float64[]; sc=Float64[]
        for k in groups
            group=filter(r->r.mode_k==k,ss)
            length(group)>=3 || continue
            push!(cc,RP.spearman_rank(getproperty.(group,:coupling_product_abs),getproperty.(group,:exact_diagonal_path_abs)))
            push!(sc,RP.spearman_rank(getproperty.(group,:susceptibility_abs),getproperty.(group,:exact_diagonal_path_abs)))
        end
        correlations[id]=Dict("coupling_rank_correlation"=>isempty(cc) ? nothing : mean(cc),
            "susceptibility_rank_correlation"=>isempty(sc) ? nothing : mean(sc),
            "mode_count"=>length(cc))
    end
    primarycorr=correlations["C33"]
    detune=if primarycorr["coupling_rank_correlation"]===nothing || primarycorr["susceptibility_rank_correlation"]===nothing
        "INCONCLUSIVE"
    elseif primarycorr["susceptibility_rank_correlation"]-primarycorr["coupling_rank_correlation"] >= CONFIG["susceptibility_rank_improvement_min"]
        "SUPPORTED"
    else
        "NOT_OBSERVED"
    end

    # Basis dependence status uses overlap, critical-mode mapping, and Gamma rank.
    allangles=Float64[]; gamma_ratios=Float64[]
    for id in keys(contexts)
        ctx=contexts[id]
        nlow=min(3,size(ctx.basis.Phi,2))
        angles=GB.principal_angles(ctx.basis.Phi[:,1:nlow],ctx.basis2.Phi[:,1:nlow])
        append!(allangles,angles)
        for p in ctx.poles.selected
            k=p.dominant_graph_mode; l=argmax(GB.mode_overlap(ctx.basis.Phi[:,k:k],ctx.op.M,
                ctx.basis2.Phi,ctx.op.M)[1,:])
            v=eval_graph(ctx,p.lambda+1e-4*(1+im))
            γ1=abs(IP.schur_gamma(v.T_hat,k).gamma)
            psi2=ctx.basis2.Phi'*GDS.psi(ctx.op,ctx.secondary.LG,p.lambda+1e-4*(1+im))*ctx.basis2.Phi
            γ2=abs(IP.schur_gamma(IP.graph_operator(ctx.basis2.Lambda,psi2,p.lambda+1e-4*(1+im)),l).gamma)
            push!(gamma_ratios,max(γ1,γ2)/max(min(γ1,γ2),eps(Float64)))
        end
    end
    robustness=if isempty(allangles)
        "BASIS_SENSITIVE"
    elseif maximum(allangles)<=deg2rad(CONFIG["backbone_principal_angle_max_deg"]) &&
           (isempty(gamma_ratios)||maximum(gamma_ratios)<=CONFIG["backbone_gamma_ratio_max"])
        "ROBUST"
    elseif maximum(allangles)<=deg2rad(45.0)
        "PARTIALLY_ROBUST"
    else
        "BASIS_SENSITIVE"
    end

    # Gate ledger: exact algebraic checks and empirical criteria stay separate.
    psi_errors=Float64[]
    for ctx in values(contexts), s in (0.2+0.9im,-0.15+2.2im,1e-3+0.3im)
        push!(psi_errors,GDS.generalized_reconstruction_error(ctx.op,ctx.primary.LG,s))
    end
    psi_med=median(psi_errors); psi_p95=matrix_p95(psi_errors)
    graph_p95=matrix_p95(c05.reconstruction_error)
    gamma_p95=isempty(c08.schur_residual) ? NaN : matrix_p95(c08.schur_residual)
    pathway_max=isempty(c11.reconstruction_error) ? NaN : maximum(c11.reconstruction_error)
    predictor_status=predictor_class(predrows)
    pstatusvec=filter(r->isfinite(r.rel_error),predrows)
    predmed=isempty(pstatusvec) ? NaN : median(getproperty.(pstatusvec,:rel_error))
    predmax=isempty(pstatusvec) ? NaN : maximum(getproperty.(pstatusvec,:rel_error))
    c0pass=expApass&&expBpass&&a3err<1e-8
    c1pass=psi_med<1e-11&&psi_p95<1e-9&&dmedian<1e-7&&dp95<1e-5
    c3pass=all(ctx->ctx.mass.spd&&ctx.basis.M_orthogonality_error<1e-10&&ctx.basis.L_diagonalization_error<1e-10,values(contexts))
    c4pass=graph_p95<1e-9
    c5pass=gamma_p95<1e-9
    c6pass=pathway_max<1e-10
    c9pass=maximum(values(powererrs))<1e-10&&sigma_special_error<1e-10
    c11pass=all(r->r.graph_basis_pass&&r.graph_identity_pass&&r.Gamma_identity_pass,cross)
    ready=c0pass&&c1pass&&c3pass&&c4pass&&c5pass&&c6pass&&c11pass
    expstatus=ready ? (predictor_status=="SUPPORTED" ? "PASS" : "PASS_LIMITED_PREDICTOR") : "PARTIAL"
    gates=NamedTuple[
        (gate="C0 input integrity",metric="A/B status, hashes, same-model bus-33 Ared reproduction",threshold_or_classification="PASS",observed="A=$(expApass), B=$(expBpass), Ared_rel=$(a3err)",status=c0pass ? "PASS" : "BLOCKED",notes="A/B output files were read only."),
        (gate="C1 generalized Psi exactness",metric="median / p95 reconstruction; Pi derivative median / p95",threshold_or_classification="<1e-11 / <1e-9; <1e-7 / <1e-5",observed="$(psi_med) / $(psi_p95); $(dmedian) / $(dp95)",status=c1pass ? "PASS" : "FAIL",notes="Complex sample points for every case; finite differences avoid Tcc poles."),
        (gate="C2 zero-frequency singularity",metric="nullspace, semisimplicity, low-frequency scaling and residue",threshold_or_classification="EXPLAINED / PARTIALLY_EXPLAINED / UNRESOLVED",observed="nullity=$(zero.nullity), semisimple=$(zero.semisimple), p=$(zero.observed_power)",status=zero.status,notes="Threshold sensitivity prevents a Jordan/residue claim; no Pi(0) or pseudoinverse."),
        (gate="C3 graph backbone and basis",metric="M SPD and generalized basis residuals",threshold_or_classification="<1e-10",observed="Mmin=$(contexts["C33"].mass.lambda_min), eM=$(contexts["C33"].basis.M_orthogonality_error), eL=$(contexts["C33"].basis.L_diagonalization_error)",status=c3pass ? "PASS" : "BLOCKED",notes="Primary backbone fixed by the configured projection rule."),
        (gate="C4 exact graph transform",metric="median / p95 reconstruction error",threshold_or_classification="<1e-11 / <1e-9",observed="p95=$(graph_p95)",status=c4pass ? "PASS" : "FAIL",notes="Imaginary-axis grid plus offsets near analyzed poles."),
        (gate="C5 exact Gamma Schur identity",metric="median / p95 scalar residual",threshold_or_classification="<1e-11 / <1e-9",observed="p95=$(gamma_p95)",status=c5pass ? "PASS" : "FAIL",notes="Complementary solves, no explicit inverse."),
        (gate="C6 exact pathway reconstruction",metric="max relative pathway sum error",threshold_or_classification="<1e-10",observed=string(pathway_max),status=c6pass ? "PASS" : "FAIL",notes="Per-mode exact pathway matrix."),
        (gate="C7 real predictor",metric="matched exact pole errors stratified by eta",threshold_or_classification="SUPPORTED / LIMITED / FALSIFIED / BLOCKED",observed="$(predictor_status), n=$(length(pstatusvec)), median=$(predmed)",status=predictor_status,notes="No predictor is called exact."),
        (gate="C8 real detuning mechanism",metric="Spearman ranking against exact diagonal pathway",threshold_or_classification="SUPPORTED / NOT_OBSERVED / INCONCLUSIVE",observed=string(detune),status=detune,notes="Susceptibility is approximate when complement modes are coupled."),
        (gate="C9 generalized dissipative power identity",metric="direct versus graph harmonic power; Sigma special-case",threshold_or_classification="<1e-10",observed="$(maximum(values(powererrs))); special=$(sigma_special_error)",status=c9pass ? "PASS" : "FAIL",notes="Hermitian imaginary part, not Hermitian{Psi}."),
        (gate="C10 backbone robustness",metric="principal-angle and Gamma sensitivity",threshold_or_classification="ROBUST / PARTIALLY_ROBUST / BASIS_SENSITIVE",observed=robustness,status=robustness,notes="Primary and Euclidean projected Hermitian choices."),
        (gate="C11 cross-bus reproducibility",metric="graph, Gamma, Psi identities across cases",threshold_or_classification="all feasible cases",observed="$(count(r->r.graph_identity_pass&&r.Gamma_identity_pass,cross))/$(length(cross))",status=c11pass ? "PASS" : "BLOCKED",notes="C0 baseline plus 30, 33, 35 and 37."),
        (gate="C12 ready for ExpD",metric="C0,C1,C3,C4,C5,C6,C11",threshold_or_classification="YES only if every required gate passes",observed=string(ready),status=ready ? "YES" : "NO",notes="Controller variation must still be exposed without semantic changes.")]
    csv_write("TABLE_C16_gate_summary.csv",gates)

    # Graph-space matrices and exact pathway objects are retained for audit.
    for (id,ctx) in contexts
        s=0.3+1.1im
        v=eval_graph(ctx,s)
        matrix_csv(joinpath(MATRICES,"$(id)_Psihat_sample.csv"),v.Psihat)
        matrix_csv(joinpath(MATRICES,"$(id)_T_hat_sample.csv"),v.T_hat)
    end
    if !isempty(pathwayrows)
        primaryrows=filter(r->r.case==contexts["C33"].case.label,pathwayrows)
        !isempty(primaryrows) && CSV.write(joinpath(MATRICES,"C33_primary_pathways.csv"),as_dataframe(primaryrows))
    end

    # Claim ledger is machine-backed by the computed gates.
    claims=[
        ("C-C01","Tq=s²M+LG+Psi is exactly equivalent to the ExpA retained operator.","EXACT",c1pass ? "SUPPORTED" : "BLOCKED","Algebraic construction from exact ExpA blocks."),
        ("C-C02","A regular Sigma(s) is a special case when Psi(s)/s is regular.","EXACT mathematical relationship","SUPPORTED","Tq=s²M+LG+Psi; Psi=sSigma recovers the B framework."),
        ("C-C03","The real graph transform preserves the generalized NEP.","EXACT",c4pass ? "SUPPORTED" : "BLOCKED","Phiᴴ Tq Phi=s²I+Lambda+Psihat."),
        ("C-C04","Gamma_k is the exact intermodal Schur self-energy in real graph coordinates.","EXACT",c5pass ? "SUPPORTED" : "BLOCKED","Exact complementary-mode solve and scalar Schur identity."),
        ("C-C05","Exact pathway entries sum to Gamma_k.","EXACT",c6pass ? "SUPPORTED" : "BLOCKED","Reconstruction checked at each accepted root."),
        ("C-C06","Condensed IEEE-39 dynamics produce non-negligible graph-mode mixing.","EMPIRICAL",maximum(c06.offdiag_ratio)>0.05 ? "SUPPORTED" : "NOT_OBSERVED","See TABLE_C06 and commutator traces."),
        ("C-C07","Complementary-mode detuning amplifies real intermodal self-energy.","EMPIRICAL",detune,"Rank comparison against exact diagonal pathways."),
        ("C-C08","The perturbative predictor explains diagonal-root to full-pole shift in a weak-coupling subset.","EMPIRICAL / ASYMPTOTIC",predictor_status,"Simple-root expansion; error, eta, complement conditioning, and graph concentration retained."),
        ("C-C09","DGeff=(1/omega) Im_H{Psihat(jomega)} is the generalized dissipative operator under the declared convention.","EXACT algebraic identity",c9pass ? "SUPPORTED" : "FAIL","Direct harmonic force-velocity power and synthetic Sigma special case."),
        ("C-C10","Approximate susceptibility ranks interactions better than coupling magnitude alone.","EMPIRICAL",detune,"Spearman correlations against exact diagonal pathway magnitude."),
        ("C-C11","Results are robust to reasonable non-fitted backbones.","EMPIRICAL",robustness,"Primary and secondary backbone comparison; see BACKBONE_AUDIT.md.")]
    open(joinpath(OUT,"CLAIM_LEDGER_EXP_C.md"),"w") do io
        println(io,"# Claim ledger — Experiment C\n\n| ID | Claim | Type | Evidence / derivation | Status |\n|---|---|---|---|---|")
        for (id,claim,typ,status,evidence) in claims
            println(io,"| $id | $claim | $typ | $evidence | $status |")
        end
    end

    # Report narrative uses measured values only; unknowns stay explicit.
    primary=contexts["C33"]
    critical=filter(r->occursin("spectral_abscissa",r.role),primary.poles.selected)
    critpole=isempty(critical) ? nothing : critical[argmax(real.([r.lambda for r in critical]))]
    critmode=critpole===nothing ? missing : critpole.dominant_graph_mode
    critpoleid=critpole===nothing ? "" : "eig$(critpole.index)"
    criticalpaths=filter(r->r.case==primary.case.label && r.pole_id==critpoleid &&
        r.evaluation_kind=="full_pole_offset" && r.diagonal_path,pathwayrows)
    dompaths=isempty(criticalpaths) ? NamedTuple[] :
        sort(criticalpaths;by=r->-r.contribution_abs)[1:min(10,length(criticalpaths))]
    open(joinpath(OUT,"REPORT_EXP_C.md"),"w") do io
        println(io,"# Experiment C — Real IEEE-39 Graph Dynamic Self-Energy\n")
        println(io,"## 1. Executive result\n\n**EXP_C_STATUS: $expstatus**\n")
        println(io,"This run uses the frozen bus-33 Experiment-A `A_reduced`/state map as its primary real input, audits it against a same-model reconstruction (relative discrepancy $a3err), and rebuilds the all-SG baseline plus nominal single-GFL replacements at buses 30, 35, and 37. Controller parameters are unchanged.\n")
        println(io,"## 2. Relation to Experiments A and B\n\nExpA input status: **$(expApass ? "PASS" : "FAIL")**. ExpB-pre status: **$(expBpass ? "PASS" : "FAIL")**. ExpB's tested generalized graph basis, exact Schur tools, and synthetic self-energy special case were reused; no A/B output was edited.\n")
        println(io,"## 3. Why Sigma(s) is not globally available\n\nThe condensed block at zero remains singular. ExpC does not evaluate `Pi_q(0)`, add an epsilon, or use a pseudoinverse. It works directly with the exact generalized object `Psi(s)=sD0+(L0-LG)+Pi_q(s)`.\n")
        println(io,"## 4. Generalized BND operator\n\nThe construction is `Tq(s)=s²M+LG+Psi(s)`. Across all five cases, the sample-point median/p95 algebraic reconstruction errors are $psi_med / $psi_p95. The analytic `Pi_q'` derivative check has median/p95 error $dmedian / $dp95; coordinate conversions and units are in TABLE_C20.\n")
        println(io,"## 5. Zero-frequency condensed singularity\n\nClassification: **$(zero.status)**. At relative tolerance $(zero.rank_tolerance), `nullity(A_cc)=$(zero.nullity)`; the tested `A_cc²` nullity and threshold sensitivity are in `ZERO_FREQUENCY_NULLITIES.csv`. Since nullities do not agree robustly across tolerances, no semisimplicity or residue claim is made. The fitted empirical low-frequency norm slope is $(zero.observed_power), over the predeclared conditioning-limited range. Dominant named right-null states: $(join(zero.dominant_states,", ")).\n")
        println(io,"## 6. Graph-backbone definition and canonicality\n\nPrimary `LG` is the preregistered M-normalized Hermitian part of `L0`, projected onto the relative-angle subspace. The secondary backbone is the Euclidean gauge-preserving Hermitian part of `L0`. Neither was fitted to pole data.\n")
        println(io,"## 7. Generalized graph basis\n\nFor C33: minimum eigenvalue of `M`=$(primary.mass.lambda_min), M-orthogonality residual=$(primary.basis.M_orthogonality_error), diagonalization residual=$(primary.basis.L_diagonalization_error), zero modes=$(primary.basis.zero_modes), negative synchronizing-backbone modes=$(length(primary.basis.negative_modes)). These eigenvalues are not called graph frequencies when negative.\n")
        println(io,"## 8. Exact graph-modal transformation\n\nThe frequency-grid p95 of `PhiᴴTqPhi - (s²I+Lambda+Psihat)` is $graph_p95. `Psi_hat` is reported as `s D0_hat + DeltaL_hat + Pi_hat`; source and diagonal/off-diagonal norms are in TABLE_C06.\n")
        println(io,"## 9. Generalized graph dissipative operator\n\nThe reported operator is `D_G_eff=(1/omega) Im_H{Psi_hat(j omega)}`. The direct harmonic power identity maximum error across cases is $(maximum(values(powererrs))); the synthetic ExpB Sigma-special-case recovery error is $sigma_special_error. Negative eigenvalues indicate an indefinite operator under this convention, not instability.\n")
        println(io,"## 10. Real IEEE-39 modal mapping\n\nThe full finite poles are from the actual PowerDynamics reduced matrices. The analysis set includes the non-gauge spectral-abscissa pair, four low-damping pairs in 0.05–5 Hz, two additional high-retained-participation pairs, and a PLL-state-dominant pair when available. Results use M-weighted right-eigenvector graph projection, not classical participation factors. See TABLE_C07.\n")
        c33gamma=filter(r->r.case==primary.case.label,gammarows)
        println(io,"## 11. Exact intermodal Schur self-energy Gamma\n\n`Gamma_k=-Tkr*Trr\\Trk` is computed with factored solves. The C33 maximum scalar Schur residual is $(isempty(c33gamma) ? missing : maximum(getproperty.(c33gamma,:schur_residual))).\n")
        println(io,"## 12. Exact pathway decomposition\n\nAt each successful diagonal root and selected full-system pole offset, the exact pathway matrix `G_k[l,m]=-psi_kl*(Trr\\)[l,m]*psi_mk` is saved with its evaluation point and reference ID; its entries sum to Gamma. Top C33 diagonal pathways for spectral-abscissa pole $critpoleid (s=$(critpole===nothing ? missing : critpole.lambda), offset $((1e-4)*(1+im))): $(join(["G_$(r.mode_k)[$(r.mode_l),$(r.mode_m)]:$(r.contribution_abs)" for r in dompaths],", ")). Individual pathways depend on the chosen basis, especially inside repeated eigenspaces.\n")
        println(io,"## 13. Uncoupled graph-modal roots\n\nRoots of `s²+nu_k+psi_kk(s)=0` were searched with actual pole seeds, analytic derivatives, line search, and a controller-pole conditioning guard. Only converged roots in the configured 0.05–5 Hz band are retained; failed searches are not filled in. See TABLE_C09.\n")
        println(io,"## 14. Pole-shift predictor\n\nFor simple nondegenerate roots, `Delta_s=-Gamma_k(s0)/(2s0+psi'_kk(s0))` is compared to actual full PowerDynamics poles. It remains an asymptotic/empirical approximation.\n")
        println(io,"## 15. Predictor validity regime\n\nClassification: **$predictor_status**; matched mode count=$(length(pstatusvec)), median relative error=$predmed, max=$predmax. Eta-bin counts/errors are in RESULTS_EXP_C.json and TABLE_C10. The matching threshold and eta bins were fixed in the config before the analysis.\n")
        println(io,"## 16. Real intermodal susceptibility\n\nThe approximate ranking uses `|psi_kl psi_lk|/|t_l0|`, compared with the exact diagonal pathway magnitudes. C33 mean Spearman correlations: coupling=$(primarycorr["coupling_rank_correlation"]), susceptibility=$(primarycorr["susceptibility_rank_correlation"]).\n")
        println(io,"## 17. Complementary-mode resonance / detuning\n\nGate C8: **$detune**. No resonance is manufactured. TABLE_C12 reports coupling product, dynamic detuning, approximate susceptibility, and exact diagonal-path size separately.\n")
        println(io,"## 18. Baseline SG versus bus-33 GFL\n\nThe two cases use case-specific graph bases because the retained second-order metric changes when the PLL angle-rate coordinate replaces the SG speed coordinate. Modes are aligned by normalized M-reference overlap and physical q-shape overlap; mode indices are not assumed to match. See TABLE_C13.\n")
        println(io,"## 19. Cross-bus 30/33/35/37 analysis\n\nC30, C33, C35, and C37 reuse the frozen nominal SimpleGFLDC parameters. Baseline C0 is included in TABLE_C15.\n")
        println(io,"## 20. Backbone sensitivity\n\nStatus: **$robustness**. The primary and secondary choices, principal angles, pole-mode alignment, and Gamma changes are reported in TABLE_C14 and BACKBONE_AUDIT.md.\n")
        println(io,"## 21. Psi component ablation\n\nOperator-only diagnostics compare `sD0`, `sD0+DeltaL`, and full `Psi`; these are not claimed to be physically realizable plants. Results are in TABLE_C17.\n")
        println(io,"## 22. Gate summary\n\n| Gate | Status | Observed |\n|---|---|---|\n")
        for gRow in gates
            println(io,"| $(gRow.gate) | $(gRow.status) | $(gRow.observed) |")
        end
        println(io,"\n## 23. Exact results\n\nExact: the generalized operator identity, the congruence transform, scalar Schur self-energy, exact pathway reconstruction, and the harmonic power identity under the declared convention, subject to their reported numerical gates.\n")
        println(io,"## 24. Approximate results\n\nApproximate: uncoupled diagonal roots, pairwise susceptibility, and first-order pole displacement. The predictor status is `$predictor_status`.\n")
        println(io,"## 25. Unsupported / blocked claims\n\nNo optimal Kp/Ki; no optimal rho; no H4 claim; no global transient-stability theorem; no EMT claim; no field claim; no universal threshold on chi_comm; no proven chi_G stability certificate; no basis-invariant pathway claim.\n")
        println(io,"## 26. Scientific interpretation\n\nThe real-case result is that the exact Schur reduction can be expressed through a generalized dynamic self-energy without requiring a finite zero-frequency static correction. Whether this yields useful low-order physical explanation is case- and backbone-dependent; the report retains measured predictor and detuning outcomes rather than transferring synthetic ExpB conclusions.\n")
        println(io,"## 27. Decision for Experiment D\n\nReady for ExpD: **$(ready ? "YES" : "NO")**. Required exact-input, graph, Gamma, pathway, and cross-bus gates are summarized above. This decision does not authorize optimization in ExpC.\n")
    end
    open(joinpath(OUT,"BACKBONE_AUDIT.md"),"w") do io
        println(io,"# Experiment C backbone sensitivity audit\n\nConclusion: **$robustness**\n")
        println(io,"The primary backbone was preregistered as `Mhalf * Pg * Herm(Minvhalf*L0*Minvhalf) * Pg * Mhalf`. The secondary backbone is the Euclidean gauge-preserving Hermitian part of `L0`. No pole-fitting criterion was used.\n")
        println(io,"Maximum principal angle: $(isempty(allangles) ? missing : rad2deg(maximum(allangles))) degrees. Maximum observed Gamma magnitude ratio: $(isempty(gamma_ratios) ? missing : maximum(gamma_ratios)). Results are case-specific; inspect TABLE_C14 and TABLE_C04 for spectra and modal mappings.\n")
        println(io,"Classification criteria were fixed in expC.toml: ROBUST requires maximum principal angle <= $(CONFIG["backbone_principal_angle_max_deg"]) degrees and maximum Gamma ratio <= $(CONFIG["backbone_gamma_ratio_max"]); PARTIALLY_ROBUST requires maximum angle <= 45 degrees; otherwise BASIS_SENSITIVE.\n")
    end

    readytext="$(ready)"
    results=Dict{String,Any}(
        "experiment"=>"BND_EXP_C","status"=>expstatus,
        "input"=>Dict("expA_status"=>(expApass ? "PASS" : "FAIL"),"expB_status"=>(expBpass ? "PASS" : "FAIL"),
            "powerdynamics_version"=>versions.powerdynamics,"networkdynamics_version"=>versions.networkdynamics,
            "julia_version"=>string(VERSION),"git_head"=>GIT_HEAD,"blas"=>sprint(show,BLAS.get_config()),
            "bus33_Ared_reproduction_relative_error"=>a3err),
        "generalized_operator"=>Dict("psi_exact"=>c1pass,"median_reconstruction_error"=>psi_med,"p95_reconstruction_error"=>psi_p95),
        "pi_derivative"=>Dict("median_relative_error"=>dmedian,"p95_relative_error"=>dp95),
        "zero_frequency"=>Dict("status"=>zero.status,"numerical_nullity"=>zero.nullity,
            "dominant_states"=>zero.dominant_states,"observed_scaling_power"=>zero.observed_power,
            "semisimple_zero_supported"=>zero.semisimple,"semisimple_residual"=>zero.semisimple_residual,
            "residue_norm"=>zero.residue===nothing ? nothing : norm(zero.residue),
            "residue_checks"=>zero.residue_checks,"nullity_sensitivity"=>zero.nullities),
        "backbone"=>Dict("primary_definition"=>CONFIG["primary_backbone"],
            "M_min_eigenvalue"=>primary.mass.lambda_min,"gauge_residual"=>primary.primary.gauge_residual,
            "negative_backbone_modes"=>length(primary.basis.negative_modes),
            "robustness_status"=>robustness,"primary_spectrum"=>primary.basis.Lambda,
            "secondary_spectrum"=>primary.basis2.Lambda),
        "graph_basis"=>Dict("M_orthogonality_error"=>primary.basis.M_orthogonality_error,
            "L_diagonalization_error"=>primary.basis.L_diagonalization_error,
            "zero_modes"=>primary.basis.zero_modes,"negative_modes"=>primary.basis.negative_modes),
        "graph_transform"=>Dict("median_error"=>median(c05.reconstruction_error),"p95_error"=>graph_p95),
        "Gamma"=>Dict("median_schur_error"=>isempty(c08.schur_residual) ? nothing : median(c08.schur_residual),
            "p95_schur_error"=>gamma_p95),
        "pathways"=>Dict("max_reconstruction_error"=>pathway_max,"dominant_primary_pathways"=>dompaths),
        "predictor"=>Dict("status"=>predictor_status,"mode_count"=>length(pstatusvec),
            "median_relative_error"=>predmed,"eta_bins"=>predictor_bins(predrows)),
        "susceptibility"=>Dict("coupling_rank_correlation"=>primarycorr["coupling_rank_correlation"],
            "susceptibility_rank_correlation"=>primarycorr["susceptibility_rank_correlation"],"detuning_mechanism"=>detune,
            "by_case"=>correlations),
        "cross_bus"=>Dict(id=>Dict("graph_identity"=>r.graph_identity_pass,"Gamma_identity"=>r.Gamma_identity_pass,
            "predictor_status"=>r.predictor_status,"analyzed_modes"=>r.number_analyzed_modes)
            for (id,r) in zip(("C0","C30","C33","C35","C37"),cross)),
        "power_identity"=>Dict("by_case_max_relative_error"=>powererrs,"synthetic_expB_sigma_special_case_error"=>sigma_special_error),
        "gates"=>Dict(r.gate=>r.status for r in gates),"ready_for_expD"=>ready,
        "run_timestamp_utc"=>string(now(UTC)))
    write_json(joinpath(OUT,"RESULTS_EXP_C.json"),results)
    tablefiles=sort(filter(f->endswith(f,".csv"),readdir(TABLES)))
    markdown_tables(tablefiles)
    write_reproduction_readme()
    render_figures()

    println("EXP_C_STATUS:"); println(expstatus)
    println("EXP_A_INPUT:"); println(expApass ? "PASS" : "FAIL")
    println("EXP_B_INPUT:"); println(expBpass ? "PASS" : "FAIL")
    println("POWERDYNAMICS:"); println(versions.powerdynamics)
    println("GENERALIZED_PSI:"); println(c1pass ? "EXACT" : "FAIL")
    println("PSI_RECON_MEDIAN_ERROR:"); println(psi_med)
    println("PSI_RECON_P95_ERROR:"); println(psi_p95)
    println("ZERO_FREQUENCY_STATUS:"); println(zero.status)
    println("ZERO_FREQUENCY_DOMINANT_STATES:"); println(join(zero.dominant_states,"; "))
    println("M_SPD:"); println(primary.mass.spd ? "YES" : "NO")
    println("GRAPH_BACKBONE_CLASS:"); println(isempty(primary.basis.negative_modes) ? "LAPLACIAN_LIKE" : "INDEFINITE")
    println("GRAPH_M_ORTHOGONALITY_ERROR:"); println(primary.basis.M_orthogonality_error)
    println("GRAPH_L_DIAGONALIZATION_ERROR:"); println(primary.basis.L_diagonalization_error)
    println("GRAPH_TRANSFORM_P95_ERROR:"); println(graph_p95)
    println("GAMMA_SCHUR_P95_ERROR:"); println(gamma_p95)
    println("PATHWAY_RECON_MAX_ERROR:"); println(pathway_max)
    println("REAL_PREDICTOR_STATUS:"); println(predictor_status)
    println("REAL_PREDICTOR_MODE_COUNT:"); println(length(pstatusvec))
    println("REAL_PREDICTOR_MEDIAN_REL_ERROR:"); println(predmed)
    println("REAL_DETUNING_MECHANISM:"); println(detune)
    println("COUPLING_RANK_CORRELATION:"); println(primarycorr["coupling_rank_correlation"])
    println("SUSCEPTIBILITY_RANK_CORRELATION:"); println(primarycorr["susceptibility_rank_correlation"])
    println("BACKBONE_ROBUSTNESS:"); println(robustness)
    println("CROSS_BUS:")
    println("baseline=$(cross[1].graph_identity_pass ? "PASS" : "BLOCKED")")
    println("30=$(cross[2].graph_identity_pass ? "PASS" : "BLOCKED")")
    println("33=$(cross[3].graph_identity_pass ? "PASS" : "BLOCKED")")
    println("35=$(cross[4].graph_identity_pass ? "PASS" : "BLOCKED")")
    println("37=$(cross[5].graph_identity_pass ? "PASS" : "BLOCKED")")
    println("MAIN_EXACT_RESULT:"); println("The exact ExpA retained operator has an exact generalized graph-modal representation through Psi(s) for the analyzed cases.")
    println("MAIN_PHYSICAL_FINDING:"); println("C33 zero-frequency behavior is $(zero.status); dominant named states are $(join(zero.dominant_states,", ")).")
    println("MAIN_APPROXIMATE_FINDING:"); println("The real-case pole predictor is $(predictor_status) with $(length(pstatusvec)) matched roots and median relative error $(predmed).")
    println("MAIN_LIMITATION:"); println("The decomposition and pathways depend on the declared graph backbone; low-frequency global Sigma(s) remains unavailable.")
    println("READY_FOR_EXP_D_PLL_DESIGN:"); println(ready ? "YES" : "NO")
    println("FILES:"); println("reports/experiment_C/REPORT_EXP_C.md; reports/experiment_C/RESULTS_EXP_C.json; reports/experiment_C/CLAIM_LEDGER_EXP_C.md; reports/experiment_C/TABLES_EXP_C.md; reports/experiment_C/BACKBONE_AUDIT.md")
    println("PUSH:"); println("NO")
end

function write_reproduction_readme()
    open(joinpath(OUT,"README_REPRODUCE.md"),"w") do io
        println(io,"# Experiment C reproduction\n\nFrom the repository root, run:\n\n```powershell\njulia --project=. experiments/bnd_expC/run_experiment_C.jl\n```\n")
        println(io,"The run reads the Experiment-A bus-33 `A_reduced` matrix/state partition and Experiment-A/B machine-readable status artifacts without modifying them. It rebuilds the all-SG baseline and bus 30/33/35/37 cases with the pinned PowerDynamics 5.0.0 model and frozen nominal SimpleGFLDC template on the first run. It retains those exact reduced matrices/state maps in `reports/experiment_C/matrices/`, keyed by hashes of the project, manifest, model adapter, and ExpA inputs. Later runs reuse only a matching cache; pass `--rebuild-models` to force a fresh PowerDynamics extraction. No tuning, package update, or network push occurs.\n")
        println(io,"The runner regenerates CSV tables, JSON results, Markdown report and tables, matrices, and figures under `reports/experiment_C/`. Julia and package versions, BLAS configuration, Git HEAD, source artifact hashes, timestamp, and the bus-33 reconstruction residual are recorded.\n")
        println(io,"Focused unit checks are in `test/bnd_expC/runtests.jl`; ExpB checks remain in `test/bnd_expB/runtests.jl`.\n")
    end
end

function render_figures()
    py=joinpath(@__DIR__,"render_figures.py")
    if Sys.which("python")===nothing
        @warn "Python not found; CSV/JSON outputs are complete but PNG figures were not rendered"
        return
    end
    Base.run(`python $py $TABLES $FIGURES $MATRICES`)
end

run()
