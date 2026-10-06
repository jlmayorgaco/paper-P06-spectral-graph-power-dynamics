using CSV, DataFrames, LinearAlgebra, TOML
include(joinpath(@__DIR__,"..","delay_dressed_replacement_frontier_20261002","DelayCharacteristic.jl"))
const DC=DelayCharacteristic
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const OLD=joinpath(ROOT,"experiments","latency_robust_pll_codesign_20261003")
BLAS.set_num_threads(1)

function vector(L,s,tau_ms)
    D=DC.delta_matrix(L,s,fill(tau_ms/1000,10)); S=svd(D)
    u=ComplexF64.(L.Q*S.U[:,end]);v=ComplexF64.(L.Q*S.V[:,end])
    u/=norm(u);v/=norm(v)
    p=abs.(conj.(u).*v);p/=sum(p)
    (;u,v,p,residual=S.S[end]/max(1,opnorm(D)))
end

function modes(id)
    if id in ("Z_zero_delay_tuned","N_nominal")
        label=id=="Z_zero_delay_tuned" ? "Z" : "N"
        table=CSV.read(joinpath(OLD,"L0_REPRODUCTION.csv"),DataFrame)
        r=only(eachrow(table[table.design_id.==label,:]))
        return [(;rank=1,tau_ms=Float64(r.tau_crit_local_ms),
            s=ComplexF64(r.critical_root_real+im*r.critical_root_imag))]
    end
    file=id=="previous_best" ? joinpath(OLD,"evaluations",
        "N_step1_multimode_lp_event_active_lp_followup_lp_followup_lp_followup_lp_followup_lp","ROOTS.csv") :
        joinpath(@__DIR__,"evaluations",id,"ROOTS.csv")
    if id!="previous_best"
        expanded=joinpath(@__DIR__,"evaluations",id,"ROOTS_EXPANDED.csv")
        isfile(expanded) && (file=expanded)
    end
    t=CSV.read(file,DataFrame)
    [(;rank=Int(r.rank),tau_ms=Float64(r.local_crossing_ms),
        s=ComplexF64(r.critical_real+im*r.critical_imag)) for r in eachrow(t)]
end

function main()
    ctx=DC.R.N.design_context(ROOT)
    refs=Dict{String,Any}()
    for id in ("Z_zero_delay_tuned","N_nominal","previous_best")
        d=TOML.parsefile(joinpath(@__DIR__,"designs",id*".toml"))
        L=DC.linearization(ctx,Float64.(d["rho"]),Float64.(d["Kp"]),Float64.(d["Ki"]))
        refs[id]=[(;rank=r.rank,vec=vector(L,r.s,r.tau_ms)) for r in modes(id)]
    end
    rows=NamedTuple[]
    ids=["previous_best","pending_41ms","pending_next_lp","full20_next_lp"]
    append!(ids,sort([splitext(basename(p))[1] for p in readdir(joinpath(@__DIR__,"designs");join=true)
                      if startswith(basename(p),"full20_step") && endswith(p,".toml")]))
    append!(ids,sort([splitext(basename(p))[1] for p in readdir(joinpath(@__DIR__,"designs");join=true)
                      if startswith(basename(p),"reserve") && endswith(p,".toml")]))
    for id in ids
        d=TOML.parsefile(joinpath(@__DIR__,"designs",id*".toml"))
        L=DC.linearization(ctx,Float64.(d["rho"]),Float64.(d["Kp"]),Float64.(d["Ki"]))
        for r in modes(id)
            vv=vector(L,r.s,r.tau_ms)
            z=refs["Z_zero_delay_tuned"][1].vec;n=refs["N_nominal"][1].vec
            old=refs["previous_best"]
            matches=[abs(dot(x.vec.v,vv.v))^2 for x in old]
            best=argmax(matches)
            family=matches[best]<0.1 ? "F4_new_orthogonal" : "F"*string(best)
            push!(rows,(;design_id=id,rank=r.rank,local_crossing_ms=r.tau_ms,
                frequency_hz=imag(r.s)/(2pi),closest_previous_best_family=family,
                closest_previous_right_MAC=matches[best],
                closest_previous_left_MAC=abs(dot(old[best].vec.u,vv.u))^2,
                closest_previous_participation_overlap=sum(sqrt.(old[best].vec.p.*vv.p))^2,
                right_MAC_Z=abs(dot(z.v,vv.v))^2,right_MAC_N=abs(dot(n.v,vv.v))^2,
                full_characteristic_residual=vv.residual,
                status="NUMERICAL_OVERLAP_NEAREST_REFERENCE_NOT_CERTIFIED_FAMILY_ID"))
        end
    end
    CSV.write(joinpath(@__DIR__,"TABLE_F3_ACTIVE_MODAL_FAMILIES.csv"),DataFrame(rows))
    foreach(println,rows)
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
