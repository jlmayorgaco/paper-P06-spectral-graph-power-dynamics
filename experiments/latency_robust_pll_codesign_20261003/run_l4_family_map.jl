using CSV, DataFrames, LinearAlgebra, TOML

module DC
include(joinpath(@__DIR__,"..","delay_dressed_replacement_frontier_20261002","DelayCharacteristic.jl"))
end

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
BLAS.set_num_threads(1)

function vectors(L,s,tau_ms)
    D=DC.DelayCharacteristic.delta_matrix(L,s,fill(tau_ms/1000,10))
    S=svd(D)
    u=ComplexF64.(L.Q*S.U[:,end]);v=ComplexF64.(L.Q*S.V[:,end])
    u/=norm(u);v/=norm(v)
    p=abs.(conj.(u).*v);p/=sum(p)
    (;u,v,p,residual=S.S[end]/max(1,opnorm(D)))
end

function modes_for(id)
    if id in ("Z","N")
        t=CSV.read(joinpath(@__DIR__,"L0_REPRODUCTION.csv"),DataFrame)
        r=only(eachrow(t[t.design_id.==id,:]))
        return [(;rank=1,tau_ms=r.tau_crit_local_ms,
            s=r.critical_root_real+im*r.critical_root_imag)]
    end
    if id=="mode_switch_eta"
        t=CSV.read(joinpath(@__DIR__,"L4_MODE_SWITCH_POINT.csv"),DataFrame)
        r=only(eachrow(t))
        return [(;rank=1,tau_ms=r.Z_family_crossing_ms,s=-0.05+im*2pi*r.Z_frequency_hz),
            (;rank=2,tau_ms=r.N_family_crossing_ms,s=-0.05+im*2pi*r.N_frequency_hz)]
    end
    path=joinpath(@__DIR__,"evaluations",id,"ROOTS.csv")
    isfile(path) || return NamedTuple[]
    t=CSV.read(path,DataFrame)
    [(;rank=r.rank,tau_ms=r.local_crossing_ms,
      s=r.critical_real+im*r.critical_imag) for r in eachrow(t[1:min(3,nrow(t)),:])]
end

function main()
    ctx=DC.DelayCharacteristic.R.N.design_context(ROOT)
    templates=Dict{String,Any}()
    for id in ("Z","N")
        d=TOML.parsefile(joinpath(@__DIR__,"designs",id=="Z" ?
            "Z_zero_delay_tuned.toml" : "N_nominal.toml"))
        L=DC.DelayCharacteristic.linearization(ctx,Float64.(d["rho"]),Float64.(d["Kp"]),Float64.(d["Ki"]))
        r=only(modes_for(id))
        templates[id]=vectors(L,r.s,r.tau_ms)
    end
    ids=("Z","N","mode_switch_eta","Z_step1_ball","Z_step1_multimode_lp","N_step1_ball",
        "N_step1_multimode_lp","N_step1_multimode_lp_next_full",
        "N_step1_multimode_lp_next_half","N_step1_multimode_lp_next_0p625",
        "N_step1_multimode_lp_event_active_lp",
        "N_step1_multimode_lp_event_active_lp_followup_lp",
        "N_step1_multimode_lp_event_active_lp_followup_lp_followup_lp",
        "N_step1_multimode_lp_event_active_lp_followup_lp_followup_lp_followup_lp",
        "N_step1_multimode_lp_event_active_lp_followup_lp_followup_lp_followup_lp_followup_lp")
    rows=NamedTuple[]
    for id in ids
        path=joinpath(@__DIR__,"designs",id=="Z" ? "Z_zero_delay_tuned.toml" :
            id=="N" ? "N_nominal.toml" : id*".toml")
        isfile(path) || continue
        d=TOML.parsefile(path)
        L=DC.DelayCharacteristic.linearization(ctx,Float64.(d["rho"]),Float64.(d["Kp"]),Float64.(d["Ki"]))
        for r in modes_for(id)
            vv=vectors(L,r.s,r.tau_ms)
            zr=abs(dot(templates["Z"].v,vv.v))^2
            nr=abs(dot(templates["N"].v,vv.v))^2
            zl=abs(dot(templates["Z"].u,vv.u))^2
            nl=abs(dot(templates["N"].u,vv.u))^2
            zp=sum(sqrt.(templates["Z"].p.*vv.p))^2
            np=sum(sqrt.(templates["N"].p.*vv.p))^2
            label=zr>nr ? "Z_like" : "N_like"
            push!(rows,(;candidate_id=id,rank=r.rank,local_crossing_ms=r.tau_ms,
                frequency_hz=imag(r.s)/(2pi),family_label=label,
                right_MAC_Z=zr,right_MAC_N=nr,left_MAC_Z=zl,left_MAC_N=nl,
                physical_component_participation_Bhattacharyya_Z=zp,
                physical_component_participation_Bhattacharyya_N=np,
                full_characteristic_normalized_sigma_min=vv.residual,
                status="NUMERICAL_LEFT_RIGHT_MODAL_OVERLAP_AND_PARTICIPATION_SIMILARITY"))
        end
        CSV.write(joinpath(@__DIR__,"L4_MODAL_FAMILY_MAP.csv"),DataFrame(rows))
        println("L4_MAP_DONE ",id);flush(stdout)
    end
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
