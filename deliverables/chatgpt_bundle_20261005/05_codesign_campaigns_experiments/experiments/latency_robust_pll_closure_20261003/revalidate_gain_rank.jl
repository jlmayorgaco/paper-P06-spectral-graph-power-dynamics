using CSV, DataFrames, TOML, LinearAlgebra
include(joinpath(@__DIR__,"..","delay_dressed_replacement_frontier_20261002","DelayCharacteristic.jl"))
const DC=DelayCharacteristic
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
BLAS.set_num_threads(1)

function main()
    ctx=DC.R.N.design_context(ROOT)
    d0=TOML.parsefile(joinpath(@__DIR__,"designs","N_nominal.toml"))
    rho=Float64.(d0["rho"])
    L0=DC.linearization(ctx,rho,Float64.(d0["Kp"]),Float64.(d0["Ki"]))
    rows=NamedTuple[]
    ids=["previous_best","pending_41ms","pending_next_lp","full20_next_lp"]
    append!(ids,sort([splitext(basename(p))[1] for p in readdir(joinpath(@__DIR__,"designs");join=true)
                      if startswith(basename(p),"full20_step") && endswith(p,".toml")]))
    append!(ids,sort([splitext(basename(p))[1] for p in readdir(joinpath(@__DIR__,"designs");join=true)
                      if startswith(basename(p),"reserve") && endswith(p,".toml")]))
    for id in ids
        d=TOML.parsefile(joinpath(@__DIR__,"designs",id*".toml"))
        Float64.(d["rho"])==rho || error("replacement changed")
        L=DC.linearization(ctx,rho,Float64.(d["Kp"]),Float64.(d["Ki"]))
        change=L.A-L0.A
        factor=(L.B-L0.B)*transpose(L0.C)
        sigma=svdvals(change)
        rank=sum(sigma .> 1e-8*maximum(sigma))
        e=norm(change-factor)/max(norm(change),eps())
        a0=norm(L.A0-L0.A0)/max(norm(L0.A0),eps())
        push!(rows,(;design_id=id,dimension=size(change,1),gain_variable_count=20,
            action_rank=rank,rank_upper_bound=10,relative_action_reconstruction_error=e,
            relative_A0_difference=a0,status=(rank<=10 && e<1e-9 && a0<1e-9 ?
            "EXACT_FACTORIZATION_NUMERICALLY_VALIDATED" : "FACTORIZATION_FAILURE")))
    end
    CSV.write(joinpath(@__DIR__,"TABLE_F11_GAIN_ACTION_RANK.csv"),DataFrame(rows))
    foreach(println,rows)
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
