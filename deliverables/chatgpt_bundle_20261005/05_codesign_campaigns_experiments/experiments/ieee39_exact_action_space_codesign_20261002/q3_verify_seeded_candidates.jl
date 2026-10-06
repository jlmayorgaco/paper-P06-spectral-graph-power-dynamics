using CSV,DataFrames,TOML,LinearAlgebra
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const DC=include(joinpath(ROOT,"experiments","delay_dressed_replacement_frontier_20261002","DelayCharacteristic.jl"))
const R=DC.R
function main()
    d=TOML.parsefile(joinpath(@__DIR__,"seed_uniform_875.toml"));ctx=R.N.design_context(ROOT)
    L=DC.linearization(ctx,Float64.(d["rho"]),Float64.(d["Kp"]),Float64.(d["Ki"]))
    rows=NamedTuple[]
    for tau_ms in (20.0,40.0), sr in (-0.02,0.05)
        tau=fill(tau_ms/1000,length(L.Ai));s=complex(sr,0.0);M=DC.delta_matrix(L,s,tau)
        sv=svdvals(M);push!(rows,(;tau_ms,s_real=sr,s_imag=0.0,sigma_min=minimum(sv),
            sigma_max=maximum(sv),relative_sigma_min=minimum(sv)/maximum(sv),
            normalized_vector_residual=minimum(sv)/max(1,norm(M))))
    end
    CSV.write(joinpath(@__DIR__,"TABLE_Q03_SEEDED_CANDIDATE_SINGULAR_VALUES.csv"),DataFrame(rows));println(DataFrame(rows))
end
main()
