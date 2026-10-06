include("LocalOracle.jl")
using .LocalOracle, LinearAlgebra, CSV, DataFrames, TOML
const L=LocalOracle
d=TOML.parsefile(joinpath(L.OUT,"FEASIBLE_NUMERICAL_REPRESENTATIVE.toml"));rows=NamedTuple[]
for bus in setdiff(30:39,Int.(d["support"]))
    for ep in (1e-5,1e-7,1e-9)
        rho=copy(d["rho"]);rho[bus-29]=1-ep
        sp=L.N.spectrum(L.CTX,rho,d["Kp"],d["Ki"])
        Aq=sp.quotient'*sp.model.Ared*sp.quotient
        beta0=1/opnorm(inv(Matrix(-Aq-.05I)))
        push!(rows,(;bus,epsilon=ep,alpha=sp.alpha,beta_zero_witness_upper=beta0,
            SG_open_loop_alpha=maximum(real.(eigvals(L.CTX.net.sg[bus].J.A))),
            modal_infeasible=sp.alpha>-.05,robust_infeasible_witness=beta0<L.BETA))
    end
end
CSV.write(joinpath(L.OUT,"ADJACENT_SUPPORT_LIMITS.csv"),DataFrame(rows))
println(DataFrame(rows))
