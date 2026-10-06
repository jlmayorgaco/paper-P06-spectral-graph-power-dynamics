include("LocalOracle.jl")
using .LocalOracle, LinearAlgebra, Random, CSV, DataFrames, TOML
const L=LocalOracle
include("StaticFrequency.jl")
include("StaticBound.jl")
sm=StaticFrequency.build();rng=MersenneTwister(20261001);rows=NamedTuple[]
for j in 1:100
    e=(.03 .+.9rand(rng,10)).*Float64.(rand(rng,10).>.3)
    sum(e)>0 || (e[1]=.5)
    kp=L.CTX.kpmin+rand(rng,10).*(L.CTX.kpmax-L.CTX.kpmin)
    ki=L.CTX.kimin+rand(rng,10).*(L.CTX.kimax-L.CTX.kimin)
    d=Dict("rho"=>1 .-e,"Kp"=>kp,"Ki"=>ki);a=L.architecture(findall(e.>0).+29,d)
    b=L.matrices(a,L.encode(a.support,d));full=100real.(-b.C*(b.A\b.B)+b.D)
    static=StaticFrequency.frequency(sm,e)
    push!(rows,(;id=j,static_Hz=static,full_Hz=full[1],max_error=maximum(abs.(full.-static)),condition=cond(StaticFrequency.matrix(sm,e))))
end
CSV.write(joinpath(L.OUT,"STATIC_MAP_WITHHELD_IDENTITY.csv"),DataFrame(rows))
println("STATIC_MAP_ERROR ",maximum(r.max_error for r in rows));flush(stdout)
t=@elapsed coeff=StaticBound.coefficients(sm,L.CTX.power)
bounds=StaticBound.root_bounds(coeff)
CSV.write(joinpath(L.OUT,"BERNSTEIN_COEFFICIENTS_FLOAT64.csv"),DataFrame(D=coeff.D,N=coeff.N,cost=coeff.cost))
result=Dict("status"=>"PROVISIONAL_FLOAT64_BOUND_REQUIRES_OUTWARD_ENCLOSURES",
    "lower_MW"=>bounds.lower,"positive_D_lower"=>bounds.positive.bound,"negative_D_lower"=>bounds.negative.bound,
    "positive_multipliers"=>bounds.positive.multipliers,"negative_multipliers"=>bounds.negative.multipliers,
    "coefficient_count"=>length(coeff.D),"coefficient_scale"=>coeff.scale,"elapsed_s"=>t,
    "max_withheld_identity_error"=>maximum(r.max_error for r in rows),"global_certificate"=>false)
open(joinpath(L.OUT,"STATIC_BOUND_FLOAT64.toml"),"w") do io;TOML.print(io,result);end
println(result)
