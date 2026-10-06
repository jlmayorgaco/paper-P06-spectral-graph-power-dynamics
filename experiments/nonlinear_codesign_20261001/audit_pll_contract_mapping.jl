include("ReducedDAE.jl")
using .ReducedDAE, Random, SHA, TOML, LinearAlgebra
const R=ReducedDAE
const OUT=joinpath(R.OUT,"pll_sector_contract")
mkpath(OUT)
ctx=R.N.design_context(R.ROOT)
rng=MersenneTwister(20261001)
maxerror=0.0
cases=0
voltages=Float64[]
for bus in 30:39
    op=R.N.trim_gfl(ctx,bus)
    push!(voltages,norm(op.u))
    for sample in 1:100
        V0=0.95+0.15rand(rng)
        reference=8randn(rng)
        delta=2rand(rng)-1
        nu=0.2randn(rng)
        du=0.1randn(rng,2)
        theta=reference+delta
        c,s=cos(reference),sin(reference)
        u=[c -s;s c]*([V0,0.0]+du)
        x=copy(op.x)
        x[3]=theta;x[4]=randn(rng);x[5]=randn(rng)
        x[6:8].+=0.1randn(rng,3)
        x[9]=2.5+0.1randn(rng)
        kp=10+80rand(rng);ki=60+500rand(rng)
        full=R.gfl_rhs(x,u,op.pars,kp,ki,:physical_supply)
        eq=-V0*sin(delta)-sin(delta)*du[1]+cos(delta)*du[2]
        predicted=[x[4]-nu,(x[5]+kp*eq-x[4])/op.pars.pll_tau,ki*eq]
        observed=[full[3]-nu,full[4],full[5]]
        global maxerror=max(maxerror,norm(predicted-observed,Inf))
        global cases+=1
    end
end
@assert maxerror<1e-8
@assert all(0.95 .<=voltages.<=1.10)
source=joinpath(@__DIR__,"ReducedDAE.jl")
out=Dict("status"=>"INSTALLED_PLL_EQUATIONS_MATCH_EXACT_MOVING_FRAME_CONTRACT",
         "cases"=>cases,"max_absolute_rhs_error"=>maxerror,
         "nominal_voltage_min"=>minimum(voltages),"nominal_voltage_max"=>maximum(voltages),
         "pll_tau"=>R.N.trim_gfl(ctx,30).pars.pll_tau,
         "ReducedDAE_sha256"=>bytes2hex(sha256(read(source))),
         "script_sha256"=>bytes2hex(sha256(read(@__FILE__))),
         "scope"=>"PLL equations only; current and DC states were varied but are not certified")
open(joinpath(OUT,"julia_mapping_audit.toml"),"w") do io
    TOML.print(io,out)
end
TOML.print(stdout,out)
