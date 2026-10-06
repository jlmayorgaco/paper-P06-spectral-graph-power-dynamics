include(joinpath(@__DIR__,"DelayedEvents.jl"))
using .DelayedEvents, LinearAlgebra, CSV, DataFrames, TOML, Random
const R=DelayedEvents.R
ctx=R.N.design_context(R.ROOT)
d=TOML.parsefile(joinpath(@__DIR__,"baseline.toml"))
p=vcat(d["rho"],log.(d["Kp"]),log.(d["Ki"]))
rng=MersenneTwister(20261003)
direction=randn(rng,30);direction./=maximum(abs.(direction))
direction.*=vcat(fill(.01,10),fill(.08,20))
CSV.write(joinpath(@__DIR__,"FD_DIRECTION.csv"),DataFrame(parameter=1:30,direction=direction))
h=.01;fdstep=1e-3;values=[]
for sign in (-1,1)
    pp=p+sign*fdstep*direction
    m=R.model(ctx,pp[1:10],exp.(pp[11:20]),exp.(pp[21:30]);bus=16,delta=100.,dc_convention=:physical_supply)
    r=DelayedEvents.tangent(m;h,sensitivities=false)
    push!(values,r.vals)
    println("FD_DONE sign=",sign," values=",r.vals);flush(stdout)
end
gradid=isempty(ARGS) ? "baseline" : ARGS[1]
G=Matrix(CSV.read(joinpath(@__DIR__,"grad",gradid,"gradients.csv"),DataFrame))[6:10,:]
analytic=G*direction;fd=(values[2]-values[1])/(2fdstep)
absolute=abs.(analytic-fd);relative=absolute./max.(abs.(analytic),1e-10)
rows=DataFrame(metric=["F","R","Vmin","Vmax","slack"],analytic=analytic,finite_difference=fd,absolute_error=absolute,relative_error=relative)
CSV.write(joinpath(@__DIR__,"TABLE_06B_EVENT_DERIVATIVE_VALIDATION.csv"),rows)
println(rows)
all((relative .<= .01) .| (absolute .<= 1e-7)) || error("DELAYED_TANGENT_FD_FAILED")
