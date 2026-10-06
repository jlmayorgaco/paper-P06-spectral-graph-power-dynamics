using CSV, DataFrames, LinearAlgebra

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(ROOT, "src", "bnd_design", "AnalyticSG.jl"))
include(joinpath(ROOT, "src", "bnd_design", "AnalyticGFLPLL.jl"))

function matrix_csv(name)
    df = CSV.read(joinpath(ROOT, "reports", "experiment_A", "matrices", "bus33_" * name * ".csv"), DataFrame)
    return Matrix{Float64}(df[:, 2:end])
end

states = CSV.read(joinpath(ROOT,"reports","experiment_A","tables","TABLE_A03_state_partition.csv"),DataFrame)
diffrows = filter(r -> r.differential_or_algebraic == "differential", states)
algrows = filter(r -> r.differential_or_algebraic == "algebraic", states)
Fx,Fy,Gx,Gy = (matrix_csv(n) for n in ("F_x","F_y","G_x","G_y"))
println("dims ", size(Fx), " ", size(Fy), " ", size(Gx), " ", size(Gy))

for bus in (30,31,39)
    op = AnalyticSG.frozen_bus_operating_point(ROOT,bus)
    J = AnalyticSG.jacobians(op.x,op.u,op.parameters)
    dr = filter(r -> Int(r.bus) == bus, diffrows)
    ar = filter(r -> Int(r.bus) == bus, algrows)
    xi = Int.(dr.reduced_state_index)
    yi = [findfirst(==(Int(r.state_index)),Int.(algrows.state_index)) for r in eachrow(ar)]
    p = collect(1:length(xi)); length(p) == 12 && (p[1:2] = [2,1])
    println("bus ",bus," x=",xi," y=",yi," Aerr=",norm(Fx[xi,xi]-J.A[p,p],Inf),
        " Berr=",norm(Fy[xi,yi]-J.B[p,:],Inf))
    println("  C signed errors ",norm(Gx[yi,xi]-J.C[:,p],Inf)," ",norm(Gx[yi,xi]+J.C[:,p],Inf))
    dA=Fx[xi,xi]-J.A[p,p]; dB=Fy[xi,yi]-J.B[p,:]
    ia=argmax(abs.(dA)); ib=argmax(abs.(dB))
    println("  Amax ",Tuple(ia)," frozen=",Fx[xi,xi][ia]," local=",J.A[p,p][ia])
    println("  Bmax ",Tuple(ib)," frozen=",Fy[xi,yi][ib]," local=",J.B[p,:][ib])
    println("  C frozen=",Gx[yi,xi]," local=",J.C[:,p])
    println("  rhs residual ",norm(length(op.x)==12 ? AnalyticSG.rhs(op.x,op.u,op.parameters) : AnalyticSG.rhs_uncontrolled(op.x,op.u,op.parameters),Inf))
end

op33 = AnalyticGFLPLL.frozen_bus33_operating_point(ROOT)
J33 = AnalyticGFLPLL.jacobians(op33.x,op33.u,op33.parameters)
dr = filter(r -> Int(r.bus) == 33, diffrows)
ar = filter(r -> Int(r.bus) == 33, algrows)
xi = Int.(dr.reduced_state_index)
yi = [findfirst(==(Int(r.state_index)),Int.(algrows.state_index)) for r in eachrow(ar)]
println("bus33 x=",xi," y=",yi," Aerr=",norm(Fx[xi,xi]-J33.A,Inf),
    " Berr=",norm(Fy[xi,yi]-J33.B,Inf))
println("  C signed errors ",norm(Gx[yi,xi]-J33.C,Inf)," ",norm(Gx[yi,xi]+J33.C,Inf))
