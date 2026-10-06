using CSV, DataFrames, DelimitedFiles, ForwardDiff, LinearAlgebra, NetworkDynamics, PowerDynamics
include(joinpath(@__DIR__,"..","..","src","bnd_model_audit_m","PDCasesM.jl"))
using .PDCasesM

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const OUT=joinpath(ROOT,"reports","experiment_M","matrices","ExpK_nominal")

function main()
    println("LOAD_INPUT_BUILD");flush(stdout)
    nw=PDCasesM.build_case("ExpK_nominal")
    eq=PDCasesM.PD39.initialize_equilibrium(nw;sparse=false,check=:error)
    input=[VIndex(16,:ZIPLoad₊Pset),VIndex(16,:ZIPLoad₊Qset)]
    output=[VIndex(38,:ctrld_gen₊machine₊ω),VIndex(39,:machine₊ω)]
    println("LOAD_INPUT_LINEARIZE");flush(stdout)
    sys=linearize_network(eq.state;in=input,out=output)
    for (name,m) in (("M",sys.M),("A",sys.A),("B",sys.B),("C",sys.C),("D",sys.D))
        mat=m isa UniformScaling ? Matrix{Float64}(I,size(sys.A,1),size(sys.A,1)) : Matrix(m)
        writedlm(joinpath(OUT,"PD_load_input_"*name*".csv"),mat,',')
    end
    pset=Float64(eq.state[input[1]]);qset=Float64(eq.state[input[2]])
    h38=Float64(eq.state[VIndex(38,:ctrld_gen₊machine₊H)])
    sn38=Float64(eq.state[VIndex(38,:ctrld_gen₊machine₊Sn)])
    h39=Float64(eq.state[VIndex(39,:machine₊H)])
    sn39=Float64(eq.state[VIndex(39,:machine₊Sn)])
    CSV.write(joinpath(OUT,"PD_load_input_metadata.csv"),DataFrame(
        Pset_pu=[pset],Qset_pu=[qset],H38_Sn38=[h38*sn38],H39_Sn39=[h39*sn39],
        equilibrium_residual=[PDCasesM.PD39._eq_residual(eq.state)]))
    x0=uflat(eq.state);p0=pflat(eq.state);t0=eq.state.t
    f=x->begin
        dx=similar(x)
        nw(dx,x,p0,t0)
        dx
    end
    ad=ForwardDiff.jacobian(f,x0)
    delta=ad-sys.A
    CSV.write(joinpath(OUT,"PD_direct_AD_check.csv"),DataFrame(
        max_abs_A_error=[maximum(abs.(delta))],relative_A_error=[norm(delta)/norm(sys.A)]))
    println("LOAD_INPUT_DONE B=",size(sys.B)," AD_error=",maximum(abs.(delta)));flush(stdout)
end

main()
