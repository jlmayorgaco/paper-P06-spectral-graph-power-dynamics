using CSV, DataFrames, LinearAlgebra, NetworkDynamics, PowerDynamics,
    OrdinaryDiffEqRosenbrock, SciMLBase
include(joinpath(@__DIR__,"..","..","src","bnd_model_audit_m","PDCasesM.jl"))
using .PDCasesM

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const OUT=joinpath(ROOT,"reports","experiment_M","tables")
const MODEL=PDCasesM.PD39.PD39Model

function pulse_network(nw,amplitude)
    vertices,edges=MODEL.copy_network_components(nw)
    d=get_defaults_dict(vertices[16])
    ps=only([s for s in keys(d) if occursin("Pset",string(s))])
    qs=only([s for s in keys(d) if occursin("Qset",string(s))])
    p0=d[ps];q0=d[qs]
    affect=(u,p,ctx)->begin
        factor=ctx.t<1.1 ? 1+amplitude : 1.0
        p[ps]=p0*factor;p[qs]=q0*factor
    end
    set_callback!(vertices[16],PresetTimeComponentCallback([1.0,1.1],
        ComponentAffect(affect,(),(ps,qs))))
    out=Network(vertices,edges)
    set_jac_prototype!(out)
    return out
end

function main()
    println("NONLINEAR_SCALING_BUILD");flush(stdout)
    nw=PDCasesM.build_case("ExpK_nominal")
    eq=PDCasesM.PD39.initialize_equilibrium(nw;sparse=false,check=:error)
    h38=Float64(eq.state[VIndex(38,:ctrld_gen₊machine₊H)])
    sn38=Float64(eq.state[VIndex(38,:ctrld_gen₊machine₊Sn)])
    h39=Float64(eq.state[VIndex(39,:machine₊H)])
    sn39=Float64(eq.state[VIndex(39,:machine₊Sn)])
    weights=[h38*sn38,h39*sn39]
    refs=[Float64(eq.state[VIndex(38,:ctrld_gen₊machine₊ω)]),
          Float64(eq.state[VIndex(39,:machine₊ω)])]
    rows=NamedTuple[]
    for amplitude in (0.001,0.0001,0.00001)
        pn=pulse_network(nw,amplitude)
        prob=SciMLBase.ODEProblem(pn,eq.state,(0.0,12.0))
        sol=SciMLBase.solve(prob,OrdinaryDiffEqRosenbrock.Rodas5P();
            callback=get_callbacks(pn),initializealg=SciMLBase.NoInit(),
            saveat=0.02,abstol=1e-12,reltol=1e-12)
        SciMLBase.successful_retcode(sol.retcode) || error("TDS failed at $amplitude")
        for t in Float64.(sol.t)
            state=NetworkDynamics.NWState(sol,t)
            omega=[Float64(state[VIndex(38,:ctrld_gen₊machine₊ω)]),
                   Float64(state[VIndex(39,:machine₊ω)])]
            f=60dot(weights,omega.-refs)/sum(weights)
            push!(rows,(case="ExpK_nominal",event_bus=16,pulse_fraction=amplitude,
                time_s=t,COI_frequency_Hz=f))
        end
        CSV.write(joinpath(OUT,"TABLE_M20_nonlinear_traces.csv"),DataFrame(rows))
        println("NONLINEAR_SCALING_DONE amplitude=",amplitude," points=",length(sol.t));flush(stdout)
    end
end

main()
