using CSV, DataFrames, LinearAlgebra, NetworkDynamics, PowerDynamics
include(joinpath(@__DIR__,"..","..","src","pd39","PD39.jl"))
using .PD39

const MODEL=PD39.PD39Model
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const OUT=joinpath(ROOT,"reports","experiment_M","tables")
mkpath(OUT)

function observed(state,bus,name)
    Float64(state[VIndex(bus,Symbol(name))])
end

function local_pole_groups(state,bus)
    sys=linearize_component(state,VIndex(bus))
    M=sys.M isa UniformScaling ? Matrix{Float64}(I,size(sys.A,1),size(sys.A,1)) : Matrix(sys.M)
    F=eigen(Matrix(sys.A),M)
    sg=ComplexF64[];gfl=ComplexF64[]
    syms=string.(sys.sym)
    for k in eachindex(F.values)
        isfinite(F.values[k]) || continue
        e=abs2.(F.vectors[:,k])
        esg=sum((e[j] for j in eachindex(e) if occursin("machine",syms[j]) ||
            occursin("avr",syms[j]) || occursin("gov",syms[j]));init=0.0)
        egfl=sum((e[j] for j in eachindex(e) if occursin("gfl",syms[j]));init=0.0)
        if esg>=egfl;push!(sg,F.values[k]) else push!(gfl,F.values[k]) end
    end
    fmt(v)=join(["$(real(z))$(imag(z)<0 ? "" : "+")$(imag(z))im" for z in sort(v;by=x->(real(x),imag(x)))],";")
    return fmt(sg),fmt(gfl)
end

function main()
    rows=NamedTuple[]
    data=MODEL.ieee39_data()
    for bus in (36,38)
        baseline=PD39.baseline_network()
        middle=MODEL.weighted_replacement_network(baseline,bus,0.5,1.0)
        mrow=data.machine[findfirst(==(bus),data.machine.bus),:]
        for rho in (0.0,0.01,0.5,0.99,1.0)
            nw=rho==0 ? baseline : rho==1 ? PD39.replace_bus(baseline,bus;
                template=PD39.simple_gfldc_template()) : rho==0.5 ? middle :
                MODEL.set_weighted_replacement_parameters(middle,bus,rho,1.0)
            try
                eq=PD39.initialize_equilibrium(nw;sparse=false,check=:error)
                prefix="ctrld_gen₊machine₊"
                sg_present=rho<1
                sn=sg_present ? observed(eq.state,bus,prefix*"Sn") : 0.0
                h=sg_present ? observed(eq.state,bus,prefix*"H") : 0.0
                torque=sg_present ? observed(eq.state,bus,prefix*"τ_m") : NaN
                sgp=sg_present ? sn*observed(eq.state,bus,prefix*"P") : 0.0
                sgq=sg_present ? sn*observed(eq.state,bus,prefix*"Q") : 0.0
                pnet=observed(eq.state,bus,"busbar₊P_MW")
                qnet=observed(eq.state,bus,"busbar₊Q_MVAr")
                scale=rho==0 ? 0.0 : rho==1 ? 1.0 : observed(eq.state,bus,"gfl₊filter₊PortScale")
                sgpoles,gflpoles=local_pole_groups(eq.state,bus)
                push!(rows,(bus,rho,epsilon=1-rho,status="EVALUATED",sg_present,
                    gfl_present=rho>0,SG_Sn_MVA=sn,SG_H_s=h,SG_H_times_Sn_MVA_s=h*sn,
                    SG_torque_pu=torque,SG_P_MW=sgp,SG_Q_Mvar=sgq,
                    GFL_device_rating_MVA="NOT_DEFINED",GFL_port_scale=scale,
                    GFL_P_MW=pnet-sgp,GFL_Q_Mvar=qnet-sgq,network_P_MW=pnet,
                    expected_frozen_SG_P_MW=(1-rho)*Float64(data.bus[findfirst(==(bus),data.bus.bus),:P])*100,
                    bus_Vbase_kV=observed(eq.state,bus,"busbar₊Vbase"),
                    SG_internal_poles=sgpoles,GFL_internal_poles=gflpoles,error=""))
            catch e
                push!(rows,(bus,rho,epsilon=1-rho,status="ERROR",sg_present=rho<1,
                    gfl_present=rho>0,SG_Sn_MVA=NaN,SG_H_s=NaN,SG_H_times_Sn_MVA_s=NaN,
                    SG_torque_pu=NaN,SG_P_MW=NaN,SG_Q_Mvar=NaN,
                    GFL_device_rating_MVA="NOT_DEFINED",GFL_port_scale=NaN,
                    GFL_P_MW=NaN,GFL_Q_Mvar=NaN,network_P_MW=NaN,
                    expected_frozen_SG_P_MW=(1-rho)*Float64(data.bus[findfirst(==(bus),data.bus.bus),:P])*100,
                    bus_Vbase_kV=NaN,SG_internal_poles="",GFL_internal_poles="",
                    error=sprint(showerror,e)))
            end
            CSV.write(joinpath(OUT,"TABLE_M14_fractional_rho_semantics.csv"),DataFrame(rows))
            println("RHO_AUDIT bus=",bus," rho=",rho," status=",rows[end].status,
                " SG_P=",rows[end].SG_P_MW);flush(stdout)
        end
    end
end

main()
