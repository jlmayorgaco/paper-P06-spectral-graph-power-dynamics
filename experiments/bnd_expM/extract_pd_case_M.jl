using CSV, DataFrames, DelimitedFiles, Graphs, LinearAlgebra, NetworkDynamics, PowerDynamics
include(joinpath(@__DIR__, "..", "..", "src", "bnd_model_audit_m", "PDCasesM.jl"))
using .PDCasesM

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const OUT = joinpath(ROOT, "reports", "experiment_M")
const CASE = length(ARGS) == 1 ? ARGS[1] : error("usage: extract_pd_case_M.jl case")
CASE in ("all_SG", "ExpG_candidate", "ExpK_nominal") || error("invalid case")
const CASEDIR = joinpath(OUT,"matrices",CASE)
mkpath(CASEDIR)

function matrixcsv(prefix, sys)
    n=size(sys.A,1)
    m=sys.M isa UniformScaling ? Matrix{Float64}(I,n,n) : Matrix(sys.M)
    for (name,value) in (("M",m),("A",sys.A),("B",sys.B),("C",sys.C),("D",sys.D))
        value===nothing && continue
        writedlm(prefix*"_"*name*".csv",Matrix(value),',')
    end
end

function defaults_table(component, case, category, id)
    [(;case,category,id,parameter=string(k),value=string(v),source="compiled_defaults")
     for (k,v) in get_defaults_dict(component)]
end

function main()
    println("BUILD_CASE ",CASE);flush(stdout)
    nw=PDCasesM.build_case(CASE)
    println("INITIALIZE ",CASE);flush(stdout)
    eq=PDCasesM.PD39.initialize_equilibrium(nw;sparse=false,check=:error)
    println("LINEARIZE ",CASE);flush(stdout)
    sys=linearize_network(eq.state)
    matrixcsv(joinpath(CASEDIR,"PD"),sys)
    mass=sys.M isa UniformScaling ? ones(size(sys.A,1)) : diag(sys.M)
    smap=DataFrame(index=1:length(sys.sym),state_name=string.(sys.sym),
        differential=mass.==1.0,mass_diagonal=mass,
        equilibrium_value=Float64.(uflat(eq.state)))
    CSV.write(joinpath(CASEDIR,"PD_state_map.csv"),smap)
    parrows=NamedTuple[];dimrows=NamedTuple[]
    for i in 1:39
        component=nw[VIndex(i)]
        append!(parrows,defaults_table(component,CASE,"bus",i))
        try
            li=linearize_component(eq.state,VIndex(i))
            matrixcsv(joinpath(CASEDIR,"PD_bus$(i)"),li)
            push!(dimrows,(case=CASE,component="bus$(i)",kind="bus",status="OK",
                states=size(li.A,1),inputs=size(li.B,2),outputs=size(li.C,1),error=""))
        catch e
            push!(dimrows,(case=CASE,component="bus$(i)",kind="bus",status="ERROR",
                states=0,inputs=0,outputs=0,error=sprint(showerror,e)))
        end
    end
    for i in 1:Graphs.ne(nw)
        component=nw[EIndex(i)]
        append!(parrows,defaults_table(component,CASE,"line",i))
        try
            li=linearize_component(eq.state,EIndex(i))
            matrixcsv(joinpath(CASEDIR,"PD_line$(i)"),li)
            push!(dimrows,(case=CASE,component="line$(i)",kind="line",status="OK",
                states=size(li.A,1),inputs=size(li.B,2),outputs=size(li.C,1),error=""))
        catch e
            push!(dimrows,(case=CASE,component="line$(i)",kind="line",status="ERROR",
                states=0,inputs=0,outputs=0,error=sprint(showerror,e)))
        end
    end
    CSV.write(joinpath(CASEDIR,"component_parameters.csv"),DataFrame(parrows))
    CSV.write(joinpath(CASEDIR,"component_dimensions.csv"),DataFrame(dimrows))
    psyms=NetworkDynamics.SII.parameter_symbols(nw)
    pvals=pflat(eq.state)
    length(psyms)==length(pvals) || error("parameter symbol/value length mismatch")
    CSV.write(joinpath(CASEDIR,"PD_parameter_values.csv"),
        DataFrame(parameter=string.(psyms),value=Float64.(pvals)))
    λ=jacobian_eigenvals(sys)
    CSV.write(joinpath(CASEDIR,"PD_poles.csv"),DataFrame(real=real.(λ),imag=imag.(λ)))
    println("CASE_DONE ",CASE," n=",size(sys.A,1)," alpha=",maximum(real.(λ)),
        " fixed_point=",eq.fixed_point);flush(stdout)
    for bus in (31,39)
        d=get_defaults_dict(nw[VIndex(bus)])
        println("BUS ",bus," ZIP_DEFAULTS ",filter(p->occursin("ZIPLoad",string(first(p))),collect(d)))
        for name in ("busbar₊P_MW","busbar₊Q_Mvar","ZIPLoad₊P","ZIPLoad₊Q",
                     "ZIPLoad₊Vset","ctrld_gen₊machine₊P","machine₊P")
            try
                println("BUS ",bus," ",name,"=",eq.state[VIndex(bus,Symbol(name))])
            catch
            end
        end
    end
    pqrows=NamedTuple[]
    for bus in (31,39)
        u=complex(eq.state[VIndex(bus,:busbar₊u_r)],eq.state[VIndex(bus,:busbar₊u_i)])
        pnet=Float64(eq.state[VIndex(bus,:busbar₊P_MW)])
        qnet=Float64(eq.state[VIndex(bus,:busbar₊Q_MVAr)])
        pload=100Float64(eq.state[VIndex(bus,:ZIPLoad₊P)])
        qload=100Float64(eq.state[VIndex(bus,:ZIPLoad₊Q)])
        vset=Float64(eq.state[VIndex(bus,:ZIPLoad₊Vset)])
        # Busbar current is the sum of all injector terminal currents. For a
        # mixed bus this sum is the separately observed net less ZIP power.
        pgen=pnet-pload; qgen=qnet-qload
        push!(pqrows,(case=CASE,bus,p_gen_MW=pgen,q_gen_Mvar=qgen,
            p_zip_MW=pload,q_zip_Mvar=qload,p_net_MW=pnet,q_net_Mvar=qnet,
            p_balance_MW=pgen+pload-pnet,q_balance_Mvar=qgen+qload-qnet,
            ZIP_Vset_pu=vset,voltage_pu=abs(u),
            generator_source="busbar_observable_minus_ZIP_observable"))
    end
    CSV.write(joinpath(CASEDIR,"component_PQ_sharing.csv"),DataFrame(pqrows))
    shares=NamedTuple[]
    for bus in 30:39
        hasload=bus in (31,39)
        pnet=Float64(eq.state[VIndex(bus,:busbar₊P_MW)])
        qnet=Float64(eq.state[VIndex(bus,:busbar₊Q_MVAr)])
        pload=hasload ? 100Float64(eq.state[VIndex(bus,:ZIPLoad₊P)]) : 0.0
        qload=hasload ? 100Float64(eq.state[VIndex(bus,:ZIPLoad₊Q)]) : 0.0
        sgp=0.0;sgq=0.0;sn=0.0
        machine_prefix=bus==39 ? "machine₊" : "ctrld_gen₊machine₊"
        try
            sn=Float64(eq.state[VIndex(bus,Symbol(machine_prefix*"Sn"))])
            sgp=sn*Float64(eq.state[VIndex(bus,Symbol(machine_prefix*"P"))])
            sgq=sn*Float64(eq.state[VIndex(bus,Symbol(machine_prefix*"Q"))])
        catch e
            if CASE=="all_SG" ||
               (CASE=="ExpG_candidate" && bus in (30,31,32,33,34,35,36,39)) ||
               (CASE=="ExpK_nominal" && bus in (38,39))
                error("SG observable missing at bus $bus: $(sprint(showerror,e))")
            end
        end
        expected=CASE=="all_SG" ? NaN : begin
            c=PDCasesM.load_frozen_candidate(CASE)
            CASE=="ExpG_candidate" ? Float64(c["generator"][bus-29]["retained_SG_MW"]) :
                Float64(c["P_initial_MW"][bus-29])*Float64(c["epsilon"][bus-29])
        end
        push!(shares,(case=CASE,bus,sg_rating_mva=sn,sg_p_MW=sgp,sg_q_Mvar=sgq,
            gfl_p_MW=pnet-pload-sgp,gfl_q_Mvar=qnet-qload-sgq,
            zip_p_MW=pload,zip_q_Mvar=qload,net_p_MW=pnet,net_q_Mvar=qnet,
            expected_analytic_sg_p_MW=expected,
            sg_p_difference_MW=sgp-expected))
    end
    CSV.write(joinpath(CASEDIR,"PD_share_audit.csv"),DataFrame(shares))
    try
        ol=open_loop_linearization(eq.state)
        matrixcsv(joinpath(CASEDIR,"PD_Zbus"),ol.Zbus)
        matrixcsv(joinpath(CASEDIR,"PD_Ynw"),ol.Ynw)
        if ol.Yinj !== nothing
            matrixcsv(joinpath(CASEDIR,"PD_Yinj"),ol.Yinj)
        end
        # The installed NetworkDynamics 1.3.0 tests use positive feedback for
        # this sign convention even though the prose example omits the flag.
        cl=NetworkDynamics.feedback(ol.Zbus,
            ol.Yinj===nothing ? ol.Ynw : ol.Ynw+ol.Yinj;pos=true)
        λcl=jacobian_eigenvals(cl)
        CSV.write(joinpath(CASEDIR,"PD_closure_poles.csv"),
            DataFrame(real=real.(λcl),imag=imag.(λcl)))
        println("CLOSURE_DONE ",CASE," count=",length(λcl))
    catch e
        write(joinpath(CASEDIR,"PD_closure_error.txt"),sprint(showerror,e))
        println("CLOSURE_ERROR ",CASE," ",sprint(showerror,e))
    end
end

main()
