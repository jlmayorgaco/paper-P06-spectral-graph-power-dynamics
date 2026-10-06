using CSV, DataFrames, DelimitedFiles, LinearAlgebra, NetworkDynamics, PowerDynamics
using Random, TOML

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const OUT=joinpath(ROOT,"reports","experiment_N")
const CASES=joinpath(OUT,"validation_cases")
include(joinpath(ROOT,"src","bnd_model_expN","PDReferenceN.jl"))
include(joinpath(ROOT,"src","bnd_model_expN","PDExactDesignN.jl"))

function save_matrix(path, matrix)
    writedlm(path,Matrix(matrix),',')
end

function save_system(prefix,sys)
    n=size(sys.A,1)
    m=sys.M isa UniformScaling ? Matrix{Float64}(I,n,n) : Matrix(sys.M)
    for (name,val) in (("M",m),("A",sys.A),("B",sys.B),("C",sys.C),("D",sys.D))
        val===nothing && continue
        save_matrix(prefix*"_"*name*".csv",val)
    end
end

function export_case(label,nw,base,ctx,rho,kp,ki; development=false)
    path=joinpath(CASES,label);mkpath(path)
    s=PDReferenceN.trim_state(nw,base,rho,kp,ki)
    res=PDReferenceN.residual_audit(nw,s)
    power=PDReferenceN.direct_power_audit(s,base,rho)
    CSV.write(joinpath(path,"PD_power_audit.csv"),power)
    pd=linearize_network(s)
    save_system(joinpath(path,"PD"),pd)
    mass=pd.M isa UniformScaling ? ones(size(pd.A,1)) : diag(pd.M)
    CSV.write(joinpath(path,"PD_state_map.csv"),DataFrame(
        index=1:length(pd.sym),state_name=string.(pd.sym),differential=mass.==1))
    λpd=jacobian_eigenvals(pd)
    CSV.write(joinpath(path,"PD_poles.csv"),DataFrame(real=real.(λpd),imag=imag.(λpd)))
    analytic=PDExactDesignN.spectrum(ctx,rho,kp,ki)
    m=analytic.model
    for (name,val) in (("Ared",m.Ared),("A_local",m.A),("B_port",m.B),
                       ("C_port",m.C),("D_port",m.D))
        save_matrix(joinpath(path,"AN_"*name*".csv"),val)
    end
    CSV.write(joinpath(path,"AN_state_map.csv"),m.state_inventory)
    CSV.write(joinpath(path,"AN_poles.csv"),
        DataFrame(real=real.(analytic.lambda),imag=imag.(analytic.lambda)))
    active_buses=[b for b in 30:39 if rho[b-29]>0]
    for bus in active_buses
        local_port=linearize_component(s,VIndex(bus))
        save_system(joinpath(path,"PD_bus$(bus)"),local_port)
    end
    CSV.write(joinpath(path,"case_definition.csv"),DataFrame(
        bus=collect(30:39),rho=Float64.(rho),Kp=Float64.(kp),Ki=Float64.(ki)))
    max_p=maximum(Float64.(power.max_P_error_pu))
    max_q=maximum(Float64.(power.max_Q_error_pu))
    bounds_ok=all(power.bounds_status .== "PASS")
    row=(;case=label,development,PD_dynamic=count(==(1.0),mass),
        AN_dynamic=m.n_dynamic,PD_algebraic=count(==(0.0),mass),
        residual_max=res.maximum,residual_norm=res.norm,
        max_P_error_pu=max_p,max_Q_error_pu=max_q,
        bounds_ok,PD_raw_alpha=maximum(real.(λpd)),AN_alpha=analytic.alpha,
        gauge_residual=analytic.gauge_residual,active_GFL_buses=length(active_buses))
    println("N_VALIDATED_CASE ",label," residual=",res.maximum,
        " P_error=",max_p," Q_error=",max_q," AN_alpha=",analytic.alpha)
    flush(stdout)
    row
end

function candidate_vectors(name)
    k0=PDExactDesignN.K0P;i0=PDExactDesignN.K0I
    if name=="ExpE_corrected"
        t=CSV.read(joinpath(ROOT,"reports","experiment_E","tables",
            "TABLE_E14_provisional_per_generator_design.csv"),DataFrame)
        sort!(t,:bus)
        return Float64.(t.rho),Float64.(t.Kp),Float64.(t.Ki)
    elseif name=="ExpG_corrected"
        t=TOML.parsefile(joinpath(ROOT,"reports","experiment_G","Z_G_FINAL.toml"))
        gens=sort(t["generator"],by=x->Int(x["bus"]))
        return [Float64(g["rho"]) for g in gens],
            [Float64(g["Kp"]) for g in gens],
            [Float64(g["Ki"]) for g in gens]
    elseif name=="ExpK_corrected"
        t=TOML.parsefile(joinpath(ROOT,"reports","experiment_K",
            "Z_K_NOMINAL_FINAL.toml"))
        return Float64.(t["rho"]),Float64.(t["Kp"]),Float64.(t["Ki"])
    else
        error("unknown candidate $name")
    end
end

function main()
    mkpath(CASES)
    println("N_VALIDATION_BASELINE");flush(stdout)
    base=PDReferenceN.frozen_baseline()
    ctx=PDExactDesignN.design_context(ROOT)
    k0=PDExactDesignN.K0P;i0=PDExactDesignN.K0I
    rows=NamedTuple[]
    push!(rows,export_case("all_SG",base.nw,base,ctx,zeros(10),fill(k0,10),fill(i0,10)))
    for bus in 30:39
        seed=zeros(10);seed[bus-29]=0.5
        println("N_VALIDATION_BUILD_SINGLE ",bus);flush(stdout)
        nw=PDReferenceN.build_architecture(base,seed)
        for (r,km,im,label) in ((0.10,k0,i0,"r010_nom"),
                                (0.50,k0,i0,"r050_nom"),
                                (0.90,k0,i0,"r090_nom"),
                                (0.99,k0,i0,"r099_nom"),
                                (0.10,0.25k0,4i0,"r010_Kpmin_Kimax"),
                                (0.90,4k0,0.25i0,"r090_Kpmax_Kimin"))
            rho=zeros(10);rho[bus-29]=r
            kp=fill(k0,10);ki=fill(i0,10)
            kp[bus-29]=km;ki[bus-29]=im
            label_full="bus$(bus)_"*label
            push!(rows,export_case(label_full,nw,base,ctx,rho,kp,ki;
                development=(bus==38 && r==0.50 && km==k0 && im==i0)))
        end
    end
    for name in ("ExpE_corrected","ExpG_corrected","ExpK_corrected")
        rho,kp,ki=candidate_vectors(name)
        println("N_VALIDATION_BUILD_MULTI ",name);flush(stdout)
        nw=PDReferenceN.build_architecture(base,rho,kp,ki)
        push!(rows,export_case(name,nw,base,ctx,rho,kp,ki))
    end
    rng=MersenneTwister(20260930)
    seed=fill(0.5,10)
    println("N_VALIDATION_BUILD_RANDOM_ARCHITECTURE");flush(stdout)
    nw=PDReferenceN.build_architecture(base,seed)
    for j in 1:5
        rho=0.1 .+ 0.8 .* rand(rng,10)
        kp=k0 .* (0.25 .+ 3.75 .* rand(rng,10))
        ki=i0 .* (0.25 .+ 3.75 .* rand(rng,10))
        push!(rows,export_case("random_$(j)",nw,base,ctx,rho,kp,ki))
    end
    CSV.write(joinpath(OUT,"TABLE_N05_validation_case_export.csv"),DataFrame(rows))
    println("N_VALIDATION_EXPORT_DONE cases=",length(rows))
end

main()
