module PDCasesM

using SHA, TOML, NetworkDynamics, PowerDynamics
include(joinpath(@__DIR__, "..", "pd39", "PD39.jl"))
using .PD39

export load_frozen_candidate, build_case

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const MODEL = PD39.PD39Model
include(joinpath(ROOT,"src","bnd_design_k","BNDDesignK.jl"))
const K_CONTEXT = Ref{Any}(nothing)
function k_context()
    if K_CONTEXT[] === nothing
        K_CONTEXT[] = BNDDesignK.design_context(ROOT)
    end
    return K_CONTEXT[]
end

function load_frozen_candidate(case::AbstractString)
    path = case == "ExpG_candidate" ? joinpath(ROOT, "reports", "experiment_G", "Z_G_FINAL.toml") :
        case == "ExpK_nominal" ? joinpath(ROOT, "reports", "experiment_K", "Z_K_NOMINAL_FINAL.toml") :
        error("not a frozen candidate: $case")
    bytes2hex(sha256(read(path))) == strip(read(path * ".sha256", String)) ||
        error("candidate SHA-256 differs from frozen sidecar: $path")
    return TOML.parsefile(path)
end

# Reproduce the builder in validate_powerdynamics_K.jl verbatim in a new audit
# namespace. This is intentionally the pre-diagnosis model, including its Vset.
function original_k_load_bus(nw, bus, rho, kp, ki)
    data = MODEL.ieee39_data()
    row = data.bus[findfirst(==(bus), data.bus.bus), :]
    row.has_load || error("expected colocalized ZIP load at $bus")
    MODEL.set_Sbase!(MODEL.BASE_MVA)
    MODEL.set_fbase!(MODEL.BASE_FREQ)
    xf=0.03; vdc=2.5; cdc=1.25; ftau=300.0; fi=600.0
    gfl=MODEL.WeightedSimpleGFLDC(
        Xf=xf, Rf=0.01, PLL_Kp=kp, PLL_Ki=ki, PLL_τ_lpf=1/(ftau*2pi),
        CC1_KP=(xf/(2pi*MODEL.BASE_FREQ))*(fi*2pi),
        CC1_KI=(xf/(2pi*MODEL.BASE_FREQ))*(fi*2pi)^2/4,
        CC1_F=0, CC1_Fcoupl=0, C_dc=cdc, V_dc=vdc,
        kp_v_dc=vdc*cdc*(5*2pi), ki_v_dc=vdc*cdc*(5*2pi)^2/4,
        port_scale=rho, name=:gfl)
    load=deepcopy(MODEL.OfficialIEEE39.load)
    parts = if rho < 1
        sg=deepcopy(bus==39 ? MODEL.OfficialIEEE39.uncontrolled_machine :
            MODEL.OfficialIEEE39.controlled_machine)
        MODEL.MTKBus(sg,gfl,load)
    else
        MODEL.MTKBus(gfl,load)
    end
    vertices, edges = MODEL.copy_network_components(nw)
    pf=MODEL.get_pfmodel(vertices[bus])
    replacement=MODEL.compile_bus(parts; name=Symbol("bus$(bus)"), vidx=bus, pf)
    MODEL.set_Sbase!(); MODEL.set_fbase!()
    MODEL.set_default!(replacement, Symbol("busbar₊Vbase"), Float64(row.base_kv))
    MODEL.OfficialIEEE39.apply_csv_params!(replacement, data.load, bus)
    ctx=k_context()
    lr=only(eachrow(ctx.net.load_audit[ctx.net.load_audit.bus .== bus,:]))
    vset=abs(ctx.net.voltage[bus]) * sqrt(lr.CSV_setpoint_load_MW/lr.initialized_load_MW)
    MODEL.set_default!(replacement, r"ZIPLoad₊Vset$", vset)
    if rho < 1
        MODEL.OfficialIEEE39.apply_csv_params!(replacement,data.machine,bus)
        Bool(row.has_avr) && MODEL.OfficialIEEE39.apply_csv_params!(replacement,data.avr,bus)
        Bool(row.has_gov) && MODEL.OfficialIEEE39.apply_csv_params!(replacement,data.gov,bus)
        mrow=data.machine[findfirst(==(bus),data.machine.bus),:]
        MODEL.set_default!(replacement,r"Sn$",(1-rho)*Float64(mrow.Sn))
    end
    vertices[bus]=replacement
    out=MODEL.Network(vertices,edges)
    MODEL.set_jac_prototype!(out)
    return out
end

function build_case(case::AbstractString)
    nw=PD39.baseline_network()
    case=="all_SG" && return nw
    c=load_frozen_candidate(case)
    for bus in 30:39
        i=bus-29
        rho=case=="ExpG_candidate" ? Float64(c["generator"][i]["rho"]) : Float64(c["rho"][i])
        kp=case=="ExpG_candidate" ? Float64(c["generator"][i]["Kp"]) : Float64(c["Kp"][i])
        ki=case=="ExpG_candidate" ? Float64(c["generator"][i]["Ki"]) : Float64(c["Ki"][i])
        rho==0 && continue
        if case=="ExpK_nominal" && bus in (31,39)
            nw=original_k_load_bus(nw,bus,rho,kp,ki)
        elseif rho==1
            nw=PD39.replace_bus(nw,bus;template=PD39.simple_gfldc_template())
        else
            nw=MODEL.weighted_replacement_network(nw,bus,rho,1.0)
        end
    end
    case=="ExpG_candidate" && return nw
    vertices, edges=MODEL.copy_network_components(nw)
    for bus in 30:39
        i=bus-29
        rho=Float64(c["rho"][i]); rho>0 || continue
        MODEL.set_default!(vertices[bus],r"PLL_Kp$",Float64(c["Kp"][i]))
        MODEL.set_default!(vertices[bus],r"PLL_Ki$",Float64(c["Ki"][i]))
    end
    out=MODEL.Network(vertices,edges)
    MODEL.set_jac_prototype!(out)
    return out
end

end
