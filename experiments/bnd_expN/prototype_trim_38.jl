# First compiled-network test of the explicit ExpN trim, bus 38 at rho=0.5.
using LinearAlgebra, NetworkDynamics, PowerDynamics
const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(ROOT, "src", "pd39", "PD39.jl"))
using .PD39

function main()
    bus=38; rho=0.5
    println("N_PROTOTYPE_BASELINE"); flush(stdout)
    nw0=PD39.baseline_network()
    eq=PD39.initialize_equilibrium(nw0; sparse=false,check=:error)
    s0=eq.state
    println("N_PROTOTYPE_MIXED_BUILD"); flush(stdout)
    nw=PD39.weighted_replacement_network(nw0,bus,rho,1.0)
    s=NWState(nw)
    mixedvars=Dict(string(x)=>x for x in NetworkDynamics.SII.variable_symbols(nw))
    mixedpars=Dict(string(x)=>x for x in NetworkDynamics.SII.parameter_symbols(nw))
    copied_states=0; copied_params=0
    for (sym,val) in zip(NetworkDynamics.SII.variable_symbols(nw0),uflat(s0))
        k=string(sym)
        if haskey(mixedvars,k)
            s[mixedvars[k]]=val
            copied_states+=1
        end
    end
    for (sym,val) in zip(NetworkDynamics.SII.parameter_symbols(nw0),pflat(s0))
        k=string(sym)
        if haskey(mixedpars,k) && k!="VIndex($bus, :ctrld_gen₊machine₊Sn)"
            s[mixedpars[k]]=val
            copied_params+=1
        end
    end
    ur=Float64(s0[VIndex(bus,:busbar₊u_r)])
    ui=Float64(s0[VIndex(bus,:busbar₊u_i)])
    sn0=Float64(s0[VIndex(bus,:ctrld_gen₊machine₊Sn)])
    p=sn0*Float64(s0[VIndex(bus,:ctrld_gen₊machine₊P)])/100
    q=sn0*Float64(s0[VIndex(bus,:ctrld_gen₊machine₊Q)])/100
    v2=ur^2+ui^2
    ir=(p*ur+q*ui)/v2
    ii=(p*ui-q*ur)/v2
    vm=sqrt(v2); theta=atan(ui,ur)
    id=p/vm; iq=-q/vm
    rf=0.01; xf=0.03
    vip_d=vm+rf*id-xf*iq
    vip_q=rf*iq+xf*id
    cki=(xf/(2pi*60))*(600*2pi)^2/4
    trimmed_states=Dict(
        "gfl₊filter₊i_f_r"=>ir, "gfl₊filter₊i_f_i"=>ii,
        "gfl₊pll₊θ"=>theta, "gfl₊pll₊Δω_rad_s"=>0.0,
        "gfl₊pll₊Δω_i_rad_s"=>0.0,
        "gfl₊cc1₊γ_d"=>vip_d/cki, "gfl₊cc1₊γ_q"=>vip_q/cki,
        "gfl₊v_dc_state"=>2.5, "gfl₊v_dc_i"=>id)
    for (k,v) in trimmed_states
        s[VIndex(bus,Symbol(k))]=v
    end
    s[VIndex(bus,:gfl₊iset_q)]=iq
    s[VIndex(bus,:gfl₊P_dc)]=vip_d*id+vip_q*iq
    println("N_PROTOTYPE_RESIDUAL"); flush(stdout)
    du=similar(uflat(s)); nw(du,uflat(s),pflat(s),0.0)
    sgp=Float64(s[VIndex(bus,:ctrld_gen₊machine₊Sn)])*
        Float64(s[VIndex(bus,:ctrld_gen₊machine₊P)])
    sgq=Float64(s[VIndex(bus,:ctrld_gen₊machine₊Sn)])*
        Float64(s[VIndex(bus,:ctrld_gen₊machine₊Q)])
    gflp=100*rho*(ur*ir+ui*ii)
    gflq=100*rho*(ui*ir-ur*ii)
    println("N_PROTOTYPE_DONE copied_states=",copied_states," copied_params=",copied_params,
        " residual_max=",maximum(abs,du)," residual_norm=",norm(du),
        " SG_P=",sgp," SG_Q=",sgq," GFL_P=",gflp," GFL_Q=",gflq,
        " target_SG_P=",(1-rho)*100p," target_GFL_P=",rho*100p,
        " machine_tau=",s[VIndex(bus,:ctrld_gen₊machine₊τ_m)])
end

main()
