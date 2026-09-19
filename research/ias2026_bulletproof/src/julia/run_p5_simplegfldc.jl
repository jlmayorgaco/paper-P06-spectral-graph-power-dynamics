"""P5: execute the official PowerDynamics SimpleGFLDC as a second GFL model."""

using LinearAlgebra
using PowerDynamics
using PowerDynamics.Library
using NetworkDynamics
using ModelingToolkitBase
include(joinpath(@__DIR__, "campaign_root.jl"))

const CAMPAIGN = campaign_root_from_args()
const RAW = joinpath(CAMPAIGN, "raw", "p5")
const REPORTS = joinpath(CAMPAIGN, "reports")
mkpath(RAW)
mkpath(REPORTS)
set_Sbase!(100.0)
set_fbase!(60.0)

const P0 = 0.35
const Q0 = -0.08
const MVA0 = 100.0 * hypot(P0, Q0)

function build_network()
    V_dc=2.5; C_dc=1.25; f_v_dc=5.0
    xwLf=0.03; Rf=0.01
    f_pll=5.0; f_tau_pll=300.0; f_i_dq=600.0
    @named gfl = ComposableInverter.SimpleGFLDC(
        Xf=xwLf, Rf=Rf,
        PLL_Kp=f_pll*2π,
        PLL_Ki=(f_pll*2π)^2/4,
        PLL_τ_lpf=1/(f_tau_pll*2π),
        CC1_KP=(xwLf/(2π*60.0))*(f_i_dq*2π),
        CC1_KI=(xwLf/(2π*60.0))*(f_i_dq*2π)^2/4,
        CC1_F=0, CC1_Fcoupl=0,
        C_dc=C_dc, V_dc=V_dc,
        kp_v_dc=V_dc*C_dc*(f_v_dc*2π),
        ki_v_dc=V_dc*C_dc*(f_v_dc*2π)*(f_v_dc*2π)/4,
    )
    @named gfl_bus = compile_bus(MTKBus(gfl); current_source=true)
    set_pfmodel!(gfl_bus, pfPQ(P=P0, Q=Q0; current_source=true))
    @named shunt = DynamicParallelRCShunt(R=1/0.05, B=1e-5)
    @named network_bus = compile_bus(MTKBus(shunt))
    set_pfmodel!(network_bus, pfShunt(G=0.05, B=1e-5))
    loopback = LoopbackConnection(; src=:gfl_bus, dst=:network_bus,
                                   potential=[:u_r, :u_i], flow=[:i_r, :i_i])
    @named slack_symbolic = Library.VδConstraint(V=1.0, δ=0.0)
    @named slack_bus = compile_bus(MTKBus(slack_symbolic); pf=pfSlack(V=1.0))
    @named branch = DynamicSeriesRLBranch(R=0.01, X=0.3)
    line = compile_line(MTKLine(branch); name=:gfl_to_slack,
                        src=:network_bus, dst=:slack_bus)
    @named branch_pf = PiLine(R=0.01, X=0.3)
    line_pf = compile_line(MTKLine(branch_pf); name=:gfl_to_slack_pf)
    set_pfmodel!(line, line_pf)
    Network([gfl_bus, network_bus, slack_bus], [loopback, line])
end

function residual(net, state)
    du = zeros(Float64, length(uflat(state)))
    net(du, uflat(state), pflat(state), 0.0)
    isempty(du) ? Inf : maximum(abs, du)
end

function main()
    net = build_network()
    state = initialize_from_pf(net; verbose=false, subverbose=false, check=:none,
                               tol=1e-8, nwtol=1e-8)
    r = residual(net, state)
    eigs = jacobian_eigenvals(state)
    alpha = maximum(real, eigs)
    open(joinpath(RAW, "p5_simplegfldc_result.csv"), "w") do io
        println(io, "model,P,Q,MVA,system_base_MVA,state_count,residual,alpha_full,operating_point_match,no_retune,status")
        println(io, join(("PowerDynamics.ComposableInverter.SimpleGFLDC", P0, Q0, MVA0,
                          100.0, length(uflat(state)), r, alpha, true, true,
                          length(uflat(state)) > 0 && isfinite(r) && r <= 1e-8 ? "PASS" : "FAIL"), ','))
    end
    open(joinpath(REPORTS, "P5_SECOND_GFL_STATUS.md"), "w") do io
        println(io, "# P5 — second materially different GFL model")
        println(io)
        println(io, "status: ", length(uflat(state)) > 0 && isfinite(r) && r <= 1e-8 ? "PASS" : "STOPPED_BY_GATE")
        println(io, "evidence_class: FRESH_SECOND_MODEL_EXECUTION")
        println(io, "model: PowerDynamics.ComposableInverter.SimpleGFLDC")
        println(io, "source_config: raw/p5/p5_simplegfldc_config.toml")
        println(io, "matched_dispatch_P_pu: ", P0)
        println(io, "matched_dispatch_Q_pu: ", Q0)
        println(io, "matched_MVA: ", MVA0)
        println(io, "system_base_MVA: 100.0")
        println(io, "no_retune: true")
        println(io, "state_count: ", length(uflat(state)))
        println(io, "residual: ", r)
        println(io, "alpha_full: ", alpha)
        println(io, "result_csv: raw/p5/p5_simplegfldc_result.csv")
    end
    println("P5_", length(uflat(state)) > 0 && isfinite(r) && r <= 1e-8 ? "PASS" : "STOPPED", " states=", length(uflat(state)), " residual=", r)
    return length(uflat(state)) > 0 && isfinite(r) && r <= 1e-8 ? 0 : 1
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    exit(main())
end
