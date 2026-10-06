using Test
using LinearAlgebra
using DataFrames

include(joinpath(@__DIR__,"..","..","src","bnd_graphpll","PhysicalGraph.jl"))
include(joinpath(@__DIR__,"..","..","src","bnd_graphpll","IdealDamping.jl"))
include(joinpath(@__DIR__,"..","..","src","bnd_graphpll","ModalGFL.jl"))
include(joinpath(@__DIR__,"..","..","src","bnd_graphpll","ModalGainLaw.jl"))
include(joinpath(@__DIR__,"..","..","src","bnd_graphpll","GraphFilterRealization.jl"))
include(joinpath(@__DIR__,"..","..","src","bnd_graphpll","NodeVariantGraphLaw.jl"))
include(joinpath(@__DIR__,"..","..","src","bnd_graphpll","SelfEnergyCorrection.jl"))
include(joinpath(@__DIR__,"..","..","src","bnd_graphpll","GraphPLLContinuation.jl"))
include(joinpath(@__DIR__,"..","..","src","bnd_graphpll","SafetyMetrics.jl"))
include(joinpath(@__DIR__,"..","..","src","bnd_validation","GraphPLLPowerDynamicsValidation.jl"))
using .PhysicalGraph
using .IdealDamping
using .ModalGFL

@testset "BND graph PLL analytical foundations" begin
    branches=DataFrame(src_bus=[1],dst_bus=[2],R=[0.02],X=[0.15],r_src=[1.0],
        G_src=[0.0],G_dst=[0.0],B_src=[0.0],B_dst=[0.0])
    Y=assemble_ybus(branches,2)
    red=kron_reduce(Y,[1])
    audit=kron_identity_audit(Y,red;trials=3)
    graph=passive_port_graph(Y,[1])
    @test audit.relative_residual < 1e-12
    @test spectral_summary(graph.Lc).psd
    @test graph.hermitian_residual < 1e-12

    M=[2.0 0.3;0.3 1.4]
    K=[3.0 -0.4;-0.4 1.7]
    d=modal_critical_damping(M,K)
    @test minimum(d.nu)>0
    @test norm(d.Minvhalf*d.D*d.Minvhalf - 2*d.U*Diagonal(sqrt.(d.nu))*d.U',Inf)<1e-10
    for ν in d.nu
        dc=2sqrt(ν)
        @test isapprox(isolated_decay_rate(ν,dc),sqrt(ν);rtol=1e-12)
        @test isolated_decay_rate(ν,dc)>=isolated_decay_rate(ν,0.7*dc)-1e-12
        @test isolated_decay_rate(ν,dc)>=isolated_decay_rate(ν,1.8*dc)-1e-12
    end

    Yp=ComplexF64[2-im  -1+0im; -1+0im  2-im]
    Lc=real.(Yp)
    modal=modal_structure(Lc,Yp)
    @test modal.classification=="EXACT_SCALAR_MODAL"
    @test size(realify_admittance(Yp))==(4,4)
end

@testset "Downstream stop gates and exact correction utility" begin
    @test ModalGainLaw.gate_status().status=="BLOCKED_UPSTREAM_GATE"
    @test GraphFilterRealization.gate_status().status=="BLOCKED_UPSTREAM_GATE"
    @test NodeVariantGraphLaw.gate_status().status=="BLOCKED_UPSTREAM_GATE"
    @test SelfEnergyCorrection.gate_status().status=="DIAGNOSTIC_ONLY"
    @test GraphPLLContinuation.gate_status().status=="BLOCKED_UPSTREAM_GATE"
    @test SafetyMetrics.hard_current_limit_status()=="NOT_MODELED: SimpleGFLDC has no hard current limiter or validated Imax parameter."
    J=[1.0 0.0;0.0 1.0]; Rinv=Matrix{Float64}(I,2,2); r=[2.0,-3.0]
    @test norm(SelfEnergyCorrection.minimum_metric_correction(J,r,Rinv)+r)<1e-12
end
