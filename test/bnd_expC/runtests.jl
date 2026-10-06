using Test
using LinearAlgebra
using Statistics

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(ROOT,"src","bnd_graph","BNDGraph.jl"))
using .BNDGraph
const GDS=BNDGraph.GeneralizedDynamicSelfEnergy
const GB=BNDGraph.GraphBackbone
const IP=BNDGraph.IntermodalPathways
const RP=BNDGraph.RealPoleAnalysis
const RSA=BNDGraph.RealSystemAdapter

@testset "Experiment C generalized self-energy" begin
    M=Matrix(Diagonal([0.5,0.8,1.1]))
    D0=[0.2 0.01 0.0; -0.02 0.3 0.01; 0.0 -0.01 0.25]
    L0=[2.0 -1.0 -1.0; -0.8 1.7 -0.9; -1.2 -0.7 1.9]
    Avc=[0.2 0.1; -0.1 0.25; 0.05 -0.2]
    Acc=[-2.0 0.1; -0.2 -3.0]
    Acq=[0.1 0.2 0.0; -0.1 0.0 0.15]
    Acv=[0.0 0.1 0.0; 0.2 0.0 -0.1]
    op=GDS.AngleOperator(M,D0,L0,Avc,Acc,Acq,Acv)
    g=ones(3)
    bb=GB.primary_backbone(M,L0,g)
    basis=GB.generalized_basis(M,bb.LG)
    s=0.3+1.7im

    @test GDS.generalized_reconstruction_error(op,bb.LG,s) < 1e-12
    @test GDS.derivative_check(op,s) < 1e-7
    @test GB.mass_audit(M).spd
    @test !GB.mass_audit(Diagonal([1.0,-0.1,2.0])).spd
    @test norm(bb.LG*g) < 1e-12
    @test basis.M_orthogonality_error < 1e-12
    @test basis.L_diagonalization_error < 1e-12

    v=GDS.psi(op,bb.LG,s)
    phat=basis.Phi'*v*basis.Phi
    That=IP.graph_operator(basis.Lambda,phat,s)
    Tdirect=basis.Phi'*GDS.retained_operator(op,s)*basis.Phi
    @test norm(That-Tdirect)/norm(Tdirect) < 1e-12

    k=2
    gm=IP.schur_gamma(That,k)
    @test gm.residual < 1e-12
    pw=IP.pathway_matrix(phat,basis.Lambda,s,k)
    @test pw.reconstruction_error < 1e-12
    @test abs(pw.gamma-gm.gamma) < 1e-10

    ω=imag(s)
    z=ComplexF64[1+0.2im,-0.3+0.8im,0.5-0.1im]
    Psi=GDS.psi(op,bb.LG,im*ω)
    q=basis.Phi*z
    vel=im*ω*q
    direct=real(dot(vel,Psi*q))/2
    DG=IP.effective_dissipative_operator(basis.Phi'*Psi*basis.Phi,ω)
    graph=ω^2*real(dot(z,DG*z))/2
    @test abs(direct-graph)/max(1,abs(direct),abs(graph)) < 1e-12

    # The Sigma framework is recovered as the exact special case Psi=s*Sigma.
    Sigma=basis.Phi'*[0.4 0.1 0; -0.05 0.6 0.02; 0 -0.03 0.5]*basis.Phi
    Shat=Sigma
    DGspecial=IP.effective_dissipative_operator(im*ω*Shat,ω)
    @test norm(DGspecial-(Shat+Shat')/2) < 1e-12
end

@testset "Predictor sign and modal alignment" begin
    ω1=2.0
    ω2=2.15
    ε=0.02
    s0=im*ω1
    T(s)=[s^2+ω1^2 ε; ε s^2+ω2^2]
    gamma=IP.schur_gamma(T(s0),1).gamma
    shift=IP.pole_shift(gamma,2s0,s0)
    coeff=[1.0,ω1^2+ω2^2,ω1^2*ω2^2-ε^2]
    roots=ComplexF64[]
    for x in ((-coeff[2]+sqrt(complex(coeff[2]^2-4coeff[1]*coeff[3])))/(2coeff[1]),
              (-coeff[2]-sqrt(complex(coeff[2]^2-4coeff[1]*coeff[3])))/(2coeff[1]))
        push!(roots,sqrt(x),-sqrt(x))
    end
    actual=roots[argmin(abs.(roots .- s0))]
    @test abs(actual-shift.predicted) < 5e-4

    M=Matrix(Diagonal([1.0,2.0]))
    L=Matrix(Diagonal([0.0,3.0]))
    b=GB.generalized_basis(M,L)
    overlap=GB.mode_overlap(b.Phi,M,b.Phi,M)
    @test minimum(diag(overlap)) > 1-1e-12
    @test length(GB.degenerate_clusters([0.0,1.0,1.0+1e-10])) == 2
    @test all(GB.principal_angles(b.Phi[:,1:1],b.Phi[:,1:1]) .< 1e-8)
    @test isfinite(RP.spearman_rank([1.0,2,3],[2.0,4,7]))
    clusters=[[1],[2,3]]
    Shat=ComplexF64[0.2 0.03 -0.04; -0.02 0.4 0.01; 0.05 -0.03 0.35]
    block=IP.block_pathway_matrix(Shat,[0.0,1.0,1.0+1e-10],0.2+1.8im,1,clusters)
    @test block.reconstruction_error < 1e-12
    blockdeg=IP.block_pathway_matrix(Shat,[0.0,1.0,1.0+1e-10],0.2+1.8im,2,clusters)
    @test size(blockdeg.gamma)==(2,2)
    @test blockdeg.reconstruction_error < 1e-12
end

@testset "The real ExpA singularity is inspected without Pi(0)" begin
    input=RSA.load_expA_primary(ROOT)
    σ=svdvals(input.matrices.Acc)
    @test minimum(σ) <= 1e-8*maximum(σ)
    runner=read(joinpath(ROOT,"experiments","bnd_expC","run_experiment_C.jl"),String)
    @test !occursin(r"pi_q\s*\(\s*op\s*,\s*0(?:\.0*)?\s*\)",runner)
    @test !occursin("pinv(",runner)
end
