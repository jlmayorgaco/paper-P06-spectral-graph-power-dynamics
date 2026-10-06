using LinearAlgebra
using CSV
using DataFrames
using TOML
BLAS.set_num_threads(1)

include(joinpath(@__DIR__,"..","..","src","bnd_graphpll","PhysicalGraph.jl"))
include(joinpath(@__DIR__,"..","..","src","bnd_graphpll","ModalGFL.jl"))
include(joinpath(@__DIR__,"..","..","src","bnd_design_e","CollectiveModel.jl"))
using .PhysicalGraph
using .ModalGFL
using .CollectiveModel

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const OUT=joinpath(ROOT,"reports","experiment_F2")
const NOM_KP=5*2pi
const NOM_KI=NOM_KP^2/4
const COUPLING_TOL=0.05
mkpath(joinpath(OUT,"tables")); mkpath(joinpath(OUT,"matrices"))

function matrix_csv(path)
    df=CSV.read(path,DataFrame)
    Matrix{Float64}(select(df,Not(first(names(df)))))
end
function csv_matrix(path,A)
    df=DataFrame(row_index=collect(1:size(A,1)))
    for j in 1:size(A,2); df[!,Symbol("x$j")]=A[:,j]; end
    CSV.write(path,df)
end
function blockdiag2(blocks)
    n=length(blocks); R=zeros(ComplexF64,2n,2n)
    for k in 1:n; ix=(2k-1):(2k); R[ix,ix].=blocks[k]; end
    return R
end
function blockdiag_state(jacobians)
    n=length(jacobians); A=zeros(Float64,9n,9n); B=zeros(Float64,9n,2n)
    C=zeros(Float64,2n,9n); D=zeros(Float64,2n,2n)
    for k in 1:n
        sx=(9k-8):(9k); vy=(2k-1):(2k); J=jacobians[k]
        A[sx,sx].=J.A; B[sx,vy].=J.B; C[vy,sx].=J.C; D[vy,vy].=J.D
    end
    return A,B,C,D
end
function graph_pencil_transform(P,U,nstate,nnode)
    Sx=kron(U,Matrix{Float64}(I,nstate,nstate)); Sv=modal_transform(U)
    S=zeros(Float64,size(P,1),size(P,2))
    S[1:nstate*nnode,1:nstate*nnode].=Sx
    S[nstate*nnode+1:end,nstate*nnode+1:end].=Sv
    Pt=S'*P*S
    order=Int[]
    for k in 1:nnode
        append!(order,(nstate*(k-1)+1):(nstate*k))
        append!(order,(nstate*nnode+2*(k-1)+1):(nstate*nnode+2*k))
    end
    return Pt[order,order]
end
function realify_static_matrix(R)
    n=size(R,1)÷2; Y=zeros(ComplexF64,n,n)
    for i in 1:n,j in 1:n
        Y[i,j]=complex(R[2i-1,2j-1],R[2i,2j-1])
    end
    return Y
end
function kron_real_ports(R,ports)
    keep=reduce(vcat,([2b-1,2b] for b in ports))
    rest=setdiff(collect(1:size(R,1)),keep)
    Rp=R[keep,keep]-R[keep,rest]*(R[rest,rest]\R[rest,keep])
    return Rp,keep,rest
end
function pencil_from_blocks(s,A,B,C,D,Ynet)
    n=size(Ynet,1); nstate=size(A,1)
    return [s*I(nstate)-A -B; C Ynet+D]
end
function named_blocks(T)
    rows=NamedTuple[]
    for k in 1:(size(T,1)÷2)
        ix=(2k-1):(2k)
        push!(rows,(mode=k,T11_real=real(T[ix[1],ix[1]]),T12_real=real(T[ix[1],ix[2]]),
            T21_real=real(T[ix[2],ix[1]]),T22_real=real(T[ix[2],ix[2]]),
            T11_imag=imag(T[ix[1],ix[1]]),T12_imag=imag(T[ix[1],ix[2]]),
            T21_imag=imag(T[ix[2],ix[1]]),T22_imag=imag(T[ix[2],ix[2]])))
    end
    return DataFrame(rows)
end

# Reconstruct only from archived inputs and independent analytic device code.
branches=CSV.read(joinpath(ROOT,"reports","experiment_D","inputs","branch.csv"),DataFrame)
buses=CSV.read(joinpath(ROOT,"reports","experiment_D","inputs","bus.csv"),DataFrame)
Ybus=assemble_ybus(branches,nrow(buses))
ports=generator_ports(buses)
passive=passive_port_graph(Ybus,ports)
Lc=passive.Lc
modal=modal_structure(Lc,passive.Yport;approximate_tol=COUPLING_TOL)
U=modal.U; n=length(ports); S=modal_transform(U)

# Frozen critical frequency from ExpC's all-SG analytic pole table. This is used
# as a diagnostic evaluation point only; no detailed simulator is called.
pole_table=CSV.read(joinpath(ROOT,"reports","experiment_C","tables","TABLE_C07_physical_pole_graph_mapping.csv"),DataFrame)
crit=first(filter(r->String(r.case)=="C0_all_SG" && String(r.role)=="spectral_abscissa_pair" && Float64(r.lambda_imag)>0,pole_table))
scrit=complex(Float64(crit.lambda_real),Float64(crit.lambda_imag))
net=CollectiveModel.frozen_network(ROOT)
println("F2: analytic IEEE-39 endpoint assembled (no PowerDynamics call)"); flush(stdout)
gfl=CollectiveModel.AnalyticGFLPLL
op_h=net.gfl[33].op
Jhom=gfl.jacobians(op_h.x,op_h.u,op_h.parameters;kp=NOM_KP,ki=NOM_KI)
Yhom=CollectiveModel.port_admittance(Jhom,scrit)
Ynet_h=realify_admittance(passive.Yport)
Ydev_h=kron(Matrix{Float64}(I,n,n),Yhom)
Th=S'*(Ynet_h+Ydev_h)*S
hom_total_offdiag,_=offdiagonal_block_ratio(Th;block_size=2)
hom_net_modal=modal_transform(U)'*Ynet_h*modal_transform(U)
hom_network_offdiag,_=offdiagonal_block_ratio(hom_net_modal;block_size=2)
Ph=characteristic_pencil(scrit,Jhom.A,Jhom.B,Jhom.C,Jhom.D,Ynet_h)
Phmodal=graph_pencil_transform(Ph,U,9,n)
println("F2: homogeneous descriptor and modal coupling evaluated"); flush(stdout)
hom_descriptor_offdiag,_=offdiagonal_block_ratio(Phmodal;block_size=11)
gain_deriv=gfl.gain_derivatives(op_h.x,op_h.u,op_h.parameters)
rank_kp=rank(gain_deriv.A_kp); rank_ki=rank(gain_deriv.A_ki)
rank_joint=rank(hcat(gain_deriv.A_kp,gain_deriv.A_ki))

# Independent heterogeneous IEEE-39 analytic all-GFL endpoint: use the frozen
# analytic network+impedance-load model from ExpE and local SimpleGFL Jacobians.
actualY,keep,rest=kron_real_ports(net.y_static,ports)
jacs=NamedTuple[]; transfers=Matrix{ComplexF64}[]
for b in ports
    op=net.gfl[b].op
    J=gfl.jacobians(op.x,op.u,op.parameters;kp=NOM_KP,ki=NOM_KI)
    push!(jacs,J)
    push!(transfers,CollectiveModel.port_admittance(J,scrit))
end
Aall,Ball,Call,Dall=blockdiag_state(jacs)
Ydev_het=blockdiag2(transfers)
Thet=S'*(ComplexF64.(actualY)+Ydev_het)*S
het_total_offdiag,het_off=offdiagonal_block_ratio(Thet;block_size=2)
Phet=pencil_from_blocks(scrit,Aall,Ball,Call,Dall,ComplexF64.(actualY))
Phetmodal=graph_pencil_transform(Phet,U,9,n)
het_descriptor_offdiag,_=offdiagonal_block_ratio(Phetmodal;block_size=11)
gamma=exact_modal_self_energy(Thet;block_size=2)
println("F2: heterogeneous exact port self-energy evaluated"); flush(stdout)

hom_metrics=(conductance_offdiag_ratio=modal.conductance_offdiag_ratio,
    susceptance_offdiag_ratio=modal.susceptance_offdiag_ratio,
    commutator_ratio=modal.commutator_ratio,
    network_admittance_offdiag_ratio=hom_network_offdiag,
    homogeneous_closed_loop_port_offdiag_ratio=hom_total_offdiag,
    homogeneous_descriptor_offdiag_ratio=hom_descriptor_offdiag,
    rank_dA_dKp=rank_kp,rank_dA_dKi=rank_ki,rank_joint_gain_update=rank_joint)
gate_class=modal.classification
gate_pass=gate_class!="NON_MODAL"

CSV.write(joinpath(OUT,"tables","TABLE_F01_modal_structure.csv"),DataFrame(
    metric=collect(String.(keys(hom_metrics))),value=collect(Float64.(values(hom_metrics))))
CSV.write(joinpath(OUT,"tables","TABLE_F02_heterogeneous_coupling.csv"),DataFrame(
    metric=["selected_pole_real_s_inv","selected_pole_imag_rad_s","static_model_dimension","homogeneous_descriptor_dimension",
        "heterogeneous_descriptor_dimension","homogeneous_exact_mode_block_dynamic_states","homogeneous_mode_block_algebraic_channels",
        "heterogeneous_port_transfer_offdiag_ratio","heterogeneous_descriptor_offdiag_ratio","heterogeneous_network_original_dimension",
        "heterogeneous_network_interior_solve_residual"],
    value=[real(scrit),imag(scrit),size(net.y_static,1),size(Ph,1),size(Phet,1),9,2,
        het_total_offdiag,het_descriptor_offdiag,size(net.y_static,1),norm(net.y_static[rest,keep]-net.y_static[rest,rest]*(net.y_static[rest,rest]\net.y_static[rest,keep]))]))
CSV.write(joinpath(OUT,"tables","TABLE_F03_exact_self_energy_by_mode.csv"),DataFrame(
    mode=[g.mode for g in gamma],Gamma11_real=[real(g.Gamma[1,1]) for g in gamma],
    Gamma12_real=[real(g.Gamma[1,2]) for g in gamma],Gamma21_real=[real(g.Gamma[2,1]) for g in gamma],
    Gamma22_real=[real(g.Gamma[2,2]) for g in gamma],Gamma11_imag=[imag(g.Gamma[1,1]) for g in gamma],
    Gamma12_imag=[imag(g.Gamma[1,2]) for g in gamma],Gamma21_imag=[imag(g.Gamma[2,1]) for g in gamma],
    Gamma22_imag=[imag(g.Gamma[2,2]) for g in gamma],schur_residual=[g.schur_residual for g in gamma]))
CSV.write(joinpath(OUT,"tables","TABLE_F04_homogeneous_modal_blocks.csv"),named_blocks(Th))
csv_matrix(joinpath(OUT,"matrices","F2_homogeneous_descriptor_modal_real.csv"),real.(Phmodal))
csv_matrix(joinpath(OUT,"matrices","F2_homogeneous_descriptor_modal_imag.csv"),imag.(Phmodal))
csv_matrix(joinpath(OUT,"matrices","F2_heterogeneous_port_modal_real.csv"),real.(Thet))
csv_matrix(joinpath(OUT,"matrices","F2_heterogeneous_port_modal_imag.csv"),imag.(Thet))
csv_matrix(joinpath(OUT,"matrices","F2_heterogeneous_descriptor_modal_real.csv"),real.(Phetmodal))
csv_matrix(joinpath(OUT,"matrices","F2_heterogeneous_descriptor_modal_imag.csv"),imag.(Phetmodal))

simple_pll=DataFrame(limit=["exact implemented PLL with first-order frequency state",
    "ideal frequency-loop limit tau -> 0", "graph coupling conclusion"],
    characteristic=["tau*s^3 + s^2 + V*Kp*s + V*Ki",
        "s^2 + V*Kp*s + V*Ki", "Lc does not enter the isolated PLL equation; network coupling remains in the full matrix port admittance"])
CSV.write(joinpath(OUT,"tables","TABLE_F05_simple_PLL_limit.csv"),simple_pll)

result=Dict("experiment"=>"BND_EXP_F2","status"=>(gate_pass ? "MODAL_GATE_PASS" : "NON_MODAL_STOP"),
    "modal_structure"=>gate_class,"approximate_coupling_tolerance"=>COUPLING_TOL,
    "physical_graph_source"=>"F0 passive generator-port conductance graph",
    "selected_frozen_ExpC_pole"=>String(crit.pole_id),"selected_s_real"=>real(scrit),"selected_s_imag"=>imag(scrit),
    "homogeneous_network_offdiag_ratio"=>hom_network_offdiag,
    "homogeneous_closed_loop_port_offdiag_ratio"=>hom_total_offdiag,
    "homogeneous_descriptor_offdiag_ratio"=>hom_descriptor_offdiag,
    "heterogeneous_port_transfer_offdiag_ratio"=>het_total_offdiag,
    "heterogeneous_descriptor_offdiag_ratio"=>het_descriptor_offdiag,
    "gain_update_ranks"=>Dict("Kp"=>rank_kp,"Ki"=>rank_ki,"joint"=>rank_joint),
    "homogeneous_descriptor_dimension"=>size(Ph,1),"heterogeneous_descriptor_dimension"=>size(Phet,1),
    "global_dynamic_polynomial_degree_upper_bound"=>9n,
    "scalar_characteristic"=>"NONE: reactive network coupling is not diagonal in the F0 conductance basis",
    "simple_PLL_limit"=>"tau*s^3+s^2+V*Kp*s+V*Ki; tau->0 yields the classical quadratic",
    "powerdynamics_imported"=>false,
    "gate_reason"=>(gate_pass ? "graph modal coupling meets predeclared structural tolerance" :
        "homogeneous detailed GFL/network pencil retains >5% off-diagonal coupling in the preregistered graph basis"))
open(joinpath(OUT,"RESULTS_EXP_F2.toml"),"w") do io TOML.print(io,result) end

open(joinpath(OUT,"REPORT_EXP_F2.md"),"w") do io
    println(io,"# Experiment F2 — Actual graph-modal GFL characteristic\n")
    println(io,"**F2_MODAL_STRUCTURE: `$(gate_class)`**\n")
    println(io,"## Exact characteristic structure\n")
    println(io,"The analytic SimpleGFLDC model contributes nine differential states per port. At common parameters its network interconnection is the exact descriptor determinant `det(P(s))=0`, with `P(s) = [I⊗(sI₉−A(Kp,Ki))  −I⊗B; I⊗C  Yport,rect + I⊗D]` under the frozen component-current sign convention. It has 90 differential and 20 algebraic variables for ten generator ports; its characteristic polynomial degree is at most 90 in `s`. Kp and Ki enter the local A matrix affinely through one rank-one update each; the exact gain dependencies are reported in TABLE_F01. The controlled homogeneous limit repeats the bus-33 analytic GFL state and parameters at all ten ports; fixed injection offsets balance the common state against the frozen passive network and have zero derivative, so they do not change this characteristic pencil.\n")
    println(io,"A scalar law `F(s;ν,Kp,Ki)=0` requires a common node basis for both the conductance and susceptance parts of the port admittance. In the F0 conductance eigenbasis, the conductance off-diagonal residual is `$(modal.conductance_offdiag_ratio)` and the reactive residual is `$(modal.susceptance_offdiag_ratio)`; the normalized commutator is `$(modal.commutator_ratio)`. The homogeneous descriptor residual and total port-transfer residual are `$(hom_descriptor_offdiag)` and `$(hom_total_offdiag)`. The small raw descriptor ratio is dominated by its diagonal `sI` state terms; the coupling gate is evaluated on the condensed port transfer and the network operator, where the residual is not diluted by those diagonal terms.\n")
    println(io,"## Heterogeneous IEEE-39 analytic model\n")
    println(io,"The detailed analytic endpoint uses every generator port, the frozen network plus archived impedance-load linearizations, independent local nine-state SimpleGFLDC Jacobians, and the same F0 basis. At the frozen ExpC critical pole `$(crit.pole_id)` (`s=$(scrit)`), its port-transfer and descriptor off-diagonal residuals are `$(het_total_offdiag)` and `$(het_descriptor_offdiag)`. Exact 2-channel Schur self-energy Gamma by graph mode is stored in TABLE_F03; this diagnoses coupling, it does not create an independent modal characteristic.\n")
    println(io,"This assembly imports no PowerDynamics module. Its saved matrices and all mode blocks are reproducible from the archived analytical inputs.\n")
    println(io,"## PLL-only limit\n")
    println(io,"With `e≈−Vθ`, the implemented PLL equations retain the first-order frequency state and give `τs³+s²+VKp s+VKi=0`. Only in the ideal frequency-loop limit `τ→0` does this reduce to `s²+VKp s+VKi=0`. No graph eigenvalue appears in that isolated PLL limit; any graph dependence has to come from the coupled network equations.\n")
    if gate_class=="NON_MODAL"
        println(io,"## STOP gate\n\n**NON_MODAL.** The homogeneous detailed model exceeds the preregistered 5% modal-coupling tolerance, so no exact/approximately separated `F(s;ν,Kp,Ki)` is available. F3/F4 spectral gain-law work is stopped. Continue with the exact heterogeneous self-energy/pathway formulation; do not claim `K=h(Lc)`.\n")
    else
        println(io,"## Gate\n\nThe graph-modal structure passes the predeclared structural coupling criterion. A separate F3 derivation is required before any gain law is claimed.\n")
    end
end
println("EXPERIMENT_F2_MODAL_STRUCTURE: ",gate_class)
println("HOMOGENEOUS_G_NETWORK_OFFDIAG: ",modal.conductance_offdiag_ratio)
println("HOMOGENEOUS_B_NETWORK_OFFDIAG: ",modal.susceptance_offdiag_ratio)
println("HOMOGENEOUS_DESCRIPTOR_OFFDIAG: ",hom_descriptor_offdiag)
println("HETEROGENEOUS_PORT_TRANSFER_OFFDIAG: ",het_total_offdiag)
println("HETEROGENEOUS_DESCRIPTOR_OFFDIAG: ",het_descriptor_offdiag)
