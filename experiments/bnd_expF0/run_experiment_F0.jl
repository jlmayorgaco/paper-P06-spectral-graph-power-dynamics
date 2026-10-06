using LinearAlgebra
using CSV
using DataFrames
using SHA
using TOML
BLAS.set_num_threads(1)

include(joinpath(@__DIR__, "..", "..", "src", "bnd_graphpll", "PhysicalGraph.jl"))
using .PhysicalGraph

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const INPUT = joinpath(ROOT, "reports", "experiment_D", "inputs")
const OUT = joinpath(ROOT, "reports", "experiment_F0")
mkpath(joinpath(OUT, "tables")); mkpath(joinpath(OUT, "matrices"))

sha256_file(path) = bytes2hex(sha256(read(path)))
function matrix_csv(path)
    df = CSV.read(path, DataFrame)
    Matrix{Float64}(select(df, Not(first(names(df)))))
end
function csv_matrix(path, A)
    n,m=size(A)
    df=DataFrame(row_index=collect(1:n))
    for j in 1:m
        df[!,Symbol("x$j")]=A[:,j]
    end
    CSV.write(path,df)
end

branch_path=joinpath(INPUT,"branch.csv")
bus_path=joinpath(INPUT,"bus.csv")
branches=CSV.read(branch_path,DataFrame)
buses=CSV.read(bus_path,DataFrame)
nbus=nrow(buses)
ports=generator_ports(buses)
port_buses=Int.(buses.bus[ports])
Ybus=assemble_ybus(branches,nbus)
passive=passive_port_graph(Ybus,ports)
Lc=passive.Lc
graph=eigen(Hermitian(Lc))
perm=sortperm(graph.values)
λc=Float64.(graph.values[perm]); U=Matrix{Float64}(graph.vectors[:,perm])

# Diagnostic candidates come from the same frozen passive network; neither is
# repaired or selected using controller performance.
Lreactive=-imag.(passive.Yport)
reactive=spectral_summary(Lreactive)
Gdirect=real.(Ybus)
Gdirect_summary=spectral_summary(Gdirect)
Lsync=nothing
sync_status="not available in frozen ExpC matrices"
sync_path=joinpath(ROOT,"reports","experiment_C","matrices","C0_LG_primary.csv")
if isfile(sync_path)
    Lsync=matrix_csv(sync_path)
    sync_status="classified from frozen ExpC C0 primary synchronizing backbone"
end

kron_audit=kron_identity_audit(Ybus,passive.reduction;trials=7)
sym_transpose=norm(Ybus-transpose(Ybus),Inf)/max(norm(Ybus,Inf),eps(Float64))
sym_port=norm(passive.Yport-transpose(passive.Yport),Inf)/max(norm(passive.Yport,Inf),eps(Float64))
full_power_min=Inf
port_power_min=Inf
for k in 1:13
    v=ComplexF64[complex(sin(0.17*i+0.39*k),cos(0.31*i-0.27*k)) for i in 1:nbus]
    vp=v[ports]
    global full_power_min=min(full_power_min,real(dot(v,Ybus*v)))
    global port_power_min=min(port_power_min,real(dot(vp,Lc*vp)))
end

Lc_summary=spectral_summary(Lc)
sync_summary=Lsync===nothing ? nothing : spectral_summary(Lsync)
gate=passive.hermitian_residual<=1e-12 && passive.imaginary_leakage<=1e-12 &&
    Lc_summary.psd && kron_audit.relative_residual<=1e-11 &&
    full_power_min>=-1e-10 && port_power_min>=-1e-10

csv_matrix(joinpath(OUT,"matrices","Ybus_real.csv"),real.(Ybus))
csv_matrix(joinpath(OUT,"matrices","Ybus_imag.csv"),imag.(Ybus))
csv_matrix(joinpath(OUT,"matrices","Yport_real.csv"),real.(passive.Yport))
csv_matrix(joinpath(OUT,"matrices","Yport_imag.csv"),imag.(passive.Yport))
csv_matrix(joinpath(OUT,"matrices","Lc_conductance.csv"),Lc)
csv_matrix(joinpath(OUT,"matrices","Lc_eigenvectors.csv"),U)
CSV.write(joinpath(OUT,"tables","TABLE_F02_graph_spectrum.csv"),DataFrame(
    mode=collect(1:length(λc)),eigenvalue_pu_conductance=λc,
    zero_mode=abs.(λc).<=Lc_summary.tolerance))

candidate_rows=NamedTuple[]
push!(candidate_rows,(candidate="Kron passive-port conductance Re(Yport)",
    source="Hermitian part of complex branch-only Kron Ybus at buses $(join(port_buses, ','))",
    lambda_min=Lc_summary.lambda_min,lambda_max=Lc_summary.lambda_max,
    zero_mode_count=Lc_summary.zero_mode_count,condition_number=Lc_summary.condition_number,
    psd=Lc_summary.psd,spd=Lc_summary.spd,symmetry_residual=Lc_summary.symmetry_residual,
    selected=true,decision="Primary: nonnegative passive dissipation operator"))
push!(candidate_rows,(candidate="Kron susceptance -Im(Yport)",source="Reactive part of the same complex port admittance",
    lambda_min=reactive.lambda_min,lambda_max=reactive.lambda_max,
    zero_mode_count=reactive.zero_mode_count,condition_number=reactive.condition_number,
    psd=reactive.psd,spd=reactive.spd,symmetry_residual=reactive.symmetry_residual,
    selected=false,decision=reactive.psd ? "PSD candidate; not primary by preregistered criterion" : "Indefinite; no eigenvalue clipping"))
if sync_summary!==nothing
    push!(candidate_rows,(candidate="ExpC C0 primary synchronizing backbone L_G",
        source=sync_status,lambda_min=sync_summary.lambda_min,lambda_max=sync_summary.lambda_max,
        zero_mode_count=sync_summary.zero_mode_count,condition_number=sync_summary.condition_number,
        psd=sync_summary.psd,spd=sync_summary.spd,symmetry_residual=sync_summary.symmetry_residual,
        selected=false,decision=sync_summary.psd ? "PSD in this case; no controller-based selection" : "Indefinite; not passed to a real square root"))
end
sync33_path=joinpath(ROOT,"reports","experiment_C","matrices","C33_LG_primary.csv")
sync33_summary=nothing
if isfile(sync33_path)
    sync33_summary=spectral_summary(matrix_csv(sync33_path))
    push!(candidate_rows,(candidate="ExpC C33 primary synchronizing backbone L_G",
        source="classified from frozen ExpC bus-33 GFL primary synchronizing backbone",
        lambda_min=sync33_summary.lambda_min,lambda_max=sync33_summary.lambda_max,
        zero_mode_count=sync33_summary.zero_mode_count,condition_number=sync33_summary.condition_number,
        psd=sync33_summary.psd,spd=sync33_summary.spd,symmetry_residual=sync33_summary.symmetry_residual,
        selected=false,decision=sync33_summary.psd ? "PSD in this case; not selected for physical graph" : "Indefinite; no real square root"))
end
CSV.write(joinpath(OUT,"tables","TABLE_F01_graph_provenance.csv"),DataFrame(
    field=["input_bus_csv","input_branch_csv","bus_order","generator_port_order","network_operator","tap_convention","interior_bus_treatment","grounding_or_reference","units","Ybus_transpose_symmetry_residual","Yport_transpose_symmetry_residual"],
    value=["reports/experiment_D/inputs/bus.csv","reports/experiment_D/inputs/branch.csv",
        join(1:nbus,","),join(port_buses,","),"branch-only conventional passive Ybus; loads and device shunts excluded",
        "PiLine_fault source ratio r_src, destination ratio 1; verified against frozen analytic assembly",
        "Kron eliminate all buses not in generator-port list using complex Ybus Schur complement",
        "No artificial ground; transformer tap ratios are retained. Gauge/zero modes reported as found.",
        "per-unit admittance / conductance on archived system base",string(sym_transpose),string(sym_port)]))
CSV.write(joinpath(OUT,"tables","TABLE_F03_kron_reduction.csv"),DataFrame(
    metric=["retained_ports","eliminated_buses","interior_block_condition_number","direct_vs_kron_relative_current_residual","direct_vs_kron_max_abs_current_residual","Ybus_transpose_symmetry_residual","Yport_transpose_symmetry_residual","full_network_min_real_power_test","port_operator_min_quadratic_test"],
    value=[length(ports),nbus-length(ports),passive.reduction.interior_condition,
        kron_audit.relative_residual,kron_audit.max_abs_residual,sym_transpose,sym_port,full_power_min,port_power_min]))
CSV.write(joinpath(OUT,"tables","TABLE_F04_candidate_classification.csv"),DataFrame(candidate_rows))

# ExpC supplied frozen critical poles and q-coordinate states. This is a
# post-preregistration diagnostic only; it cannot alter the selected graph.
mode_diag_status="separate read-only diagnostic: run experiments/bnd_expF0/run_mode_overlap_F0.py"
c0mat=joinpath(ROOT,"reports","experiment_C","matrices")
c0table=joinpath(ROOT,"reports","experiment_C","tables","TABLE_C07_physical_pole_graph_mapping.csv")
if get(ENV,"BND_F0_RUN_MODE_DIAG","0")=="1" && all(isfile.([joinpath(c0mat,"C0_M.csv"),joinpath(c0mat,"C0_Ared.csv"),
                joinpath(c0mat,"C0_states.csv"),c0table]))
    println("F0_MODE_DIAGNOSTIC: loading frozen ExpC q coordinates")
    M=matrix_csv(joinpath(c0mat,"C0_M.csv"))
    Ared=matrix_csv(joinpath(c0mat,"C0_Ared.csv"))
    states=CSV.read(joinpath(c0mat,"C0_states.csv"),DataFrame)
    state_bus=Dict{Int,Int}()
    for row in eachrow(states)
        name=String(row.state_name)
        if occursin("machine",name) && occursin("δ",name)
            m=match(r"VIndex\((\d+),",name)
            m===nothing || (state_bus[parse(Int,m.captures[1])]=Int(row.state_index))
        end
    end
    qbus=sort(collect(keys(state_bus)))
    if length(qbus)==length(port_buses)==size(M,1) && qbus==port_buses
        qidx=[state_bus[b] for b in qbus]
        ME=eigen(Symmetric((M+M')/2))
        if minimum(ME.values)>0
            Mh=ME.vectors*Diagonal(sqrt.(ME.values))*ME.vectors'
            Mih=ME.vectors*Diagonal(1 ./ sqrt.(ME.values))*ME.vectors'
            GE=eigen(Hermitian((Mih*Lc*Mih + (Mih*Lc*Mih)')/2))
            gp=sortperm(GE.values); Ug=Matrix(GE.vectors[:,gp])
            pole_rows=CSV.read(c0table,DataFrame)
            filter!(r -> String(r.case)=="C0_all_SG" &&
                String(r.role) in ("spectral_abscissa_pair","lowest_damping_band_pair","high_retained_participation") &&
                Float64(r.lambda_imag)>0, pole_rows)
            unique!(pole_rows,:pole_id)
            Efull=eigen(Ared)
            println("F0_MODE_DIAGNOSTIC: matching frozen critical poles")
            mode_rows=NamedTuple[]; qcols=Vector{Vector{ComplexF64}}()
            for row in eachrow(pole_rows)
                target=complex(Float64(row.lambda_real),Float64(row.lambda_imag))
                j=argmin(abs.(Efull.values .- target))
                abs(Efull.values[j]-target)<=1e-6 || continue
                qmw=ComplexF64.(Mh*Efull.vectors[qidx,j])
                nq=norm(qmw); nq>eps(Float64) || continue
                qmw./=nq; push!(qcols,qmw)
                energy=abs2.(Ug'*qmw)
                order=sortperm(energy;rev=true)
                push!(mode_rows,(pole_id=String(row.pole_id),lambda_real=real(target),
                    lambda_imag=imag(target),frequency_hz=Float64(row.frequency_hz),
                    q_norm=Float64(row.q_norm),top_graph_mode=order[1],
                    top_mode_energy=energy[order[1]],top3_graph_mode_energy=sum(energy[order[1:min(3,end)]])))
            end
            if !isempty(mode_rows)
                CSV.write(joinpath(OUT,"tables","TABLE_F05_graph_vs_real_modes.csv"),DataFrame(mode_rows))
                Q=hcat(qcols...); SE=svd(Q); r=count(x->x>1e-10*maximum(SE.S),SE.S)
                Qorth=Matrix(SE.U[:,1:r])
                bands=[("low_modes_1_3",1:min(3,length(λc))),
                    ("middle_modes_4_7",4:min(7,length(λc))),
                    ("high_modes_8_10",max(1,length(λc)-2):length(λc))]
                angle_rows=NamedTuple[]
                for (label,ix) in bands
                    ix=collect(ix); isempty(ix) && continue
                    σ=clamp.(svdvals(Ug[:,ix]'*Qorth),0.0,1.0)
                    for (j,v) in enumerate(σ)
                        push!(angle_rows,(graph_band=label,principal_angle_index=j,
                            cosine_overlap=v,angle_degrees=rad2deg(acos(v)),
                            critical_q_subspace_rank=r))
                    end
                end
                CSV.write(joinpath(OUT,"tables","TABLE_F06_principal_angles.csv"),DataFrame(angle_rows))
                mode_diag_status="computed from frozen ExpC C0 critical poles with M-weighted q projection"
            end
        else
            mode_diag_status="not computed: ExpC q mass matrix is not positive definite"
        end
    else
        mode_diag_status="not computed: ExpC angle-state ordering does not match generator ports"
    end
end

result=(experiment="BND_EXP_F0",status=gate ? "PASS" : "FAIL",
    primary_graph="Lc = Hermitian(Yport) = Re(Yport), complex passive Kron reduction then generator-port selection",
    input_bus_sha256=sha256_file(bus_path),input_branch_sha256=sha256_file(branch_path),
    nbus=nbus,branch_count=nrow(branches),generator_ports=port_buses,
    Ybus_transpose_symmetry_residual=sym_transpose,Yport_transpose_symmetry_residual=sym_port,
    interior_block_condition_number=passive.reduction.interior_condition,
    Lc_lambda_min=Lc_summary.lambda_min,Lc_lambda_max=Lc_summary.lambda_max,
    Lc_condition_number=Lc_summary.condition_number,Lc_zero_mode_count=Lc_summary.zero_mode_count,
    Lc_psd=Lc_summary.psd,Lc_spd=Lc_summary.spd,Lc_hermitian_residual=passive.hermitian_residual,
    Lc_imaginary_leakage=passive.imaginary_leakage,
    Kron_relative_current_residual=kron_audit.relative_residual,
    Kron_max_abs_current_residual=kron_audit.max_abs_residual,
    full_network_min_real_power_test=full_power_min,port_operator_min_quadratic_test=port_power_min,
    reactive_candidate_lambda_min=reactive.lambda_min,reactive_candidate_lambda_max=reactive.lambda_max,
    reactive_candidate_psd=reactive.psd,
    synchronizing_backbone_available=Lsync!==nothing,
    synchronizing_backbone_lambda_min=sync_summary===nothing ? "" : sync_summary.lambda_min,
    synchronizing_backbone_lambda_max=sync_summary===nothing ? "" : sync_summary.lambda_max,
    synchronizing_backbone_psd=sync_summary===nothing ? "" : sync_summary.psd,
    synchronizing_backbone_C33_lambda_min=sync33_summary===nothing ? "" : sync33_summary.lambda_min,
    synchronizing_backbone_C33_lambda_max=sync33_summary===nothing ? "" : sync33_summary.lambda_max,
    synchronizing_backbone_C33_psd=sync33_summary===nothing ? "" : sync33_summary.psd,
    mode_comparison_status=mode_diag_status,
    gate_reason=gate ? "PSD passive network-derived graph and port-current Kron identity verified" : "one or more passivity/Kron/PSD checks failed")
open(joinpath(OUT,"RESULTS_EXP_F0.toml"),"w") do io
    TOML.print(io,Dict(string(k)=>v for (k,v) in pairs(result)))
end

open(joinpath(OUT,"REPORT_EXP_F0.md"),"w") do io
    println(io,"# Experiment F0 — Physical IEEE-39 generator-port graph\n")
    println(io,"**F0_STATUS: ",gate ? "PASS" : "FAIL","**\n")
    println(io,"## Preregistered primary operator\n")
    println(io,"Before any controller analysis, the primary graph was fixed as the real Hermitian part of the passive complex branch admittance after Kron reduction onto generator buses: `Lc = Re(Yport)`. The branch-only network excludes generator, load, and converter device admittances. Generator-port order is `$(join(port_buses, ", "))`; all other buses are eliminated. Values are per-unit conductance on the archived system base. The source-side transformer ratio is retained exactly. No artificial ground, clipping, or absolute-eigenvalue operation is used.\n")
    println(io,"## Physical and numerical checks\n")
    println(io,"- `Ybus` and its Kron port matrix are complex-symmetric to relative residuals `$(sym_transpose)` and `$(sym_port)`.\n")
    println(io,"- `Lc` eigenvalue range: `$(Lc_summary.lambda_min)` to `$(Lc_summary.lambda_max)`; zero modes: `$(Lc_summary.zero_mode_count)`; condition number: `$(Lc_summary.condition_number)`.\n")
    println(io,"- The complex Kron map reproduces port currents for seven deterministic port-voltage trials with relative residual `$(kron_audit.relative_residual)` (maximum absolute residual `$(kron_audit.max_abs_residual)`).\n")
    println(io,"- Minimum tested real power over the full passive network: `$(full_power_min)`; minimum port quadratic form: `$(port_power_min)`.\n")
    println(io,"- The candidate table classifies reactive susceptance and ExpC C0/C33 synchronizing backbones separately. The C33 backbone is tested directly; no indefinite operator is modified to make it PSD.\n")
    println(io,"- ExpC critical modes are compared with the graph basis only as a diagnostic: `$(mode_diag_status)`. The M-weighted mode-energy and principal-angle tables are retained separately.\n")
    println(io,"## Gate\n")
    println(io,gate ? "**F0_PASS.** A network-derived PSD operator exists and the port-current relation is validated. This establishes a mathematically admissible graph spectrum; it does not establish that the detailed GFL model is modal in this basis." : "**F0_FAIL.** The primary physical/passivity/Kron gate did not pass; spectral gain-law work is blocked.")
    println(io,"\nInputs are the frozen CSV copies under `reports/experiment_D/inputs`; SHA-256 values are in `RESULTS_EXP_F0.toml`. The experiment code imports no PowerDynamics module.\n")
end

println("EXPERIMENT_F0_STATUS: ",gate ? "PASS" : "FAIL")
println("PRIMARY_GRAPH: Re(Yport), ports=",join(port_buses,","))
println("LC_SPECTRUM: [",Lc_summary.lambda_min,", ",Lc_summary.lambda_max,"] zeros=",Lc_summary.zero_mode_count,
    " cond=",Lc_summary.condition_number)
println("KRON_RELATIVE_CURRENT_RESIDUAL: ",kron_audit.relative_residual)
println("REACTIVE_CANDIDATE_PSD: ",reactive.psd," spectrum=[",reactive.lambda_min,", ",reactive.lambda_max,"]")
if sync_summary!==nothing
    println("EXPC_SYNCHRONIZING_BACKBONE_PSD: ",sync_summary.psd," spectrum=[",sync_summary.lambda_min,", ",sync_summary.lambda_max,"]")
end
if sync33_summary!==nothing
    println("EXPC33_SYNCHRONIZING_BACKBONE_PSD: ",sync33_summary.psd," spectrum=[",sync33_summary.lambda_min,", ",sync33_summary.lambda_max,"]")
end
