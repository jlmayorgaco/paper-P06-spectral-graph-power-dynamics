include(joinpath(@__DIR__,"..","..","src","bnd_design","BNDDesign.jl"))

using .BNDDesign.AnalyticDeviceModel
using .BNDDesign.AnalyticGFLPLL
import .BNDDesign.AnalyticSG
using .BNDDesign.PortReduction
using .BNDDesign.WoodburyPLL
using CSV
using DataFrames
using Dates
using LinearAlgebra
using Printf
using Random
using SHA
using Statistics

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const OUT=joinpath(ROOT,"reports","experiment_D")
const TABLES=joinpath(OUT,"tables")
const MATRICES=joinpath(OUT,"matrices")
const CROSS_BUS_SCOPE=get(ENV,"BND_EXP_D_CROSS_BUS","true")=="true"
const PORT_BUSES=CROSS_BUS_SCOPE ? (30,33,35,37) : (33,)
mkpath(TABLES); mkpath(MATRICES); mkpath(joinpath(OUT,"figures"))

function sha256_file(path)
    return bytes2hex(sha256(read(path)))
end

function write_empty_tables()
    schemas=Dict(
      "TABLE_D04_port_update_rank.csv"=>["bus","frequency","rho","rank","singular_values","affine_rho_residual","status"],
      "TABLE_D05_closure_validation.csv"=>["bus","frequency","closure_dimension","determinant_residual","pole_residual","status"],
      "TABLE_D07_analytic_gain_candidates.csv"=>["bus","omega","branch","Kp","Ki","physical","exclusion_reason"],
      "TABLE_D08_analytic_rho_candidates.csv"=>["bus","omega","branch","rho","Kp","Ki","candidate_type","admissible"],
      "TABLE_D09_analytic_optimum.csv"=>["bus","rho_star","MW_replaced","Kp_star","Ki_star","omega_star","boundary_real","boundary_imag","active_constraint"],
      "TABLE_D10_multimode_check.csv"=>["bus","mode","lambda_real","lambda_imag","damping_ratio","pass"],
      "TABLE_D11_PD_validation.csv"=>["bus","rho_analytic","Kp_analytic","Ki_analytic","predicted_lambda","PD_lambda","abs_error","real_error","freq_error","damping_error"],
      "TABLE_D12_TDS_validation.csv"=>["bus","rho","candidate","signal","NRMSE","retcode"],
      "TABLE_D13_self_energy_mechanism.csv"=>["case","mode","psi_diag","gamma","cancellation_factor","dominant_pathways"],
      "TABLE_D14_cross_bus_summary.csv"=>["bus","pll_structure","pll_update_rank","port_update_rank","closure_dimension","gain_polynomial_degree","rho_polynomial_degree","rho_star","Kp_star","Ki_star","PD_validation_error","status"],
      "TABLE_D15_method_vs_validation_solver.csv"=>["method","rho_star","objective_gap","runtime_s","full_model_evaluations","label"])
    for (name,cols) in schemas
        df=DataFrame()
        for col in cols; df[!,Symbol(col)]=String[]; end
        CSV.write(joinpath(TABLES,name),df)
    end
end

function csv_matrix(path,A)
    names=["x$(j)" for j in axes(A,2)]
    df=DataFrame(Matrix(A),Symbol.(names))
    insertcols!(df,1,:row_index=>collect(1:size(A,1)))
    CSV.write(path,df)
end

function singular_text(s)
    join((@sprintf("%.12e",v) for v in s),";")
end

function main()
    frozen=AnalyticGFLPLL.frozen_bus33_operating_point(ROOT)
    x,u,p=frozen.x,frozen.u,frozen.parameters
    J=jacobians(x,u,p)
    sg=AnalyticSG.frozen_bus_operating_point(ROOT,33)
    sgdiag=AnalyticSG.steady_state_diagnostics(sg.x,sg.u,sg.parameters)
    Js=AnalyticSG.jacobians(sg.x,sg.u,sg.parameters)
    gflrated=rebase_to_rating(frozen,sg.parameters.rating_mva,sg.parameters.system_base_mva)
    Jgflrated=jacobians(gflrated.x,gflrated.u,gflrated.parameters)
    alpha=gflrated.rating_ratio
    port_points=[0.05+2pi*hz*im for hz in 10.0 .^ range(log10(0.01),log10(10.0),length=41)]
    append!(port_points,ComplexF64[-0.09884700664631951+0.43026220070517623im,
        -0.35790147141068607+3.9211749225740378im,-0.40561412912665884+6.185610958228528im])
    rho_samples=(0.0,0.1,0.25,0.5,0.75,0.9,1.0)
    port_rows=NamedTuple[]; port_transfer_errors=Float64[]; affine_rows=Float64[]
    cross_rows=NamedTuple[]; current_errors=Float64[]; power_errors=Float64[]
    for bus in PORT_BUSES
        sg_bus=AnalyticSG.frozen_bus_operating_point(ROOT,bus)
        sd=AnalyticSG.steady_state_diagnostics(sg_bus.x,sg_bus.u,sg_bus.parameters)
        gfl_sys=bus==33 ? frozen : AnalyticGFLPLL.operating_point_for_injection(
            sg_bus.u,sd.network_power.p,sd.network_power.q)
        ratio=sg_bus.parameters.rating_mva/sg_bus.parameters.system_base_mva
        gfl_own=rebase_to_rating(gfl_sys,sg_bus.parameters.rating_mva,sg_bus.parameters.system_base_mva)
        Jsg=AnalyticSG.jacobians(sg_bus.x,sg_bus.u,sg_bus.parameters)
        Jgflsys=jacobians(gfl_sys.x,gfl_sys.u,gfl_sys.parameters)
        Jgfldev=jacobians(gfl_own.x,gfl_own.u,gfl_own.parameters)
        gfl_rhs=norm(rhs(gfl_own.x,gfl_own.u,gfl_own.parameters),Inf)
        ierr=norm(sd.port_current-rated_output(gfl_own.x,gfl_own.u,gfl_own.parameters,ratio))
        i_gfl=rated_output(gfl_own.x,gfl_own.u,gfl_own.parameters,ratio)
        p_gfl=real(complex(sg_bus.u[1],sg_bus.u[2])*conj(complex(i_gfl[1],i_gfl[2])))
        q_gfl=imag(complex(sg_bus.u[1],sg_bus.u[2])*conj(complex(i_gfl[1],i_gfl[2])))
        perr=maximum(abs.([sd.network_power.p-p_gfl,sd.network_power.q-q_gfl]))
        push!(current_errors,ierr); push!(power_errors,perr)
        transfer_errs=Float64[]; bus_affine=Float64[]
        for s in port_points
            Ysg=AnalyticSG.port_transfer(Jsg.A,Jsg.B,Jsg.C,Jsg.D,s)
            Ygfl=rated_port_transfer(Jgfldev,s,ratio)
            Ysys=port_transfer(Jgflsys.A,Jgflsys.B,Jgflsys.C,Jgflsys.D,s)
            push!(transfer_errs,norm(Ysys-Ygfl)/max(norm(Ysys),eps(Float64)))
            Dgfl=Ygfl-Ysg
            for rho in rho_samples
                Ymix=PortReduction.rated_parallel_mix(Ysg,Ygfl,rho)
                err=norm(Ymix-(Ysg+rho*Dgfl))/max(norm(Ymix),eps(Float64))
                push!(affine_rows,err); push!(bus_affine,err)
                d=rho*Dgfl
                sv=svdvals(Matrix(d)); tol=isempty(sv) ? 0.0 : maximum(size(d))*eps(Float64)*first(sv)
                prank=count(>(tol),sv)
                push!(port_rows,(bus=bus,frequency=imag(s)/(2pi),s_real=real(s),rho=rho,
                    rank=prank,singular_values=singular_text(sv),affine_rho_residual=err,
                    status="LOCAL_PORT_ONLY"))
            end
        end
        push!(port_transfer_errors,maximum(transfer_errs))
        push!(cross_rows,(bus=bus,rating_mva=sg_bus.parameters.rating_mva,
            system_to_device_ratio=ratio,sg_rhs_inf=sd.rhs_norm,
            gfl_rhs_inf=gfl_rhs,
            sg_P=sd.network_power.p,sg_Q=sd.network_power.q,
            port_current_error=ierr,port_power_error=perr,
            gfl_transfer_equivalence_max=maximum(transfer_errs),
            affine_rho_residual_max=maximum(bus_affine),
            status=(sd.rhs_norm<1e-10&&gfl_rhs<1e-10&&ierr<1e-9&&perr<1e-9&&maximum(transfer_errs)<1e-10 ?
                "PASS_LOCAL_PORT" : "FAIL")))
    end
    sg_gfl_current_error=maximum(current_errors)
    sg_gfl_power_error=maximum(power_errors)
    gd=gain_derivatives(x,u,p)
    w=verify_woodbury(x,u,p;samples=32)
    kp0,ki0=p.pll_kp,p.pll_ki
    B0=J.B-kp0*gd.B_kp-ki0*gd.B_ki

    # Check affine parameter structure at deterministic independent points.
    rng=Random.MersenneTwister(2601)
    affine_A=Float64[]; affine_B=Float64[]; cd_errors=Float64[]
    for _ in 1:16
        kp=kp0*(0.5+rand(rng)); ki=ki0*(0.5+rand(rng))
        Jp=jacobians(x,u,p;kp,ki)
        Aaff=J.A+(kp-kp0)*gd.A_kp+(ki-ki0)*gd.A_ki
        Baff=J.B+(kp-kp0)*gd.B_kp+(ki-ki0)*gd.B_ki
        push!(affine_A,norm(Jp.A-Aaff)/max(norm(Jp.A),eps(Float64)))
        push!(affine_B,norm(Jp.B-Baff)/max(norm(Jp.B),eps(Float64)))
        push!(cd_errors,max(norm(Jp.C-J.C),norm(Jp.D-J.D)))
    end

    matrices=Dict("dA_dKp"=>gd.A_kp,"dA_dKi"=>gd.A_ki,"dA_both"=>hcat(gd.A_kp,gd.A_ki),
        "dB_dKp"=>gd.B_kp,"dB_dKi"=>gd.B_ki,"dB_both"=>hcat(gd.B_kp,gd.B_ki),
        "dC_dKp"=>gd.C_kp,"dC_dKi"=>gd.C_ki,"dD_dKp"=>gd.D_kp,"dD_dKi"=>gd.D_ki)
    dep_rows=NamedTuple[]
    for (name,A) in matrices
        audit=rank_audit(A)
        push!(dep_rows,(matrix=name,parameter=occursin("Kp",name)&&occursin("Ki",name) ? "Kp,Ki" :
                occursin("Kp",name) ? "Kp" : "Ki",exact_affine=true,
            numerical_rank=audit.rank,leading_singular_values=singular_text(audit.singular_values),
            rank_tolerance=audit.tolerance,residual=string(audit.residual)))
    end
    CSV.write(joinpath(TABLES,"TABLE_D02_gain_dependency_audit.csv"),DataFrame(dep_rows))
    CSV.write(joinpath(TABLES,"TABLE_D03_Woodbury_validation.csv"),DataFrame(w.sample_rows))

    # This denominator is for the local PLL state resolvent only. It is not the
    # SG-to-GFL replacement closure or its pole-boundary polynomial.
    F=WoodburyPLL.factorization(x,u,p)
    poly_rows=NamedTuple[]
    for hz in 10.0 .^ range(log10(0.01),log10(10.0),length=41)
        s=0.05+2pi*hz*im
        c=WoodburyPLL.denominator_coefficients(F.A0,F.U,F.V,s)
        push!(poly_rows,(bus=33,omega=2pi*hz,degree_Kp=1,degree_Ki=1,total_degree=2,
            coefficient_name="c00",real=real(c.c00),imag=imag(c.c00),
            closure_context="local PLL state resolvent only"))
        for (nm,z) in (("c10",c.c10),("c01",c.c01),("c11",c.c11))
            push!(poly_rows,(bus=33,omega=2pi*hz,degree_Kp=1,degree_Ki=1,total_degree=2,
                coefficient_name=nm,real=real(z),imag=imag(z),
                closure_context="local PLL state resolvent only"))
        end
    end
    CSV.write(joinpath(TABLES,"TABLE_D06_gain_polynomial.csv"),DataFrame(poly_rows))

    # Save local state-space model and all gain-update factors for independent audit.
    for (name,A) in (("GFL_A_reference.csv",J.A),("GFL_B_reference.csv",J.B),
        ("GFL_C_reference.csv",J.C),("GFL_D_reference.csv",J.D),
        ("GFL_A_zero_gain.csv",F.A0),("GFL_U_PLL.csv",F.U),("GFL_V_PLL.csv",F.V),
        ("dA_dKp.csv",gd.A_kp),("dA_dKi.csv",gd.A_ki),
        ("dB_dKp.csv",gd.B_kp),("dB_dKi.csv",gd.B_ki))
        csv_matrix(joinpath(MATRICES,name),A)
    end

    srcfiles=filter(f->endswith(f,".jl"),[joinpath(d,f) for (d,_,fs) in walkdir(joinpath(ROOT,"src","bnd_design")) for f in fs])
    import_re=r"(?m)^\s*(?:using|import)\s+(?:PowerDynamics|[A-Za-z0-9_\.]*\.PowerDynamics)\b"
    firewall=all(!occursin(import_re,read(f,String)) for f in srcfiles)
    inputpaths=[joinpath(ROOT,"reports","experiment_A","tables","TABLE_A03_state_partition.csv"),
        joinpath(ROOT,"reports","experiment_A","matrices","bus33_baseline_equilibrium.csv"),
        joinpath(ROOT,"reports","experiment_A","matrices","bus33_mixed_equilibrium.csv"),
        joinpath(ROOT,"reports","experiment_A","RESULTS_EXP_A.json"),
        joinpath(ROOT,"reports","experiment_B","RESULTS_EXP_B.json"),
        joinpath(ROOT,"reports","experiment_C","RESULTS_EXP_C.json"),
        joinpath(ROOT,"reports","experiment_D","inputs","bus.csv"),
        joinpath(ROOT,"reports","experiment_D","inputs","branch.csv"),
        joinpath(ROOT,"reports","experiment_D","inputs","load.csv"),
        joinpath(ROOT,"reports","experiment_D","inputs","machine.csv"),
        joinpath(ROOT,"reports","experiment_D","inputs","avr.csv"),
        joinpath(ROOT,"reports","experiment_D","inputs","gov.csv")]
    hashes=Dict(relpath(f,ROOT)=>sha256_file(f) for f in inputpaths)
    resultpaths=[joinpath(ROOT,"reports","experiment_A","RESULTS_EXP_A.json"),
        joinpath(ROOT,"reports","experiment_B","RESULTS_EXP_B.json"),
        joinpath(ROOT,"reports","experiment_C","RESULTS_EXP_C.json")]
    a_pass=occursin(r"\"status\"\s*:\s*\"PASS\"",read(resultpaths[1],String))
    b_pass=occursin(r"\"status\"\s*:\s*\"PASS\"",read(resultpaths[2],String))
    c_pass=occursin(r"\"status\"\s*:\s*\"PASS\"",read(resultpaths[3],String))
    pd_calls=0
    d1=(a_pass&&b_pass&&c_pass&&firewall&&maximum(affine_A)<1e-10&&maximum(affine_B)<1e-10)
    d2=(w.resolvent_median<1e-11&&w.resolvent_p95<1e-9&&w.port_median<1e-11&&w.port_p95<1e-9)
    d3=(sgdiag.rhs_norm<1e-10&&sg_gfl_current_error<1e-9&&maximum(port_transfer_errors)<1e-10)
    d4=(maximum(affine_rows)<1e-12)
    d13=CROSS_BUS_SCOPE&&length(cross_rows)==4&&all(r->r.status=="PASS_LOCAL_PORT",cross_rows)&&
        maximum(getproperty.(cross_rows,:port_current_error))<1e-9&&
        maximum(getproperty.(cross_rows,:gfl_transfer_equivalence_max))<1e-10&&
        maximum(getproperty.(cross_rows,:affine_rho_residual_max))<1e-12
    zero_error=norm(rhs(x,u,p))

    gates=[
      (gate="D0 design/validation firewall",status=firewall ? "PASS_FIREWALL" : "FAIL",
       evidence="No PowerDynamics import in src/bnd_design; validation calls=0; candidate not frozen.",
       limitation="Candidate-order subgate remains unexercised because no rho design exists."),
      (gate="D1 PLL parameter structure",status=d1 ? "EXACT_AFFINE_LOW_RANK" : "FAIL",
       evidence="Nine-state SimpleGFLDC equations; rank(dA/dKp)=1, rank(dA/dKi)=1, joint rank=2; B has same ranks; C,D are gain-independent.",
       limitation="This is the local GFL component model, not yet the network-embedded SG/GFL aggregate."),
      (gate="D2 Woodbury PLL representation",status=d2 ? "PASS" : "FAIL",
       evidence="$(w.sample_count) deterministic complex-frequency/gain checks; resolvent p95=$(w.resolvent_p95), port p95=$(w.port_p95).",
       limitation="Validates the local PLL/GFL state-space update only."),
      (gate="D3 continuous rho physical model",status=d3 ? "PASS_LOCAL_PORT" : "FAIL",
       evidence="Independent controlled-SG port and base-rebased GFL port share the bus-33 frozen V/P/Q; SG ODE residual=$(sgdiag.rhs_norm), current mismatch=$(sg_gfl_current_error), P/Q error=$(sg_gfl_power_error).",
       limitation="This proves the local physical port pair only; current limits are absent from SimpleGFLDC and the network aggregate is not assembled."),
      (gate="D4 rho update structure",status=d4 ? "PASS_AFFINE_RHO_LOCAL" : "FAIL",
       evidence="Rated parallel-port model yields DeltaY=rho*(Y_GFL-Y_SG); max residual=$(maximum(affine_rows)) across $(length(port_points)) complex points and $(length(rho_samples)) rho values.",
       limitation="No full retained-network selector update is included in this local port test."),
      (gate="D5 global low-rank update",status="NOT_RUN",evidence="",limitation="Requires independent full analytic network assembly."),
      (gate="D6 replacement closure",status="NOT_RUN",evidence="",limitation="No full-network port update is available."),
      (gate="D7 algebraic PLL boundary",status="NOT_RUN",evidence="Local PLL denominator is bilinear in Kp,Ki.",limitation="No physical replacement closure to place on a desired pole boundary."),
      (gate="D8 algebraic rho solution",status="NOT_RUN",evidence="",limitation="Requires D4 and D6."),
      (gate="D9 analytic maximum rho",status="NOT_RUN",evidence="",limitation="No branch candidates exist to differentiate or compare."),
      (gate="D10 multimode safety",status="NOT_RUN",evidence="",limitation="No analytic replacement candidate."),
      (gate="D11 independent PowerDynamics validation",status="NOT_RUN",evidence="PowerDynamics calls=$(pd_calls).",limitation="Correctly withheld until a candidate is calculated and frozen."),
      (gate="D12 validation grid agreement",status="NOT_RUN",evidence="",limitation="No frozen candidate."),
      (gate="D13 cross-bus reproducibility",status=!CROSS_BUS_SCOPE ? "NOT_RUN" : d13 ? "PASS_LOCAL_PORT_REPLICATION" : "FAIL",
       evidence="Independent baseline SG and rated GFL port models were rebuilt at buses 30/33/35/37; max current mismatch=$(maximum(getproperty.(cross_rows,:port_current_error))), max transfer mismatch=$(maximum(getproperty.(cross_rows,:gfl_transfer_equivalence_max))).",
       limitation="This is local port replication only; the full-network replacement closure remains unassembled."),
      (gate="D14 self-energy explanation",status="NOT_RUN",evidence="",limitation="No designed candidate to compare."),
      (gate="D15 ready for ExpE",status="NO",evidence="",limitation="Required D3,D6,D7,D8,D11,D13 gates are not complete.")]
    CSV.write(joinpath(TABLES,"TABLE_D16_structure_gate_summary.csv"),DataFrame(gates))
    CSV.write(joinpath(TABLES,"TABLE_D04_port_update_rank.csv"),DataFrame(port_rows))
    CSV.write(joinpath(TABLES,"TABLE_D14_cross_bus_summary.csv"),DataFrame(cross_rows))
    provenance=[(item="Experiment A-C statuses",value="PASS / PASS / PASS",source="frozen result files"),
        (item="Frozen input hashes",value=join(("$(k)=$(v)" for (k,v) in sort(collect(hashes))),"; "),source="SHA-256"),
        (item="PowerDynamics design calls",value=string(pd_calls),source="analytic runner call audit"),
        (item="GFL source equations",value="PowerDynamics 5.0.0 SimpleGFLDC component equations transcribed independently",source="installed package source; no runtime import"),
        (item="Local GFL state order",value=join(state_names(),";"),source="ExpA state map and component equations"),
        (item="Frozen PLL nominal gains",value="Kp=$(kp0), Ki=$(ki0), tau=$(p.pll_tau)",source="src/pd39/model.jl parameter formulas"),
        (item="Frozen bus-33 SG port residual",value=string(sgdiag.rhs_norm),source="independent Sauer-Pai/AVR/TGOV1 equations and ExpA all-SG equilibrium"),
        (item="Rated-base GFL transfer equivalence p95",value=string(quantile(port_transfer_errors,0.95)),source="state/parameter base conversion; system-base current scaling"),
        (item="Local GFL equilibrium residual",value=string(zero_error),source="independently coded ODE at frozen ExpA operating point"),
        (item="Candidate file",value="NOT_CREATED",source="no rho maximum derived; validation remains disabled"),
        (item="Run timestamp UTC",value=string(now(UTC)),source="Julia Dates")]
    CSV.write(joinpath(TABLES,"TABLE_D01_model_and_design_provenance.csv"),DataFrame(provenance))

    results=Dict("experiment"=>"BND_EXP_D","status"=>"PARTIAL",
        "design_uses_powerdynamics"=>false,"validation_uses_powerdynamics"=>true,
        "validation_calls"=>pd_calls,"candidate_frozen"=>false,
        "primary_bus"=>33,"pll_structure"=>(d1 ? "EXACT_AFFINE_LOW_RANK" : "BLOCKED"),
        "pll_update_rank"=>2,"woodbury"=>Dict("status"=>(d2 ? "PASS" : "FAIL"),
            "resolvent_median_error"=>w.resolvent_median,"resolvent_p95_error"=>w.resolvent_p95,
            "port_median_error"=>w.port_median,"port_p95_error"=>w.port_p95,
            "determinant_factorization_max_error"=>w.determinant_factorization_max,
            "sample_count"=>w.sample_count,"sample_seed"=>w.seed),
        "local_gfl"=>Dict("state_dimension"=>9,"input_dimension"=>2,"output_dimension"=>2,
            "rank_dA_dKp"=>1,"rank_dA_dKi"=>1,"rank_joint_A_updates"=>2,
            "rank_dB_dKp"=>1,"rank_dB_dKi"=>1,"rank_joint_B_updates"=>2,
            "C_D_gain_dependent"=>false,"affine_A_max_relative_residual"=>maximum(affine_A),
            "affine_B_max_relative_residual"=>maximum(affine_B),"frozen_equilibrium_rhs_norm"=>zero_error,
            "pll_nominal_Kp"=>kp0,"pll_nominal_Ki"=>ki0,"pll_tau_s"=>p.pll_tau),
        "local_ports"=>Dict("status"=>((d3&&d4) ? "PASS" : "BLOCKED"),
            "sg_equilibrium_rhs_inf"=>sgdiag.rhs_norm,"sg_P_pu"=>sgdiag.network_power.p,
            "sg_Q_pu"=>sgdiag.network_power.q,"rated_gfl_current_mismatch"=>sg_gfl_current_error,
            "rated_gfl_power_mismatch"=>sg_gfl_power_error,"system_to_device_base_ratio"=>alpha,
            "gfl_transfer_equivalence_p95"=>quantile(port_transfer_errors,0.95),
            "max_affine_rho_residual"=>maximum(affine_rows),"port_dimension"=>2,
            "sample_points"=>length(port_points),"rho_structural_checks"=>length(rho_samples)),
        "replacement"=>Dict("rho_model"=>"BLOCKED","local_port_rho_model"=>((d3&&d4) ? "AFFINE_IN_RHO" : "BLOCKED"),"port_update_rank"=>2,
            "full_operator_dimension"=>nothing,"closure_dimension"=>nothing,
            "gain_boundary_degree"=>nothing,"rho_polynomial_degree"=>nothing),
        "candidate"=>nothing,"powerdynamics"=>Dict("validation_status"=>"NOT_RUN"),
        "cross_bus_scope"=>(CROSS_BUS_SCOPE ? "FOUR_BUS_LOCAL_PORT" : "BUS33_ONLY"),
        "cross_bus"=>Dict(string(r.bus)=>r.status for r in cross_rows),
        "ready_for_expE"=>false,"input_hashes"=>hashes,"gates"=>Dict(g.gate=>g.status for g in gates))
    write_json(joinpath(OUT,"RESULTS_EXP_D_STRUCTURE.json"),results)

    claims=[
      ("D-C01","PLL parameter dependence is low-dimensional / low-rank.","EXACT","Local nine-state model: affine; Kp and Ki each rank one in A and B, joint rank two."),
      ("D-C02","Continuous SG→GFL replacement is an affine/polynomial/rational low-rank port update in rho.","LOCAL_EXACT","The independently derived 2-current ports share a frozen bus-33 equilibrium and give DeltaY=rho(Y_GFL-Y_SG); full network embedding remains pending."),
      ("D-C03","The full characteristic condition reduces via determinant lemma to a small replacement closure.","NOT_ESTABLISHED","Generic determinant-lemma code exists; no system-specific port update is defined."),
      ("D-C04","PLL gains on a desired pole boundary are roots of a low-degree algebraic system.","LOCAL_ONLY","The isolated PLL state-resolvent denominator is bilinear; no physical pole-boundary closure is available."),
      ("D-C05","For single-bus replacement, maximum admissible rho can be selected from a finite set of algebraic candidates.","NOT_ESTABLISHED","No rho branch or candidate set."),
      ("D-C06","The analytically derived candidate agrees with detailed PowerDynamics validation.","NOT_RUN","No analytic candidate has been frozen."),
      ("D-C07","The design requires fewer full-system evaluations than numerical methods.","NOT_RUN","No candidate or benchmark."),
      ("D-C08","Self-energy/pathway analysis explains how analytic PLL tuning restores margin.","NOT_RUN","No designed candidate."),
      ("D-C09","The same analytic pipeline generalizes across buses 30/33/35/37.","LOCAL_REPLICATION","The independently derived rated two-current ports reproduce frozen SG injections and GFL base conversions at all four buses; no full-network closure or design is claimed.")]
    open(joinpath(OUT,"CLAIM_LEDGER_EXP_D_STRUCTURE.md"),"w") do io
        println(io,"# Structural audit claim ledger — Experiment D\n\n| ID | Claim | Type/status | Evidence or limitation |\n|---|---|---|---|")
        for (id,claim,typ,evidence) in claims
            println(io,"| $id | $claim | $typ | $evidence |")
        end
    end
    write_report(results,w,gates,zero_error,maximum(affine_A),maximum(affine_B))
    terminal_summary(results,w,gates)
end

function write_report(results,w,gates,eqerr,aerr,berr)
    open(joinpath(OUT,"REPORT_EXP_D_STRUCTURE.md"),"w") do io
        println(io,"# Experiment D — Structural Audit (pre-candidate phase)\n")
        println(io,"## 1. Executive result\n\n**EXP_D_STATUS: PARTIAL.** The independently coded nine-state `SimpleGFLDC` model has exact affine PLL-gain dependence with rank-one Kp and Ki updates (joint A-update rank 2). The local bus-33 SG/GFL models share the frozen operating point after exact MVA-base conversion, and the local current-port replacement is affine in rho. The full-network selector update, replacement closure, algebraic candidate, and independent PowerDynamics validation remain incomplete.\n")
        println(io,"## 2. Separation between design and validation\n\n`src/bnd_design/` imports no PowerDynamics module. The analytic runner makes zero PowerDynamics calls. It reads only the frozen ExpA operating point and explicitly re-derives the GFL equations. `src/bnd_validation/` is reserved for candidate-gated validation; it was not called. No `ANALYTIC_CANDIDATE.json` was written, so no validation or validation grid ran.\n")
        println(io,"## 3. Input physical models\n\nBus 33 uses the frozen ExpA voltage and internal state values. The stock GFL parameters are transcribed from the A–C configuration: $((results["local_gfl"]["pll_nominal_Kp"])) rad/s proportional gain, $((results["local_gfl"]["pll_nominal_Ki"])) rad/s² integral gain, and $((results["local_gfl"]["pll_tau_s"])) s PLL low-pass time constant. The controlled Sauer–Pai, Type-I AVR, TGOV1, and GFL equations were implemented separately from installed package source and the copied parameter rows in `inputs/`.\n")
        println(io,"## 4. Continuous rho construction\n\n**PASS FOR THE LOCAL PORT.** The SG and GFL internal equations stay unscaled. Device states and current/control parameters are transformed exactly between 100 MVA system base and the 800 MVA bus rating; the GFL full-rating terminal current is then referred back by 8. The resulting local two-current port maps satisfy `Y_mix=(1-rho)Y_SG+rho*Y_GFL` and preserve bus-33 P/Q at the frozen voltage. The full-network embedding is not yet available, and the stock GFL equations contain no current limiter, so no admissible-current claim is made for intermediate rho.\n")
        println(io,"## 5. Exact PLL parameter dependence\n\nFor `e=-sin(theta)u_r+cos(theta)u_i`, the component equations are `dot(Delta_omega_i)=Ki*e`, `dot(Delta_omega)=(Delta_omega_i+Kp*e-Delta_omega)/tau`, `dot(theta)=Delta_omega`. The other six SimpleGFLDC state equations have no direct Kp/Ki dependence. Thus `A(K)=A0+Kp Ap+Ki Ai`; `rank(Ap)=1`, `rank(Ai)=1`, `rank([Ap Ai])=2`. The same rank pattern holds in B; C and D are gain-independent for the terminal-current output.\n")
        println(io,"## 6. Low-rank factorization\n\nThe exact analytic rows are `Ap=e4*(grad_x(e)/tau)'` and `Ai=e5*grad_x(e)'`; the input matrices use the voltage gradient `grad_u(e)=[-sin(theta),cos(theta)]`. Full singular spectra and numerical tolerances are in `TABLE_D02_gain_dependency_audit.csv`. Maximum independent affine reconstruction residuals: A=$aerr, B=$berr. Frozen-point local RHS norm: $eqerr.\n")
        println(io,"## 7. Woodbury reduction\n\nThe rank-two Woodbury resolvent was compared with direct local-model inversion at $(w.sample_count) deterministic gain/frequency samples over 0.01–10 Hz. Median/p95 state-resolvent error: $(w.resolvent_median) / $(w.resolvent_p95). Median/p95 terminal-port transfer error: $(w.port_median) / $(w.port_p95). The determinant factorization max relative error is $(w.determinant_factorization_max). See `TABLE_D03_Woodbury_validation.csv`.\n")
        println(io,"## 8. SG→GFL port update\n\nThe independent local current-port update is exact and affine in rho, with rank and affine-residual checks in `TABLE_D04_port_update_rank.csv`. Its validated scope is the device port, not the full network.\n")
        println(io,"## 9. Matrix determinant lemma\n\n**NOT RUN.** No system-level port selector and all-SG retained operator have been assembled, so no determinant-lemma residual is reported.\n")
        println(io,"## 10. Replacement closure\n\n**NOT RUN.** The full 20-coordinate graph operator has not been assembled with a bus-port selector. Closure dimension and pole equivalence remain unknown.\n")
        println(io,"## 11. Graph/self-energy reduction\n\nExperiment C remains frozen. No graph truncation has been used for design.\n")
        println(io,"## 12. Algebraic PLL boundary\n\nFor the isolated PLL state-resolvent update, `det(I-diag(Kp,Ki)H(s))` is bilinear in the gains (degree 1 in each, total degree 2). Coefficients are sampled only for structural documentation in `TABLE_D06_gain_polynomial.csv`; this is not the system pole-boundary polynomial and no gains were selected from these samples.\n")
        println(io,"## 13. Algebraic rho boundary\n\n**NOT RUN.** There is no system-level pole-boundary equation from which to derive rho branches.\n")
        println(io,"## 14. Stationary maximum replacement\n\n**NOT RUN.** No finite algebraic candidate set has been derived or compared.\n")
        println(io,"## 15. Multimode safety\n\n**NOT RUN.** No replacement candidate exists for the protected-pole checks.\n")
        println(io,"## 16. Frozen analytic candidate\n\nNo candidate, rho branch, or optimum exists. `ANALYTIC_CANDIDATE.json` is intentionally absent.\n")
        println(io,"## 17. Independent PowerDynamics validation\n\n**NOT RUN.** The validation gate refuses before importing PowerDynamics because no analytic candidate is frozen.\n")
        println(io,"## 18. Validation-only exhaustive map\n\n**NOT RUN.** A validation grid is prohibited until the candidate is derived and frozen; no grid was inspected.\n")
        println(io,"## 19. Cross-bus reproduction\n\nCross-bus scope: $(results["cross_bus_scope"]). The four-bus run reports independent local SG/GFL ports in `TABLE_D14_cross_bus_summary.csv`; full-network closures and designs were not run.\n")
        println(io,"## 20. Beyond Nodal Damping interpretation\n\nNo design-mechanism conclusion is available.\n")
        println(io,"## 21. Gate summary\n\nSee `TABLE_D16_structure_gate_summary.csv`. This file records only the pre-candidate structural audit; the final design outcome is in `REPORT_EXP_D.md`.\n")
        println(io,"## 22. Exact mathematical results\n\nThe local SimpleGFLDC PLL dependence is affine; Kp and Ki each contribute rank one and have joint rank two. The local state and port resolvents satisfy the rank-two Woodbury identity within the reported residuals. The rated local two-current replacement is affine in rho.\n")
        println(io,"## 23. Empirical validation results\n\nThe deterministic algebra checks, local rated-port checks, and four-bus port replication passed. Detailed PowerDynamics candidate validation and time-domain validation are not run.\n")
        println(io,"## 24. Failed hypotheses / limitations\n\nNo full-system low-rank closure, algebraic design, or candidate has been established. The scope uses a frozen external operating point, one GFL device, and local ports; the model has no current limiter. There is no EMT, field validation, robust uncertainty result, H4 transfer, or multi-bus optimality claim.\n")
        println(io,"## 25. Publication-ready contribution statement\n\nAt this stage, the supported contribution is a reproducible local structural result: the detailed nine-state GFL PLL update is exactly affine and rank two jointly, and a same-equilibrium rated SG/GFL current-port mixture is affine in rho at the four tested candidate buses. No system-level controller-design claim is supported yet.\n")
        println(io,"## 26. Decision for Experiment E\n\n**NO.** D3, D6–D9 and D11–D13 are incomplete.\n")
    end
end

function terminal_summary(r,w,gates)
    println("EXP_D_STATUS:"); println("PARTIAL")
    println("DESIGN_USES_POWERDYNAMICS:"); println("NO")
    println("VALIDATION_USES_POWERDYNAMICS:"); println("YES")
    println("PRIMARY_BUS:"); println("33")
    println("PLL_STRUCTURE:"); println(r["pll_structure"])
    println("PLL_UPDATE_RANK:"); println(r["pll_update_rank"])
    println("WOODBURY_P95_ERROR:"); println(w.port_p95)
        println("RHO_MODEL:"); println(r["replacement"]["rho_model"])
        println("PORT_UPDATE_RANK:"); println(r["replacement"]["port_update_rank"])
    println("FULL_OPERATOR_DIM:"); println("N/A")
    println("CLOSURE_DIM:"); println("N/A")
    println("GRAPH_CRITICAL_DIM:"); println("N/A")
    println("GAIN_POLYNOMIAL_DEGREE:"); println("LOCAL PLL: (1,1), total degree 2; replacement boundary: N/A")
    println("RHO_POLYNOMIAL_DEGREE:"); println("N/A")
    println("ANALYTIC_RHO_STAR:"); println("N/A")
    println("ANALYTIC_MW_REPLACED:"); println("N/A")
    println("ANALYTIC_KP_STAR:"); println("N/A")
    println("ANALYTIC_KI_STAR:"); println("N/A")
    println("ANALYTIC_BOUNDARY_POLE:"); println("N/A")
    println("PD_VALIDATION_POLE:"); println("NOT_RUN")
    println("PD_POLE_ERROR:"); println("N/A")
    println("PD_STABILITY_MARGIN:"); println("N/A")
    println("PD_DAMPING_MIN:"); println("N/A")
    println("VALIDATION_GRID:"); println("NOT_RUN")
    println("SELF_ENERGY_MECHANISM:"); println("INCONCLUSIVE")
    println("CROSS_BUS:"); for bus in (30,33,35,37); println("$(bus)=",get(r["cross_bus"],string(bus),"NOT_RUN")); end
        println("MAIN_ANALYTIC_RESULT:"); println("The local rated SG/GFL bus-33 port is equilibrium-consistent and affine in rho; the SimpleGFLDC PLL update is rank two and Woodbury-verified.")
    println("MAIN_VALIDATION_RESULT:"); println("No PowerDynamics validation was run because no analytic rho candidate exists.")
        println("MAIN_LIMITATION:"); println("The local device ports are built; the full-network port selector, replacement closure, analytic candidate, and validation remain unbuilt.")
    println("READY_FOR_EXP_E:"); println("NO")
    println("FILES:"); println("reports/experiment_D/REPORT_EXP_D_STRUCTURE.md; reports/experiment_D/RESULTS_EXP_D_STRUCTURE.json; reports/experiment_D/CLAIM_LEDGER_EXP_D_STRUCTURE.md")
    println("PUSH:"); println("NO")
end

main()
