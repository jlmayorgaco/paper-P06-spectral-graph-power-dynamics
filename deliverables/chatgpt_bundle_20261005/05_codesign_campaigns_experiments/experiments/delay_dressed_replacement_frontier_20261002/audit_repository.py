from pathlib import Path
import hashlib, json, subprocess, tomllib, sys, platform
ROOT=Path(__file__).resolve().parents[2]; OUT=Path(__file__).resolve().parent
sources=[
"Project.toml","Manifest.toml",
"src/pd39/PD39.jl","src/pd39/model.jl","src/pd39/equilibrium.jl","src/pd39/stability.jl","src/pd39/tds.jl",
"experiments/analytic_iteration_20261001/run_experiment.jl","experiments/analytic_iteration_20261001/validate_pd.jl",
"experiments/analytic_iteration_20261001/report.py","experiments/nonlinear_codesign_20261001/ReducedDAE.jl",
"experiments/nonlinear_codesign_20261001/PDPhysicalReference.jl","src/bnd_model_expN/PDExactDesignN.jl",
"src/bnd_design_e/CollectiveModel.jl",
"reports/analytic_iteration_20261001/refined/protocol.toml",
"reports/analytic_iteration_20261001/refined/baseline.toml",
"reports/analytic_iteration_20261001/refined/joint_final.toml",
"reports/analytic_iteration_20261001/refined/fixed_gains_final.toml",
"reports/analytic_iteration_20261001/refined/RESULT.json",
"reports/analytic_iteration_20261001/refined/RESULTADOS_ES.txt",
"reports/analytic_iteration_20261001/refined/source_hashes_before.toml",
"reports/analytic_iteration_20261001/refined/source_hashes_after.toml",
"reports/experiment_N/MODEL_FREEZE.json","reports/experiment_N/FINAL_SUMMARY_EXP_N.md",
"reports/experiment_P/FINAL_SUMMARY_EXP_P.md","reports/experiment_P/README_REPRODUCE.md","reports/experiment_P/PROVENANCE.toml",
"src/bnd_graphpll/PhysicalGraph.jl","src/bnd_graphpll/ModalGainLaw.jl",
"src/bnd_graphpll/NodeVariantGraphLaw.jl","src/bnd_graphpll/ModalGFL.jl",
"src/bnd_graphpll/SelfEnergyCorrection.jl","src/bnd_graphpll/GraphPLLContinuation.jl",
"src/bnd_graphpll/GraphFilterRealization.jl","src/bnd_graphpll/IdealDamping.jl","src/bnd_graphpll/SafetyMetrics.jl"]
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
hashes={r:sha(ROOT/r) for r in sources if (ROOT/r).is_file()}
missing=[r for r in sources if not (ROOT/r).is_file()]
head=subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip()
branch=subprocess.check_output(["git","branch","--show-current"],cwd=ROOT,text=True).strip()
status=subprocess.check_output(["git","status","--porcelain=v1"],cwd=ROOT,text=True,errors="replace")
dirty=[x[3:].strip().replace("\\","/") for x in status.splitlines() if not x[3:].strip().replace("\\","/").startswith("experiments/delay_dressed_replacement_frontier_20261002/")]
man=tomllib.loads((ROOT/"Manifest.toml").read_text(encoding="utf-8")); deps=man.get("deps",{})
pkgs={}
for name in ["PowerDynamics","NetworkDynamics","SciMLBase","OrdinaryDiffEqRosenbrock","OrdinaryDiffEqNonlinearSolve","ModelingToolkitBase","ForwardDiff","DelayDiffEq","NonlinearEigenproblems"]:
 e=deps.get(name,[]);e=[e] if isinstance(e,dict) else e;pkgs[name]=[x.get("version","stdlib/path") for x in e if isinstance(x,dict)]
events=[{"id":"design_bus8_plus100","bus":8,"delta_load_MW":100},{"id":"holdout_bus8_minus100","bus":8,"delta_load_MW":-100},{"id":"holdout_bus16_plus100","bus":16,"delta_load_MW":100},{"id":"holdout_bus16_minus100","bus":16,"delta_load_MW":-100},{"id":"holdout_bus29_plus100","bus":29,"delta_load_MW":100},{"id":"holdout_bus29_minus100","bus":29,"delta_load_MW":-100}]
data={
"snapshot":{"git_head":head,"branch":branch,"preexisting_dirty_tree":True,"preexisting_dirty_path_count":len(dirty),"preexisting_dirty_paths":dirty},
"environment":{"os":platform.platform(),"python":sys.version,"julia_version":"1.11.9","manifest_versions":pkgs,"verified_loaded_versions":{"PowerDynamics":"5.0.0","NetworkDynamics":"1.3.0","SciMLBase":"3.53.2"},"delay_solver_in_manifest":False,"nonlinear_eigen_solver_in_manifest":False},
"source_sha256":hashes,"missing_expected_sources":missing,
"model":{"physical":"PowerDynamics official IEEE-39, 100 MVA/60 Hz; ten SG/GFL candidate ports 30:39; mixed devices share the bus voltage.",
"reduction":"ReducedDAE.jl eliminates passive algebraic network buses with a Schur/Kron solve and retains device dynamic states for fixed support.",
"rho":"ExpN frozen contract: P/Q share at frozen operating voltage; SG Sn scales by epsilon=1-rho with H fixed; GFL terminal current port weighted by rho; use initialized device dispatch, not net bus injection.",
"pll":"Physical PLL_LPF error e_ang=-sin(theta)u_r+cos(theta)u_i; d(xi)/dt=Ki*e_ang; d(DeltaOmega)/dt=(xi+Kp*e_ang-DeltaOmega)/tau_lpf; d(theta)/dt=DeltaOmega. Existing tau_lpf=1/(2*pi*300) s is a frequency-output low-pass, not measurement delay.",
"delay":"No pure measurement/computation delay is implemented in canonical models.",
"solvers":"Equilibrium PD39.initialize_equilibrium; mixed trim PDPhysicalReference.trim_state/residual_audit; ODE eigenspectrum from ReducedDAE/PowerDynamics linearization; event ODE uses Rodas5P. No DDE solver in frozen manifest.",
"graph":"src/bnd_graphpll contains passive port graph, modal gain law, self-energy correction, graph continuation/filter realization, and safety metrics. Registered F2 physical modal separation failed (off-diagonal reactive ratio 0.2914 > 0.05)."},
"baseline_contract":{"prior_campaign":"reports/analytic_iteration_20261001/refined","initial_rho":"uniform 0.875, buses 30:39","nominal_Kp":"2*pi*5 rad/s","nominal_Ki":"(2*pi*5)^2/4 rad/s^2","design_event":events[0],"frozen_events":events[1:],"limits":{"modal_decay_margin_s_inv":0.05,"frequency_peak_Hz":0.5,"rocof_Hz_s":0.5,"window_s":0.5,"voltage_pu":[0.9,1.1],"minimum_actuator_fraction_slack":0.002,"hard_current_dc_limits_certified":False},"reproduction_status":"NOT_RUN_YET; prior numbers are reference only","prior_joint_holdouts":"4/5 fail at least one declared frequency/actuator constraint"}}
(OUT/"BASELINE_MANIFEST.json").write_text(json.dumps(data,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
doc=["# Model provenance and audit","",f"- Repository HEAD {head}, branch {branch}.","- The pre-existing working tree was dirty in {len(dirty)} paths before this experiment folder was created. Exact path list is in BASELINE_MANIFEST.json; this experiment will not modify those files.","- Julia 1.11.9; PowerDynamics 5.0.0; NetworkDynamics 1.3.0; SciMLBase 3.53.2. Source hashes are in BASELINE_MANIFEST.json.","- Canonical model sources: experiments/nonlinear_codesign_20261001/PDPhysicalReference.jl, experiments/nonlinear_codesign_20261001/ReducedDAE.jl, src/pd39/model.jl, src/bnd_model_expN/PDExactDesignN.jl, and src/bnd_design_e/CollectiveModel.jl.","- The physical PLL uses e_ang=-sin(theta)u_r+cos(theta)u_i. Kp and Ki both act on this error, but the frequency-output channel also has the existing 300-Hz low-pass. No pure time delay is present.","- The frozen event set is design bus 8 +100 MW plus holdouts bus 8 -100 MW, bus 16 +/-100 MW, and bus 29 +/-100 MW. Metrics use a 0.5-s window. Declared exploratory limits: frequency 0.5 Hz, RoCoF 0.5 Hz/s, voltage 0.9-1.1 pu, modal decay 0.05 s^-1, actuator slack 0.002. Current/DC hard limits are not certified.","- ExpN uses dispatch-weighted component semantics: SG rating scales with epsilon=1-rho while H is fixed; GFL terminal current is scaled by rho. Preserve initialized generator P/Q.","- Equilibrium: PowerDynamics initialize_equilibrium and PDPhysicalReference trim. Spectrum: fixed-support reduced ODE and full PowerDynamics Jacobian. Nonlinear integration is ODE Rodas5P; the manifest contains no DDE or nonlinear-eigen solver.","- The prior joint candidate fails 4/5 frozen external events; prior ExpN/ExpP summaries also do not establish globality. These remain warnings, not new reproduced results.","","## Hashed sources",""]
doc += [f"- {r}: SHA-256 {h}" for r,h in hashes.items()]
if missing: doc += ["","Missing expected paths:"]+[f"- {r}" for r in missing]
(OUT/"MODEL_PROVENANCE.md").write_text("\n".join(doc)+"\n",encoding="utf-8")
print("audit artifacts written",len(hashes),"source hashes;",len(dirty),"pre-existing dirty paths;",len(missing),"missing sources")
