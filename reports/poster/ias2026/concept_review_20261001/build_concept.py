"""Organize existing evidence only; no new model solves or optimizations."""
from pathlib import Path
import csv
import hashlib
import html
import json
import shutil

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[3]
SRC = ROOT / "reports/experiment_Q2B/CERTIFIED_SEARCH"
design = pd.read_csv(SRC / "FINAL_DESIGN_BY_BUS.csv")
marginal = pd.read_csv(SRC / "MARGINAL_SECURITY_VALUES.csv")
events = pd.read_csv(SRC / "PD/EVENTS.csv")
summary = json.loads((SRC / "TERMINAL_SUMMARY.json").read_text())
retained = design[design.retained_MW > 0].copy()
assert abs(design.retained_MW.sum() - summary["retained_MW"]) < 1e-8
assert np.max(np.abs(marginal.total - 1)) < 3e-6

navy, blue, gold, red, muted = "#15364a", "#147da0", "#b86d18", "#b5333b", "#596d78"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "axes.titleweight": "bold", "axes.labelcolor": navy,
                     "text.color": navy, "figure.facecolor": "white"})
fig, axes = plt.subplots(1, 2, figsize=(13, 4.6), gridspec_kw={"width_ratios": [1, 1.2]})
y = np.arange(len(retained))
axes[0].barh(y, retained.retained_MW, color=blue, height=.62)
axes[0].set_yticks(y, [f"Bus {b}" for b in retained.bus])
axes[0].invert_yaxis()
axes[0].set_xlim(0, retained.retained_MW.max() * 1.30)
axes[0].set_xlabel("Retained synchronous dispatch [MW]")
axes[0].set_title("Where synchronous generation remains", loc="left", pad=16)
for j, val in enumerate(retained.retained_MW):
    axes[0].text(val + 2, j, f"{val:.2f}", va="center", fontsize=10)
positions = np.arange(len(marginal))
axes[1].bar(positions - .19, marginal.peak_frequency, .36, label="Frequency peaks", color=blue)
axes[1].bar(positions + .19, marginal.robust, .36, label="Normalized model robustness", color=gold)
axes[1].axhline(0, color=muted, linewidth=.7)
axes[1].set_xticks(positions, marginal.bus)
axes[1].set_xlabel("Retained SG bus")
axes[1].set_ylabel("Multiplier-weighted marginal contribution")
axes[1].set_ylim(-.07, 1.18)
axes[1].set_title("What limits further local removal", loc="left", pad=16)
axes[1].legend(loc="upper center", fontsize=9, frameon=False)
fig.text(.5, .015, "Local KKT interpretation, not a physical partition of MW. Signed contributions sum to one; modal and RoCoF multipliers are zero.",
         ha="center", fontsize=9, color=muted)
fig.tight_layout(rect=[0, .055, 1, 1])
for ext in ("png", "svg"):
    fig.savefig(OUT / f"01_retention_and_constraints.{ext}", dpi=180)
plt.close(fig)

chosen = events.set_index("event").loc[["bus16_100MW", "bus8_100MW", "bus29_100MW"]]
ratios = chosen[["Fpeak", "Rpeak"]].to_numpy() / .5
fig, ax = plt.subplots(figsize=(10.5, 4))
for i in range(3):
    for j in range(2):
        fail = ratios[i, j] > 1
        face = "#fae9e7" if fail else "#e7f2f1"
        ax.add_patch(plt.Rectangle((j, i), 1, 1, facecolor=face, edgecolor="white", linewidth=6))
        actual = chosen.iloc[i][["Fpeak", "Rpeak"][j]]
        unit = "Hz" if j == 0 else "Hz/s"
        label = f"{actual:.3f} {unit}  |  {ratios[i,j]:.3f} × limit"
        ax.text(j+.5, i+.43, label, ha="center", va="center", fontsize=12, color=red if fail else navy)
        ax.text(j+.5, i+.71, "EXCEEDS LIMIT" if fail else "WITHIN LIMIT", ha="center", fontsize=9,
                fontweight="bold", color=red if fail else blue)
ax.set_xlim(0,2)
ax.set_ylim(3,0)
ax.set_xticks([.5,1.5], ["Maximum frequency deviation", "Maximum windowed RoCoF"])
ax.xaxis.tick_top()
ax.set_yticks([.5,1.5,2.5], ["Bus 16 · design event", "Bus 8 · external check", "Bus 29 · external check"])
ax.tick_params(length=0, pad=10)
for sp in ax.spines.values():
    sp.set_visible(False)
fig.text(.5,.02,"PowerDynamics nonlinear trajectories; +100 MW ZIP-setpoint changes; T = 0.5 s; ten generator-bus outputs.",
         ha="center",fontsize=9,color=muted)
fig.tight_layout(rect=[0,.055,1,1])
for ext in ("png","svg"):
    fig.savefig(OUT / f"02_event_scope.{ext}",dpi=180)
plt.close(fig)

shutil.copyfile(SRC / "figures/FIG_C03_LINEAR_VS_PD.png", OUT / "03_original_trajectory_comparison.png")
matrix = [
    ["H4 / IAS26-110","Specific four-site nominal blocker","Collective instability is prior art; M1 remains blocked","Context / backup","reports/poster/ias2026/deliverables/IAS26-120/IAS26-120_ADVERSARIAL_AUDIT.md"],
    ["TX3","Reproducible finite externalities and reconstruction","Material harm on limiting mode not established","Supplement / separate paper","reports/papers/tx3_connected_intervention_calculus/TX3_FINAL_CLAIM_FREEZE.md"],
    ["CDW + IEEE-68","Branch-sensitivity ranking transfers across benchmarks/models","Context reversals do not transfer; materiality caveats","Alternative engineering story","reports/poster/ias2026/research/contextual_dynamic_weakness/ieee68_replication/docs/CDW68_NOVELTY_ASSESSMENT.md"],
    ["B/C + F0-F7","Exact graph-modal diagnostics and conditional predictor","Physical graph scalar separation fails F2","Technical background","reports/experiment_F_program/FINAL_SUMMARY.md"],
    ["Q2B CERTIFIED_SEARCH","Joint co-design; numerical local optimum; nonlinear design-event validation","Fixed support/event; normalized uncertainty; all gain bounds active; no current limiter","Recommended poster core","reports/experiment_Q2B/CERTIFIED_SEARCH/RESULTADO_LOCAL_Y_BRECHA_GLOBAL.md"],
    ["Nonlinear comparison certificate","Regional guarantee for all inputs in declared bounded class","Microscopic normalized input envelope; design frozen","Future theory","reports/nonlinear_codesign_20261001/coupled_dq_contract/CONTROL_MODAL_COMPARACION_Y_CODESIGN_ES.txt"],
]
with (OUT/"EVIDENCE_MATRIX.csv").open("w",newline="",encoding="utf-8-sig") as f:
    writer=csv.writer(f)
    writer.writerow(["line","supported_result","main_limit","editorial_use","source"])
    writer.writerows(matrix)

source_paths = ["FINAL_DESIGN_BY_BUS.csv","MARGINAL_SECURITY_VALUES.csv","PD/EVENTS.csv",
                "TERMINAL_SUMMARY.json","Z_LOCAL_SECURE_FINAL.toml",
                "RESULTADO_LOCAL_Y_BRECHA_GLOBAL.md","figures/FIG_C03_LINEAR_VS_PD.png"]
manifest={"purpose":"Poster content organization from existing results; no new simulations",
          "sources":{str((SRC/p).relative_to(ROOT)):hashlib.sha256((SRC/p).read_bytes()).hexdigest()
                     for p in source_paths},
          "retained_MW_from_csv":float(design.retained_MW.sum()),
          "gfl_fraction_from_summary":summary["GFL_fraction"]}
(OUT/"SOURCE_MANIFEST.json").write_text(json.dumps(manifest,indent=2),encoding="utf-8")

rows="".join(f"<tr><td>{int(r.bus)}</td><td>{r.retained_MW:.2f}</td><td>{r.Kp:.3f}</td><td>{r.Ki:.3f}</td></tr>"
             for r in design.itertuples())
document = """<!doctype html><html lang="en"><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Beyond Nodal Damping — Where Should Synchronous Generation Remain?</title>
<style>
*{box-sizing:border-box}body{margin:0;background:#edf1f2;color:#173547;font-family:Arial,Helvetica,sans-serif;line-height:1.48}
main{max-width:1280px;margin:30px auto;background:white;padding:42px;box-shadow:0 2px 20px #16364a12}
.eyebrow{text-transform:uppercase;letter-spacing:.12em;font-size:12px;font-weight:700;color:#97601a}
h1{font-size:44px;line-height:1.07;max-width:920px;margin:15px 0}h2{font-size:22px;line-height:1.2;margin-top:0}
.lede{max-width:950px;font-size:19px}.scope{background:#edf4f5;padding:15px 18px;border-left:4px solid #147da0;font-size:14px}
.stats{display:grid;grid-template-columns:repeat(3,1fr);gap:20px;margin:25px 0}.stat{padding:20px;background:#15364a;color:white}
.stat strong{display:block;font-size:34px}.stat span{font-size:13px;display:block}
.grid{display:grid;grid-template-columns:1fr 1fr;gap:26px}.panel{border-top:3px solid #b86d18;padding-top:20px;margin-top:25px}
.wide{grid-column:1/-1}img{width:100%;height:auto;display:block}
p{margin:10px 0}.small,figcaption{font-size:13px;color:#4c6472}.formula{font-family:Georgia,serif;font-size:23px;padding:15px;background:#f6f7f7}
.note{background:#fff6e8;padding:16px;border-left:4px solid #b86d18}a{color:#126785}
table{border-collapse:collapse;width:100%;font-size:13px}th,td{border-bottom:1px solid #dce5e9;text-align:right;padding:6px}
th:first-child,td:first-child{text-align:left}footer{margin-top:30px;border-top:1px solid #c9d5dc;padding-top:20px}
@media(max-width:850px){main{margin:0;padding:22px}h1{font-size:34px}.grid{display:block}.stats{gap:8px}.stat{padding:12px}.stat strong{font-size:26px}}
@media print{body{background:white}main{box-shadow:none;margin:0;padding:15px}.panel{break-inside:avoid}}
</style><main>
<div class="eyebrow">Beyond Nodal Damping · Working concept · 1 October 2026</div>
<h1>Where Should Synchronous Generation Remain?</h1>
<p class="lede">Joint SG–GFL and PLL design with numerical local optimality and nonlinear event validation.</p>
<p class="scope"><b>Central finding.</b> In this IEEE-39 local design, frequency peaks and normalized model robustness limit further removal of synchronous generation. The spectral-margin constraint is slack. The remaining generators serve different marginal roles.</p>
<div class="stats"><div class="stat"><strong>__GFL__%</strong><span>of original dispatch assigned to GFL</span></div>
<div class="stat"><strong>__SG__ MW</strong><span>synchronous dispatch retained at five buses</span></div>
<div class="stat"><strong>25 variables</strong><span>five retentions + twenty local PLL gains</span></div></div>
<p class="small">A model-specific research result. Fixed design event: +100 MW ZIP-load setpoint at bus 16.
Causal phase measurements use T = 0.5 s at ten generator buses. This is a content concept, not the final competition poster.</p>
<div class="grid">
<section class="panel"><h2>1. Co-design a feasible replacement</h2>
<div class="formula">min<sub>ε,Kp,Ki</sub> Σ P<sub>i</sub><sup>0</sup> ε<sub>i</sub></div>
<p>Enforce the full physical spectrum, a normalized robustness requirement, and frequency / RoCoF peaks of the analytic model. Preserve competing time-domain peaks as separate active constraints.</p>
<p class="small">The numerical local certificate applies to the analytic security problem on the declared support. Nonlinear validation checks the selected design; it does not establish optimality for nonlinear trajectory constraints.</p>
</section>
<section class="panel"><h2>2. Explain each retained location</h2>
<p>Frequency-peak multipliers dominate at buses 30 and 39. The normalized model-robustness constraint dominates at buses 33, 35 and 37. The modal-margin and RoCoF multipliers are zero at this point.</p>
<p class="small">These are local KKT contributions, not a new law, causal proof, physical shares of reserve, or prices transferable to another grid. All twenty gains are at their declared bounds.</p>
</section>
<section class="panel wide"><img src="01_retention_and_constraints.png" alt="Retained synchronous MW by bus and signed normalized marginal contributions">
<figcaption>Left: the retained dispatch. Right: multiplier-weighted contributions to the local stationarity condition. Negative contributions are retained.</figcaption></section>
<section class="panel wide"><h2>3. Predict, then validate the nonlinear response</h2>
<img src="03_original_trajectory_comparison.png" alt="Original analytic and nonlinear PowerDynamics frequency and RoCoF trajectories">
<figcaption>Unmodified source figure from Q2B/CERTIFIED_SEARCH. The horizontal coordinate is time after the event.
The displayed load step is a ZIP-setpoint change; effective demand varies with voltage.
PowerDynamics peaks: |Δf| = 0.468755 Hz and |RoCoF| = 0.171771 Hz/s; both declared limits are 0.5.</figcaption></section>
<section class="panel wide"><h2>4. Test the design outside its target event</h2>
<img src="02_event_scope.png" alt="Design event passes while two external locations each violate one requirement">
<figcaption>The two external failures remain visible. An all-location or all-input guarantee has not been demonstrated for this design.</figcaption></section>
<section class="panel"><h2>Numerical evidence behind the design</h2>
<ul><li>KKT, LICQ and positive reduced Hessian at the local point.</li>
<li>Riccati / bounded-real inequality checked at 256-bit precision across all frequencies.</li>
<li>Physical matrices, trim and 143 finite poles reconciled with PowerDynamics.</li>
<li>Selected nonlinear event validated after freezing the design.</li></ul>
<p class="small">Numerical checks use declared tolerances, not outward-rounded formal proofs.
The robustness norm acts on a frozen-coordinate reduced operator; it is not a load-disturbance envelope.
The global lower bound is still zero.</p></section>
<section class="panel"><h2>Reproducible design values</h2>
<table><thead><tr><th>Bus</th><th>SG MW</th><th>Kp</th><th>Ki</th></tr></thead><tbody>__ROWS__</tbody></table>
<p class="small">PLL gains are in the frozen model convention. They do not add a sustained active-power resource. Hard converter-current protection is not modeled.</p></section>
</div>
<footer>
<div class="note"><b>Required before a stronger performance claim:</b> optimize a matched common-gain baseline under the same constraints; test gain-bound and physically meaningful parameter sensitivity; expand the declared disturbance set if claiming spatial generality. These experiments have not been run as part of this review.</div>
<p><a href="ELECCION_PARA_COMPETENCIA_ES.txt">Competition strategy and scientific defense (Spanish)</a> ·
<a href="DECISION_Y_GUION_ES.txt">Decision and evidence review (Spanish)</a> ·
<a href="EVIDENCE_MATRIX.csv">Research-line comparison</a> ·
<a href="SOURCE_MANIFEST.json">Source hashes</a> ·
<a href="../../../experiment_Q2B/CERTIFIED_SEARCH/RESULTADO_LOCAL_Y_BRECHA_GLOBAL.md">Full result and scope</a></p>
<p class="small"><b>Prior art acknowledged:</b>
<a href="https://doi.org/10.1109/ISGT-Europe47291.2020.9248794">Zhao &amp; Flynn: high GFL shares and control</a>;
<a href="https://arxiv.org/abs/2412.15446">Chatterjee &amp; Geng: GFM/GFL allocation and tuning</a>;
<a href="https://arxiv.org/abs/2604.17603">Wang &amp; Geng: stability-constrained OPF and nodal shadow prices</a>.
No first-in-literature or global-optimum claim is made.</p>
</footer></main></html>"""
document=document.replace("__GFL__",f'{100*summary["GFL_fraction"]:.2f}').replace("__SG__",f'{summary["retained_MW"]:.2f}').replace("__ROWS__",rows)
(OUT/"POSTER_CONCEPT.html").write_text(document,encoding="utf-8")
print(json.dumps({"status":"generated_from_existing_evidence","output":str(OUT),
                  "retained_MW":float(design.retained_MW.sum()),"source_files_hashed":len(source_paths)}))
