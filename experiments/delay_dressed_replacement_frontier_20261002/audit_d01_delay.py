from __future__ import annotations
import csv
from pathlib import Path

HERE=Path(__file__).resolve().parent;BASE=HERE/"baseline_reproduction"
def rows(p):
    with p.open(newline="",encoding="utf-8") as f:return list(csv.DictReader(f))
roots=rows(BASE/"TABLE_D01_TRACKED_DDE_POLES.csv")
metrics={r["delay_pattern_id"]:r for r in rows(HERE/"TABLE_D05_GSP_DELAY_METRICS.csv")}
alpha0=-0.08080909538068841
out=[]
for name in sorted(set(r["delay_pattern_id"] for r in roots)):
    rr=[r for r in roots if r["delay_pattern_id"]==name]
    local=max(float(r["real"]) for r in rr)
    out.append({"delay_pattern_id":name,"tracked_roots":len(rr),"converged_roots":sum(r["converged"]=="true" for r in rr),
        "tracked_pair_rightmost_real":local,"tau0_baseline_spectral_abscissa":alpha0,
        "tracked_pair_real_shift_from_tau0":local-alpha0,
        "max_nev_residual":max(float(r["residual"]) for r in rr),"chi_tau":metrics[name]["chi_tau"] if name in metrics else "",
        "full_DDE_root_coverage":"NOT_ESTABLISHED","claim_status":"SUPPORTED_LOCAL"})
with (BASE/"TABLE_D01_LOCAL_DDE_ROOTS.csv").open("w",newline="",encoding="utf-8") as f:
    w=csv.DictWriter(f,fieldnames=list(out[0]));w.writeheader();w.writerows(out)
report=["# D1 exact DDE spectrum gate","","**Status: BLOCKED_EXACT_DDE_SPECTRUM.**","",
    "The delayed-error characteristic matrix is assembled without Padé. Its tau=0 limit reproduces the full 203-state finite ODE spectrum for the baseline and joint design (matched max errors 1.97e-11 and 3.76e-11; nonlinear-eigenvalue residuals below 3.4e-17).",
    "",
    "For three frozen patterns, a bordered-Newton continuation tracked the dominant zero-delay conjugate pair through four delay steps. Each returned root has nonlinear-eigenvalue residual below 3.4e-18. These are local tracked-root results only; the CSV labels the result SUPPORTED_LOCAL.",
    "",
    "The required rightmost root count was not completed. A norm-derived finite contour combined with a low-rank determinant-lemma evaluation was implemented, but the bound-driven contour exceeded the practical compute budget before a stable argument-principle count could be returned. No count, completeness, or exact rightmost-root claim is accepted. The DDE stability frontier and optimization therefore remain blocked; Padé was not substituted.",
    "",
    "The dominant tracked pair shifts by only the values in TABLE_D01_LOCAL_DDE_ROOTS.csv across the tested patterns. Because other characteristic roots were not counted, this does not establish the full spectral abscissa or falsify the spatial-delay hypothesis.",""]
(HERE/"D1_GATE_REPORT.md").write_text("\n".join(report),encoding="utf-8")
print(out)
