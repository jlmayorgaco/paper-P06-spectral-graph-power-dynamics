# ruff: noqa: E501  -- evidence labels, docstrings and verbatim source quotes kept on one line
"""Assemble the final cross-tool validation deliverables from the run outputs.

Reads only files written by the pack scripts, the supplementary checks and Phase E, plus
the frozen inputs (holdout_lines.json, the frozen F2c predictions). It computes no new
physics except the band eigenvalues for the figure (the same live internal solve as
supp_dynamic_checks.py).

Writes (research root):
  results/CROSS_TOOL_EVIDENCE_TABLE.csv
  results/PANDAPOWER_STATIC_PARITY.csv
  results/ANDES_DYNAMIC_PARITY.csv
  results/ANDES_LINE_SENSITIVITY_HOLDOUT.csv
  results/CROSS_TOOL_VALIDATION_FIGURE.pdf / .png (+ _source.csv)
  validation_runs/cross_tool_final/phaseE/holdout_metrics.json
Usage (tx3-analysis venv, research root): python assemble_final.py .
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

repo = Path(sys.argv[1]).resolve()
HERE = Path(__file__).resolve().parent
RES = repo / "results"
HANDOFF = repo / "validation_inputs/cross_tool_handoff/claude_cross_tool_handoff"
FOREST = HANDOFF / "dynamic_forest_ieee39_validation/dynamic_forest_ieee39_validation"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"


def jload(p: Path) -> dict:
    return json.loads(p.read_text())


py_ref = jload(HERE / "pack/results/python_reference.json")
pp_as = jload(HERE / "pack/results/pandapower_parity.json")
pp_fix = jload(HERE / "pack_tapside/results/pandapower_parity.json")
st = jload(HERE / "supp/static_supp.json")
an_st = jload(HERE / "supp/andes_static_supp.json")
reg = jload(HERE / "supp/andes_regularization_diag.json")
an_dyn = jload(HERE / "pack/results/andes_dynamic_summary.json")
dyn = jload(HERE / "supp/dynamic_supp.json")
dyn_df = pd.read_csv(HERE / "supp/dynamic_supp.csv")
e1_int_info = jload(HERE / "phaseE/internal_E1.json")
e1_an_info = jload(HERE / "phaseE/andes_E1.json")
pid = jload(HERE / "phaseE/port_identity_diag.json")
holdout = jload(FOREST / "holdout_lines.json")
f2c = pd.read_csv(FOREST / "data/F2c_holdout_reequilibrated_port_line_sensitivity.csv")
e1i = pd.read_csv(HERE / "phaseE/internal_E1.csv")
e1a = pd.read_csv(HERE / "phaseE/andes_E1.csv")
assert (
    list(f2c.line_index)
    == list(e1i.line_index)
    == list(e1a.line_index)
    == holdout["lines"]
)

# ------------------------------------------------------------------ evidence rows
rows: list[dict] = []


def add(phase, check, comparison, measured, gate, ok, note=""):
    verdict = ok if isinstance(ok, str) else ("PASS" if ok else "FAIL")
    rows.append(
        {
            "phase": phase,
            "check": check,
            "comparison": comparison,
            "measured": measured,
            "gate": gate,
            "verdict": verdict,
            "note": note,
        }
    )


# Phase A
add(
    "A",
    "PF convergence mismatch",
    "Python",
    py_ref["pf_max_mismatch"],
    "<1e-10",
    py_ref["pf_max_mismatch"] < 1e-10,
)
add(
    "A",
    "pytest ieee39 network+baseline",
    "Python",
    py_ref["pytest_stdout"].splitlines()[-1],
    "all pass",
    py_ref["pytest_returncode"] == 0,
)
add(
    "A",
    "Ybus max abs error (pu)",
    "Python vs frozen ANDES npy + restored shunts",
    py_ref["ybus_max_abs_error"],
    "<=1e-9",
    py_ref["ybus_max_abs_error"] <= 1e-9,
    "pack method; stale npy omits bus-4/5 shunts, pack adds them",
)
add(
    "A",
    "Vm max abs error (pu)",
    "Python vs frozen ANDES PF",
    py_ref["vm_max_abs_error_pu"],
    "<=1e-5",
    py_ref["vm_max_abs_error_pu"] <= 1e-5,
)
add(
    "A",
    "Va gauge-aligned max error (rad)",
    "Python vs frozen ANDES PF",
    py_ref["va_max_gauge_aligned_error_rad"],
    "<=1e-5",
    py_ref["va_max_gauge_aligned_error_rad"] <= 1e-5,
)
a = st["A_supp"]
add(
    "A",
    "PV generator P max error (pu)",
    "Python vs frozen ANDES PF",
    a["pv_p"],
    "<=1e-8",
    a["pv_p"] <= 1e-8,
    "supplementary",
)
add(
    "A",
    "slack generator P error (pu)",
    "Python vs frozen ANDES PF",
    a["slack_p"],
    "<=1e-8",
    "FAIL (literal); EXPLAINED",
    "ANDES line equations add +1e-8 to r and x (line.py:200); mirrored in Python -> 1.5e-12 (DIAGNOSTIC)",
)
add(
    "A",
    "generator Q max error (pu)",
    "Python vs frozen ANDES PF",
    max(a["pv_q"], a["slack_q"]),
    "<=1e-5",
    max(a["pv_q"], a["slack_q"]) <= 1e-5,
    "supplementary",
)

# Phase B
for tag, pk, sp in (
    ("as supplied", pp_as, st["B_as_supplied"]),
    ("translation-corrected", pp_fix, st["B_corrected"]),
):
    cmp = f"pandapower 3.4.0 ({tag}) vs Python"
    add("B", "PF converged", cmp, pk["converged"], "True", bool(pk["converged"]))
    add(
        "B",
        "Ybus max abs error (pu)",
        cmp,
        sp["ybus_max_abs_error"],
        "<=1e-9",
        sp["ybus_max_abs_error"] <= 1e-9,
        "supplementary; as supplied, from_ppc puts the tap on the HV winding of branches 35, 37, 38"
        if tag == "as supplied"
        else "supplementary; branches 35, 37, 38 reoriented HV->LV (tap 1/t, z*t^2)",
    )
    add(
        "B",
        "Vm max abs error (pu)",
        cmp,
        pk["vm_max_abs_error_pu"],
        "<=1e-5",
        pk["vm_max_abs_error_pu"] <= 1e-5,
        "pack",
    )
    add(
        "B",
        "Va gauge-aligned max error (rad)",
        cmp,
        pk["va_max_gauge_aligned_error_rad"],
        "<=1e-5",
        pk["va_max_gauge_aligned_error_rad"] <= 1e-5,
        "pack",
    )
    add(
        "B",
        "generator P max error incl. slack (pu)",
        cmp,
        max(sp["gen_P_max_abs_error_pu"], sp["slack_P_abs_error_pu"]),
        "<=1e-8",
        max(sp["gen_P_max_abs_error_pu"], sp["slack_P_abs_error_pu"]) <= 1e-8,
        "supplementary",
    )
    add(
        "B",
        "generator Q max error incl. slack (pu)",
        cmp,
        max(sp["gen_Q_max_abs_error_pu"], sp["slack_Q_abs_error_pu"]),
        "<=1e-5",
        max(sp["gen_Q_max_abs_error_pu"], sp["slack_Q_abs_error_pu"]) <= 1e-5,
        "supplementary",
    )
    add(
        "B",
        "branch terminal P/Q max error (MVA), pack formula",
        cmp,
        pk["canonical_branch_terminal_flow_max_error_MVA"],
        "<=1e-3",
        pk["canonical_branch_terminal_flow_max_error_MVA"] <= 1e-3,
        "pack",
    )
    add(
        "B",
        "branch terminal P/Q max error (MVA), pandapower native",
        cmp,
        sp["native_branch_terminal_flow_max_error_MVA"],
        "<=1e-3",
        sp["native_branch_terminal_flow_max_error_MVA"] <= 1e-3,
        "supplementary",
    )

# Phase C
add(
    "C",
    "case file SHA256 == canonical source",
    "ANDES 2.0.0 live",
    an_st["case_sha256"],
    an_st["canonical_source_sha256"],
    an_st["case_hash_matches_canonical_source"],
)
add(
    "C",
    "live PF == frozen ANDES reference (Vm, Va, P, Q)",
    "ANDES live vs frozen",
    max(an_st["vs_frozen_reference"].values()),
    "==0 (bit-identical)",
    max(an_st["vs_frozen_reference"].values()) == 0.0,
)
add(
    "C",
    "live Ybus (Line+Shunt) vs Python (pu)",
    "ANDES live vs Python",
    an_st["ybus_live_andes_vs_python_max_abs_error"],
    "<=1e-9",
    an_st["ybus_live_andes_vs_python_max_abs_error"] <= 1e-9,
    "stale npy not used for any verdict",
)
add(
    "C",
    "shunts and line order/taps identical to canonical",
    "ANDES live vs canonical JSON",
    an_st["andes_shunts"]
    == [
        {"bus": int(s["bus"]), "g": s["g"], "b": s["b"]}
        for s in an_st["canonical_shunts"]
    ]
    and an_st["line_order_bus1_bus2_tap_identical_to_canonical"]
    and an_st["line_b1_b2_g1_g2_all_zero"],
    "True",
    True,
)
v = an_st["vs_python_pf"]
add(
    "C",
    "Vm max abs error (pu)",
    "ANDES live vs Python",
    v["vm_max_abs_error_pu"],
    "<=1e-5",
    v["vm_max_abs_error_pu"] <= 1e-5,
)
add(
    "C",
    "Va gauge-aligned max error (rad)",
    "ANDES live vs Python",
    v["va_max_gauge_aligned_error_rad"],
    "<=1e-5",
    v["va_max_gauge_aligned_error_rad"] <= 1e-5,
)
add(
    "C",
    "PV generator P max error (pu)",
    "ANDES live vs Python",
    v["pv_p_max_abs_error_pu"],
    "<=1e-8",
    v["pv_p_max_abs_error_pu"] <= 1e-8,
)
add(
    "C",
    "slack generator P error (pu)",
    "ANDES live vs Python",
    v["slack_p_abs_error_pu"],
    "<=1e-8",
    "FAIL (literal); EXPLAINED",
    "unchanged at ANDES tol 1e-12 -> not convergence; ANDES +1e-8 impedance regularization",
)
add(
    "C",
    "generator Q max error (pu)",
    "ANDES live vs Python",
    max(v["pv_q_max_abs_error_pu"], v["slack_q_abs_error_pu"]),
    "<=1e-5",
    max(v["pv_q_max_abs_error_pu"], v["slack_q_abs_error_pu"]) <= 1e-5,
)
r2 = reg["python_with_andes_1e-8_regularization_vs_andes"]
add(
    "C",
    "DIAGNOSTIC all static errors with ANDES regularization mirrored",
    "ANDES (tol 1e-12) vs Python(+1e-8)",
    max(r2[k] for k in r2 if k != "python_pf_mismatch"),
    "<=1e-8 (diagnostic)",
    max(r2[k] for k in r2 if k != "python_pf_mismatch") <= 1e-8,
    "attribution test only; frozen model unchanged",
)

# Phase D
r3 = dyn["R3 first-order AVR"]
add(
    "D",
    "R3 critical-band alpha max error (1/s)",
    "ANDES fresh vs internal (pack)",
    an_dyn["R3_max_alpha_error"],
    "<=1e-4",
    an_dyn["R3_max_alpha_error"] <= 1e-4,
)
add(
    "D",
    "R3 critical-band frequency max error (Hz)",
    "ANDES fresh vs internal (pack)",
    an_dyn["R3_max_frequency_error_hz"],
    "<=1e-3",
    an_dyn["R3_max_frequency_error_hz"] <= 1e-3,
)
add(
    "D",
    "R3 RHP count agrees, every subset",
    "ANDES fresh vs internal (pack)",
    an_dyn["R3_all_RHP_counts_agree"],
    "True",
    bool(an_dyn["R3_all_RHP_counts_agree"]),
)
for stage, s in dyn.items():
    add(
        "D",
        f"{stage[:2]} band (0.3-1.5 Hz) nearest-mode mismatch, both directions",
        "ANDES fresh vs internal live",
        s["max_band_mismatch"],
        "<=1e-4",
        s["max_band_mismatch"] <= 1e-4,
        f"{s['n_subsets']} subsets",
    )
    add(
        "D",
        f"{stage[:2]} RHP count agrees, every subset (live)",
        "ANDES fresh vs internal live",
        s["rhp_agree_all"],
        "True",
        s["rhp_agree_all"],
    )
    add(
        "D",
        f"{stage[:2]} whole dynamic spectrum worst match (1/s)",
        "ANDES fresh vs internal live",
        s["max_whole_spectrum_match"],
        "informative",
        "INFO",
        "only unmatched internal eigenvalues: 12 zero-gain PSS filter states (-1/4.2, -1/0.75), column coupling exactly 0",
    )
    add(
        "D",
        f"{stage[:2]} repeatability: fresh vs frozen ANDES spectra",
        "ANDES",
        s["max_fresh_vs_frozen_andes"],
        "==0",
        s["max_fresh_vs_frozen_andes"] == 0.0,
    )


# ------------------------------------------------------------------ Phase E metrics
def metrics(pred: np.ndarray, ref: np.ndarray) -> dict:
    """pred/ref complex arrays aligned by holdout line; ref is the independent result."""
    pr, rr = pred.real, ref.real
    signs = int(np.sum(np.sign(pr) == np.sign(rr)))
    rho = float(spearmanr(pr, rr).statistic)
    top_p = set(np.argsort(pr)[:5])
    top_r = set(np.argsort(rr)[:5])
    rel = np.abs(pred - ref) / np.abs(ref)
    strong = np.abs(rr) > np.median(np.abs(rr))
    strong_opp = int(np.sum(strong & (np.sign(pr) != np.sign(rr))))
    worst = int(np.argmax(rel))
    return {
        "sign_agreement": signs,
        "n": int(len(pr)),
        "spearman_re": rho,
        "top5_stabilizing_overlap": len(top_p & top_r),
        "median_relative_complex_error": float(np.median(rel)),
        "max_relative_complex_error": float(rel.max()),
        "worst_line": int(holdout["lines"][worst]),
        "strong_branches_opposite_sign": strong_opp,
        "sign_mismatch_lines": [
            int(holdout["lines"][k])
            for k in range(len(pr))
            if np.sign(pr[k]) != np.sign(rr[k])
        ],
        "PASS_frozen_strong": bool(
            rho >= 0.90
            and signs >= 11
            and len(top_p & top_r) >= 4
            and np.median(rel) <= 0.10
            and strong_opp == 0
        ),
        "PREFERRED_12_12_and_rho_0p95": bool(signs == 12 and rho >= 0.95),
    }


c = lambda re, im: re.to_numpy() + 1j * im.to_numpy()  # noqa: E731
andes = c(e1a.andes_dlambda_real, e1a.andes_dlambda_imag)
int_dae = c(e1i.internal_dae_dlambda_real, e1i.internal_dae_dlambda_imag)
int_port = c(e1i.internal_port_ds_real, e1i.internal_port_ds_imag)
f2_port = c(f2c.port_reeq_ds_real, f2c.port_reeq_ds_imag)
f2_dae = c(f2c.dae_dlambda_real, f2c.dae_dlambda_imag)
M = {
    "E1a internal DAE vs ANDES (equation-equivalent R3, same cases)": metrics(
        int_dae, andes
    ),
    "E1b internal port vs ANDES (equation-equivalent R3, same cases)": metrics(
        int_port, andes
    ),
    "E1c internal port vs internal DAE (R3)": metrics(int_port, int_dae),
    "E2a frozen F2c GFL port vs ANDES R3 (CROSS-MODEL)": metrics(f2_port, andes),
    "E2b frozen F2c GFL DAE vs ANDES R3 (CROSS-MODEL)": metrics(f2_dae, andes),
    "E2c frozen F2c GFL DAE vs internal R3 DAE (same tool, converter model changed)": metrics(
        f2_dae, int_dae
    ),
    "E0 frozen F2c GFL port vs frozen F2c GFL DAE (frozen record)": metrics(
        f2_port, f2_dae
    ),
}
(HERE / "phaseE/holdout_metrics.json").write_text(json.dumps(M, indent=2))
for name, m in M.items():
    gate = "rho>=0.90, signs>=11/12, top5>=4/5, median rel<=10%, no strong opposite"
    measured = (
        f"signs {m['sign_agreement']}/12; rho {m['spearman_re']:.4f}; top5 {m['top5_stabilizing_overlap']}/5; "
        f"median rel {m['median_relative_complex_error']:.2e}; max rel {m['max_relative_complex_error']:.2e} "
        f"(L{m['worst_line']}); strong opposite {m['strong_branches_opposite_sign']}"
    )
    note = "preferred 12/12 & rho>=0.95: " + (
        "yes" if m["PREFERRED_12_12_and_rho_0p95"] else "no"
    )
    if m["sign_mismatch_lines"]:
        note += f"; sign mismatch on lines {m['sign_mismatch_lines']}"
    add(
        "E",
        "12-line branch-sensitivity holdout",
        name,
        measured,
        gate,
        m["PASS_frozen_strong"],
        note,
    )
add(
    "E",
    "port characteristic zero |mu(lambda0)+1| (R3 static)",
    "internal",
    e1_int_info["port_mu_at_lambda0_plus_1_abs"],
    "<=1e-8",
    "FAIL (literal); EXPLAINED",
    "base/flagship equilibria differ by 2.5e-7 in this static-injection mode (tol-independent); residual off the port rows",
)
add(
    "E",
    "det identity det(T_S)/det(T_0)=det(I+M) (R3 static)",
    "internal",
    e1_int_info["det_identity_rel_error_at_j2pi0p7"],
    "<=1e-8",
    "FAIL (literal); EXPLAINED",
    "same cause; frozen GFL boundaries: <=1.4e-13 over 160 checks (FC18)",
)
add(
    "E",
    "gamma=1 critical mode alpha",
    "internal vs ANDES",
    abs(e1_int_info["lambda0_real"] - e1_an_info["lambda0_real"]),
    "<=1e-4",
    abs(e1_int_info["lambda0_real"] - e1_an_info["lambda0_real"]) <= 1e-4,
    f"internal {e1_int_info['lambda0_real']:+.6f}, ANDES {e1_an_info['lambda0_real']:+.6f} at {e1_an_info['lambda0_freq_hz']:.4f} Hz",
)

ev = pd.DataFrame(rows)
ev.to_csv(RES / "CROSS_TOOL_EVIDENCE_TABLE.csv", index=False)

# ------------------------------------------------------------------ per-element CSVs
pps = pd.read_csv(HERE / "supp/pandapower_static_parity_long.csv")
pps = pps[pps.comparison.str.startswith("B_")].rename(
    columns={"comparison": "variant", "reference": "python", "other": "pandapower"}
)
pps["variant"] = pps.variant.map(
    {"B_as_supplied": "as_supplied", "B_corrected": "translation_corrected"}
)
pps.to_csv(RES / "PANDAPOWER_STATIC_PARITY.csv", index=False)
dyn_df.to_csv(RES / "ANDES_DYNAMIC_PARITY.csv", index=False)
hold = pd.DataFrame(
    {
        "line_index": e1a.line_index,
        "line": e1a.line,
        "andes_R3_dlambda_real": andes.real,
        "andes_R3_dlambda_imag": andes.imag,
        "internal_R3_dae_dlambda_real": int_dae.real,
        "internal_R3_dae_dlambda_imag": int_dae.imag,
        "internal_R3_port_ds_real": int_port.real,
        "internal_R3_port_ds_imag": int_port.imag,
        "frozen_F2c_GFL_port_ds_real": f2_port.real,
        "frozen_F2c_GFL_port_ds_imag": f2_port.imag,
        "frozen_F2c_GFL_dae_dlambda_real": f2_dae.real,
        "frozen_F2c_GFL_dae_dlambda_imag": f2_dae.imag,
    }
)
hold["rel_err_internal_dae_vs_andes"] = np.abs(int_dae - andes) / np.abs(andes)
hold["rel_err_internal_port_vs_andes"] = np.abs(int_port - andes) / np.abs(andes)
hold["rel_err_frozen_gfl_port_vs_andes_CROSS_MODEL"] = np.abs(f2_port - andes) / np.abs(
    andes
)
hold["rank_andes_R3"] = pd.Series(andes.real).rank().astype(int).to_numpy()
hold["rank_internal_R3_port"] = pd.Series(int_port.real).rank().astype(int).to_numpy()
hold["rank_frozen_F2c_GFL_port"] = pd.Series(f2_port.real).rank().astype(int).to_numpy()
hold.to_csv(RES / "ANDES_LINE_SENSITIVITY_HOLDOUT.csv", index=False)

# ------------------------------------------------------------------ figure
sys.path[:0] = [str(repo / "src"), str(repo / "experiments")]
import F1_reconciliation as f1  # noqa: E402
from ibr_cycles.models.ieee39_case import ReplacementPlan, solve_case  # noqa: E402

spectra = np.load(HERE / "pack/results/F1_andes_spectra.npz")
band_int, band_an = [], []
for stage, services in f1.STAGES.items():
    for members in f1.SUBSETS:
        label = "+".join(map(str, members)) or "BASE"
        plan = (
            ReplacementPlan.of({b: 1.0 for b in members}, device="static_power")
            if members
            else ReplacementPlan.of({})
        )
        lam = np.linalg.eigvals(solve_case(plan, machine_services=services).system.A)
        an = spectra[f"{stage}|{label}"]
        for arr, sink in ((lam, band_int), (an, band_an)):
            f = arr.imag / (2 * np.pi)
            sel = arr[(f >= 0.3) & (f <= 1.5)]
            sink.extend((stage[:2], label, z.real, z.imag / (2 * np.pi)) for z in sel)
bi = pd.DataFrame(band_int, columns=["stage", "subset", "alpha", "freq_hz"])
ba = pd.DataFrame(band_an, columns=["stage", "subset", "alpha", "freq_hz"])

pyf = pd.read_csv(HERE / "pack/results/python_branch_flows.csv")
ppf_fix = pd.read_csv(HERE / "pack_tapside/results/pandapower_branch_flows.csv")
ppf_as = pd.read_csv(HERE / "pack/results/pandapower_branch_flows.csv")
cols = ["P_from_MW", "Q_from_Mvar", "P_to_MW", "Q_to_Mvar"]

plt.rcParams.update(
    {
        "font.size": 8,
        "axes.edgecolor": INK2,
        "axes.labelcolor": INK,
        "xtick.color": INK2,
        "ytick.color": INK2,
        "axes.linewidth": 0.6,
        "font.family": "DejaVu Sans",
    }
)
fig, ax = plt.subplots(1, 3, figsize=(10.5, 3.5), constrained_layout=True)
# (A)
x = pyf[cols].to_numpy().ravel()
ax[0].axline((0, 0), slope=1, color=GRID, lw=1.0, zorder=0)
ax[0].scatter(
    x,
    ppf_as[cols].to_numpy().ravel(),
    s=16,
    facecolors="none",
    edgecolors=ORANGE,
    lw=0.8,
    label=f"as supplied (max err {pp_as['canonical_branch_terminal_flow_max_error_MVA']:.0f} MVA)",
)
ax[0].scatter(
    x,
    ppf_fix[cols].to_numpy().ravel(),
    s=9,
    color=BLUE,
    edgecolors="white",
    lw=0.3,
    label=f"translation-corrected (max err {pp_fix['canonical_branch_terminal_flow_max_error_MVA']:.1e} MVA)",
)
ax[0].set_xlabel("Python branch terminal P/Q (MW, Mvar)")
ax[0].set_ylabel("pandapower 3.4.0 (MW, Mvar)")
ax[0].set_title("(A) Static branch-flow parity", loc="left", color=INK, fontsize=9)
ax[0].legend(frameon=False, fontsize=6.5, loc="upper left")
# (B)
ax[1].scatter(
    ba.alpha,
    ba.freq_hz,
    s=34,
    facecolors="none",
    edgecolors=ORANGE,
    lw=0.9,
    label="ANDES 2.0.0",
)
ax[1].scatter(
    bi.alpha,
    bi.freq_hz,
    s=8,
    color=BLUE,
    edgecolors="white",
    lw=0.3,
    label="Python (project DAE)",
)
ax[1].axvline(0, color=GRID, lw=1.0, zorder=0)
ax[1].set_xlabel("Re λ (1/s)")
ax[1].set_ylabel("frequency (Hz)")
worst_band = max(s["max_band_mismatch"] for s in dyn.values())
ax[1].set_title(
    "(B) Electromechanical band, R2 + R3", loc="left", color=INK, fontsize=9
)
ax[1].text(
    0.98,
    0.03,
    f"12 cases · worst nearest-mode\nmismatch {worst_band:.1e} 1/s",
    transform=ax[1].transAxes,
    ha="right",
    va="bottom",
    fontsize=6.5,
    color=INK2,
)
ax[1].legend(frameon=False, fontsize=6.5, loc="center left")
# (C)
xa = andes.real
ax[2].axline((0, 0), slope=1, color=GRID, lw=1.0, zorder=0)
ax[2].axhline(0, color=GRID, lw=0.6, zorder=0)
ax[2].axvline(0, color=GRID, lw=0.6, zorder=0)
m1, m2 = (
    M["E1b internal port vs ANDES (equation-equivalent R3, same cases)"],
    M["E2a frozen F2c GFL port vs ANDES R3 (CROSS-MODEL)"],
)
ax[2].scatter(
    xa,
    f2_port.real,
    s=22,
    marker="D",
    facecolors="none",
    edgecolors=ORANGE,
    lw=0.8,
    label=f"frozen GFL port (cross-model)\nρ={m2['spearman_re']:.2f}, signs {m2['sign_agreement']}/12",
)
ax[2].scatter(
    xa,
    int_port.real,
    s=14,
    color=BLUE,
    edgecolors="white",
    lw=0.3,
    label=f"port, same R3 case\nρ={m1['spearman_re']:.2f}, signs {m1['sign_agreement']}/12",
)
for k, li in enumerate(holdout["lines"]):
    if (
        abs(xa[k]) > 0.03 or xa[k] > 0.01
    ):  # selective direct labels: strongest branches only
        ax[2].annotate(
            f"L{li}",
            (xa[k], int_port.real[k]),
            xytext=(4, -3),
            textcoords="offset points",
            fontsize=6,
            color=INK2,
        )
ax[2].set_xlabel("ANDES Re dλ/dγ_e (R3 equation-equivalent)")
ax[2].set_ylabel("predicted Re dλ/dγ_e (1/s)")
ax[2].set_title("(C) 12-line holdout sensitivity", loc="left", color=INK, fontsize=9)
ax[2].legend(frameon=False, fontsize=6.0, loc="lower right")
for a_ in ax:
    a_.grid(color=GRID, lw=0.4)
    a_.set_axisbelow(True)
    for sp_ in ("top", "right"):
        a_.spines[sp_].set_visible(False)
fig.savefig(RES / "CROSS_TOOL_VALIDATION_FIGURE.pdf")
fig.savefig(RES / "CROSS_TOOL_VALIDATION_FIGURE.png", dpi=200)
src = pd.concat(
    [
        pd.DataFrame(
            {
                "panel": "A",
                "series": "as_supplied",
                "x": x,
                "y": ppf_as[cols].to_numpy().ravel(),
            }
        ),
        pd.DataFrame(
            {
                "panel": "A",
                "series": "translation_corrected",
                "x": x,
                "y": ppf_fix[cols].to_numpy().ravel(),
            }
        ),
        pd.DataFrame(
            {
                "panel": "B",
                "series": "ANDES",
                "x": ba.alpha,
                "y": ba.freq_hz,
                "label": ba.stage + "|" + ba.subset,
            }
        ),
        pd.DataFrame(
            {
                "panel": "B",
                "series": "Python",
                "x": bi.alpha,
                "y": bi.freq_hz,
                "label": bi.stage + "|" + bi.subset,
            }
        ),
        pd.DataFrame(
            {
                "panel": "C",
                "series": "internal_R3_port",
                "x": xa,
                "y": int_port.real,
                "label": holdout["lines"],
            }
        ),
        pd.DataFrame(
            {
                "panel": "C",
                "series": "frozen_F2c_GFL_port",
                "x": xa,
                "y": f2_port.real,
                "label": holdout["lines"],
            }
        ),
    ]
)
src.to_csv(RES / "CROSS_TOOL_VALIDATION_FIGURE_source.csv", index=False)

print(ev.to_string(index=False, max_colwidth=70))
print(json.dumps(M, indent=1))
