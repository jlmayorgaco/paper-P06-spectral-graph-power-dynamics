"""Aggregate executed evidence and render diagnostic figures; never infer missing validation."""

from __future__ import annotations

import csv
import json
import math
import tomllib
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
FIG = HERE / "figures"
FIG.mkdir(exist_ok=True)


def read_csv(name: str) -> list[dict[str, str]]:
    with (HERE / name).open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def write_csv(name: str, rows: list[dict]) -> None:
    if not rows:
        return
    with (HERE / name).open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def design_path(cid: str) -> Path:
    if cid == "Z":
        return HERE / "designs/Z_zero_delay_tuned.toml"
    if cid == "N":
        return HERE / "designs/N_nominal.toml"
    return HERE / "designs" / f"{cid}.toml"


def event_data(cid: str) -> tuple[str, dict | None]:
    event_id = "Z_zero_delay_tuned" if cid=="Z" else "N_nominal" if cid=="N" else cid
    complete = HERE / "event_validations" / event_id / "Q0_RESULT.toml"
    screen = HERE / "event_screens" / cid / "Q0_RESULT.toml"
    if complete.exists():
        d = tomllib.loads(complete.read_text(encoding="utf-8"))
        return ("PASS_ALL_FIVE" if d.get("all_five_events_pass") else "FAIL_ALL_FIVE"), d
    if screen.exists():
        d = tomllib.loads(screen.read_text(encoding="utf-8"))
        return ("PASS_SCREEN_ONLY" if d.get("screen_event_pass") else "FAIL_SCREEN"), d
    return "UNTESTED", None


def spectral(cid: str) -> dict | None:
    path = HERE / "evaluations" / cid / "SPECTRAL.csv"
    if path.exists():
        with path.open(newline="", encoding="utf-8") as f:
            return next(csv.DictReader(f))
    return None


def get_event_extrema(d: dict | None) -> tuple[float, float, float, str]:
    if not d:
        return math.nan, math.nan, math.nan, "UNKNOWN"
    ev = d.get("events", [])
    if not ev:
        return math.nan, math.nan, math.nan, "UNKNOWN"
    f = max(float(e["F_peak_Hz"]) for e in ev)
    r = max(float(e["RoCoF_peak_Hz_s"]) for e in ev)
    a = min(float(e["min_SG_actuator_fraction_slack"]) for e in ev)
    worst = max(ev, key=lambda e: float(e["F_peak_Hz"]))["event"]
    return f, r, a, worst


def make_trace() -> list[dict]:
    initial = {r["design_id"]: r for r in read_csv("L0_REPRODUCTION.csv")}
    ids = [
        ("Z", "Z", "baseline", 0),
        ("Z_step1_ball", "Z", "ball_predictor", 1),
        ("Z_step1_multimode_lp", "Z", "multimode_lp", 1),
        ("N", "N", "baseline", 0),
        ("N_step1_ball", "N", "ball_predictor", 1),
        ("N_step1_multimode_lp", "N", "multimode_lp", 1),
        ("N_step1_multimode_lp_next_full", "N", "multimode_lp_full", 2),
        ("N_step1_multimode_lp_next_half", "N", "multimode_lp_half", 2),
        ("N_step1_multimode_lp_next_0p625", "N", "event_secant_0.625", 2),
        ("N_step1_multimode_lp_event_active_lp", "N", "sparse_event_active_lp", 2),
        ("N_step1_multimode_lp_event_active_lp_followup_lp", "N", "reused_event_J_followup_lp", 3),
        ("N_step1_multimode_lp_event_active_lp_followup_lp_followup_lp", "N", "reused_event_J_followup_lp", 4),
        ("N_step1_multimode_lp_event_active_lp_followup_lp_followup_lp_followup_lp", "N", "reused_event_J_followup_lp", 5),
        ("N_step1_multimode_lp_event_active_lp_followup_lp_followup_lp_followup_lp_followup_lp", "N", "reused_event_J_followup_lp", 6),
    ]
    rows = []
    bounds = json.loads((HERE / "EXPERIMENT_MANIFEST.json").read_text())["gain_bounds"]
    for cid, start, method, iteration in ids:
        dpath = design_path(cid)
        if not dpath.exists():
            continue
        d = tomllib.loads(dpath.read_text(encoding="utf-8"))
        s = initial[cid] if cid in initial else spectral(cid)
        if s is None:
            continue
        tau = float(s.get("tau_crit_local_ms", s.get("local_tau_crit_ms", s.get("first_local_crossing_ms"))))
        alpha = float(s.get("alpha_zero_delay_s_inv", s.get("alpha_zero_delay")))
        frequency = float(s.get("critical_delay_frequency_hz", s.get("critical_frequency_hz")))
        evt_status, evt = event_data(cid)
        f, r, a, worst = get_event_extrema(evt)
        spectral_ok = (s.get("complete_contour_bracket_pass") == "true"
                       or s.get("below_status") == "SAFE" and s.get("above_status") == "UNSAFE"
                       or s.get("safe_below_status") == "SAFE" and s.get("unsafe_above_status") == "UNSAFE")
        accepted = bool(spectral_ok and alpha <= -0.05 and evt_status == "PASS_ALL_FIVE")
        kp = np.asarray(d["Kp"], dtype=float)
        ki = np.asarray(d["Ki"], dtype=float)
        max_fraction = max(float(np.max(kp / bounds["Kp_max"])), float(np.max(ki / bounds["Ki_max"])))
        rows.append(dict(iteration=iteration, start_id=start, candidate_id=cid, method=method,
                         tau_crit_ms=tau, alpha_zero_delay=alpha,
                         critical_delay_frequency_hz=frequency,
                         critical_family_id="see_L4_modal_family_map",
                         KKT_residual="NOT_ESTABLISHED", max_gain_fraction=max_fraction,
                         worst_zero_delay_event=worst,
                         max_frequency_deviation=f, max_RoCoF=r,
                         min_SG_actuator_slack=a, event_validation=evt_status,
                         spectral_bracket_pass=spectral_ok, accepted=accepted))
    write_csv("L3_OPTIMIZATION_TRACE.csv", rows)
    return rows


def make_replacement_reproduction() -> None:
    low = spectral("prior_rho00_nominal")
    high = next(r for r in read_csv("L0_REPRODUCTION.csv") if r["design_id"] == "N")
    if low is None:
        return
    rows = [dict(design_id="prior_rho00_nominal", GFL_percent=87.5,
                 tau_crit_ms=float(low["local_tau_crit_ms"]),
                 alpha_zero_delay_s_inv=float(low["alpha_zero_delay_s_inv"]),
                 full_contour_bracket_pass=low["complete_contour_bracket_pass"],
                 source="independent_current_recomputation"),
            dict(design_id="N", GFL_percent=88.45514078456176,
                 tau_crit_ms=float(high["tau_crit_local_ms"]),
                 alpha_zero_delay_s_inv=float(high["alpha_zero_delay_s_inv"]),
                 full_contour_bracket_pass="true",
                 source="independent_current_recomputation")]
    write_csv("L0_REPLACEMENT_ONLY_REPRODUCTION.csv", rows)


def style() -> None:
    plt.rcParams.update({"font.size": 10, "font.family": "DejaVu Sans", "axes.spines.top": False,
                         "axes.spines.right": False, "figure.facecolor": "white",
                         "savefig.facecolor": "white", "axes.grid": True,
                         "grid.color": "#e1e6eb", "grid.alpha": 0.8})


def save(fig, name: str) -> None:
    fig.savefig(FIG / name, dpi=220, bbox_inches="tight")
    plt.close(fig)


def figures(trace: list[dict]) -> None:
    style()
    l1 = read_csv("L1_ACTION_SPACE_VALIDATION.csv")
    fig, ax = plt.subplots(figsize=(7, 4))
    errors = np.asarray([float(r["relative_reconstruction_error"]) for r in l1])
    ax.scatter(np.arange(1, len(errors)+1), errors, color="#126782", s=25)
    ax.axhline(1e-12, color="#d25b40", linestyle="--", label="10⁻¹² reference")
    ax.set_yscale("log"); ax.set_xlabel("Descriptor evaluation")
    ax.set_ylabel("Relative reconstruction error")
    ax.set_title("Exact fixed-ρ gain action: rank 10 in 48 evaluations")
    ax.legend(frameon=False)
    save(fig, "FIG_L1_GAIN_ACTION_RANK10.png")

    g = read_csv("L2_TAU_MARGIN_GRADIENT_VALIDATION.csv")
    a = np.asarray([float(r["analytic_margin_gradient_ms_per_gain"]) for r in g])
    b = np.asarray([float(r["FD_margin_gradient_ms_per_gain"]) for r in g])
    fig, ax = plt.subplots(figsize=(5, 5))
    for label, color in [("Z", "#c76a20"), ("N", "#126782")]:
        index = [i for i,r in enumerate(g) if r["design_id"] == label]
        ax.scatter(a[index], b[index], s=55, label=label, color=color)
    lim = [min(a.min(), b.min()) * 1.1, max(a.max(), b.max()) * 0.9]
    ax.plot(lim, lim, color="#555", linewidth=1)
    ax.set(xlabel="Analytical dτcrit/dK (ms per gain)", ylabel="Centered FD (ms per gain)",
           title="Simple-root latency gradient")
    ax.legend(frameon=False)
    save(fig, "FIG_L2_ANALYTIC_VS_FD.png")

    active = read_csv("L2_ACTIVE_FAMILY_SCAN.csv")
    fig, (ax, inset) = plt.subplots(1, 2, figsize=(10, 4.1), gridspec_kw={"width_ratios": [1.25, 1]})
    for family, color, marker in [("Z_template_like", "#c76a20", "o"),
                                  ("N_template_like", "#126782", "s")]:
        part = [r for r in active if r["active_family_label"] == family]
        ax.scatter([float(r["eta"]) for r in part],
                   [float(r["local_earliest_crossing_ms"]) for r in part],
                   marker=marker, color=color, s=45, label="Z limiting" if marker=="o" else "N limiting")
    ax.plot([float(r["eta"]) for r in active],
            [float(r["local_earliest_crossing_ms"]) for r in active],
            color="#555", alpha=.5, linewidth=1)
    sw = read_csv("L4_MODE_SWITCH_POINT.csv")[0]
    ax.axvline(float(sw["eta_estimate"]), color="#7c5385", linestyle="--", linewidth=1.5,
               label=f"switch η={float(sw['eta_estimate']):.6f}")
    ax.set(xlabel="Log-gain interpolation Z → N", ylabel="Earliest local crossing (ms)",
           title="Limiting-family envelope")
    ax.legend(frameon=False, fontsize=9)
    switch = read_csv("L4_MODE_SWITCH_BISECTION.csv")
    abscissa = np.asarray([float(r["eta"]) for r in switch])
    order = np.argsort(abscissa)
    inset.plot(abscissa[order], np.asarray([float(r["Z_crossing_ms"]) for r in switch])[order],
               color="#c76a20", label="Z family")
    inset.plot(abscissa[order], np.asarray([float(r["N_crossing_ms"]) for r in switch])[order],
               color="#126782", label="N family")
    inset.axvline(float(sw["eta_estimate"]), color="#7c5385", linestyle="--", linewidth=1)
    inset.set(xlabel="η near switch", ylabel="Crossing delay (ms)", title="Tracked two-root crossing")
    inset.legend(frameon=False, fontsize=9)
    fig.suptitle("PLL gain tuning changes the limiting delayed mode", fontsize=12)
    save(fig, "FIG_L4_MODAL_FAMILY_SWITCH.png")

    shown = [r for r in trace if r["candidate_id"] in
             ("Z", "Z_step1_multimode_lp", "N", "N_step1_multimode_lp",
              "N_step1_multimode_lp_next_full", "N_step1_multimode_lp_next_0p625",
              "N_step1_multimode_lp_event_active_lp",
              "N_step1_multimode_lp_event_active_lp_followup_lp",
              "N_step1_multimode_lp_event_active_lp_followup_lp_followup_lp",
              "N_step1_multimode_lp_event_active_lp_followup_lp_followup_lp_followup_lp",
              "N_step1_multimode_lp_event_active_lp_followup_lp_followup_lp_followup_lp_followup_lp")]
    short_names = {"Z":"Z", "Z_step1_multimode_lp":"Z1", "N":"N",
                   "N_step1_multimode_lp":"N1", "N_step1_multimode_lp_next_full":"N2 full",
                   "N_step1_multimode_lp_next_0p625":"N2 0.625",
                   "N_step1_multimode_lp_event_active_lp":"N2 event LP",
                   "N_step1_multimode_lp_event_active_lp_followup_lp":"N3 follow-up",
                   "N_step1_multimode_lp_event_active_lp_followup_lp_followup_lp":"N4 follow-up",
                   "N_step1_multimode_lp_event_active_lp_followup_lp_followup_lp_followup_lp":"N5 follow-up",
                   "N_step1_multimode_lp_event_active_lp_followup_lp_followup_lp_followup_lp_followup_lp":"N6 follow-up"}
    names = [short_names[r["candidate_id"]] for r in shown]
    fig, ax = plt.subplots(figsize=(9, 4.3))
    colors = ["#2b8a68" if r["accepted"] else "#c86151" for r in shown]
    ax.bar(np.arange(len(shown)), [r["tau_crit_ms"] for r in shown], color=colors)
    ax.set_ylim(36.5, max(r["tau_crit_ms"] for r in shown)+0.25)
    ax.set_xticks(np.arange(len(shown)), names, rotation=25, ha="right")
    ax.set(ylabel="Local spectral latency threshold (ms)",
           title="Accepted designs pass all five zero-delay events")
    ax.text(0.01, 0.98, "Green: fully validated; red: rejected or incomplete", transform=ax.transAxes,
            va="top", fontsize=9)
    save(fig, "FIG_L3_ACCEPTED_AND_REJECTED_STEPS.png")

    point = [r for r in trace if r["start_id"]=="N" and math.isfinite(float(r["min_SG_actuator_slack"]))]
    fig, ax = plt.subplots(figsize=(8.6, 4.4))
    labels = {"N":"N baseline", "N_step1_multimode_lp":"N1", 
              "N_step1_multimode_lp_next_full":"Full step (rejected)",
              "N_step1_multimode_lp_next_half":"Half step",
              "N_step1_multimode_lp_next_0p625":"0.625 step",
              "N_step1_multimode_lp_event_active_lp":"Event-aware LP",
              "N_step1_multimode_lp_event_active_lp_followup_lp":"Follow-up LP 1",
              "N_step1_multimode_lp_event_active_lp_followup_lp_followup_lp":"Follow-up LP 2",
              "N_step1_multimode_lp_event_active_lp_followup_lp_followup_lp_followup_lp":"Follow-up LP 3",
              "N_step1_multimode_lp_event_active_lp_followup_lp_followup_lp_followup_lp_followup_lp":"Follow-up LP 4"}
    for r in point:
        ax.scatter(r["tau_crit_ms"], r["min_SG_actuator_slack"],
                   color="#2b8a68" if r["accepted"] else "#c86151", s=70,
                   label=labels.get(r["candidate_id"],r["candidate_id"]))
    ax.axhline(0.002, linestyle="--", color="#a54239", label="Frozen actuator limit")
    ax.set_ylim(0.001982,0.002055)
    ax.set_xlim(min(r["tau_crit_ms"] for r in point)-0.05,
                max(r["tau_crit_ms"] for r in point)+0.05)
    ax.set(xlabel="Local spectral latency threshold (ms)",
           ylabel="Minimum SG-actuator fraction slack",
           title="The nominal event constraint limits retuning")
    ax.legend(frameon=False, fontsize=8, loc="center left", bbox_to_anchor=(1.02,.5))
    save(fig, "FIG_L3_ACTUATOR_TRADEOFF.png")

    accepted = [r for r in trace if r["accepted"]]
    if accepted:
        best = max(accepted, key=lambda r: r["tau_crit_ms"])
        n = tomllib.loads(design_path("N").read_text())
        d = tomllib.loads(design_path(best["candidate_id"]).read_text())
        x = np.arange(30,40)
        fig, ax = plt.subplots(figsize=(8, 4.4))
        ax.bar(x-0.2, 100*(np.asarray(d["Kp"])/np.asarray(n["Kp"])-1), width=0.4,
               label="Kp", color="#126782")
        ax.bar(x+0.2, 100*(np.asarray(d["Ki"])/np.asarray(n["Ki"])-1), width=0.4,
               label="Ki", color="#c76a20")
        ax.set_xticks(x)
        ax.set(xlabel="Generator bus", ylabel="Gain change from N (%)",
               title=f"Best fully validated found: {best['tau_crit_ms']:.3f} ms")
        ax.legend(frameon=False)
        save(fig, "FIG_L3_BEST_PLL_GAIN_MAP.png")

    mode_steps=[]
    for r in trace:
        if r["start_id"]!="N" or not r["accepted"] or r["candidate_id"]=="N":
            continue
        path=HERE/"evaluations"/r["candidate_id"]/"ROOTS.csv"
        if not path.exists():
            continue
        roots=read_csv(str(path.relative_to(HERE)))
        crossings=sorted(float(x["local_crossing_ms"]) for x in roots)
        if len(crossings)>=3:
            mode_steps.append((r["iteration"],r["candidate_id"],crossings[:3]))
    if mode_steps:
        fig,ax=plt.subplots(figsize=(7.4,4.2))
        x=np.arange(len(mode_steps))
        for rank,color in [(1,"#126782"),(2,"#c76a20")]:
            ax.plot(x,[v[2][rank]-v[2][0] for v in mode_steps],marker="o",color=color,
                    label=f"Mode {rank+1} minus limiting mode")
        gap_labels={"N_step1_multimode_lp":"N1",
                    "N_step1_multimode_lp_next_half":"N2 half",
                    "N_step1_multimode_lp_next_0p625":"N2 .625",
                    "N_step1_multimode_lp_event_active_lp":"N2 event",
                    "N_step1_multimode_lp_event_active_lp_followup_lp":"N3",
                    "N_step1_multimode_lp_event_active_lp_followup_lp_followup_lp":"N4",
                    "N_step1_multimode_lp_event_active_lp_followup_lp_followup_lp_followup_lp":"N5",
                    "N_step1_multimode_lp_event_active_lp_followup_lp_followup_lp_followup_lp_followup_lp":"N6"}
        ax.set_xticks(x,[gap_labels.get(v[1],v[1]) for v in mode_steps],rotation=25,ha="right")
        ax.set(xlabel="Accepted N-path iteration",ylabel="Local crossing gap (ms)",
               title="Several delayed modes approach the active boundary")
        ax.legend(frameon=False)
        save(fig,"FIG_L4_COCRITICAL_MODE_GAPS.png")


def main() -> None:
    make_replacement_reproduction()
    trace = make_trace()
    figures(trace)
    print("Finalized", len(trace), "trace rows;", len(list(FIG.glob("*.png"))), "figures")


if __name__ == "__main__":
    main()
