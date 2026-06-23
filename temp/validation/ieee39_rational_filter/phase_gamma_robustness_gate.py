"""Gate Gamma robustness: are resolvent findings artifacts of identity weights?

This script repeats the equal-cost coordination part of Gate 1 on exactly the
same generated networks, actions, seeds, costs, and flow-feasibility rules, but
changes the protection weights used in the protected resolvent margin.

The goal is not to make the result look good.  It asks whether two Gate-1
observations survive reasonable Gamma choices:

1. Modal damping and protected-resolvent certificates rank levers differently.
2. Equal-cost coordinated packages have CV < 0 under the protected resolvent.

The dynamic model is the same explicit second-order model used by Gate 1.  IEEE
cases are topology/flow templates, not full ANDES IBR validations.
"""

from __future__ import annotations

import csv
import json
import math
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import scipy
from scipy.linalg import eigvals, svdvals
from scipy.stats import pearsonr, spearmanr

import phase_coordination_lever_gate as gate1


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs" / "phase_gamma_robustness_gate"
RNG_SEED = gate1.RNG_SEED
RANDOM_GAMMA_COUNT = 8
RANDOM_GAMMA_SIGMA = 0.50


@dataclass(frozen=True)
class GammaProfile:
    name: str
    kind: str
    Wy: np.ndarray
    Wu: np.ndarray
    description: str


def write_json(path: Path, data: Any) -> None:
    def default(obj: Any) -> Any:
        if isinstance(obj, np.generic):
            return obj.item()
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return str(obj)

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, default=default)


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def git_commit_id() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            cwd=ROOT,
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except Exception:
        return "unavailable"


def version_manifest() -> dict[str, str]:
    manifest = {
        "python": sys.version.replace("\n", " "),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
    }
    try:
        import andes  # type: ignore

        manifest["andes"] = getattr(andes, "__version__", "available_unknown_version")
    except Exception:
        manifest["andes"] = "not_imported"
    return manifest


def normalize_diag(weights: np.ndarray) -> np.ndarray:
    weights = np.asarray(weights, dtype=float)
    med = float(np.median(weights))
    if med <= 0:
        return weights
    return weights / med


def gamma_profiles(case: gate1.NetworkCase, rng: np.random.Generator) -> list[GammaProfile]:
    n = case.n
    profiles: list[GammaProfile] = []

    eye = np.eye(n)
    profiles.append(
        GammaProfile(
            name="identity",
            kind="physical",
            Wy=eye,
            Wu=eye,
            description="Baseline: identity input and frequency-output protection weights.",
        )
    )

    profiles.append(
        GammaProfile(
            name="inertia_output",
            kind="physical",
            Wy=np.diag(normalize_diag(1.0 / np.sqrt(case.M))),
            Wu=eye,
            description="Output frequency deviations at low-inertia buses receive larger weight.",
        )
    )

    type_weights = np.where(case.M <= np.median(case.M), 2.0, 1.0)
    profiles.append(
        GammaProfile(
            name="type_input_low_inertia",
            kind="physical",
            Wy=eye,
            Wu=np.diag(normalize_diag(type_weights)),
            description="Inputs at low-inertia buses are weighted as renewable/IBR-like disturbance locations.",
        )
    )

    profiles.append(
        GammaProfile(
            name="frequency_output",
            kind="physical",
            Wy=eye,
            Wu=eye,
            description=(
                "Frequency-focused control. In this model the resolvent output already selects omega, "
                "so this is an explicit equivalence/control check against identity."
            ),
        )
    )

    for k in range(RANDOM_GAMMA_COUNT):
        wy = normalize_diag(np.exp(rng.normal(0.0, RANDOM_GAMMA_SIGMA, n)))
        wu = normalize_diag(np.exp(rng.normal(0.0, RANDOM_GAMMA_SIGMA, n)))
        profiles.append(
            GammaProfile(
                name=f"random_diag_{k:02d}",
                kind="random",
                Wy=np.diag(wy),
                Wu=np.diag(wu),
                description="Positive diagonal random input/output weights, median-normalized.",
            )
        )
    return profiles


def protected_resolvent_score(
    case: gate1.NetworkCase,
    w: np.ndarray,
    D: np.ndarray,
    M: np.ndarray,
    profile: GammaProfile,
) -> float:
    A = gate1.state_matrix(case, w, D, M)
    if np.max(eigvals(A).real) > 1e-7:
        return gate1.BAD_MARGIN
    n2 = A.shape[0]
    n = case.n
    B = np.vstack([np.zeros((n, n)), np.diag(1.0 / M)])
    C = np.hstack([np.zeros((n, n)), np.eye(n)])
    eye = np.eye(n2)
    g = 0.0
    for om in gate1.RESOLVENT_FREQS:
        H = profile.Wy @ C @ np.linalg.solve(1j * om * eye - A, B) @ profile.Wu
        g = max(g, float(svdvals(H)[0]))
    return float(math.log(max(g, 1e-12)))


def score(
    case: gate1.NetworkCase,
    actions: tuple[gate1.Action, ...],
    metric_name: str,
    base_score: float,
    profile: GammaProfile | None,
) -> dict[str, Any]:
    w, D, M, feasible = gate1.apply_actions(case, actions)
    if not feasible:
        return {"actions": actions, "score": gate1.BAD_MARGIN, "delta": gate1.BAD_MARGIN, "feasible": False}
    if metric_name == "S_zeta":
        value = gate1.metric_zeta(case, w, D, M)
    elif metric_name == "S_g":
        if profile is None:
            raise ValueError("S_g requires a Gamma profile")
        value = protected_resolvent_score(case, w, D, M, profile)
    else:
        raise ValueError(metric_name)
    return {"actions": actions, "score": value, "delta": float(value - base_score), "feasible": value < gate1.BAD_MARGIN / 2}


def best_single(
    case: gate1.NetworkCase,
    lever: str,
    eps: float,
    metric_name: str,
    base_score: float,
    profile: GammaProfile | None,
) -> dict[str, Any]:
    rows = [
        score(case, (action,), metric_name, base_score, profile)
        for action in gate1.candidate_actions(case, lever, eps)
    ]
    feasible = [r for r in rows if r["feasible"]]
    return min(feasible, key=lambda r: r["delta"]) if feasible else {"feasible": False, "delta": gate1.BAD_MARGIN}


def ranked_singles(
    case: gate1.NetworkCase,
    lever: str,
    eps: float,
    metric_name: str,
    base_score: float,
    profile: GammaProfile | None,
    topk: int,
) -> list[dict[str, Any]]:
    rows = [
        score(case, (action,), metric_name, base_score, profile)
        for action in gate1.candidate_actions(case, lever, eps)
    ]
    feasible = [r for r in rows if r["feasible"]]
    return sorted(feasible, key=lambda r: r["delta"])[:topk]


def best_combined_screened(
    case: gate1.NetworkCase,
    lever_a: str,
    lever_b: str,
    eps_a: float,
    eps_b: float,
    metric_name: str,
    base_score: float,
    profile: GammaProfile | None,
) -> dict[str, Any]:
    cand_a = ranked_singles(case, lever_a, eps_a, metric_name, base_score, profile, gate1.TOPK_COMBINED)
    cand_b = ranked_singles(case, lever_b, eps_b, metric_name, base_score, profile, gate1.TOPK_COMBINED)
    rows = []
    for ra in cand_a:
        for rb in cand_b:
            rows.append(score(case, tuple(ra["actions"]) + tuple(rb["actions"]), metric_name, base_score, profile))
    feasible = [r for r in rows if r["feasible"]]
    return min(feasible, key=lambda r: r["delta"]) if feasible else {"feasible": False, "delta": gate1.BAD_MARGIN}


def lever_ranking(single: dict[str, dict[str, Any]]) -> str:
    return ">".join(
        lever
        for lever, _ in sorted(
            [(lever, single[lever]["delta"]) for lever in ["w", "D", "M"] if single[lever]["feasible"]],
            key=lambda item: item[1],
        )
    )


def descriptors(case: gate1.NetworkCase) -> dict[str, float]:
    A = gate1.state_matrix(case, case.w, case.D, case.M)
    poles = eigvals(A)
    departure = math.sqrt(max(float(np.linalg.norm(A, "fro") ** 2 - np.sum(np.abs(poles) ** 2)), 0.0))
    L = gate1.laplacian(case.n, case.edges, case.w)
    inv_sqrt_m = np.diag(1.0 / np.sqrt(case.M))
    Lm = inv_sqrt_m @ L @ inv_sqrt_m
    Dm = inv_sqrt_m @ np.diag(case.D) @ inv_sqrt_m
    denom = np.linalg.norm(Lm, "fro") * np.linalg.norm(Dm, "fro")
    chi = float(np.linalg.norm(Lm @ Dm - Dm @ Lm, "fro") / max(denom, 1e-14))
    return {
        "henrici_departure": float(departure),
        "commutator_chi": chi,
    }


def run_case_profile(case: gate1.NetworkCase, profile: GammaProfile) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    base_z = gate1.metric_zeta(case, case.w, case.D, case.M)
    base_g = protected_resolvent_score(case, case.w, case.D, case.M, profile)
    if base_z >= gate1.BAD_MARGIN / 2 or base_g >= gate1.BAD_MARGIN / 2:
        return [], {"case": case.name, "gamma": profile.name, "base_ok": False}

    single_eps = {lever: gate1.EQUAL_COST_J / gate1.COSTS[lever] for lever in ["w", "D", "M"]}
    half_eps = {lever: gate1.EQUAL_COST_J / (2.0 * gate1.COSTS[lever]) for lever in ["w", "D", "M"]}

    single_z = {
        lever: best_single(case, lever, single_eps[lever], "S_zeta", base_z, None)
        for lever in ["w", "D", "M"]
    }
    single_g = {
        lever: best_single(case, lever, single_eps[lever], "S_g", base_g, profile)
        for lever in ["w", "D", "M"]
    }
    rank_z = lever_ranking(single_z)
    rank_g = lever_ranking(single_g)

    rows: list[dict[str, Any]] = []
    for la, lb in [("w", "D"), ("w", "M"), ("D", "M")]:
        if not single_g[la]["feasible"] or not single_g[lb]["feasible"]:
            continue
        combined_g = best_combined_screened(case, la, lb, half_eps[la], half_eps[lb], "S_g", base_g, profile)
        if not combined_g["feasible"]:
            continue
        best_single_g = min(float(single_g[la]["delta"]), float(single_g[lb]["delta"]))
        cv_g = float(combined_g["delta"] - best_single_g)

        if not single_z[la]["feasible"] or not single_z[lb]["feasible"]:
            cv_z = float("nan")
        else:
            combined_z = best_combined_screened(case, la, lb, half_eps[la], half_eps[lb], "S_zeta", base_z, None)
            cv_z = (
                float(combined_z["delta"] - min(float(single_z[la]["delta"]), float(single_z[lb]["delta"])))
                if combined_z["feasible"]
                else float("nan")
            )
        rows.append(
            {
                "case": case.name,
                "family": case.family,
                "n": case.n,
                "gamma": profile.name,
                "gamma_kind": profile.kind,
                "pair": f"{la}+{lb}",
                "S_g_cv": cv_g,
                "S_g_cv_negative": bool(cv_g < 0.0),
                "S_zeta_cv": cv_z,
                "S_zeta_cv_negative": bool(cv_z < 0.0) if not math.isnan(cv_z) else "",
                "same_cv_sign": bool((cv_g < 0.0) == (cv_z < 0.0)) if not math.isnan(cv_z) else "",
                "S_zeta_single_ranking": rank_z,
                "S_g_single_ranking": rank_g,
                "same_single_ranking": bool(rank_z == rank_g),
                "S_g_single_a": float(single_g[la]["delta"]),
                "S_g_single_b": float(single_g[lb]["delta"]),
                "S_g_combined": float(combined_g["delta"]),
                "combined_action": gate1.action_label(combined_g["actions"]),
            }
        )
    desc = descriptors(case)
    summary = {
        "case": case.name,
        "family": case.family,
        "n": case.n,
        "gamma": profile.name,
        "gamma_kind": profile.kind,
        "base_ok": True,
        "S_zeta_single_ranking": rank_z,
        "S_g_single_ranking": rank_g,
        "same_single_ranking": bool(rank_z == rank_g),
        "cv_sign_agreement_fraction": float(np.mean([r["same_cv_sign"] for r in rows])) if rows else float("nan"),
        "cv_disagreement_fraction": float(1.0 - np.mean([r["same_cv_sign"] for r in rows])) if rows else float("nan"),
        **desc,
    }
    return rows, summary


def corr(x: list[float], y: list[float]) -> dict[str, float]:
    mask = np.isfinite(x) & np.isfinite(y)
    xx = np.asarray(x, dtype=float)[mask]
    yy = np.asarray(y, dtype=float)[mask]
    if len(xx) < 3 or np.std(xx) <= 1e-14 or np.std(yy) <= 1e-14:
        return {"n": int(len(xx)), "pearson_r": float("nan"), "pearson_p": float("nan"), "spearman_r": float("nan"), "spearman_p": float("nan")}
    pr = pearsonr(xx, yy)
    sr = spearmanr(xx, yy)
    return {
        "n": int(len(xx)),
        "pearson_r": float(pr.statistic),
        "pearson_p": float(pr.pvalue),
        "spearman_r": float(sr.statistic),
        "spearman_p": float(sr.pvalue),
    }


def aggregate(rows: list[dict[str, Any]], summaries: list[dict[str, Any]]) -> dict[str, Any]:
    df = pd.DataFrame(rows)
    sdf = pd.DataFrame(summaries)
    gamma_rows = []
    for gamma in sorted(df["gamma"].unique()):
        gsub = df[df["gamma"] == gamma]
        ssub = sdf[sdf["gamma"] == gamma]
        kind = str(gsub["gamma_kind"].iloc[0])
        row = {
            "gamma": gamma,
            "gamma_kind": kind,
            "n_pair_cases": int(len(gsub)),
            "single_ranking_agreement_pct": float(ssub["same_single_ranking"].mean()),
            "cv_sign_agreement_pct": float(gsub["same_cv_sign"].mean()),
        }
        for pair in ["w+D", "w+M", "D+M"]:
            psub = gsub[gsub["pair"] == pair]
            row[f"{pair}_cv_lt_0_pct"] = float(psub["S_g_cv_negative"].mean()) if len(psub) else float("nan")
            row[f"{pair}_median_cv"] = float(psub["S_g_cv"].median()) if len(psub) else float("nan")
        gamma_rows.append(row)

    corr_rows = []
    for gamma in sorted(sdf["gamma"].unique()):
        sub = sdf[sdf["gamma"] == gamma]
        rank_disagree = [0.0 if x else 1.0 for x in sub["same_single_ranking"].tolist()]
        cv_disagree = [float(x) for x in sub["cv_disagreement_fraction"].tolist()]
        for descriptor_name in ["henrici_departure", "commutator_chi"]:
            descriptor_values = [float(x) for x in sub[descriptor_name].tolist()]
            for target_name, target_values in [
                ("ranking_disagreement", rank_disagree),
                ("cv_disagreement_fraction", cv_disagree),
            ]:
                c = corr(descriptor_values, target_values)
                corr_rows.append(
                    {
                        "gamma": gamma,
                        "descriptor": descriptor_name,
                        "target": target_name,
                        **c,
                    }
                )

    physical = [r for r in gamma_rows if r["gamma_kind"] == "physical"]
    random_rows = [r for r in gamma_rows if r["gamma_kind"] == "random"]
    discrepancy_robust = all(r["single_ranking_agreement_pct"] < 0.50 for r in physical)
    coordination_robust = all(
        r["w+D_cv_lt_0_pct"] > 0.50 and r["w+M_cv_lt_0_pct"] > 0.50 for r in physical
    )
    # Interpretability is judged across physical gamma profiles by the strongest
    # absolute Spearman relation with either descriptor and ranking disagreement.
    physical_corr = [
        r
        for r in corr_rows
        if r["gamma"] in {g["gamma"] for g in physical} and r["target"] == "ranking_disagreement"
    ]
    max_abs_spearman = max(
        [abs(r["spearman_r"]) for r in physical_corr if np.isfinite(r["spearman_r"])],
        default=float("nan"),
    )
    interpretable = bool(np.isfinite(max_abs_spearman) and max_abs_spearman >= 0.35)

    if discrepancy_robust and interpretable:
        status = "DISCREPANCY_ROBUST_INTERPRETABLE"
        implication = "Modal and protected-resolvent certificates diverge across physical Gamma choices, with descriptor correlation."
    elif discrepancy_robust:
        status = "DISCREPANCY_ROBUST_WEAK_INTERPRETABILITY"
        implication = "The certificate discrepancy survives physical Gamma choices, but descriptor correlations are weak."
    else:
        status = "DISCREPANCY_NOT_ROBUST"
        implication = "The certificate discrepancy can be removed or weakened by a reasonable Gamma choice."

    return {
        "gamma_rows": gamma_rows,
        "correlation_rows": corr_rows,
        "n_cases": int(len(sdf["case"].unique())),
        "n_profiles": int(len(sdf["gamma"].unique())),
        "physical_profiles": [r["gamma"] for r in physical],
        "random_profiles": [r["gamma"] for r in random_rows],
        "discrepancy_robust": discrepancy_robust,
        "coordination_robust": coordination_robust,
        "max_abs_spearman_descriptor_to_ranking_disagreement": max_abs_spearman,
        "status": status,
        "implication": implication,
    }


def write_report(status: dict[str, Any], agg: dict[str, Any]) -> None:
    lines = [
        "# Gate Gamma Robustness",
        "",
        "This gate repeats the equal-cost coordination test using the same networks,",
        "actions, costs, and seeds as Gate 1, but changes the protected-resolvent",
        "input/output weights.",
        "",
        f"- Commit: `{status['git_commit']}`",
        f"- Random seed: `{status['random_seed']}`",
        f"- Cases: `{agg['n_cases']}`",
        f"- Gamma profiles: `{agg['n_profiles']}`",
        f"- Random diagonal Gamma samples: `{RANDOM_GAMMA_COUNT}`",
        "",
        "## Gamma Robustness Table",
        "",
        "| Gamma | kind | ranking agreement S_zeta vs S_g | CV sign agreement | % CV<0 w+D | % CV<0 w+M | % CV<0 D+M |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for r in agg["gamma_rows"]:
        lines.append(
            f"| {r['gamma']} | {r['gamma_kind']} | {100*r['single_ranking_agreement_pct']:.1f}% | "
            f"{100*r['cv_sign_agreement_pct']:.1f}% | {100*r['w+D_cv_lt_0_pct']:.1f}% | "
            f"{100*r['w+M_cv_lt_0_pct']:.1f}% | {100*r['D+M_cv_lt_0_pct']:.1f}% |"
        )
    lines.extend(
        [
            "",
            "## Descriptor Correlations",
            "",
            "Correlations use case-level descriptor values and certificate-disagreement scores.",
            "",
            "| Gamma | descriptor | target | n | Pearson r | Spearman r |",
            "|---|---|---|---:|---:|---:|",
        ]
    )
    for r in agg["correlation_rows"]:
        if r["gamma_kind"] if "gamma_kind" in r else False:
            pass
        if r["target"] != "ranking_disagreement":
            continue
        lines.append(
            f"| {r['gamma']} | {r['descriptor']} | {r['target']} | {r['n']} | "
            f"{r['pearson_r']:.3g} | {r['spearman_r']:.3g} |"
        )
    lines.extend(
        [
            "",
            "## Verdict",
            "",
            f"- Status: **{agg['status']}**",
            f"- Implication: {agg['implication']}",
            f"- Discrepancy robust over physical Gamma choices: `{agg['discrepancy_robust']}`",
            f"- Coordination robust over physical Gamma choices for w+D and w+M: `{agg['coordination_robust']}`",
            f"- Max absolute Spearman descriptor/ranking-disagreement correlation: `{agg['max_abs_spearman_descriptor_to_ranking_disagreement']:.3g}`",
            "",
            "## Scope",
            "",
            "- Gamma matrices are fixed at the base case for each network; actions change the system, not the protection standard.",
            "- The frequency-output profile is intentionally equivalent to identity because Gate 1 already measured omega output.",
            "- This is still a reduced dynamic-model gate, not a full ANDES IBR validation.",
        ]
    )
    (OUT / "phase_gamma_robustness_gate_report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir(parents=True, exist_ok=True)
    rng_cases = np.random.default_rng(RNG_SEED)
    cases = gate1.synthetic_cases(rng_cases) + gate1.ieee_cases(rng_cases)

    all_rows: list[dict[str, Any]] = []
    summaries: list[dict[str, Any]] = []
    for ci, case in enumerate(cases, start=1):
        rng_gamma = np.random.default_rng(RNG_SEED + 1000 + ci)
        profiles = gamma_profiles(case, rng_gamma)
        print(f"[gamma] {ci}/{len(cases)} {case.name} profiles={len(profiles)}")
        for profile in profiles:
            rows, summary = run_case_profile(case, profile)
            all_rows.extend(rows)
            summaries.append(summary)

    agg = aggregate(all_rows, summaries)
    status = {
        "gate": "gamma_robustness",
        "overall_status": agg["status"],
        "git_commit": git_commit_id(),
        "random_seed": RNG_SEED,
        "random_gamma_count": RANDOM_GAMMA_COUNT,
        "random_gamma_sigma": RANDOM_GAMMA_SIGMA,
        "version_manifest": version_manifest(),
        "aggregate": agg,
    }
    write_csv(OUT / "gamma_cv_rows.csv", all_rows)
    write_csv(OUT / "gamma_case_summaries.csv", summaries)
    write_csv(OUT / "gamma_aggregate_table.csv", agg["gamma_rows"])
    write_csv(OUT / "gamma_descriptor_correlations.csv", agg["correlation_rows"])
    write_json(OUT / "phase_gamma_robustness_gate_status.json", status)
    write_report(status, agg)
    print(json.dumps(status, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
