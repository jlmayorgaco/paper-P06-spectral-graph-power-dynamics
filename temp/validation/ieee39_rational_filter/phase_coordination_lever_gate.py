"""Gate 1: value of coordinated planning levers under flow feasibility.

This gate tests whether coordinated packages of topology, damping, and inertia
actions outperform the best single lever at equal normalized cost.  It is
deliberately conservative about claims:

* The ground truth is the complete eigensolve of the tested second-order dynamic
  model, not a modal screen.
* Topology actions are accepted only if the post-action DC flow remains feasible.
* Both a modal damping margin and an input-output resolvent margin are reported.

The test bed is synthetic connected Watts-Strogatz graphs plus two IEEE-derived
topologies loaded from local ANDES cases (WSCC9 and IEEE14).  The IEEE cases are
used for topology/flow structure only; the dynamic model in this gate is the
explicit low-inertia swing model controlled by the parameters w, D, and M.
"""

from __future__ import annotations

import csv
import json
import math
import shutil
import subprocess
import sys
from dataclasses import dataclass
from itertools import product
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
import scipy
from scipy.linalg import eigvals, svdvals


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs" / "phase_coordination_lever_gate"
RNG_SEED = 20260609

ZETA_REQ = 0.05
COSTS = {"w": 1.0, "D": 2.0, "M": 8.0}
EPS_GRID = [0.05, 0.10, 0.20]
ETA_SYN = 1e-12
EQUAL_COST_J = 1.0
TOPK_COMBINED = 6
N_SYNTHETIC_PER_SIZE = 12
SYNTHETIC_SIZES = [5, 6, 7, 8, 9]
FLOW_MARGIN = 1.25
FLOW_ABS_MIN = 0.15
RESOLVENT_FREQS = np.linspace(0.05, 20.0, 80)
BAD_MARGIN = 1e6


@dataclass(frozen=True)
class NetworkCase:
    name: str
    family: str
    n: int
    edges: tuple[tuple[int, int], ...]
    w: np.ndarray
    M: np.ndarray
    D: np.ndarray
    P: np.ndarray
    fmax: np.ndarray


@dataclass(frozen=True)
class Action:
    lever: str
    target: int
    eps: float


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


def laplacian(n: int, edges: Iterable[tuple[int, int]], w: np.ndarray) -> np.ndarray:
    L = np.zeros((n, n), dtype=float)
    for k, (i, j) in enumerate(edges):
        wk = float(w[k])
        L[i, i] += wk
        L[j, j] += wk
        L[i, j] -= wk
        L[j, i] -= wk
    return L


def connected(n: int, edges: tuple[tuple[int, int], ...], w: np.ndarray) -> bool:
    vals = np.linalg.eigvalsh(laplacian(n, edges, w))
    return bool(vals[1] > 1e-9) if n > 1 else True


def dc_theta(n: int, edges: tuple[tuple[int, int], ...], w: np.ndarray, P: np.ndarray) -> np.ndarray:
    L = laplacian(n, edges, w)
    # Remove the slack/reference coordinate and solve L theta = P.
    keep = np.arange(1, n)
    theta = np.zeros(n)
    theta[keep] = np.linalg.solve(L[np.ix_(keep, keep)], P[keep])
    return theta


def dc_flows(n: int, edges: tuple[tuple[int, int], ...], w: np.ndarray, P: np.ndarray) -> np.ndarray:
    theta = dc_theta(n, edges, w, P)
    return np.array([w[k] * (theta[i] - theta[j]) for k, (i, j) in enumerate(edges)], dtype=float)


def flow_feasible(case: NetworkCase, w: np.ndarray) -> bool:
    if not connected(case.n, case.edges, w):
        return False
    try:
        flows = dc_flows(case.n, case.edges, w, case.P)
    except np.linalg.LinAlgError:
        return False
    return bool(np.all(np.abs(flows) <= case.fmax + 1e-10))


def state_matrix(case: NetworkCase, w: np.ndarray, D: np.ndarray, M: np.ndarray) -> np.ndarray:
    n = case.n
    L = laplacian(n, case.edges, w)
    invM = np.diag(1.0 / M)
    return np.block(
        [
            [np.zeros((n, n)), np.eye(n)],
            [-invM @ L, -invM @ np.diag(D)],
        ]
    )


def zeta_min_from_A(A: np.ndarray) -> float:
    poles = eigvals(A)
    zetas = []
    for s in poles:
        if abs(s.imag) < 1e-7 or abs(s) < 1e-10:
            continue
        zetas.append(float(-s.real / abs(s)))
    if not zetas:
        return float("inf")
    return min(zetas)


def metric_zeta(case: NetworkCase, w: np.ndarray, D: np.ndarray, M: np.ndarray) -> float:
    A = state_matrix(case, w, D, M)
    if np.max(eigvals(A).real) > 1e-7:
        return BAD_MARGIN
    zeta = max(zeta_min_from_A(A), 1e-12)
    return float(math.log(ZETA_REQ / zeta))


def metric_resolvent(case: NetworkCase, w: np.ndarray, D: np.ndarray, M: np.ndarray) -> float:
    A = state_matrix(case, w, D, M)
    if np.max(eigvals(A).real) > 1e-7:
        return BAD_MARGIN
    n2 = A.shape[0]
    n = case.n
    B = np.vstack([np.zeros((n, n)), np.diag(1.0 / M)])
    C = np.hstack([np.zeros((n, n)), np.eye(n)])
    eye = np.eye(n2)
    g = 0.0
    for om in RESOLVENT_FREQS:
        H = C @ np.linalg.solve(1j * om * eye - A, B)
        g = max(g, float(svdvals(H)[0]))
    return float(math.log(max(g, 1e-12)))


def metric(case: NetworkCase, w: np.ndarray, D: np.ndarray, M: np.ndarray, metric_name: str) -> float:
    if metric_name == "S_zeta":
        return metric_zeta(case, w, D, M)
    if metric_name == "S_g":
        return metric_resolvent(case, w, D, M)
    raise ValueError(metric_name)


def apply_actions(
    case: NetworkCase,
    actions: Iterable[Action],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, bool]:
    w = case.w.copy()
    D = case.D.copy()
    M = case.M.copy()
    for action in actions:
        if action.lever == "w":
            w[action.target] *= 1.0 + action.eps
        elif action.lever == "D":
            D[action.target] *= 1.0 + action.eps
        elif action.lever == "M":
            M[action.target] *= 1.0 + action.eps
        else:
            raise ValueError(action.lever)
    feasible = flow_feasible(case, w)
    return w, D, M, feasible


def candidate_actions(case: NetworkCase, lever: str, eps: float) -> list[Action]:
    if lever == "w":
        return [Action("w", k, eps) for k in range(len(case.edges))]
    if lever in {"D", "M"}:
        return [Action(lever, i, eps) for i in range(case.n)]
    raise ValueError(lever)


def evaluate_action(
    case: NetworkCase,
    actions: Iterable[Action],
    metric_name: str,
    base_score: float,
) -> dict[str, Any]:
    action_tuple = tuple(actions)
    w, D, M, feasible = apply_actions(case, action_tuple)
    if not feasible:
        return {
            "actions": action_tuple,
            "score": BAD_MARGIN,
            "delta": BAD_MARGIN,
            "feasible": False,
        }
    score = metric(case, w, D, M, metric_name)
    return {
        "actions": action_tuple,
        "score": score,
        "delta": float(score - base_score),
        "feasible": bool(score < BAD_MARGIN / 2),
    }


def best_single(
    case: NetworkCase,
    lever: str,
    eps: float,
    metric_name: str,
    base_score: float,
) -> dict[str, Any]:
    rows = [
        evaluate_action(case, [action], metric_name, base_score)
        for action in candidate_actions(case, lever, eps)
    ]
    feasible = [r for r in rows if r["feasible"]]
    return min(feasible, key=lambda r: r["delta"]) if feasible else {"feasible": False, "delta": BAD_MARGIN}


def ranked_singles(
    case: NetworkCase,
    lever: str,
    eps: float,
    metric_name: str,
    base_score: float,
    topk: int,
) -> list[dict[str, Any]]:
    rows = [
        evaluate_action(case, [action], metric_name, base_score)
        for action in candidate_actions(case, lever, eps)
    ]
    feasible = [r for r in rows if r["feasible"]]
    return sorted(feasible, key=lambda r: r["delta"])[:topk]


def action_label(actions: Iterable[Action]) -> str:
    return "+".join(f"{a.lever}{a.target}@{a.eps:.4g}" for a in actions)


def best_combined_screened(
    case: NetworkCase,
    lever_a: str,
    lever_b: str,
    eps_a: float,
    eps_b: float,
    metric_name: str,
    base_score: float,
) -> dict[str, Any]:
    cand_a = ranked_singles(case, lever_a, eps_a, metric_name, base_score, TOPK_COMBINED)
    cand_b = ranked_singles(case, lever_b, eps_b, metric_name, base_score, TOPK_COMBINED)
    rows = []
    for ra, rb in product(cand_a, cand_b):
        rows.append(evaluate_action(case, list(ra["actions"]) + list(rb["actions"]), metric_name, base_score))
    feasible = [r for r in rows if r["feasible"]]
    return min(feasible, key=lambda r: r["delta"]) if feasible else {"feasible": False, "delta": BAD_MARGIN}


def regular_ring_edges(n: int, k: int) -> set[tuple[int, int]]:
    edges: set[tuple[int, int]] = set()
    half = max(1, k // 2)
    for i in range(n):
        for d in range(1, half + 1):
            j = (i + d) % n
            edges.add((min(i, j), max(i, j)))
    return edges


def watts_strogatz_edges(n: int, k: int, beta: float, rng: np.random.Generator) -> tuple[tuple[int, int], ...]:
    edges = regular_ring_edges(n, k)
    edge_list = list(edges)
    for edge in edge_list:
        if rng.random() >= beta:
            continue
        i, j = edge
        edges.discard(edge)
        choices = [v for v in range(n) if v != i and (min(i, v), max(i, v)) not in edges]
        if not choices:
            edges.add(edge)
            continue
        new_j = int(rng.choice(choices))
        edges.add((min(i, new_j), max(i, new_j)))
    return tuple(sorted(edges))


def make_case(name: str, family: str, n: int, edges: tuple[tuple[int, int], ...], rng: np.random.Generator) -> NetworkCase:
    if not edges:
        raise ValueError("empty edge set")
    w = rng.uniform(0.5, 1.5, len(edges))
    M = rng.uniform(0.35, 1.2, n)
    sg_count = max(1, n // 4)
    sg_idx = rng.choice(n, size=sg_count, replace=False)
    M[sg_idx] *= rng.uniform(6.0, 12.0, sg_count)
    D = rng.uniform(0.03, 0.12, n)
    P = rng.normal(0.0, 1.0, n)
    P -= np.mean(P)
    # Scale dispatch so the base DC flow is feasible with a nontrivial margin.
    flows0 = dc_flows(n, edges, w, P)
    scale = 0.55 / max(float(np.max(np.abs(flows0))), 1e-8)
    P *= scale
    flows = dc_flows(n, edges, w, P)
    fmax = np.maximum(FLOW_MARGIN * np.abs(flows), FLOW_ABS_MIN)
    return NetworkCase(name=name, family=family, n=n, edges=edges, w=w, M=M, D=D, P=P, fmax=fmax)


def synthetic_cases(rng: np.random.Generator) -> list[NetworkCase]:
    cases = []
    for n in SYNTHETIC_SIZES:
        k = min(4, n - 1)
        made = 0
        attempts = 0
        while made < N_SYNTHETIC_PER_SIZE and attempts < 500:
            attempts += 1
            edges = watts_strogatz_edges(n, k, beta=0.25, rng=rng)
            w_probe = np.ones(len(edges))
            if not connected(n, edges, w_probe):
                continue
            case = make_case(f"ws_n{n}_{made:02d}", "synthetic_ws", n, edges, rng)
            if metric_zeta(case, case.w, case.D, case.M) < BAD_MARGIN / 2:
                cases.append(case)
                made += 1
    return cases


def ieee_case_from_andes(case_name: str, path: Path, rng: np.random.Generator) -> NetworkCase:
    sheets = pd.read_excel(path, sheet_name=None)
    bus = sheets["Bus"]
    line = sheets["Line"]
    active_line = line[pd.to_numeric(line.get("u", 1), errors="coerce").fillna(1).astype(float) > 0]
    active_bus = bus[pd.to_numeric(bus.get("u", 1), errors="coerce").fillna(1).astype(float) > 0]
    bus_ids = [int(x) for x in active_bus["idx"].tolist()]
    bus_map = {bus_id: i for i, bus_id in enumerate(bus_ids)}
    edges = []
    weights = []
    for _, row in active_line.iterrows():
        b1, b2 = int(row["bus1"]), int(row["bus2"])
        if b1 not in bus_map or b2 not in bus_map:
            continue
        x = abs(float(row.get("x", 0.0)))
        if x <= 1e-8:
            continue
        i, j = bus_map[b1], bus_map[b2]
        edges.append((min(i, j), max(i, j)))
        weights.append(1.0 / x)
    # Merge parallel branches.
    merged: dict[tuple[int, int], float] = {}
    for edge, weight in zip(edges, weights):
        merged[edge] = merged.get(edge, 0.0) + weight
    final_edges = tuple(sorted(merged.keys()))
    case = make_case(case_name, "ieee_topology", len(bus_ids), final_edges, rng)
    # Preserve relative IEEE susceptance weights but normalize their scale.
    raw_w = np.array([merged[e] for e in final_edges], dtype=float)
    raw_w /= np.median(raw_w)
    return NetworkCase(
        name=case.name,
        family=case.family,
        n=case.n,
        edges=case.edges,
        w=raw_w,
        M=case.M,
        D=case.D,
        P=case.P,
        fmax=np.maximum(FLOW_MARGIN * np.abs(dc_flows(case.n, case.edges, raw_w, case.P)), FLOW_ABS_MIN),
    )


def ieee_cases(rng: np.random.Generator) -> list[NetworkCase]:
    paths = {
        "wscc9_topology": Path(r"C:\Users\walla\anaconda3\Lib\site-packages\andes\cases\wscc9\wscc9.xlsx"),
        "ieee14_topology": Path(r"C:\Users\walla\anaconda3\Lib\site-packages\andes\cases\ieee14\ieee14_full.xlsx"),
    }
    cases = []
    for name, path in paths.items():
        if path.exists():
            cases.append(ieee_case_from_andes(name, path, rng))
    return cases


def run_case(case: NetworkCase) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    rows_syn: list[dict[str, Any]] = []
    rows_cv: list[dict[str, Any]] = []
    base = {
        "S_zeta": metric(case, case.w, case.D, case.M, "S_zeta"),
        "S_g": metric(case, case.w, case.D, case.M, "S_g"),
    }
    if any(v >= BAD_MARGIN / 2 for v in base.values()):
        return [], [], {"case": case.name, "base_ok": False}

    feasible_topology_actions = 0
    for eps in EPS_GRID:
        feasible_topology_actions += sum(
            1
            for action in candidate_actions(case, "w", eps)
            if evaluate_action(case, [action], "S_zeta", base["S_zeta"])["feasible"]
        )

    lever_pairs = [("w", "D"), ("w", "M"), ("D", "M")]
    for metric_name in ["S_zeta", "S_g"]:
        base_score = base[metric_name]
        for eps in EPS_GRID:
            best = {
                lever: best_single(case, lever, eps, metric_name, base_score)
                for lever in ["w", "D", "M"]
            }
            for la, lb in lever_pairs:
                if not best[la]["feasible"] or not best[lb]["feasible"]:
                    continue
                combined = evaluate_action(
                    case,
                    list(best[la]["actions"]) + list(best[lb]["actions"]),
                    metric_name,
                    base_score,
                )
                if not combined["feasible"]:
                    continue
                da = float(best[la]["delta"])
                db = float(best[lb]["delta"])
                dab = float(combined["delta"])
                syn = dab - da - db
                syn_norm = syn / (abs(da) + abs(db) + ETA_SYN)
                rows_syn.append(
                    {
                        "case": case.name,
                        "family": case.family,
                        "n": case.n,
                        "metric": metric_name,
                        "pair": f"{la}+{lb}",
                        "eps": eps,
                        "delta_a": da,
                        "delta_b": db,
                        "delta_ab": dab,
                        "syn": syn,
                        "syn_norm": syn_norm,
                        "action_a": action_label(best[la]["actions"]),
                        "action_b": action_label(best[lb]["actions"]),
                        "action_ab": action_label(combined["actions"]),
                    }
                )

        # Equal-cost value of coordination.
        single_eps = {lever: EQUAL_COST_J / COSTS[lever] for lever in ["w", "D", "M"]}
        half_eps = {lever: EQUAL_COST_J / (2.0 * COSTS[lever]) for lever in ["w", "D", "M"]}
        single = {
            lever: best_single(case, lever, single_eps[lever], metric_name, base_score)
            for lever in ["w", "D", "M"]
        }
        ranking = sorted(
            [(lever, single[lever]["delta"]) for lever in ["w", "D", "M"] if single[lever]["feasible"]],
            key=lambda item: item[1],
        )
        for la, lb in lever_pairs:
            if not single[la]["feasible"] or not single[lb]["feasible"]:
                continue
            combined = best_combined_screened(
                case,
                la,
                lb,
                half_eps[la],
                half_eps[lb],
                metric_name,
                base_score,
            )
            if not combined["feasible"]:
                continue
            best_single_delta = min(float(single[la]["delta"]), float(single[lb]["delta"]))
            cv = float(combined["delta"] - best_single_delta)
            rows_cv.append(
                {
                    "case": case.name,
                    "family": case.family,
                    "n": case.n,
                    "metric": metric_name,
                    "pair": f"{la}+{lb}",
                    "delta_single_a": float(single[la]["delta"]),
                    "delta_single_b": float(single[lb]["delta"]),
                    "delta_best_single": best_single_delta,
                    "delta_combined": float(combined["delta"]),
                    "cv": cv,
                    "cv_negative": bool(cv < 0.0),
                    "single_a": action_label(single[la]["actions"]),
                    "single_b": action_label(single[lb]["actions"]),
                    "combined": action_label(combined["actions"]),
                    "single_ranking": ">".join(lever for lever, _ in ranking),
                }
            )
    summary = {
        "case": case.name,
        "family": case.family,
        "n": case.n,
        "n_edges": len(case.edges),
        "base_S_zeta": base["S_zeta"],
        "base_S_g": base["S_g"],
        "base_flow_max_abs": float(np.max(np.abs(dc_flows(case.n, case.edges, case.w, case.P)))),
        "topology_actions_feasible_over_eps_grid": feasible_topology_actions,
        "topology_actions_tested_over_eps_grid": len(case.edges) * len(EPS_GRID),
        "base_ok": True,
    }
    return rows_syn, rows_cv, summary


def aggregate(rows_syn: list[dict[str, Any]], rows_cv: list[dict[str, Any]], case_summaries: list[dict[str, Any]]) -> dict[str, Any]:
    syn_df = pd.DataFrame(rows_syn)
    cv_df = pd.DataFrame(rows_cv)
    agg_rows = []
    for metric_name in ["S_zeta", "S_g"]:
        for pair in ["w+D", "w+M", "D+M"]:
            syn_sub = syn_df[(syn_df["metric"] == metric_name) & (syn_df["pair"] == pair)]
            cv_sub = cv_df[(cv_df["metric"] == metric_name) & (cv_df["pair"] == pair)]
            for eps in EPS_GRID:
                syn_eps = syn_sub[syn_sub["eps"] == eps]
                agg_rows.append(
                    {
                        "metric": metric_name,
                        "pair": pair,
                        "eps": eps,
                        "n_syn": int(len(syn_eps)),
                        "median_syn_norm": float(syn_eps["syn_norm"].median()) if len(syn_eps) else float("nan"),
                        "pct_syn_norm_lt_minus_0p10": float((syn_eps["syn_norm"] < -0.10).mean()) if len(syn_eps) else float("nan"),
                        "n_cv": int(len(cv_sub)),
                        "median_cv": float(cv_sub["cv"].median()) if len(cv_sub) else float("nan"),
                        "pct_cv_lt_0": float((cv_sub["cv"] < 0.0).mean()) if len(cv_sub) else float("nan"),
                    }
                )

    disagreement = []
    for case in sorted(set(cv_df["case"].tolist())) if len(cv_df) else []:
        for pair in ["w+D", "w+M", "D+M"]:
            z = cv_df[(cv_df["case"] == case) & (cv_df["pair"] == pair) & (cv_df["metric"] == "S_zeta")]
            g = cv_df[(cv_df["case"] == case) & (cv_df["pair"] == pair) & (cv_df["metric"] == "S_g")]
            if len(z) and len(g):
                disagreement.append(
                    {
                        "case": case,
                        "pair": pair,
                        "zeta_cv_negative": bool(z.iloc[0]["cv"] < 0),
                        "resolvent_cv_negative": bool(g.iloc[0]["cv"] < 0),
                        "zeta_single_ranking": str(z.iloc[0]["single_ranking"]),
                        "resolvent_single_ranking": str(g.iloc[0]["single_ranking"]),
                        "same_cv_sign": bool((z.iloc[0]["cv"] < 0) == (g.iloc[0]["cv"] < 0)),
                        "same_single_ranking": bool(str(z.iloc[0]["single_ranking"]) == str(g.iloc[0]["single_ranking"])),
                    }
                )
    dis_df = pd.DataFrame(disagreement)
    feasible_ratios = [
        s["topology_actions_feasible_over_eps_grid"] / max(s["topology_actions_tested_over_eps_grid"], 1)
        for s in case_summaries
        if s.get("base_ok")
    ]
    verdict = verdict_from_agg(agg_rows)
    return {
        "aggregate_rows": agg_rows,
        "n_cases": len([s for s in case_summaries if s.get("base_ok")]),
        "n_synthetic_cases": len([s for s in case_summaries if s.get("base_ok") and s.get("family") == "synthetic_ws"]),
        "n_ieee_topology_cases": len([s for s in case_summaries if s.get("base_ok") and s.get("family") == "ieee_topology"]),
        "median_topology_feasible_ratio": float(np.median(feasible_ratios)) if feasible_ratios else float("nan"),
        "metric_disagreement": {
            "n_pair_cases": int(len(dis_df)) if len(dis_df) else 0,
            "pct_same_cv_sign": float(dis_df["same_cv_sign"].mean()) if len(dis_df) else float("nan"),
            "pct_same_single_ranking": float(dis_df["same_single_ranking"].mean()) if len(dis_df) else float("nan"),
        },
        "verdict": verdict,
    }


def verdict_from_agg(agg_rows: list[dict[str, Any]]) -> dict[str, Any]:
    cv_support = [
        r
        for r in agg_rows
        if r["eps"] == EPS_GRID[0] and not math.isnan(r["pct_cv_lt_0"]) and r["pct_cv_lt_0"] > 0.50
    ]
    syn_support = [
        r
        for r in agg_rows
        if not math.isnan(r["pct_syn_norm_lt_minus_0p10"]) and r["pct_syn_norm_lt_minus_0p10"] > 0.50
    ]
    if cv_support:
        status = "COORDINATION_VALUE_OBSERVED"
        implication = "Equal-cost mixed packages beat the best single lever in a majority of cases for at least one metric/pair."
    else:
        status = "NO_CONSISTENT_COORDINATION_VALUE"
        implication = "Do not claim coordinated retrofitting value. In this gate, the best single lever usually matches or beats mixed packages at equal cost."
    return {
        "status": status,
        "implication": implication,
        "n_cv_majority_rows": len(cv_support),
        "n_syn_majority_rows": len(syn_support),
    }


def write_report(status: dict[str, Any], agg: dict[str, Any]) -> None:
    lines = [
        "# Gate 1: Coordination Value of Planning Levers",
        "",
        "This gate tests whether coordinated packages of topology, damping, and inertia",
        "outperform the best single lever under equal normalized cost and DC flow",
        "feasibility.  The ground truth is the complete eigensolve of the tested",
        "second-order dynamic model.",
        "",
        f"- Commit: `{status['git_commit']}`",
        f"- Random seed: `{status['random_seed']}`",
        f"- Costs: `c_w={COSTS['w']}`, `c_D={COSTS['D']}`, `c_M={COSTS['M']}`",
        f"- Equal-cost budget: `{EQUAL_COST_J}`",
        f"- Combined-package search: top-{TOPK_COMBINED} screened candidates per lever",
        f"- Cases passing base filters: `{agg['n_cases']}` "
        f"({agg['n_synthetic_cases']} synthetic, {agg['n_ieee_topology_cases']} IEEE-topology)",
        f"- Median topology-action feasibility ratio: `{agg['median_topology_feasible_ratio']:.3f}`",
        "",
        "## Aggregate Table",
        "",
        "| metric | pair | eps | n | median SynNorm | % SynNorm < -0.10 | median CV | % CV < 0 |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for r in agg["aggregate_rows"]:
        lines.append(
            f"| {r['metric']} | {r['pair']} | {r['eps']:.2f} | {r['n_syn']} | "
            f"{r['median_syn_norm']:.4g} | {100*r['pct_syn_norm_lt_minus_0p10']:.1f}% | "
            f"{r['median_cv']:.4g} | {100*r['pct_cv_lt_0']:.1f}% |"
        )
    md = agg["metric_disagreement"]
    lines.extend(
        [
            "",
            "## Metric Agreement",
            "",
            f"- Pair-case comparisons: `{md['n_pair_cases']}`",
            f"- Same CV sign for modal and resolvent margins: `{100*md['pct_same_cv_sign']:.1f}%`",
            f"- Same single-lever ranking for modal and resolvent margins: `{100*md['pct_same_single_ranking']:.1f}%`",
            "",
            "## Verdict",
            "",
            f"- Status: **{agg['verdict']['status']}**",
            f"- Implication: {agg['verdict']['implication']}",
            f"- Rows with CV majority support: `{agg['verdict']['n_cv_majority_rows']}`",
            f"- Rows with nonlinear-synergy majority support: `{agg['verdict']['n_syn_majority_rows']}`",
            "",
            "## Scope",
            "",
            "- This is not a full ANDES IBR validation. IEEE cases are used as topology/flow templates.",
            "- The resolvent margin uses identity protection scalings because no protection matrices were provided.",
            "- The equal-cost package search is favorable but screened: it evaluates combinations among the top single-lever candidates.",
        ]
    )
    (OUT / "phase_coordination_lever_gate_report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(RNG_SEED)
    cases = synthetic_cases(rng) + ieee_cases(rng)

    all_syn: list[dict[str, Any]] = []
    all_cv: list[dict[str, Any]] = []
    summaries: list[dict[str, Any]] = []
    for idx, case in enumerate(cases, start=1):
        print(f"[gate1] {idx}/{len(cases)} {case.name} n={case.n} e={len(case.edges)}")
        rows_syn, rows_cv, summary = run_case(case)
        all_syn.extend(rows_syn)
        all_cv.extend(rows_cv)
        summaries.append(summary)

    agg = aggregate(all_syn, all_cv, summaries)
    status = {
        "gate": "coordination_lever_value",
        "overall_status": agg["verdict"]["status"],
        "git_commit": git_commit_id(),
        "random_seed": RNG_SEED,
        "costs": COSTS,
        "equal_cost_budget": EQUAL_COST_J,
        "eps_grid": EPS_GRID,
        "topk_combined": TOPK_COMBINED,
        "resolvent_freq_grid_rad_s": {
            "min": float(RESOLVENT_FREQS[0]),
            "max": float(RESOLVENT_FREQS[-1]),
            "count": int(len(RESOLVENT_FREQS)),
        },
        "version_manifest": version_manifest(),
        "aggregate": agg,
    }

    write_csv(OUT / "coordination_syn_rows.csv", all_syn)
    write_csv(OUT / "coordination_cv_rows.csv", all_cv)
    write_csv(OUT / "coordination_case_summaries.csv", summaries)
    write_csv(OUT / "coordination_aggregate_table.csv", agg["aggregate_rows"])
    write_json(OUT / "phase_coordination_lever_gate_status.json", status)
    write_report(status, agg)
    print(json.dumps(status, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
