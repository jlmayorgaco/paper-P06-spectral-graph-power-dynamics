"""Recompute the registered local/collective eta diagnostic on a small exact stratum."""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results"
sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_tx4_robustness as exact  # noqa: E402
from audit_existing_robustness import PARAMS  # noqa: E402
from ibr_cycles.models.ieee39_case import ReplacementPlan, solve_case  # noqa: E402
from ibr_cycles.models.ieee39_devices import ConverterParameters  # noqa: E402
from ibr_cycles.models.port_admittance import build_action_space  # noqa: E402


H4 = (30, 33, 35, 37)


def make_case(params: dict[str, float], portfolio: tuple[int, ...]):
    network = exact.network_for_epsilon(params["epsilon"])
    converter = ConverterParameters(voltage_control=True, voltage_gain=params["g"], voltage_leak=0.05)
    converters = None
    if 37 in portfolio:
        converters = {37: ConverterParameters(voltage_control=True, voltage_gain=params["g"] * params["h"], voltage_leak=0.05)}
    plan = ReplacementPlan.of(
        {bus: 1.0 for bus in portfolio},
        q_policy="matched",
        device="gfl",
        machine_damping=params["damping"],
    )
    return solve_case(
        plan,
        network=network,
        converter=converter,
        converters=converters,
        machine_scaling={"ka": params["k"], "ta": params["t"]},
        machine_services={"inertia": params["inertia"]},
        tol=1e-9,
    )


def h4_mode(case) -> complex:
    eig = np.linalg.eigvals(case.system.A)
    transverse = eig[np.abs(eig) >= exact.GAUGE_TOL]
    frequency = np.abs(transverse.imag) / (2.0 * np.pi)
    em = transverse[(frequency >= exact.MODE_LOW_HZ) & (frequency <= exact.MODE_HIGH_HZ)]
    if not em.size:
        raise RuntimeError("No H4 EM-band mode available")
    return complex(em[int(np.argmax(em.real))])


def one_row(campaign: str, source_row: pd.Series) -> dict[str, object]:
    params = {name: float(source_row[name]) for name in PARAMS}
    base_case = make_case(params, ())
    h4_case = make_case(params, H4)
    s_h = h4_mode(h4_case)
    values: list[float] = []
    local_values: list[float] = []
    for omitted in H4:
        subset = tuple(bus for bus in H4 if bus != omitted)
        proper_case = make_case(params, subset)
        space = build_action_space(base_case, proper_case, subset)
        split = space.split(s_h)
        values.append(float(split["sigma_min_i_plus_q"]))
        m = space.m(s_h)
        total = np.eye(space.dimension, dtype=np.complex128) + m
        for index in range(space.order):
            block = total[2 * index : 2 * index + 2, 2 * index : 2 * index + 2]
            local_values.append(float(np.linalg.svd(block, compute_uv=False)[-1]))
    return {
        "campaign": campaign,
        "condition_id": str(source_row.condition_id),
        **params,
        "s_H_real": float(s_h.real),
        "s_H_imag": float(s_h.imag),
        "s_H_frequency_hz": float(abs(s_h.imag) / (2.0 * np.pi)),
        "eta_corrected": float(max(0.0, min(values))),
        "collective_sigma_min_min": float(max(0.0, min(values))),
        "local_sigma_min_min": float(max(0.0, min(local_values))),
        "n_proper_portfolios": len(values),
        "status": "OK",
        "eta_definition": "min_i sigma_min(I+Q_{H4\\i}(s_H))",
        "source": "EXACT_DAE_STRATIFIED_RECONSTRUCTION",
    }


def main() -> None:
    start = time.perf_counter()
    rows: list[dict[str, object]] = []
    failures: list[dict[str, object]] = []
    for campaign in ("QMC", "MC"):
        truth = pd.read_csv(RESULTS / f"TX4_{campaign}_EXACT_BLOCKER_TRUTH.csv")
        selected = truth[truth.EXACT_H4].sort_values("condition_id").head(64)
        for _, source_row in selected.iterrows():
            try:
                rows.append(one_row(campaign, source_row))
            except Exception as exc:
                failures.append({"campaign": campaign, "condition_id": source_row.condition_id, "status": "FAIL", "error": f"{type(exc).__name__}: {exc}"})
    result = pd.DataFrame(rows)
    result.to_csv(RESULTS / "TX4_ETA_LOCAL_COLLECTIVE_EXACT.csv", index=False)
    pd.DataFrame(failures).to_csv(RESULTS / "TX4_ETA_LOCAL_COLLECTIVE_EXACT_FAILURES.csv", index=False)
    summary = {
        "n_requested_max": 128,
        "n_selected": int(len(rows) + len(failures)),
        "n_computed": int(len(rows)),
        "n_failures": int(len(failures)),
        "eta_nonnegative": bool((result.eta_corrected >= 0).all()) if len(result) else False,
        "runtime_s": time.perf_counter() - start,
    }
    (RESULTS / "TX4_ETA_LOCAL_COLLECTIVE_EXACT_SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
