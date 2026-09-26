# ruff: noqa: E501  -- table labels kept on one line
"""CC04 - Phase 4: does the connected cumulant add information over existing diagnostics?

Preregistered (configs/ias2026/connected_cumulants_prereg_v1.yaml, de09db3b). Data: the
frozen FC04 Pg-matched set at the E12 census policy (10 failing, 25 stable size-4
portfolios) plus the ten CC03 boundaries as a descriptive group. Orientations and the
decision rule are fixed in the preregistration; AUC is threshold-free.
"""

from __future__ import annotations

import json
import sys
import time
from itertools import permutations
from multiprocessing import Pool

import numpy as np
import pandas as pd
from _cc import (
    E12_CENSUS,
    E13_PREDICTORS,
    FC04_TABLE,
    ROOT,
    Realization,
    band_critical,
    out_dir,
    write_json,
)

from _overnight import pin_blas_threads
from ibr_cycles.cycles.connected import (
    boolean_mobius,
    characteristic_values,
    connected_share,
    cumulants,
    subsets,
)

OUT = out_dir("CC04")
SEED = 20260913
PERMUTATIONS = 10000
ORIENT = {  # +1: higher = unstable; -1: lower = unstable (preregistered)
    "abs_chi_axis": +1,
    "abs_mu_axis": +1,
    "cycle_score_axis": +1,
    "closure_distance_axis": -1,
    "nu_at_mode": +1,
    "pg_mw": +1,
    "sn_mva": +1,
    "gscr": -1,
    "min_scr": -1,
}
OLD = ("abs_mu_axis", "cycle_score_axis", "pg_mw", "sn_mva", "gscr", "min_scr")


def portfolio_pg(members) -> float:
    payload = json.loads((ROOT / "configs/ias2026/ieee39_network.json").read_text())
    p0 = {int(r["bus"]): float(r["p0"]) for r in payload["pv"]}
    return 100.0 * sum(p0[b] for b in members)


def evaluate(task):
    pin_blas_threads()
    role, name = task
    members = tuple(sorted(int(b) for b in name.split("+")))
    r = Realization(members, None)
    lam = band_critical(r.vertex_spectrum(members))
    out = {
        "role": role,
        "portfolio": name,
        "lambda_re": lam.real,
        "lambda_hz": lam.imag / (2 * np.pi),
    }
    blocks = tuple(range(len(members)))
    for tag, s in (("axis", complex(0.0, lam.imag)), ("mode", lam)):
        q, d, k, _ = r.q_matrix(s)
        f = characteristic_values(q, [2] * len(members))
        chi = cumulants(f, blocks)
        full = frozenset(blocks)
        if tag == "axis":
            mu = boolean_mobius(f, blocks)
            out["abs_chi_axis"] = abs(chi[full])
            out["abs_mu_axis"] = abs(mu[full])
            out["closure_distance_axis"] = float(
                np.min(np.abs(np.linalg.eigvals(q) + 1.0))
            )
            m = d @ k
            best = 0.0
            for chosen in subsets(blocks):
                if len(chosen) < 2:
                    continue
                head, rest = chosen[0], chosen[1:]
                for tail in permutations(rest):
                    cyc = (head, *tail)
                    prod = np.eye(2, dtype=complex)
                    for i, a in enumerate(cyc):
                        b = cyc[(i + 1) % len(cyc)]
                        prod = prod @ m[2 * a : 2 * a + 2, 2 * b : 2 * b + 2]
                    best = max(best, float(np.abs(np.linalg.eigvals(prod)).max()))
            out["cycle_score_axis"] = best
        else:
            out["abs_F_at_mode"] = abs(f[full])
            out["nu_at_mode"] = connected_share(chi, blocks)
            out["abs_chi_at_mode"] = abs(chi[full])
    out["pg_mw_network"] = portfolio_pg(members)
    return out


def auc(scores, labels):
    """Mann-Whitney AUC: P(score_unstable > score_stable), ties counted 1/2."""

    pos, neg = scores[labels], scores[~labels]
    greater = (pos[:, None] > neg[None, :]).sum()
    ties = (pos[:, None] == neg[None, :]).sum()
    return float((greater + 0.5 * ties) / (len(pos) * len(neg)))


def main(argv) -> int:
    started = time.time()
    fc04 = pd.read_csv(FC04_TABLE)
    tasks = [(r.role, r.portfolio) for r in fc04.itertuples()]
    with Pool(16) as pool:
        rows = pool.map(evaluate, tasks, chunksize=1)
    df = pd.DataFrame(rows)
    census = pd.read_csv(E12_CENSUS).set_index("members")
    e13 = pd.read_csv(E13_PREDICTORS).set_index("members")
    df["pg_mw"] = (
        fc04.set_index("portfolio").loc[df.portfolio, "replaced_pg_mw"].to_numpy()
    )
    df["sn_mva"] = census.loc[df.portfolio, "replaced_mw"].to_numpy()
    df["gscr"] = e13.loc[df.portfolio, "gscr"].to_numpy()
    df["min_scr"] = e13.loc[df.portfolio, "min_scr"].to_numpy()
    df["pg_check_abs_diff"] = (df.pg_mw - df.pg_mw_network).abs()
    df.to_csv(OUT / "CC04_portfolios.csv", index=False)

    labels = (df.role == "failing").to_numpy()
    rng = np.random.default_rng(SEED)
    perms = [rng.permutation(labels) for _ in range(PERMUTATIONS)]
    table = []
    for name, sign in ORIENT.items():
        x = sign * df[name].to_numpy(dtype=float)
        a = auc(x, labels)
        null = np.array([auc(x, p) for p in perms])
        table.append(
            {
                "diagnostic": name,
                "orientation": "higher = unstable" if sign > 0 else "lower = unstable",
                "auc": a,
                "permutation_p": float((np.sum(null >= a) + 1) / (PERMUTATIONS + 1)),
                "failing_median": float(df.loc[labels, name].median()),
                "stable_median": float(df.loc[~labels, name].median()),
            }
        )
    auc_table = pd.DataFrame(table)
    auc_table.to_csv(OUT / "CC04_auc.csv", index=False)

    a = auc_table.set_index("diagnostic")
    best_old = float(a.loc[list(OLD), "auc"].max())
    promote = bool(
        a.loc["abs_chi_axis", "auc"] >= 0.80
        and a.loc["abs_chi_axis", "auc"] >= best_old + 0.05
        and a.loc["abs_chi_axis", "permutation_p"] <= 0.01
    )
    cc03 = json.loads((OUT.parent / "CC03" / "CC03_boundaries.json").read_text())[
        "events"
    ]
    big = [e for e in cc03 if e["size"] >= 3]
    rem_major = sum(e["remainder_share"] >= 0.25 for e in big) > len(big) / 2
    clarify = bool(rem_major or a.loc["nu_at_mode", "auc"] >= 0.80)
    verdict = "PROMOTE" if promote else ("CLARIFY" if clarify else "NEGATIVE")
    boundaries = pd.DataFrame(
        [
            {
                "event": e["event"],
                "H": e["H"],
                "abs_chi_H": e["abs_chi_H"],
                "abs_mu_H": e["abs_mu_H"],
                "cycle_score_M": e["cycle_score_M"],
                "closure_distance": e["closure_distance"],
                "nu_H": e["nu_H"],
                "remainder_share": e["remainder_share"],
                "sn_mva": float(census.loc[e["H"], "replaced_mw"]),
                "pg_mw": portfolio_pg([int(b) for b in e["H"].split("+")]),
                "gscr": float(e13.loc[e["H"], "gscr"])
                if e["H"] in e13.index
                else np.nan,
                "min_scr": float(e13.loc[e["H"], "min_scr"])
                if e["H"] in e13.index
                else np.nan,
            }
            for e in cc03
        ]
    )
    boundaries.to_csv(OUT / "CC04_boundary_group.csv", index=False)
    summary = {
        "n_failing": int(labels.sum()),
        "n_stable": int((~labels).sum()),
        "best_old_auc": best_old,
        "auc_abs_chi": float(a.loc["abs_chi_axis", "auc"]),
        "auc_closure_distance": float(a.loc["closure_distance_axis", "auc"]),
        "auc_nu": float(a.loc["nu_at_mode", "auc"]),
        "remainder_share_majority_ge_0.25": bool(rem_major),
        "remainder_shares": {e["event"]: e["remainder_share"] for e in big},
        "verdict": verdict,
        "pg_check_max_abs_diff": float(df.pg_check_abs_diff.max()),
        "elapsed_s": round(time.time() - started, 1),
    }
    write_json(OUT / "CC04_summary.json", summary)
    pd.set_option("display.width", 250)
    print(auc_table.to_string(index=False))
    print(boundaries.to_string(index=False))
    print(json.dumps(summary, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
