"""CDW00: generate the frozen preregistration inputs (no model evaluation).

- holdout policy points H01..H24: Latin hypercube, seed 20260920, over u in [0,1]
  (g = u^2, as the TX4 P-suites), k in [0.5, 2.3], t in [0.5, 3.0], h in [0, 2];
- CDW holdout envelope draws: the TX4 PCV05 envelope definitions (EM-f, EM-u, EC,
  EMC) with a NEW seed 20260921 (+ envelope index), N = 40 per envelope; converter
  groups are drawn for all nine census candidates (TX4 drew them for the core only).
Writes results/prereg_inputs/{holdout_policies.json, cdw_envelope_draws.json}.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.stats import qmc

OUT = Path(__file__).resolve().parents[2] / "results" / "prereg_inputs"
SEED_POLICY = 20260920
SEED_ENV = 20260921
N_POLICY = 24
N_ENV = 40
SG_BUSES = (30, 31, 32, 33, 34, 35, 36, 37, 38, 39)
CONV_BUSES = (30, 31, 32, 33, 34, 35, 36, 37, 38)
MACHINE_GROUPS = {  # identical to PCV05_robustness.py (TX4, frozen)
    "M": (("m",), 0.80, 1.20),
    "XP": (("xd1", "xq1"), 0.90, 1.10),
    "KA": (("ka",), 0.85, 1.15),
    "TE": (("ta",), 0.85, 1.15),
    "PSS_K": (("pss_gain",), 0.80, 1.20),
    "PSS_T": (("pss_washout", "pss_wash_lag", "pss_lag"), 0.80, 1.20),
}
CONV_GROUPS = {
    "PLL": (("kp_pll", "ki_pll"), 0.80, 1.20),
    "OUTER": (("kp_p", "ki_p", "kp_q", "ki_q"), 0.80, 1.20),
    "CURRENT": (("kp_i", "ki_i"), 0.80, 1.20),
    "TAU_P": (("tau_p",), 0.80, 1.20),
    "XF": (("xf",), 0.90, 1.10),
}


def lhs(dim, seed, n):
    return qmc.LatinHypercube(d=dim, rng=np.random.default_rng(seed)).random(n)


def sc(u, lo, hi):
    return round(float(lo + u * (hi - lo)), 10)


def policies():
    rows = []
    for i, (u, k, t, h) in enumerate(lhs(4, SEED_POLICY, N_POLICY)):
        rows.append(
            {
                "id": f"H{i + 1:02d}",
                "g": round(float(u) ** 2, 10),
                "k": sc(k, 0.5, 2.3),
                "t": sc(t, 0.5, 3.0),
                "h": sc(h, 0.0, 2.0),
            }
        )
    return rows


def draws():
    mg, cg = list(MACHINE_GROUPS), list(CONV_GROUPS)
    nu, nc = len(SG_BUSES) * len(mg), len(CONV_BUSES) * len(cg)
    unit = lambda r: {  # noqa: E731
        str(b): {g: sc(r[i * len(mg) + j], *MACHINE_GROUPS[g][1:]) for j, g in enumerate(mg)}
        for i, b in enumerate(SG_BUSES)
    }
    conv = lambda r: {  # noqa: E731
        str(b): {g: sc(r[i * len(cg) + j], *CONV_GROUPS[g][1:]) for j, g in enumerate(cg)}
        for i, b in enumerate(CONV_BUSES)
    }
    out = []
    for i, r in enumerate(lhs(len(mg), SEED_ENV + 0, N_ENV)):
        out.append({"envelope": "EM-f", "draw": i, "fleet": {g: sc(r[j], *MACHINE_GROUPS[g][1:]) for j, g in enumerate(mg)}, "unit": {}, "conv": {}})
    for i, r in enumerate(lhs(nu, SEED_ENV + 1, N_ENV)):
        out.append({"envelope": "EM-u", "draw": i, "fleet": {}, "unit": unit(r), "conv": {}})
    for i, r in enumerate(lhs(nc, SEED_ENV + 2, N_ENV)):
        out.append({"envelope": "EC", "draw": i, "fleet": {}, "unit": {}, "conv": conv(r)})
    for i, r in enumerate(lhs(nu + nc, SEED_ENV + 3, N_ENV)):
        out.append({"envelope": "EMC", "draw": i, "fleet": {}, "unit": unit(r[:nu]), "conv": conv(r[nu:])})
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for name, payload in (("holdout_policies.json", policies()), ("cdw_envelope_draws.json", draws())):
        p = OUT / name
        p.write_text(json.dumps(payload, indent=1), encoding="utf-8")
        print(name, hashlib.sha256(p.read_bytes()).hexdigest())


if __name__ == "__main__":
    main()
