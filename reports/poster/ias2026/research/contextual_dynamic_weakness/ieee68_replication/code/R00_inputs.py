# ruff: noqa: E501
"""R0/R3: frozen design inputs of CDW68, generated from fixed seeds BEFORE any portfolio stability result.

Writes inputs/cdw68_design_v1.json and inputs/cdw68_inputs_manifest.json (sha256 of every frozen input).
"""

from __future__ import annotations

import _r68 as R  # noqa: I001

import hashlib
import json

import numpy as np

SEED_POLICY = 20260913
SEED_DRAW_A5, SEED_DRAW_A10 = 20260914, 20260915
SEED_DRAW_B5, SEED_DRAW_B10 = 20260916, 20260917
SEED_R11 = 20260918
N_CAND = 2000
U_RANGE, K_RANGE = (0.0, 1.0), (0.5, 2.0)


def lhs(rng, n, d):
    out = np.empty((n, d))
    for j in range(d):
        out[:, j] = (rng.permutation(n) + rng.random(n)) / n
    return out


def maximin(rng, n, d, n_cand):
    best, best_d = None, -1.0
    for _ in range(n_cand):
        x = lhs(rng, n, d)
        diff = x[:, None, :] - x[None, :, :]
        dist = np.sqrt((diff ** 2).sum(-1))
        dmin = dist[np.triu_indices(n, 1)].min()
        if dmin > best_d:
            best, best_d = x, dmin
    return best, float(best_d)


def draws(seed, n, names, half_width):
    rng = np.random.default_rng(seed)
    x = lhs(rng, n, len(names))
    return [{nm: float(1.0 + half_width * (2.0 * x[i, j] - 1.0)) for j, nm in enumerate(names)} for i in range(n)]


def main():
    rng = np.random.default_rng(SEED_POLICY)
    x, dmin = maximin(rng, 16, 2, N_CAND)
    pol = []
    for i in range(16):
        u = U_RANGE[0] + x[i, 0] * (U_RANGE[1] - U_RANGE[0])
        k = K_RANGE[0] + x[i, 1] * (K_RANGE[1] - K_RANGE[0])
        pol.append({"id": f"P68_{i + 1:02d}", "u": float(u), "g": float(u * u), "k": float(k),
                    "split": "discovery" if i < 4 else "holdout"})
    mach = ["h", "xd1q1", "ka", "pss_k"]
    conv = ["pll", "current", "outer"]
    dA5 = draws(SEED_DRAW_A5, 20, mach + conv, 0.05)
    dA10 = draws(SEED_DRAW_A10, 20, mach + conv, 0.10)
    dB5 = draws(SEED_DRAW_B5, 10, mach, 0.05)
    dB10 = draws(SEED_DRAW_B10, 10, mach, 0.10)

    def pack(lst, env, model):
        out = []
        for i, d in enumerate(lst):
            out.append({"id": f"{model}_{env}_{i:02d}", "envelope": env, "machine": {k: d[k] for k in mach},
                        "converter": {k: d[k] for k in conv if k in d} if model == "A" else {}})
        return out

    elig = [b["e"] for b in R.branches() if b["u"] != 0.0]
    r11 = sorted(np.random.default_rng(SEED_R11).choice(elig, size=20, replace=False).tolist())
    design = {
        "version": "cdw68_design_v1",
        "candidates_V68": list(R.V68),
        "candidate_rule": "the six physical plants (G1-G12) with the largest documented scheduled active power: G12 13.50, G11 10.00, G9 8.00, G6 7.00, G3 6.50, G4 6.32 pu; area equivalents G13-G16 (incl. slack G16) excluded; extends the repository's Gate 3 rule (G9, G6, G3, G4); NETS and NYPS both represented; no tie occurs",
        "policies_A": pol,
        "policy_design": {"domain": {"u": U_RANGE, "k": K_RANGE, "g": "u^2"}, "seed": SEED_POLICY, "n_candidates": N_CAND, "min_pairwise_distance_unit_square": dmin,
                          "source": "the validated IEEE-68 policy window of Gate 3 (configs/ieee68/G3_preregistration.yaml): g in [0,1] (square-root spacing, here g = u^2), k in [0.5, 2]"},
        "policies_B": [{"id": p["id"].replace("P68", "B68"), "k": p["k"], "from": p["id"], "split": p["split"]} for p in pol],
        "draws_A": pack(dA5, "E05", "A") + pack(dA10, "E10", "A"),
        "draws_B": pack(dB5, "E05", "B") + pack(dB10, "E10", "B"),
        "draw_seeds": {"A_E05": SEED_DRAW_A5, "A_E10": SEED_DRAW_A10, "B_E05": SEED_DRAW_B5, "B_E10": SEED_DRAW_B10},
        "draws_ranking_subset_A": [f"A_E05_{i:02d}" for i in range(5)] + [f"A_E10_{i:02d}" for i in range(5)],
        "eligible_branches": elig,
        "r11_branches": r11,
        "r11_seed": SEED_R11,
        "r11_policies": [p["id"] for p in pol if p["split"] == "discovery"],
        "ablation_policies": [p["id"] for p in pol[:6]],
        "ranking_target": list(R.V68),
        "gammas": [1.10, 1.25, 1.50],
        "tds": {"disturbance": "50 MVAr (0.5 pu on 100 MVA) shunt reactor at bus 3, on at t = 1 s, off at t = 11 s (Canizares et al. 2017, IEEE TPWRS 32(1), Benchmark 6 case (b); Singh & Pal 2013)",
                "window_s": [12.0, 30.0], "t_end_s": 31.0},
    }
    R.write_json(R.INPUTS / "cdw68_design_v1.json", design)
    files = sorted(p for p in R.INPUTS.rglob("*") if p.is_file() and p.name != "cdw68_inputs_manifest.json")
    man = {str(p.relative_to(R.P68)).replace("\\", "/"): hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    man["configs/ieee68/ieee68_network.json"] = hashlib.sha256(R.NETWORK.read_bytes()).hexdigest()
    R.write_json(R.INPUTS / "cdw68_inputs_manifest.json", man)
    print(json.dumps({"policies": [(p["id"], round(p["u"], 4), round(p["g"], 4), round(p["k"], 4), p["split"]) for p in pol], "dmin": dmin, "r11": r11, "n_elig": len(elig)}, indent=1))


if __name__ == "__main__":
    main()
