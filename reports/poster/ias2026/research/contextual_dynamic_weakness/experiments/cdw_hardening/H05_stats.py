# ruff: noqa: E501
"""H05 GOLD-A statistics (prereg H5; statistical plan) and L7 fresh-draw robustness (H01e)."""

from __future__ import annotations

import _hinfra as HI  # noqa: I001  (first: sets sys.path for the CDW modules)

import numpy as np
import pandas as pd

import _cdw as C
import _hstats as ST
import _infra as I
import E01_census as E1


def lattice_rev(recs):
    """Core-lattice (V4) stable-context reversal: global, same-mode (lvB) and EM same-mode (lvC)."""

    by = {r["S"]: r for r in recs}
    modes = {k: [E1.compact(m) for m in v.get("modes", [])] for k, v in by.items()}
    cls = {"A": {}, "B": {}, "C": {}}
    for S in C.subsets(C.V4):
        lS = C.label(S)
        if by[lS]["status"] != "STABLE":
            continue
        cm = E1.crit_mode(modes[lS])
        for i in C.V4:
            if i in S:
                continue
            lSi = C.label(S + (i,))
            d = by[lSi]["alpha"] - by[lS]["alpha"]
            c = -1 if d <= -C.TAU_MAT else (1 if d >= C.TAU_MAT else 0)
            cls["A"].setdefault(i, set()).add(c)
            if cm is not None and modes[lSi]:
                m, _ = E1.match(cm, modes[lSi])
                if m is not None and m["crit"]:
                    cls["B"].setdefault(i, set()).add(c)
                    if HI.in_band(by[lS]["lam_hz"]) and HI.in_band(by[lSi]["lam_hz"]):
                        cls["C"].setdefault(i, set()).add(c)
    return {lvl: any({-1, 1} <= s for s in d.values()) for lvl, d in cls.items()}, by["BASE"]["status"]


def l7():
    rows = []
    for r in I.Store("H_H01e").all():
        if not r.get("ok"):
            continue
        t = r["task"]
        rev, st0 = lattice_rev(r["records"])
        rows.append({"env": t["env"], "draw": t["draw"], "base_status": st0, **{f"rev_{k}": v for k, v in rev.items()}})
    df = pd.DataFrame(rows)
    df.to_csv(HI.RESULTS / "H05_L7_fresh_draws.csv", index=False)
    cov = df.groupby("env")[["rev_A", "rev_B", "rev_C"]].mean()
    return {"L7_cov_by_env_global": cov.rev_A.round(3).to_dict(), "L7_cov_by_env_samemode": cov.rev_B.round(3).to_dict(),
            "L7_cov_by_env_EM_samemode": cov.rev_C.round(3).to_dict(),
            "L7_n_envs_ge_0.5": int((cov.rev_A >= 0.5).sum()), "L7_pass": bool((cov.rev_A >= 0.5).sum() >= 3),
            "L7_EM_n_envs_ge_0.5": int((cov.rev_C >= 0.5).sum())}


def run():
    pol = pd.read_csv(HI.RESULTS / "H03_policy_summary.csv")
    rk = pd.read_csv(HI.RESULTS / "H04_ranking.csv")
    pairs = pd.read_parquet(HI.RESULTS / "H03_nested_pairs.parquet")
    out = {}
    stab = pol[pol.base_stable]
    for split in ("new", "old", "discovery"):
        s = stab[stab.split == split]
        for lvl in ("A", "B", "C", "D"):
            pt, lo, hi, _ = ST.cluster_boot(s[f"{lvl}_has"].astype(float), stat=np.mean)
            out[f"{split}_cov_{lvl}"] = [pt, lo, hi]
        pc = pairs[(pairs.level == "C") & pairs.pid.isin(s.pid)]
        mag = pc.groupby("pid").mag.median()
        out[f"{split}_C_mag_per_policy_median"] = list(ST.cluster_boot(mag.to_numpy())[:3]) if len(mag) else None
        out[f"{split}_C_mag_IQR"] = mag.quantile([0.25, 0.75]).tolist() if len(mag) else None
        for strat in ("FULL", "EM", "SAME"):
            q = rk[(rk.split == split) & rk.base_stable & (rk.stratum == strat)]
            out[f"{split}_{strat}_pstar"] = list(ST.cluster_boot(q.p_star.to_numpy())[:3]) if len(q) else None
            out[f"{split}_{strat}_pstar_IQR"] = q.p_star.quantile([0.25, 0.75]).tolist() if len(q) else None
            out[f"{split}_{strat}_regret"] = list(ST.cluster_boot(q.frac_regret.to_numpy())[:3]) if len(q) else None
            out[f"{split}_{strat}_regret_IQR"] = q.frac_regret.quantile([0.25, 0.75]).tolist() if len(q) else None
    # old -> new shift (FULL)
    fn = rk[(rk.stratum == "FULL") & rk.base_stable]
    out["shift_pstar_new_minus_old"] = ST.boot_diff_medians(fn[fn.split == "new"].p_star, fn[fn.split == "old"].p_star)
    out["shift_regret_new_minus_old"] = ST.boot_diff_medians(fn[fn.split == "new"].frac_regret, fn[fn.split == "old"].frac_regret)
    # T1 sign test (F_conf)
    k, n, p = ST.sign_test_less(fn[fn.split == "new"].p_star, 0.90)
    out["T1_sign_test"] = {"k_below_0.90": k, "n": n, "p_one_sided": p}
    out.update(l7())
    HI.write_json("H05_stats.json", out)
    return out


if __name__ == "__main__":
    import json

    print(json.dumps(run(), indent=1, default=str))
