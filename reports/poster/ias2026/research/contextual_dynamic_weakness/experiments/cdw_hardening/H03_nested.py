# ruff: noqa: E501
"""H03 nested sign reversal (levels A, B, C, D) with curvature witnesses (prereg H3; theory
theory/CDW_CONTEXTUAL_CURVATURE_THEOREM.md, Propositions C and D).

Old census (discovery, old holdout) is a RE-ANALYSIS; the new holdout is confirmatory."""

from __future__ import annotations

import _hinfra as HI  # noqa: I001  (first: sets sys.path for the CDW modules)

import numpy as np
import pandas as pd

import _cdw as C
import _hdata as HD
import E01_census as E1

V9 = C.V9


def _nested_pairs(masks_lo, masks_hi):
    """All (a, b) with a in masks_lo, b in masks_hi, a strict subset of b."""

    if not len(masks_lo) or not len(masks_hi):
        return np.empty((0, 2), int)
    a = np.asarray(masks_lo)[:, None]
    b = np.asarray(masks_hi)[None, :]
    ok = ((a & b) == a) & (a != b)
    ia, ib = np.nonzero(ok)
    return np.stack([np.asarray(masks_lo)[ia], np.asarray(masks_hi)[ib]], axis=1)


def witness(alpha, status, hz, S1, S2, bi, direction):
    """Canonical chain S1 -> S2 (added units in increasing bus order): largest curvature term with
    the required sign, its unit, the Prop. D bound check and EM-clean flags."""

    add = [b for b in range(9) if (S2 >> b & 1) and not (S1 >> b & 1)]
    m = len(add)
    Sr, terms, clean_any, clean_res = S1, [], False, False
    for b in add:
        i_, j_, ij = Sr | (1 << bi), Sr | (1 << b), Sr | (1 << bi) | (1 << b)
        d = alpha[ij] - alpha[i_] - alpha[j_] + alpha[Sr]
        terms.append((d, b))
        four = (Sr, i_, j_, ij)
        em4 = all(status[x] == "STABLE" and HD.in_em(hz[x]) for x in four)
        if em4 and ((direction == "s2d" and d > 0) or (direction == "d2s" and d < 0)):
            clean_any = True
            clean_res = clean_res or abs(d) >= C.TAU_RES
        Sr = j_
    D = (alpha[S2 | (1 << bi)] - alpha[S2]) - (alpha[S1 | (1 << bi)] - alpha[S1])
    if direction == "s2d":
        d, b = max(terms)
        bound_ok = d >= D / m - 1e-9 * max(1.0, abs(D))
    else:
        d, b = min(terms)
        bound_ok = d <= D / m + 1e-9 * max(1.0, abs(D))
    return d, V9[b], m, bound_ok, clean_any, clean_res, D


def analyse_policy(pid, cen, mg):
    d = cen[pid]
    alpha, status, hz, modes = d["alpha"], d["status"], d["hz"], d["modes"]
    g = mg[mg.pid == pid]
    pairs = []
    for lvl in ("A", "B", "C"):
        e = g[g[f"lv{lvl}"]]
        for i in V9:
            ei = e[e.i == i]
            neg = ei[ei.cls == -1]
            pos = ei[ei.cls == 1]
            dl = dict(zip(ei["mask"], ei.delta, strict=True))
            bi = HD.BIT[i]
            for direction, lo, hi in (("s2d", neg["mask"].to_numpy(), pos["mask"].to_numpy()),
                                      ("d2s", pos["mask"].to_numpy(), neg["mask"].to_numpy())):
                for S1, S2 in _nested_pairs(lo, hi):
                    S1, S2 = int(S1), int(S2)
                    wd, wj, m, bok, cany, cres, D = witness(alpha, status, hz, S1, S2, bi, direction)
                    rec = {"pid": pid, "level": lvl, "i": i, "S1": HD.label_of(S1), "S2": HD.label_of(S2), "dir": direction,
                           "d1": dl[S1], "d2": dl[S2], "mag": min(abs(dl[S1]), abs(dl[S2])), "size1": bin(S1).count("1"),
                           "size2": bin(S2).count("1"), "m": m, "D": D, "wit_d": wd, "wit_j": wj, "wit_bound_ok": bok,
                           "em_clean": cany, "em_clean_res": cres, "hz1": hz[S1], "hz2": hz[S2]}
                    if lvl == "C":
                        c1, c2 = E1.crit_mode(modes[S1]), E1.crit_mode(modes[S2])
                        ok = c1 is not None and c2 is not None and abs(c1["hz"] - c2["hz"]) <= C.DF_MAX
                        rec["mac12"] = E1.mac_c(c1["phi"], c2["phi"]) if ok else np.nan
                        rec["lvD"] = bool(ok and rec["mac12"] >= C.MAC_MIN)
                    pairs.append(rec)
    # arbitrary-context statistics (old A1 / A2 definitions) for L1 and H13 targets
    st = g[g.lvA]
    n_rev = sum(((st.i == i) & (st.cls == -1)).any() and ((st.i == i) & (st.cls == 1)).any() for i in V9)
    tr = g[g.lvB]
    n_rev_tr = sum(((tr.i == i) & (tr.cls == -1)).any() and ((tr.i == i) & (tr.cls == 1)).any() for i in V9)
    return pairs, {"n_rev_arbitrary": int(n_rev), "n_rev_arbitrary_tracked": int(n_rev_tr)}


def summarise(pairs: pd.DataFrame, extra: dict, cen) -> pd.DataFrame:
    rows = []
    for pid in sorted(cen):
        base = cen[pid]["alpha"][0]
        rec = {"pid": pid, "split": HI.split_of(pid), "base_alpha": base, "base_stable": bool(base < 0), **extra[pid]}
        p = pairs[pairs.pid == pid] if len(pairs) else pairs
        for lvl in ("A", "B", "C", "D"):
            q = p[p.level == "C"] if lvl == "D" else p[p.level == lvl]
            if lvl == "D":
                q = q[q.lvD.fillna(False).astype(bool)]
            rec[f"{lvl}_has"] = bool(len(q) > 0)
            rec[f"{lvl}_n_pairs"] = int(len(q))
            rec[f"{lvl}_n_units"] = int(q.i.nunique()) if len(q) else 0
            rec[f"{lvl}_n_s2d"] = int((q.dir == "s2d").sum()) if len(q) else 0
            rec[f"{lvl}_n_d2s"] = int((q.dir == "d2s").sum()) if len(q) else 0
            rec[f"{lvl}_both_dirs"] = bool(rec[f"{lvl}_n_s2d"] > 0 and rec[f"{lvl}_n_d2s"] > 0)
            rec[f"{lvl}_max_mag"] = float(q.mag.max()) if len(q) else np.nan
            rec[f"{lvl}_med_mag"] = float(q.mag.median()) if len(q) else np.nan
            rec[f"{lvl}_med_abs_wit"] = float(q.wit_d.abs().median()) if len(q) else np.nan
            rec[f"{lvl}_frac_em_clean"] = float(q.em_clean.mean()) if len(q) else np.nan
            rec[f"{lvl}_min_m"] = int(q.m.min()) if len(q) else np.nan
        rows.append(rec)
    return pd.DataFrame(rows)


def run():
    cen = HD.load_all_census()
    mg = HD.all_marginals(cen)
    mg.to_parquet(HI.RESULTS / "H03_marginals.parquet", index=False)
    allpairs, extra = [], {}
    for pid in sorted(cen):
        p, x = analyse_policy(pid, cen, mg)
        allpairs += p
        extra[pid] = x
    pairs = pd.DataFrame(allpairs)
    pairs.to_parquet(HI.RESULTS / "H03_nested_pairs.parquet", index=False)
    pol = summarise(pairs, extra, cen)
    pol.to_csv(HI.RESULTS / "H03_policy_summary.csv", index=False)
    gate = {"bound_check_all_ok": bool(pairs.wit_bound_ok.all()) if len(pairs) else None}
    for split in ("new", "old", "discovery"):
        s = pol[(pol.split == split) & pol.base_stable]
        gate[f"{split}_n_base_stable"] = int(len(s))
        gate[f"{split}_n_base_unstable"] = int(((pol.split == split) & ~pol.base_stable).sum())
        for lvl in ("A", "B", "C", "D"):
            gate[f"{split}_frac_{lvl}"] = float(s[f"{lvl}_has"].mean()) if len(s) else np.nan
            gate[f"{split}_frac_{lvl}_both_dirs"] = float(s[f"{lvl}_both_dirs"].mean()) if len(s) else np.nan
        gate[f"{split}_frac_arbitrary_global"] = float((s.n_rev_arbitrary > 0).mean()) if len(s) else np.nan
        gate[f"{split}_frac_arbitrary_tracked"] = float((s.n_rev_arbitrary_tracked > 0).mean()) if len(s) else np.nan
    gate["H3_primary_frac_C_new"] = gate["new_frac_C"]
    gate["H3_pass"] = bool(gate["new_frac_C"] >= 0.75)
    gate["H3_D_secondary_pass"] = bool(gate["new_frac_D"] >= 0.75)
    gate["L1_pass"] = bool(gate["new_frac_arbitrary_global"] >= 0.75)
    gate["L2_pass"] = bool(gate["new_frac_A"] >= 0.75)
    gate["L3_pass"] = bool(gate["new_frac_B"] >= 0.50)
    gate["L4_pass"] = gate["H3_pass"]
    nc = pairs[(pairs.level == "C") & pairs.pid.str.startswith("HARDENING")]
    gate["new_C_pairs"] = int(len(nc))
    gate["new_C_dir_counts"] = nc.dir.value_counts().to_dict() if len(nc) else {}
    gate["new_C_mag_quartiles"] = nc.mag.quantile([0.25, 0.5, 0.75]).round(5).tolist() if len(nc) else []
    gate["new_C_frac_em_clean_witness"] = float(nc.em_clean.mean()) if len(nc) else np.nan
    HI.write_json("H03_gate.json", gate)
    return gate


if __name__ == "__main__":
    import json

    print(json.dumps(run(), indent=1, default=str))
