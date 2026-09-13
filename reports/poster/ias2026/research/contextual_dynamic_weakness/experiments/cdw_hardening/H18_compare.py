# ruff: noqa: E501
"""H18 cross-model matrix: ALT-WECC (ANDES, TX3 library GFL chain) vs the custom GFL (prereg H18).

Custom-GFL values come from the existing census (E01 / H_H01), E4 / H06 link runs, H09 corridors
and H18c topology evaluations; ALT values from results/hardening/alt/cases."""

from __future__ import annotations

import _hinfra as HI  # noqa: I001

import json

import numpy as np
import pandas as pd

import _analysis as AN
import _cdw as C
import _hdata as HD
import _infra as I
import E01_census as E1
import E34_sens as E34
import H09_corridors as H9
import H14_topology as H14

ALT = HI.RESULTS / "alt"
POL = H14.H18_POLICIES
LINES = (3, 9, 15, 17, 22, 26, 31, 32, 35, 39, 42, 44)
CORR = ("TXother", "TXall", "K2_01")
ACTS = H14.H18_ACTIONS
VALID = ("STABLE", "UNSTABLE")


def lattice_reversal(recs):
    """recs {label: {'status','alpha','hz','modes'}} over the V4 lattice -> (global, same-mode, EM same-mode, frac EM crit)."""

    cls = {"A": {}, "B": {}, "C": {}}
    n_stab = n_em = 0
    for S in C.subsets(C.V4):
        lS = C.label(S)
        r = recs.get(lS)
        if r is None or r["status"] != "STABLE":
            continue
        n_stab += 1
        n_em += HI.in_band(r["hz"])
        cm = E1.crit_mode(r["modes"]) if r["modes"] else None
        for i in C.V4:
            if i in S:
                continue
            ri = recs.get(C.label(S + (i,)))
            if ri is None or ri["status"] not in VALID:
                continue
            d = ri["alpha"] - r["alpha"]
            c = -1 if d <= -C.TAU_MAT else (1 if d >= C.TAU_MAT else 0)
            cls["A"].setdefault(i, set()).add(c)
            if cm is not None and ri["modes"]:
                m, _ = E1.match(cm, ri["modes"])
                if m is not None and m["crit"]:
                    cls["B"].setdefault(i, set()).add(c)
                    if HI.in_band(r["hz"]) and HI.in_band(ri["hz"]):
                        cls["C"].setdefault(i, set()).add(c)
    out = {k: any({-1, 1} <= s for s in d.values()) for k, d in cls.items()}
    out["frac_stable_ctx_EM"] = n_em / n_stab if n_stab else np.nan
    return out


def load_alt():
    rows = [json.loads(p.read_text()) for p in (ALT / "cases").glob("*.json")]
    for r in rows:
        r["modes_c"] = [E1.compact(m) for m in r.get("modes", [])]
    return rows


def custom_core(cen, pid):
    d = cen[pid]
    out = {}
    for m in range(512):
        lab = HD.label_of(m)
        if all(int(u) in C.V4 for u in (lab.split("+") if lab != "BASE" else [])):
            out[lab] = {"status": d["status"][m], "alpha": d["alpha"][m], "hz": d["hz"][m], "modes": d["modes"][m]}
    return out


def custom_lines():
    old = pd.read_parquet(I.RESULTS / "CDW_E4_link_sensitivity.parquet")
    old = old[(old.semantics == "SPR") & (old.target == "H4") & old.env.isna()]
    out = {}
    for pid in POL:
        if pid.startswith("HARDENING"):
            continue
        g = old[old.pid == pid].set_index("idx")
        out[pid] = {e: -(g.loc[e, "large_x1.5_alpha"] - g.loc[e, "alpha0"]) for e in LINES}
    try:
        new = pd.read_parquet(HI.RESULTS / "H06_links_raw.parquet")
        new = new[(new.family == "link") & (new.target == "H4") & new.env.isna()]
        for pid in POL:
            if pid.startswith("HARDENING"):
                g = new[new.pid == pid].set_index("idx")
                out[pid] = {e: -(g.loc[e, "g1.50_alpha"] - g.loc[e, "alpha0"]) for e in LINES}
    except FileNotFoundError:
        pass
    return out


def run():
    alt = load_alt()
    A = {(r["pid"], r["S"], r["net"]): r for r in alt}
    cen = {**HD.load_census("E01"), **HD.load_census("H_H01")}
    cen = {k: v for k, v in cen.items() if k in POL}
    rows = []
    cl = custom_lines()
    h9 = pd.read_csv(HI.RESULTS / "H09_corridors.csv") if (HI.RESULTS / "H09_corridors.csv").exists() else None
    tc = {}
    for r in I.Store("H_H18c").all():
        if r.get("ok"):
            tc[(r["task"]["pid"], r["task"]["action"])] = r["record"]["alpha"]
    st = E34.static_link_baselines(C.V4).set_index("e")
    counts = {"alt_cases": len(alt), "alt_status": pd.Series([r["status"] for r in alt]).value_counts().to_dict(),
              "alt_init_residual_max": float(np.nanmax([r.get("init_residual", np.nan) for r in alt])),
              "alt_frac_init_ok_core": float(np.mean([(r.get("init_residual", 1) <= 1e-6) for r in alt if r["net"] == "NOMINAL"]))}
    for pid in POL:
        rec = {"pid": pid}
        # ------------------------------------------------ A / B: core lattice ---
        altcore = {}
        for lab in [C.label(s) for s in C.subsets(C.V4)]:
            r = A.get((pid, lab, "NOMINAL"))
            if r is not None:
                altcore[lab] = {"status": r["status"], "alpha": r.get("alpha", np.nan), "hz": r.get("lam_hz", np.nan), "modes": r["modes_c"]}
        ra = lattice_reversal(altcore)
        rc = lattice_reversal(custom_core(cen, pid))
        rec.update({f"alt_rev_{k}": v for k, v in ra.items()})
        rec.update({f"gfl_rev_{k}": v for k, v in rc.items()})
        rec["alt_base_status"] = altcore.get("BASE", {}).get("status")
        rec["alt_H4_status"] = altcore.get("30+33+35+37", {}).get("status")
        rec["alt_H4_alpha"] = altcore.get("30+33+35+37", {}).get("alpha")
        rec["alt_H4_hz"] = altcore.get("30+33+35+37", {}).get("hz")
        st_alt = {lab: v["status"] for lab, v in altcore.items()}
        rec["alt_n_limit_active_core"] = int(sum(v == "LIMIT_ACTIVE" for v in st_alt.values()))
        rec["alt_n_unresolved_core"] = int(sum(v == "BOUNDARY_OR_UNRESOLVED" for v in st_alt.values()))
        # ------------------------------------------------ C / C2: lines ----------
        h4 = A.get((pid, "30+33+35+37", "NOMINAL"))
        a0 = h4.get("alpha") if h4 and h4["status"] in VALID else np.nan
        alt_eff, alt_fd = {}, {}
        for e in LINES:
            r15 = A.get((pid, "30+33+35+37", f"L{e}_g1.5"))
            rp, rm = A.get((pid, "30+33+35+37", f"L{e}_g1.001")), A.get((pid, "30+33+35+37", f"L{e}_g0.999"))
            alt_eff[e] = -(r15["alpha"] - a0) if (r15 and r15["status"] in VALID) else np.nan
            alt_fd[e] = -(rp["alpha"] - rm["alpha"]) / 0.002 if (rp and rm and rp["status"] in VALID and rm["status"] in VALID) else np.nan
        ce = cl.get(pid, {})
        x = np.array([ce.get(e, np.nan) for e in LINES])
        y = np.array([alt_eff[e] for e in LINES])
        fd = np.array([alt_fd[e] for e in LINES])
        rec["C_kendall"] = AN.kendall(x, y)
        rec["C_spearman"] = AN.spearman(x, y)
        rec["C_sign_agree_material"] = float(np.mean(np.sign(x[np.abs(x) >= C.TAU_MAT]) == np.sign(y[np.abs(x) >= C.TAU_MAT]))) if np.any(np.abs(x) >= C.TAU_MAT) else np.nan
        s1 = st.loc[list(LINES), "S1_absP"].to_numpy(float)
        rec["C2_rho_altFD"] = AN.spearman(fd, y)
        rec["C2_rho_S1"] = AN.spearman(s1, y)
        rec["C2_diff"] = rec["C2_rho_altFD"] - rec["C2_rho_S1"] if np.isfinite(rec["C2_rho_altFD"]) and np.isfinite(rec["C2_rho_S1"]) else np.nan
        rec["C_n_finite"] = int((np.isfinite(x) & np.isfinite(y)).sum())
        # ------------------------------------------------ D: corridors -----------
        if h9 is not None:
            g = h9[(h9.pid == pid) & (h9.budget == "L1_0.50") & (h9["set"].isin(["old", "new"])) & h9.corridor.isin(CORR)]
            gfl_top = g.sort_values("eff", ascending=False).corridor.iloc[0] if len(g) else None
            ae = {}
            for c in CORR:
                r = A.get((pid, "30+33+35+37", f"C_{c}"))
                ae[c] = -(r["alpha"] - a0) if (r and r["status"] in VALID) else np.nan
            ae_s = pd.Series(ae).dropna()
            rec["D_gfl_top"] = gfl_top
            rec["D_alt_top"] = ae_s.idxmax() if len(ae_s) else None
            rec["D_top_agree"] = bool(gfl_top is not None and rec["D_alt_top"] == gfl_top)
        # ------------------------------------------------ E: topology ------------
        agree, n = 0, 0
        g0 = tc.get((pid, "NOMINAL"))
        for a in ACTS:
            ga = tc.get((pid, a))
            r = A.get((pid, "30+33+35+37", f"T_{a}"))
            if ga is None or g0 is None or r is None or r["status"] not in VALID or not np.isfinite(a0):
                continue
            dg, da = ga - g0, r["alpha"] - a0
            if abs(dg) >= C.TAU_MAT:
                n += 1
                agree += np.sign(dg) == np.sign(da)
            rec[f"E_{a}_gfl"] = dg
            rec[f"E_{a}_alt"] = da
        rec["E_n"] = n
        rec["E_agree"] = agree / n if n else np.nan
        rows.append(rec)
    df = pd.DataFrame(rows)
    df.to_csv(HI.RESULTS / "H18_crossmodel.csv", index=False)
    tested = df[df.alt_base_status == "STABLE"]
    gate = {**counts, "n_tested_policies": int(len(tested))}

    fa = float(tested.alt_rev_A.mean()) if len(tested) else np.nan
    gate["A_frac"] = fa
    gate["A_verdict"] = "TRANSFERS" if fa >= 0.5 else ("PARTIAL" if fa > 0 else "MODEL-SPECIFIC")
    gate["A_gfl_frac_same_policies"] = float(tested.gfl_rev_A.mean()) if len(tested) else np.nan
    em_share = float(tested.alt_rev_frac_stable_ctx_EM.mean()) if len(tested) else np.nan
    gate["B_meaningful"] = bool(em_share >= 0.5)
    gate["B_alt_mean_frac_stable_ctx_EM"] = em_share
    fb = float(tested.alt_rev_B.mean()) if len(tested) else np.nan
    gate["B_frac"] = fb
    gate["B_verdict"] = ("TRANSFERS" if fb >= 0.5 else ("PARTIAL" if fb > 0 else "MODEL-SPECIFIC")) if gate["B_meaningful"] else "NOT_MEANINGFUL"
    gate["B_EM_frac"] = float(tested.alt_rev_C.mean()) if len(tested) else np.nan
    mk = float(df.C_kendall.median())
    gate["C_median_kendall"] = mk
    gate["C_median_spearman"] = float(df.C_spearman.median())
    gate["C_verdict"] = "TRANSFERS" if mk >= 0.5 else ("PARTIAL" if mk >= 0.2 else "MODEL-SPECIFIC")
    md = float(df.C2_diff.median())
    gate["C2_median_diff"] = md
    gate["C2_median_rho_altFD"] = float(df.C2_rho_altFD.median())
    gate["C2_median_rho_S1"] = float(df.C2_rho_S1.median())
    gate["C2_verdict"] = "TRANSFERS" if md >= 0.2 else ("PARTIAL" if md >= 0 else "MODEL-SPECIFIC")
    if "D_top_agree" in df:
        fd_ = float(df.D_top_agree.mean())
        gate["D_frac_top1_agree"] = fd_
        gate["D_verdict"] = "TRANSFERS" if fd_ >= 0.75 else ("PARTIAL" if fd_ >= 0.5 else "MODEL-SPECIFIC")
    tot_n = df.E_n.sum()
    agree_all = float((df.E_agree * df.E_n).sum() / tot_n) if tot_n else np.nan
    gate["E_pooled_sign_agree"] = agree_all
    gate["E_n_pairs"] = int(tot_n)
    gate["E_verdict"] = "TRANSFERS" if agree_all >= 0.75 else ("PARTIAL" if agree_all >= 0.5 else "MODEL-SPECIFIC")
    gate["L8_pass"] = gate["A_verdict"] == "TRANSFERS"
    HI.write_json("H18_gate.json", gate)
    return gate


if __name__ == "__main__":
    print(json.dumps(run(), indent=1, default=str))
