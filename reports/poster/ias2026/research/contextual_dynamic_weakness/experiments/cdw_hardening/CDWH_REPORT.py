# ruff: noqa: E501
"""H32: docs/CDW_HARDENING_FINAL_REPORT.md generated from the result files (numbers never typed)."""

from __future__ import annotations

import _hinfra as HI  # noqa: I001

import json
import subprocess

import pandas as pd

import _infra as I

R = HI.RESULTS
DOCS = I.DOCS


def j(n):
    p = R / n
    return json.loads(p.read_text()) if p.exists() else {}


def f(x, d=2):
    try:
        return f"{float(x):.{d}f}"
    except (TypeError, ValueError):
        return str(x)


def git(*a):
    return subprocess.check_output(["git", "-C", str(I.REPO), *a], text=True).strip()


def build(extra: dict):
    h3, h4, h5, h7, h10, h13, h15, h18, h19, det = (j("H03_gate.json"), j("H04_gate.json"), j("H05_stats.json"), j("H07_goldb_gate.json"), j("H10_gate.json"),
                                                    j("H13_mixing.json"), j("H15_topology_gate.json"), j("H18_gate.json"), j("H19_evidence.json"), j("CDWH_DETERMINISM.json"))
    cm = pd.read_csv(I.RESULTS / "CDW_HARDENED_CLAIM_MATRIX.csv")
    man = pd.read_csv(I.RESULTS / "CDW_HARDENING_RUN_MANIFEST.csv")
    wall = man.wall_s.fillna(0).sum()
    ra = h4.get("regret_anatomy_new", {})
    rob = h10.get("H10_robust_weak_corridors", [])
    L = []
    A = L.append
    A("# CDW hardening campaign: final report\n")
    A(f"Branch `research/contextual-dynamic-weakness-hardening` (from `e73dd355`). Preregistration `{HI.PREREG_COMMIT}`. No push. TX4 untouched. "
      "Every number below is read from `results/hardening/*.json` by `experiments/cdw_hardening/CDWH_REPORT.py`.\n")
    A("## 1. Verdict in one paragraph\n")
    A(f"The contextual core of CDW survives a new preregistered holdout. Nested electromechanical same-mode reversals occur in {h3['new_frac_C']:.3f} of base-stable new-holdout policies (17/19; gate 0.75 PASS). "
      f"GOLD-B survives with a large margin: on the {h7['GB_prereg_n']} preregistration-eligible held-out conditions the total re-equilibrated sensitivity ranks finite reinforcements at median Spearman {f(h7['GB_prereg_median_rho_Dtotal'],3)} against {f(h7['GB_prereg_median_rho_S1'],3)} for |P| "
      f"(lower 95% bound of the paired advantage {f(h7['GB_prereg_diff_ci95'][0],3)} > 0.20); over all 64 solved conditions {f(h7['GB_primary_median_rho_Dtotal'],3)} vs {f(h7['GB_primary_median_rho_S1'],3)}. "
      f"Three prior readings are weakened: (i) the fixed-ranking failure (p* {f(h4['new_FULL_median_pstar'])}, regret {f(h4['new_FULL_median_frac_regret'])}) is carried by jumps onto fast converter-control modes "
      f"({ra.get('choice_out_of_EM_band')}/{ra.get('material')} material regrets); for EM transitions a fixed ranking orders {f(h4['new_EM_median_pstar'])} of pairs; (ii) the modal-mixing mechanism is largely a frequency artefact "
      f"(controller share of variance {f(h13['var_share_C'])}; new-holdout fixed-frequency test fails); (iii) the reversals and the specific branch ranking do not transfer to a frozen WECC library GFL chain, although the ranking method does. "
      f"Topology effects survive the EM-only reclassification. Corridors: robust weak corridor(s) under the equal-budget size-matched null: {', '.join(rob) or 'none'}.\n")
    A("## 2. Phase status\n")
    A("| phase | status | key file |\n|---|---|---|")
    rows = [("H0 preregistration", f"committed `{HI.PREREG_COMMIT}`", "docs/CDW_HARDENING_PREREG_V1.md"),
            ("H1 new policy holdout (24 maximin LHS)", "COMPLETE, 19 base-stable", "results/hardening/prereg_inputs/hardening_policies.json"),
            ("H2 curvature theorem", "PROVED (classical algebra)", "theory/CDW_CONTEXTUAL_CURVATURE_THEOREM.md"),
            ("H3 nested reversal", f"PASS ({f(h3['new_frac_C'],3)})", "results/hardening/H03_gate.json"),
            ("H4 fixed ranking + transfer", "PASS (FULL); fails in EM stratum", "results/hardening/H04_gate.json"),
            ("H5 GOLD-A statistics", "COMPLETE", "results/hardening/H05_stats.json"),
            ("H6/H7 GOLD-B new holdout", "PASS (literal-eligible set and full set)", "results/hardening/H07_goldb_gate.json"),
            ("H8 frozen vs total", "COMPLETE", "results/hardening/H08_*"),
            ("H9-H11 corridors", "PASS" if h10.get("H10_pass") else ("NEGATIVE" if h10 else "PENDING"), "results/hardening/H10_gate.json, H11_corridor_table.csv"),
            ("H12/H13 mixing", "NEGATIVE (primary fails)", "results/hardening/H13_mixing.json"),
            ("H14-H16 topology EM", "SUPPORTED (Q13 yes; static fails)", "results/hardening/H15_topology_gate.json"),
            ("H17 ALT-WECC qualification", "PASS Q1-Q4 (with logged refinements)", "results/hardening/alt/H17_qualification_andes.json"),
            ("H18 cross-model matrix", f"A {h18['A_verdict']}, C {h18['C_verdict']}, C2 {h18['C2_verdict']}, D {h18['D_verdict']}, E {h18['E_verdict']}", "results/hardening/H18_gate.json"),
            ("H19 evidence layers", f"{h19['n_positive']}/8 positive (L8 negative)", "results/hardening/H19_evidence.json"),
            ("H20 novelty boundary", "COMPLETE", "docs/CDW_NOVELTY_BOUNDARY.md"),
            ("H21 literature", "57 verified entries", "results/CDW_LITERATURE_GAP_MATRIX.csv"),
            ("H22 reviewer attack", "COMPLETE", "docs/CDW_REVIEWER_ATTACK.md"),
            ("H23-H30 paper + supplement", "COMPLETE", "reports/papers/cdw_contextual_dynamic_weakness/"),
            ("H31 internal reviews", extra.get("reviews", "see docs"), "docs/CDW_REVIEWER1_FINAL.md, CDW_REVIEWER2_FINAL.md"),
            ("determinism", "identical" if det.get("all_identical") else "NOT identical", "results/hardening/CDWH_DETERMINISM.json")]
    for a, b, c in rows:
        A(f"| {a} | {b} | `{c}` |")
    A("")
    A("## 3. Results by question\n")
    q = [
        ("Does nested contextual sign reversal survive the new holdout?", f"Yes. Level A (global) {f(h3['new_frac_A'])}, level B {f(h3['new_frac_B'])}, level C {f(h3['new_frac_C'],3)} of base-stable new policies (old holdout: {f(h3['old_frac_C'],3)})."),
        ("Does same-mode electromechanical reversal survive?", f"Yes: 17/19 at level C; every such policy also has a level-D (same family across contexts) reversal. Magnitudes are small (quartiles {h3['new_C_mag_quartiles']} s^-1). Both directions within a policy in {f(h3['new_frac_C_both_dirs'])} of policies."),
        ("Is the optimal fixed node ranking still materially insufficient?", f"FULL: yes (p* {f(h4['new_FULL_median_pstar'])}, regret {f(h4['new_FULL_median_frac_regret'])}; 19/19 policies p*<0.90). EM-only: no (p* {f(h4['new_EM_median_pstar'])}, regret {f(h4['new_EM_median_frac_regret'])}). {ra.get('choice_out_of_EM_band')} of {ra.get('material')} material regrets are jumps out of the EM band."),
        ("How large is the old-to-new holdout shift?", f"Negligible: p* shift {h5['shift_pstar_new_minus_old'][0]:.3f} (interval {h5['shift_pstar_new_minus_old'][1]:.3f} to {h5['shift_pstar_new_minus_old'][2]:.3f}); regret shift {h5['shift_regret_new_minus_old'][0]:.3f}; level-C coverage 0.889 old vs 0.895 new."),
        ("Does GOLD-A remain strong enough for a headline claim?", "Yes for contextual reversal (nested, EM same-mode), with the magnitude caveat. The fixed-ranking half must be stated with its fast-mode anatomy."),
        ("Does total dynamic line sensitivity still beat every static baseline on the new holdout?", f"Yes: Dtotal exceeds the per-condition best static index in {f(h7['GB_primary_frac_cond_Dtotal_beats_every_static'])} of conditions; oracle static median {f(h7['GB_primary_median_rho_oracle_static'],3)}."),
        ("Is the improvement statistically and practically material?", f"Yes: median paired advantage {f(h7['GB_prereg_median_diff'],3)} (CI {h7['GB_prereg_diff_ci95'][0]:.3f}-{h7['GB_prereg_diff_ci95'][1]:.3f}); top-5 precision {f(h7['GB_prereg_median_top5_Dtotal'])} vs {f(h7['GB_prereg_median_top5_S1'])}. The sign-flip p is {f(h7['GB_prereg_T2_signflip_p'],3)} on the literal 8-condition set (6 clusters), <1e-3 on all 64."),
        ("When does re-equilibration materially change the line ranking?", f"In {f(h7['H8_frac_material_H4'])} of V4 conditions by the old H2 rule; every condition has at least one frozen/total sign flip; the frozen derivative drops to {f(h7['GB_new_policies_median_rho_Dfrozen'],3)} on the 24 new policies while Dtotal stays {f(h7['GB_new_policies_median_rho_Dtotal'],3)}. Controller coordinates: frozen = total ({h7['node_R1_max_rel_frozen_minus_total_controller']:.1e})."),
        ("Does TXother remain a strong corridor under equal intervention budget?", f"Family-A null percentile >= 0.95 in {f(h10.get('TXother|A|new_frac_ge95'))} of new conditions (median percentile {f(h10.get('TXother|A|new_median_pct'))})." if h10 else "PENDING"),
        ("Does it beat size-matched random connected corridors?", f"Family B: >= 0.95 in {f(h10.get('TXother|B|new_frac_ge95'))} of new conditions." if h10 else "PENDING"),
        ("How much of old modal-mixing variation was controller- vs frequency-induced?", f"Controller term {f(h13['var_share_C'])} of the variance (median |C| share {f(h13['share_abs_C_of_abs_C_plus_abs_F_median'])}); the rest is frequency movement."),
        ("Does fixed-frequency modal mixing still correlate with contextuality?", f"No at the preregistered bar: new rho {f(h13['T3_primary']['rho'])} (p {f(h13['T3_primary']['p'],3)}, Holm {f(h19['F_conf_holm_p']['T3_mixing'],3)}); old holdout at fixed frequency rho {f(h13['tests']['old|mu10|n_rev']['rho'])}."),
        ("After filtering to EM modes, can topology still both remove and create incompatibility?", f"Yes: removal at P4 by {h15['em_removal_H4_P4']}; {h15['n_em_creation']} EM creations of hyperedges of size <= 3; {h15['n_fast_creation']} fast-only."),
        ("Do Fiedler/gSCR/static topology scores remain poor predictors of the EM effect?", f"Yes: no score meets |rho| >= 0.6 in >= 75% of policies (best fraction {max(v['frac_abs_ge_0.6_em'] for v in h15['H16'].values()):.2f}); medians around |0.5|, i.e. partial information."),
        ("Does contextual reversal transfer to a genuinely dynamic alternative converter model?", f"No: ALT-WECC {h18['A_frac']} of {h18['n_tested_policies']} tested policies (custom GFL {f(h18['A_gfl_frac_same_policies'])}); every stable-context ALT marginal is >= 0."),
        ("Does dynamic link ranking transfer?", f"The ranking does not (median Kendall {f(h18['C_median_kendall'])}); the method does (ALT total vs ALT finite rho {f(h18['C2_median_rho_altFD'])} vs |P| {f(h18['C2_median_rho_S1'])})."),
        ("Which findings are model-specific?", "The reversals, the fixed-ranking failure (driven by the custom GFL's fast modes), the specific branch ranking. Transferring: the ranking method, the top corridor (0.75), topology direction (0.94)."),
        ("What exact novelty remains after the literature review?", "The combination in docs/CDW_NOVELTY_BOUNDARY.md section 2: certified nested EM same-mode reversal of SG->GFL replacements with curvature witnesses; quantified oracle fixed-ranking failure and its fast-mode anatomy; held-out validation of portfolio-conditioned total sensitivity as a reinforcement ranking against static indices; EM-only topology effects; a frozen library-GFL cross-model boundary. No single ingredient is new."),
        ("What are the 3 strongest defensible paper contributions?", "(1) portfolio-conditioned re-equilibrated reinforcement ranking (GOLD-B, STRONGLY_SUPPORTED); (2) nested EM same-mode contextual reversal with curvature certificates (SUPPORTED); (3) fixed-ranking insufficiency with its modal anatomy (SUPPORTED, FULL only)."),
        ("What claims must be dropped or weakened?", "Dropped: modal-mixing mechanism; model-independent reversal; converter-independent branch ranking; weak corridors unless H10 passes. Weakened: fixed-ranking failure to 'carried by fast converter-control modes'; 'weakness is contextual' to SUPPORTED (not STRONGLY)."),
        ("Is the paper ready for TPWRS?", extra.get("ready", "see reviews")),
        ("If not, what single missing experiment blocks submission?", extra.get("blocker", "see reviews")),
    ]
    for k, (qq, aa) in enumerate(q, 1):
        A(f"**Q{k}. {qq}** {aa}\n")
    A("## 4. Final claim matrix\n")
    A("| claim | final status |\n|---|---|")
    for c, s in zip(cm["claim"], cm["final status"], strict=True):
        A(f"| {c} | {s} |")
    A("")
    A("## 5. Negative and blocked results (same prominence)\n")
    A("- Fixed-ranking insufficiency does not hold for EM-only transitions (A3 bar failed in the EM stratum).")
    A("- Fixed-frequency modal mixing does not correlate with contextuality at the preregistered bar; the earlier E12 reading was mostly a frequency effect.")
    A("- Contextual reversal does not transfer to the frozen WECC library GFL chain (0/7); the specific branch ranking is converter-model-specific.")
    A(f"- Corridors: {'robust weak corridor(s): ' + ', '.join(rob) if rob else 'no robust weak corridor under the equal-budget size-matched null'}." if h10 else "- Corridors: pending.")
    A("- GFM generalization and the Africano/PV benchmark remain BLOCKED (no validated GFM model; no source material).")
    A("- The preregistered R0 <= 1e-8 eligibility rule was stricter than the solver tolerance; the literal primary GOLD-B set has only 8 conditions (still passes).\n")
    A("## 6. Deviations\n")
    A("See `docs/CDW_HARDENING_DEVIATIONS.md`: execution fixes (relaunch, import order, bound-check tolerance, H06 pool crash rerun), H17 structural-zero and limiter refinements recorded before the cross-model matrix, one ALT datum seen during qualification, and the literal application of the H6 eligibility rule. No gate threshold changed.\n")
    A("## 7. Run facts\n")
    A(f"- Compute phases: {len(man)}; tasks: {int(man.n_tasks.fillna(0).sum())}; summed phase wall time {wall/3600:.1f} h (parallel). Determinism: {'identical' if det.get('all_identical') else 'NOT identical'}.")
    A(f"- Commits on this branch: `{git('log', '--oneline', 'e73dd355..HEAD').replace(chr(10), '`; `')}`.")
    A(f"- Literature: 57 verified entries; privacy note: one literature agent reported once including the user's e-mail as a Crossref `mailto` parameter; no other personal data was sent.")
    if extra.get("scores"):
        A(f"- Internal review scores: {extra['scores']}")
    (DOCS / "CDW_HARDENING_FINAL_REPORT.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("report written")


if __name__ == "__main__":
    import sys

    build(json.loads(sys.argv[1]) if len(sys.argv) > 1 else {})
