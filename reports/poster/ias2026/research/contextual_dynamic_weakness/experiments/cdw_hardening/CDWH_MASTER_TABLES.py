# ruff: noqa: E501
"""H27/H32: results/CDW_HARDENED_CLAIM_MATRIX.csv, results/CDW_HARDENING_MASTER_RESULTS.csv,
results/CDW_HARDENING_RUN_MANIFEST.csv (all values read from the result files)."""

from __future__ import annotations

import _hinfra as HI  # noqa: I001

import hashlib
import json

import pandas as pd

import _infra as I

R = HI.RESULTS
TOP = I.RESULTS


def j(n):
    p = R / n
    return json.loads(p.read_text()) if p.exists() else {}


def f(x, d=2):
    try:
        return f"{float(x):.{d}f}"
    except (TypeError, ValueError):
        return str(x)


def claim_matrix():
    h3, h4, h5, h7, h10, h13, h15, h18, h19 = (j("H03_gate.json"), j("H04_gate.json"), j("H05_stats.json"), j("H07_goldb_gate.json"),
                                               j("H10_gate.json"), j("H13_mixing.json"), j("H15_topology_gate.json"), j("H18_gate.json"), j("H19_evidence.json"))
    rob = h10.get("H10_robust_weak_corridors")
    rv, ex = j("H31_revision.json"), j("H31_explore.json")
    from scipy.stats import beta
    k_c, n_c = 17, 19
    cp = (beta.ppf(0.025, k_c, n_c - k_c + 1), beta.ppf(0.975, k_c + 1, n_c - k_c))
    ld = rv.get("levelD_base_stable", {})
    dc = rv.get("direction_counts_new", {})
    sc = rv.get("screened_ranking_median", {}).get("new", {})
    lo = rv.get("loco_fixed_branch_list", {}).get("new", {})
    dt = rv.get("dtot_vs_topology_median", {})
    fm = rv.get("fast_modes_new", {})
    mv = rv.get("mixing_variance_terms", {})
    rows = [
        dict(claim="Chain identity locates an interaction term of the matching sign for every nested reversal", prereg_hypothesis="H2 (theory)",
             old_evidence="E1b: alpha neither sub- nor supermodular at 39 policies", new_hardening_evidence="theory/CDW_CONTEXTUAL_CURVATURE_THEOREM.md; Prop. D bound verified on all 208291 nested pairs",
             discovery="n/a", old_holdout="n/a", new_holdout="bound holds 100%", uncertainty="n/a", alternate_dynamic_model="n/a", full_mode="n/a",
             electromechanical_only="94% of level-C reversals have an EM-clean witness", same_mode="n/a", effect_size="n/a", statistic="exact algebra",
             final_status="PROVED (classical algebra)", allowed_wording="by the classical chain identity each nested reversal locates an interaction term of the matching sign",
             prohibited_wording="new identity; novel curvature theory; certificate as an empirical finding (the bound holds for every nested pair by construction)"),
        dict(claim="Nested EM same-mode contextual sign reversal", prereg_hypothesis="H3 primary (level C >= 0.75 new holdout)",
             old_evidence="A1 18/18, A2 16/18 (arbitrary contexts)", new_hardening_evidence=f"level C {h3.get('new_frac_C'):.3f} new (17/19); level D 17/19 ({ld.get('pairs')} pairs, s2d {ld.get('s2d')} / d2s {ld.get('d2s')}); both directions at level D {ld.get('policies_both_dirs')}/19",
             discovery=f(h3.get("discovery_frac_C")), old_holdout=f(h3.get("old_frac_C")), new_holdout=f(h3.get("new_frac_C")),
             uncertainty=f"fresh draws EM same-mode coverage by env {h5.get('L7_cov_by_env_EM_samemode')}", alternate_dynamic_model=f"ALT-WECC {h18.get('B_verdict')} (0/7)",
             full_mode=f"level A {f(h3.get('new_frac_A'))}", electromechanical_only=f"level C {f(h3.get('new_frac_C'))}", same_mode=f"level D {f(h3.get('new_frac_D'))}",
             effect_size=f"level-D magnitude quartiles {ld.get('mag_quartiles')} s^-1 (about 0.3% damping ratio), max {ld.get('max_mag', float('nan')):.3f}; coverage 0 at tau 0.05",
             statistic=f"bootstrap interval {h5.get('new_cov_C')}; exact binomial [{cp[0]:.2f}, {cp[1]:.2f}]; gate passes on the point estimate only",
             final_status="SUPPORTED (small effects; one network; not reproduced by the library converter)",
             allowed_wording="in 17/19 base-stable new-holdout policies the same replacement moves the same tracked EM mode in opposite directions in two nested stable portfolios; effects are about 0.01 s^-1",
             prohibited_wording="universal; large effect; model-independent; planning-relevant damping change; strongly supported"),
        dict(claim="alpha_perp is neither submodular nor supermodular (both nested directions)", prereg_hypothesis="H3/Prop. C",
             old_evidence="E1b all 39 policies", new_hardening_evidence=f"level-C both directions in {f(h3.get('new_frac_C_both_dirs'))} of new policies; level A both directions {f(h3.get('new_frac_A_both_dirs'))}",
             discovery=f(h3.get("discovery_frac_A_both_dirs")), old_holdout=f(h3.get("old_frac_A_both_dirs")), new_holdout=f(h3.get("new_frac_A_both_dirs")),
             uncertainty="n/a", alternate_dynamic_model="ALT global alpha pinned by an isolated pole (uninformative)", full_mode="level A",
             electromechanical_only=f"level C both directions {f(h3.get('new_frac_C_both_dirs'))}; EM-clean squares with both signs in 18/19 new, 16/18 old",
             same_mode=f"level D both directions {dc.get('D', {}).get('both_policies')}/19", effect_size="second differences beyond 1e-3 s^-1 of both signs", statistic="one-step squares (direct)", final_status="SUPPORTED",
             allowed_wording="on the tested class alpha_perp is neither submodular nor supermodular; no algorithmic consequence follows because the stable family is not a lattice", prohibited_wording="greedy fails; new non-submodularity result"),
        dict(claim="Oracle fixed node ranking materially insufficient", prereg_hypothesis="H4-new (A3 rule)", old_evidence="old holdout p* 0.752, regret 0.369",
             new_hardening_evidence=f"new p* {f(h4.get('new_FULL_median_pstar'))}, regret {f(h4.get('new_FULL_median_frac_regret'))}",
             discovery=f"p* {f(h4.get('discovery_FULL_median_pstar'))}", old_holdout=f"p* {f(h4.get('old_FULL_median_pstar'))}", new_holdout=f"p* {f(h4.get('new_FULL_median_pstar'))}",
             uncertainty="n/a", alternate_dynamic_model="not tested", full_mode=f"regret {f(h4.get('new_FULL_median_frac_regret'))}",
             electromechanical_only=f"p* {f(h4.get('new_EM_median_pstar'))}, regret {f(h4.get('new_EM_median_frac_regret'))} (fails A3 bar)",
             same_mode=f"p* {f(h4.get('new_SAME_median_pstar'))}, regret {f(h4.get('new_SAME_median_frac_regret'))}",
             effect_size=f"{h4.get('regret_anatomy_new', {}).get('share_out_of_band', float('nan')):.2f} of material regrets are jumps out of the EM band",
             statistic=f"T1 sign test p (Holm) {h19.get('F_conf_holm_p', {}).get('T1_goldA_sign')}; regret-optimal order {rv.get('min_regret_ranking_median', {}).get('new|FULL')}; screened {f(sc.get('regret_screened_on_feasible'), 3)} on feasible contexts",
             final_status="SUPPORTED (FULL only; driven by singularity-type crossings; nearly adequate with a stability screen)",
             allowed_wording="a fixed ranking errs materially in about a third of decisions because it cannot anticipate aperiodic singularity crossings at high converter share; screened for stability it errs in about 6%",
             prohibited_wording="fixed rankings fail for electromechanical dynamics; weak buses do not exist; fast converter-control modes (unqualified)"),
        dict(claim="Fixed rankings do not transfer across policies", prereg_hypothesis="L6", old_evidence="optimal orders differ across policies",
             new_hardening_evidence=f"offdiag regret {f(h4.get('transfer_new_median_offdiag_regret_by_test'))} vs diag {f(h4.get('transfer_new_median_diag_regret'))}",
             discovery="in matrix", old_holdout="in matrix", new_holdout=f(h4.get("transfer_new_median_offdiag_regret_by_test")), uncertainty="n/a",
             alternate_dynamic_model="n/a", full_mode="yes", electromechanical_only="not computed", same_mode="not computed",
             effect_size=f"{h4.get('n_distinct_optimal_orders_full')} distinct optimal orders among {h4.get('transfer_n_policies')} policies",
             statistic=f"distance-degradation Spearman {f(h4.get('distance_vs_degradation_spearman'))} (descriptive); median Kendall between optimal orders {f(rv.get('optimal_order_distance', {}).get('median_pairwise_kendall'))}",
             final_status="NOT SUPPORTED as a separate claim (transfer adds 0.026 regret; L6 rule met)",
             allowed_wording="the failure is present in sample; transferring a ranking across policies adds little regret", prohibited_wording="rankings do not transfer"),
        dict(claim="Weakness is contextual (8-layer evidence rule)", prereg_hypothesis="H19", old_evidence="GOLD-A",
             new_hardening_evidence=f"{h19.get('n_positive')}/8 layers positive incl. L2 and L5 (L8 cross-model negative)", discovery="n/a", old_holdout="n/a",
             new_holdout=f"{h19.get('n_positive')}/8", uncertainty="L7 positive", alternate_dynamic_model="L8 negative", full_mode="yes",
             electromechanical_only="L5 fails in EM stratum", same_mode="L4 positive", effect_size="see layers", statistic="frozen rule",
             final_status="SUPPORTED (rule met; downgraded from STRONGLY because L5 does not hold in the EM stratum)",
             allowed_wording="a context-independent node ranking does not suffice on this benchmark", prohibited_wording="nodes have no intrinsic weakness anywhere"),
        dict(claim="Portfolio-conditioned total branch sensitivity ranks finite reinforcements (GOLD-B)", prereg_hypothesis="H7 primary",
             old_evidence="old holdout median rho 0.987 vs 0.50", new_hardening_evidence=f"prereg-eligible {h7.get('GB_prereg_n')} cond: {f(h7.get('GB_prereg_median_rho_Dtotal'),3)} vs {f(h7.get('GB_prereg_median_rho_S1'),3)}; all 64: {f(h7.get('GB_primary_median_rho_Dtotal'),3)} vs {f(h7.get('GB_primary_median_rho_S1'),3)}",
             discovery="selection of S1", old_holdout="0.987 (old campaign)", new_holdout=f(h7.get("GB_prereg_median_rho_Dtotal"), 3),
             uncertainty=f"fresh draws {f(h7.get('GB_fresh_draws_median_rho_Dtotal'),3)} vs {f(h7.get('GB_fresh_draws_median_rho_S1'),3)}",
             alternate_dynamic_model=f"method transfers (C2 {h18.get('C2_verdict')}); ranking model-specific (C {h18.get('C_verdict')})",
             full_mode=f(h7.get("GB_primary_median_rho_Dtotal"), 3), electromechanical_only=f(h7.get("GB_EM_median_rho_Dtotal"), 3), same_mode=f(h7.get("GB_SAME_median_rho_Dtotal"), 3),
             effect_size=f"median paired advantage {f(h7.get('GB_prereg_median_diff'),3)} (CI {h7.get('GB_prereg_diff_ci95')})",
             statistic=f"T2 sign-flip p {f(h7.get('GB_prereg_T2_signflip_p'),3)} (6 clusters, not significant); new policies: frozen {f(h7.get('GB_new_policies_median_rho_Dfrozen'),3)}, learned fixed list {f(lo.get('rho_fixed_list'),3)} (post hoc); doublings {f(dt.get('dbl', {}).get('rho_full'),3)}, outages {f(dt.get('out', {}).get('rho_full'),3)}",
             final_status="STRONGLY_SUPPORTED (preregistered gate on 8 conditions; T2 not significant)",
             allowed_wording="the re-equilibrated sensitivity of the current portfolio, a known construction, ranks finite reinforcements, doublings and outages better than the frozen derivative, a fixed list learned from other conditions, and static indices",
             prohibited_wording="novel sensitivity formula; static metrics are useless; converter-independent ranking; computationally cheaper (not shown)"),
        dict(claim="Re-equilibration matters for network/operating-point coordinates, not controller coordinates", prereg_hypothesis="H8",
             old_evidence="R1 structural check", new_hardening_evidence=f"controller rel diff {h7.get('node_R1_max_rel_frozen_minus_total_controller')}; branch material {f(h7.get('H8_frac_material_H4'))} of conditions",
             discovery="n/a", old_holdout="n/a", new_holdout=f"frozen rho on new policies {f(h7.get('GB_new_policies_median_rho_Dfrozen'),3)} vs total {f(h7.get('GB_new_policies_median_rho_Dtotal'),3)}",
             uncertainty="n/a", alternate_dynamic_model="n/a", full_mode="yes", electromechanical_only="n/a", same_mode="n/a",
             effect_size=f"load indirect dominates {h7.get('node_opoint_frac_indirect_dominates', {}).get('load')}", statistic="descriptive", final_status="SUPPORTED",
             allowed_wording="under SPR, frozen = total for controller coordinates; the operating-point term can dominate for loads, setpoints and branches", prohibited_wording="first to include re-equilibration"),
        dict(claim="Robust weak corridor under equal budget and size-matched null", prereg_hypothesis="H10 primary",
             old_evidence="TXother/TXall top-3 under unnormalized x1.5 (confounded by size)", new_hardening_evidence=f"robust set: {rob}",
             discovery="see H11 table", old_holdout="see H10 percentiles", new_holdout=f"TXother frac>=0.95: {f(h10.get('TXother|A|new_frac_ge95'))}; TXall: {f(h10.get('TXall|A|new_frac_ge95'))}",
             uncertainty="fresh draws included", alternate_dynamic_model=f"corridor top-1 {h18.get('D_verdict')}", full_mode="yes",
             electromechanical_only=f"TXother EM frac {f(h10.get('TXother|A|new_EM_frac_ge95'))}", same_mode="n/a", effect_size="percentile vs size null",
             statistic="percentile", final_status=("SUPPORTED" if h10.get("H10_pass") else ("NEGATIVE" if h10 else "PENDING")),
             allowed_wording=("<corridor> is in the top 5% of equal-budget size-matched groups in X of held-out conditions" if h10.get("H10_pass") else "the unnormalized corridor ranking reflects corridor size; no robust weak corridor under the preregistered null"),
             prohibited_wording="weak corridor without equal-budget/null qualifier"),
        dict(claim="Controller-weighted modal mixing explains contextuality", prereg_hypothesis="H13 primary", old_evidence="E12 Q2 rho 0.657 (at omega_c(theta))",
             new_hardening_evidence=f"fixed omega: old {f(h13.get('tests', {}).get('old|mu10|n_rev', {}).get('rho'))}, new {f(h13.get('T3_primary', {}).get('rho'))} (p {f(h13.get('T3_primary', {}).get('p'),3)}, Holm {f(h19.get('F_conf_holm_p', {}).get('T3_mixing'),3)})",
             discovery="n/a", old_holdout=f(h13.get("tests", {}).get("old|mu10|n_rev", {}).get("rho")), new_holdout=f(h13.get("T3_primary", {}).get("rho")),
             uncertainty="n/a", alternate_dynamic_model="n/a", full_mode="n/a", electromechanical_only="n/a", same_mode="n/a",
             effect_size=f"variance terms: controller {f(mv.get('var_C_over_var_total'))}, frequency {f(mv.get('var_F_over_var_total'))}, 2cov {f(mv.get('2cov_over_var_total'))}", statistic="permutation + Holm (low power at n = 19)", final_status="NEGATIVE (not confirmed)",
             allowed_wording="the earlier mixing correlation is largely a frequency effect; the mixing interpretation is not confirmed", prohibited_wording="mixing causes reversal"),
        dict(claim="Topology alone removes and creates EM incompatibility; static scores do not predict", prereg_hypothesis="H14-H16 (Q13)",
             old_evidence="E7 (full alpha, fast modes included)", new_hardening_evidence=f"EM removal at P4 {h15.get('em_removal_H4_P4')}; {h15.get('n_em_creation')} EM creations; 0 fast-only",
             discovery="included", old_holdout="H01-H06 included", new_holdout="HARDENING_H01-H06 included", uncertainty="n/a",
             alternate_dynamic_model=f"topology direction {h18.get('E_verdict')} ({f(h18.get('E_pooled_sign_agree'))})", full_mode="yes", electromechanical_only="yes (Q13 YES)",
             same_mode="tracked EM mode", effect_size="see H14", statistic=f"no static score predicts (max frac {max(v['frac_abs_ge_0.6_em'] for v in h15.get('H16', {}).values()) if h15 else 'n/a'})",
             final_status="SUPPORTED (supporting observation: all 45 creations are outages, 42 at discovery or old-holdout policies)", allowed_wording="single admissible branch actions both remove and create EM incompatibilities, as expected of N-1 small-signal behavior; seven static scores fail the preregistered rule",
             prohibited_wording="topology switching is novel; static scores are useless"),
        dict(claim="Contextual reversal transfers to a validated library GFL", prereg_hypothesis="H18-A", old_evidence="E23 static injection 7/7 (not dynamic)",
             new_hardening_evidence=f"ALT-WECC {h18.get('A_frac')} of {h18.get('n_tested_policies')}", discovery="n/a", old_holdout="4 policies", new_holdout="4 policies",
             uncertainty="n/a", alternate_dynamic_model=h18.get("A_verdict"), full_mode="0/7", electromechanical_only="0/7", same_mode="0/7",
             effect_size="global alpha pinned at -0.100 s^-1 by an isolated REPCA1 state (test uninformative); EM-tracked (post hoc): stabilizing at 3/7, nested reversal 1/7; QFLAG=1: 0/7", statistic="descriptive",
             final_status="NEGATIVE (preregistered test uninformative; exploratory tracked test 1/7)",
             allowed_wording="the reversals were not reproduced by the frozen WECC library chain; switching on its voltage control did not produce them either", prohibited_wording="model-independent reversal; converter voltage control enables reversals"),
        dict(claim="Fixed-ranking regret arises at singularity-induced-bifurcation-type crossings", prereg_hypothesis="none (post hoc, reviewer-requested)", old_evidence="n/a",
             new_hardening_evidence=f"{fm.get('frac_fast_real', float('nan')):.1%} of fast critical eigenvalues real; alpha quantiles {fm.get('fast_alpha_quantiles')} s^-1; six sweeps jump within one 0.025 step where sigma_min(g_z) is 0.002-0.007",
             discovery="n/a", old_holdout="n/a", new_holdout="yes", uncertainty="n/a", alternate_dynamic_model="n/a", full_mode="yes", electromechanical_only="n/a", same_mode="n/a",
             effect_size=f"participation theta_PLL {f(ex.get('fast_participation_mean_top', {}).get('gfl:theta'))}, i_q {f(ex.get('fast_participation_mean_top', {}).get('gfl:i_q'))}", statistic="descriptive",
             final_status="EXPLORATORY", allowed_wording="the events carry the signature of a singularity-induced bifurcation of the phasor DAE; their physical course needs network dynamics or EMT",
             prohibited_wording="physical instability of real converters; EMT-validated"),
        dict(claim="GFM and Africano/PV generalization", prereg_hypothesis="scope", old_evidence="E17 BLOCKED", new_hardening_evidence="no validated GFM in repository",
             discovery="n/a", old_holdout="n/a", new_holdout="n/a", uncertainty="n/a", alternate_dynamic_model="n/a", full_mode="n/a", electromechanical_only="n/a",
             same_mode="n/a", effect_size="n/a", statistic="n/a", final_status="BLOCKED", allowed_wording="open", prohibited_wording="any GFM or Africano/PV claim"),
    ]
    df = pd.DataFrame(rows)
    df.columns = [c.replace("_", " ") for c in df.columns]
    df.to_csv(TOP / "CDW_HARDENED_CLAIM_MATRIX.csv", index=False)
    return df


def master_results():
    rows = []
    for name in ("H03_gate.json", "H04_gate.json", "H05_stats.json", "H07_goldb_gate.json", "H10_gate.json", "H13_mixing.json", "H15_topology_gate.json",
                 "H18_gate.json", "H19_evidence.json", "CDWH_DETERMINISM.json"):
        d = j(name)

        def walk(prefix, v):
            if isinstance(v, dict):
                for k, x in v.items():
                    walk(f"{prefix}.{k}", x)
            elif isinstance(v, list) and v and isinstance(v[0], (dict, list)):
                rows.append({"source": name, "key": prefix, "value": json.dumps(v, default=str)[:500]})
            else:
                rows.append({"source": name, "key": prefix, "value": json.dumps(v, default=str)[:500]})
        for k, v in d.items():
            walk(k, v)
    pd.DataFrame(rows).to_csv(TOP / "CDW_HARDENING_MASTER_RESULTS.csv", index=False)
    return len(rows)


def run_manifest():
    st = json.loads(HI.STATUS.read_text())
    rows = []
    for ph, v in st["phases"].items():
        store = {"H01": "H_H01", "H01e": "H_H01e", "H06": "H_H06", "H08n": "H_H08n", "H09": "H_H09", "H12": "H_H12", "H12b": "H_H12b",
                 "H14": "H_H14", "H18c": "H_H18c", "H10": "H_H10"}.get(ph)
        n = len(list((I.RAW / store).glob("*.json"))) if store else 0
        rows.append({"phase": ph, "state": v.get("state"), "wall_s": v.get("wall_s"), "n_tasks": (v.get("tasks") or {}).get("n_tasks"),
                     "n_errors": v.get("n_errors"), "raw_files": n, "store": f"raw/{store}", "updated": v.get("updated")})
    alt = list((R / "alt" / "cases").glob("*.json"))
    rows.append({"phase": "H17/H18 ALT-WECC", "state": "COMPLETE", "wall_s": None, "n_tasks": len(alt), "n_errors": 0, "raw_files": len(alt),
                 "store": "results/hardening/alt/cases", "updated": None})
    man = pd.DataFrame(rows)
    man["prereg_commit"] = HI.PREREG_COMMIT
    man["environment"] = json.dumps({k: st.get("environment", {}).get(k) for k in ("python", "numpy", "scipy", "pandas")})
    man.to_csv(TOP / "CDW_HARDENING_RUN_MANIFEST.csv", index=False)
    ins = [{"file": p.name, "sha256": hashlib.sha256(p.read_bytes()).hexdigest()} for p in sorted((R / "prereg_inputs").glob("*.json"))]
    pd.DataFrame(ins).to_csv(R / "prereg_inputs_sha256.csv", index=False)
    return man


if __name__ == "__main__":
    print(claim_matrix()[["claim", "final status"]].to_string())
    print(master_results(), "master rows")
    print(run_manifest().to_string())
