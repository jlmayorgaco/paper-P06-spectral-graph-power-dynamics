# CDW hardening campaign: final report

Branch `research/contextual-dynamic-weakness-hardening` (from `e73dd355`). Preregistration `05b507e3`. No push. TX4 untouched. Every number below is read from `results/hardening/*.json` by `experiments/cdw_hardening/CDWH_REPORT.py`.

## 1. Verdict in one paragraph

The contextual core of CDW survives a new preregistered holdout. Nested electromechanical same-mode reversals occur in 0.895 of base-stable new-holdout policies (17/19; gate 0.75 PASS). GOLD-B survives with a large margin: on the 8 preregistration-eligible held-out conditions the total re-equilibrated sensitivity ranks finite reinforcements at median Spearman 0.985 against 0.491 for |P| (lower 95% bound of the paired advantage 0.477 > 0.20); over all 64 solved conditions 0.987 vs 0.497. Three prior readings are weakened: (i) the fixed-ranking failure (p* 0.75, regret 0.38) is carried by jumps onto fast converter-control modes (2089/2265 material regrets); for EM transitions a fixed ranking orders 0.97 of pairs; (ii) the modal-mixing mechanism is largely a frequency artefact (controller share of variance 0.19; new-holdout fixed-frequency test fails); (iii) the reversals and the specific branch ranking do not transfer to a frozen WECC library GFL chain, although the ranking method does. Topology effects survive the EM-only reclassification. Corridors: robust weak corridor(s) under the equal-budget size-matched null: none.

## 2. Phase status

| phase | status | key file |
|---|---|---|
| H0 preregistration | committed `05b507e3` | `docs/CDW_HARDENING_PREREG_V1.md` |
| H1 new policy holdout (24 maximin LHS) | COMPLETE, 19 base-stable | `results/hardening/prereg_inputs/hardening_policies.json` |
| H2 curvature theorem | PROVED (classical algebra) | `theory/CDW_CONTEXTUAL_CURVATURE_THEOREM.md` |
| H3 nested reversal | PASS (0.895) | `results/hardening/H03_gate.json` |
| H4 fixed ranking + transfer | PASS (FULL); fails in EM stratum | `results/hardening/H04_gate.json` |
| H5 GOLD-A statistics | COMPLETE | `results/hardening/H05_stats.json` |
| H6/H7 GOLD-B new holdout | PASS (literal-eligible set and full set) | `results/hardening/H07_goldb_gate.json` |
| H8 frozen vs total | COMPLETE | `results/hardening/H08_*` |
| H9-H11 corridors | PENDING | `results/hardening/H10_gate.json, H11_corridor_table.csv` |
| H12/H13 mixing | NEGATIVE (primary fails) | `results/hardening/H13_mixing.json` |
| H14-H16 topology EM | SUPPORTED (Q13 yes; static fails) | `results/hardening/H15_topology_gate.json` |
| H17 ALT-WECC qualification | PASS Q1-Q4 (with logged refinements) | `results/hardening/alt/H17_qualification_andes.json` |
| H18 cross-model matrix | A MODEL-SPECIFIC, C MODEL-SPECIFIC, C2 TRANSFERS, D TRANSFERS, E TRANSFERS | `results/hardening/H18_gate.json` |
| H19 evidence layers | 7/8 positive (L8 negative) | `results/hardening/H19_evidence.json` |
| H20 novelty boundary | COMPLETE | `docs/CDW_NOVELTY_BOUNDARY.md` |
| H21 literature | 57 verified entries | `results/CDW_LITERATURE_GAP_MATRIX.csv` |
| H22 reviewer attack | COMPLETE | `docs/CDW_REVIEWER_ATTACK.md` |
| H23-H30 paper + supplement | COMPLETE | `reports/papers/cdw_contextual_dynamic_weakness/` |
| H31 internal reviews | see docs | `docs/CDW_REVIEWER1_FINAL.md, CDW_REVIEWER2_FINAL.md` |
| determinism | identical | `results/hardening/CDWH_DETERMINISM.json` |

## 3. Results by question

**Q1. Does nested contextual sign reversal survive the new holdout?** Yes. Level A (global) 1.00, level B 0.95, level C 0.895 of base-stable new policies (old holdout: 0.889).

**Q2. Does same-mode electromechanical reversal survive?** Yes: 17/19 at level C; every such policy also has a level-D (same family across contexts) reversal. Magnitudes are small (quartiles [0.01127, 0.0127, 0.01653] s^-1). Both directions within a policy in 0.37 of policies.

**Q3. Is the optimal fixed node ranking still materially insufficient?** FULL: yes (p* 0.75, regret 0.38; 19/19 policies p*<0.90). EM-only: no (p* 0.97, regret 0.03). 2089 of 2265 material regrets are jumps out of the EM band.

**Q4. How large is the old-to-new holdout shift?** Negligible: p* shift -0.002 (interval -0.036 to 0.041); regret shift 0.007; level-C coverage 0.889 old vs 0.895 new.

**Q5. Does GOLD-A remain strong enough for a headline claim?** Yes for contextual reversal (nested, EM same-mode), with the magnitude caveat. The fixed-ranking half must be stated with its fast-mode anatomy.

**Q6. Does total dynamic line sensitivity still beat every static baseline on the new holdout?** Yes: Dtotal exceeds the per-condition best static index in 1.00 of conditions; oracle static median 0.499.

**Q7. Is the improvement statistically and practically material?** Yes: median paired advantage 0.492 (CI 0.477-0.608); top-5 precision 1.00 vs 0.20. The sign-flip p is 0.061 on the literal 8-condition set (6 clusters), <1e-3 on all 64.

**Q8. When does re-equilibration materially change the line ranking?** In 0.28 of V4 conditions by the old H2 rule; every condition has at least one frozen/total sign flip; the frozen derivative drops to 0.712 on the 24 new policies while Dtotal stays 0.996. Controller coordinates: frozen = total (2.2e-07).

**Q9. Does TXother remain a strong corridor under equal intervention budget?** PENDING

**Q10. Does it beat size-matched random connected corridors?** PENDING

**Q11. How much of old modal-mixing variation was controller- vs frequency-induced?** Controller term 0.19 of the variance (median |C| share 0.33); the rest is frequency movement.

**Q12. Does fixed-frequency modal mixing still correlate with contextuality?** No at the preregistered bar: new rho 0.55 (p 0.016, Holm 0.033); old holdout at fixed frequency rho 0.21.

**Q13. After filtering to EM modes, can topology still both remove and create incompatibility?** Yes: removal at P4 by ['dbl13', 'dbl43', 'dbl1', 'dbl0']; 45 EM creations of hyperedges of size <= 3; 0 fast-only.

**Q14. Do Fiedler/gSCR/static topology scores remain poor predictors of the EM effect?** Yes: no score meets |rho| >= 0.6 in >= 75% of policies (best fraction 0.38); medians around |0.5|, i.e. partial information.

**Q15. Does contextual reversal transfer to a genuinely dynamic alternative converter model?** No: ALT-WECC 0.0 of 7 tested policies (custom GFL 0.71); every stable-context ALT marginal is >= 0.

**Q16. Does dynamic link ranking transfer?** The ranking does not (median Kendall 0.02); the method does (ALT total vs ALT finite rho 0.99 vs |P| 0.79).

**Q17. Which findings are model-specific?** The reversals, the fixed-ranking failure (driven by the custom GFL's fast modes), the specific branch ranking. Transferring: the ranking method, the top corridor (0.75), topology direction (0.94).

**Q18. What exact novelty remains after the literature review?** The combination in docs/CDW_NOVELTY_BOUNDARY.md section 2: certified nested EM same-mode reversal of SG->GFL replacements with curvature witnesses; quantified oracle fixed-ranking failure and its fast-mode anatomy; held-out validation of portfolio-conditioned total sensitivity as a reinforcement ranking against static indices; EM-only topology effects; a frozen library-GFL cross-model boundary. No single ingredient is new.

**Q19. What are the 3 strongest defensible paper contributions?** (1) portfolio-conditioned re-equilibrated reinforcement ranking (GOLD-B, STRONGLY_SUPPORTED); (2) nested EM same-mode contextual reversal with curvature certificates (SUPPORTED); (3) fixed-ranking insufficiency with its modal anatomy (SUPPORTED, FULL only).

**Q20. What claims must be dropped or weakened?** Dropped: modal-mixing mechanism; model-independent reversal; converter-independent branch ranking; weak corridors unless H10 passes. Weakened: fixed-ranking failure to 'carried by fast converter-control modes'; 'weakness is contextual' to SUPPORTED (not STRONGLY).

**Q21. Is the paper ready for TPWRS?** see reviews

**Q22. If not, what single missing experiment blocks submission?** see reviews

## 4. Final claim matrix

| claim | final status |
|---|---|
| Chain identity certifies nested reversals by curvature of the matching sign | PROVED |
| Nested EM same-mode contextual sign reversal | SUPPORTED |
| alpha_perp is neither submodular nor supermodular (both nested directions) | SUPPORTED |
| Oracle fixed node ranking materially insufficient | SUPPORTED (FULL only; carried by fast modes) |
| Fixed rankings do not transfer across policies | SUPPORTED |
| Weakness is contextual (8-layer evidence rule) | SUPPORTED (rule met; downgraded from STRONGLY because L5 does not hold in the EM stratum) |
| Portfolio-conditioned total branch sensitivity ranks finite reinforcements (GOLD-B) | STRONGLY_SUPPORTED |
| Re-equilibration matters for network/operating-point coordinates, not controller coordinates | SUPPORTED |
| Robust weak corridor under equal budget and size-matched null | PENDING |
| Controller-weighted modal mixing explains contextuality | NEGATIVE |
| Topology alone removes and creates EM incompatibility; static scores do not predict | SUPPORTED |
| Contextual reversal transfers to a validated library GFL | NEGATIVE |
| GFM and Africano/PV generalization | BLOCKED |

## 5. Negative and blocked results (same prominence)

- Fixed-ranking insufficiency does not hold for EM-only transitions (A3 bar failed in the EM stratum).
- Fixed-frequency modal mixing does not correlate with contextuality at the preregistered bar; the earlier E12 reading was mostly a frequency effect.
- Contextual reversal does not transfer to the frozen WECC library GFL chain (0/7); the specific branch ranking is converter-model-specific.
- Corridors: pending.
- GFM generalization and the Africano/PV benchmark remain BLOCKED (no validated GFM model; no source material).
- The preregistered R0 <= 1e-8 eligibility rule was stricter than the solver tolerance; the literal primary GOLD-B set has only 8 conditions (still passes).

## 6. Deviations

See `docs/CDW_HARDENING_DEVIATIONS.md`: execution fixes (relaunch, import order, bound-check tolerance, H06 pool crash rerun), H17 structural-zero and limiter refinements recorded before the cross-model matrix, one ALT datum seen during qualification, and the literal application of the H6 eligibility rule. No gate threshold changed.

## 7. Run facts

- Compute phases: 11; tasks: 2981; summed phase wall time 0.6 h (parallel). Determinism: identical.
- Commits on this branch: `eb57d907 results(cdw-hardening): GOLD-A/B hardening, mixing, topology EM, ALT-WECC cross-model, literature`; `4de412d2 code(cdw-hardening): hardening modules, curvature theory note, deviation log (before H18 matrix)`; `05b507e3 prereg(cdw-hardening): H0 hardening preregistration, before any new numerics`.
- Literature: 57 verified entries; privacy note: one literature agent reported once including the user's e-mail as a Crossref `mailto` parameter; no other personal data was sent.
