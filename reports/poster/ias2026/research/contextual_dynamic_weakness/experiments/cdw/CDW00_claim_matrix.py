"""CDW00: write the initial claim matrix results/CDW_CLAIM_MATRIX_V1.csv (no computation)."""

import csv
from pathlib import Path

OUT = Path(__file__).resolve().parents[2] / "results" / "CDW_CLAIM_MATRIX_V1.csv"
COLS = ["claim_id", "claim", "status", "source_or_test", "tx4_link", "not_claimed"]
ROWS = [
    # inherited proved
    ("I-P01", "Exact transverse quotient; det(sI-A) = s^2 det(sI-A_perp) under rotation covariance, D=0, constant Pm", "INHERITED_PROVED", "TX4 Thm 1", "V01", "validity with governors or D != 0"),
    ("I-P02", "H is an antichain; any-order safety iff no hyperedge in T (stable base)", "INHERITED_PROVED", "TX4 Prop 1-2", "V02", ""),
    ("I-P03", "Boundary localization: H changes only at imaginary-axis crossings", "INHERITED_PROVED", "TX4 Thm 2", "V02", ""),
    ("I-P04", "Descriptor-affine local actions; network-closure factorization; boundaries are zeros of det(I+Q_SS)", "INHERITED_PROVED", "TX4 Lemma 1, Thm 3, Cor 1", "V08", "which order closes a given boundary"),
    ("I-P05", "Boundary sensitivity ds*/da (simple eigenvalue), port form = eigenvalue form", "INHERITED_PROVED", "TX4 Prop 3", "V10", "finite-step prediction"),
    ("I-P06", "Symmetry-deflated zero-frequency port", "INHERITED_PROVED", "TX4 Thm 4", "T08", ""),
    ("I-P07", "Minimum-cardinality destabilization is NP-complete", "INHERITED_PROVED", "TX4", "V26", ""),
    ("I-P08", "Connected-cumulant identities (secondary)", "INHERITED_PROVED", "TX4", "V22", "any causal or bridge meaning"),
    # inherited validated (benchmark-specific)
    ("I-V01", "P4: all proper subsets of {30,33,35,37} stable, H4 unstable (alpha +0.127)", "INHERITED_VALIDATED", "TX4 PCV02/G2", "V03", "a general weak-bus set"),
    ("I-V02", "Policy alone reorganizes H at fixed network/dispatch; static metrics cannot encode it", "INHERITED_VALIDATED", "TX4", "V06,V07", "static metrics wrong in general"),
    ("I-V03", "Port derivative ranks 12 lines' finite 1.5x effect (rho 0.937) at one condition", "INHERITED_VALIDATED", "TX4 PCV04", "V11", "ranking at other conditions"),
    ("I-V04", "Non-composability persists across 3/4 envelopes; witness identity not robust", "INHERITED_VALIDATED", "TX4 PCV05", "V15,V16", "robust contextuality"),
    ("I-V05", "Frozen-equilibrium line derivative rho 0.615 vs re-equilibrated 1.000 (12 lines, one condition)", "INHERITED_VALIDATED", "dynamic-forest report", "", "frozen/total gap elsewhere"),
    ("I-V06", "Oscillatory boundaries subcritical Hopf; small-signal verdicts confirmed by phasor TDS", "INHERITED_VALIDATED", "TX4", "V30", "EMT validity"),
    # inherited refuted
    ("I-R01", "Static quantities alone identify failing coalitions", "REFUTED", "TX4", "V27", ""),
    ("I-R02", "GFL branch ranking transfers across converter models", "REFUTED", "TX4", "V28", ""),
    ("I-R03", "kappa=4 coalition robust to primary frequency control", "REFUTED", "TX4", "V29", ""),
    ("I-R04", "Simple-cycle / holonomy / dominant path / effective rank explains failures", "REFUTED", "TX4", "N01,N05,N06", ""),
    ("I-R05", "|chi| suppression as design target", "REFUTED", "TX4", "V24", ""),
    # new hypotheses
    ("C-H1", "Contextual sign reversal recurs across holdout policies (stable contexts)", "HYPOTHESIS", "E1 gate A1", "", "universal weak buses"),
    ("C-H1t", "Reversal survives mode tracking", "HYPOTHESIS", "E1 gate A2", "", ""),
    ("C-H1b", "Every fixed node-only ranking is materially insufficient for next-replacement decisions", "HYPOTHESIS", "E1 gate A3", "", "that rankings are useless"),
    ("C-H1c", "Contextuality persists under envelopes even where the witness changes", "HYPOTHESIS", "E2 gate A4", "V15,V16", "probability"),
    ("C-HS", "alpha_perp is neither submodular nor supermodular on the tested class", "HYPOTHESIS", "E1b", "", "greedy always fails"),
    ("C-IV", "IFT total derivatives match finite re-equilibrated differences", "NOT_YET_TESTED", "E3/E4 gate IV", "", ""),
    ("C-R1", "Frozen = total for dynamic-only controls under SPR; frozen = 0 for setpoints (structural)", "NOT_YET_TESTED", "E3 structural check", "", "empirical generality"),
    ("C-H2n", "Frozen and total node sensitivities disagree materially", "HYPOTHESIS", "E3 H2-node", "", ""),
    ("C-H2l", "Frozen and total link sensitivities disagree materially", "HYPOTHESIS", "E4 H2-link", "I-V05", ""),
    ("C-GB", "Total dynamic link measure beats best static baseline on held-out finite effects", "HYPOTHESIS", "E4 GOLD-B", "V11", "universal weak lines"),
    ("C-H3", "Robust dynamic weak corridors exist", "HYPOTHESIS", "E6", "", "nominal top line = corridor"),
    ("C-H4", "Topology alone can remove/create incompatibilities; static topology scores do not predict", "HYPOTHESIS", "E7", "V13", "cost-effectiveness"),
    ("C-EX", "Local stability-equivalent exchange rates are first-order valid", "HYPOTHESIS", "E8", "", "economic exchange rate"),
    ("C-H5", "Plan-level design succeeds where single-boundary tuning leaves an unsafe subset", "HYPOTHESIS", "E9 GOLD-D", "", "global impossibility"),
    ("C-SH", "Context-averaged attribution (Shapley) hides large sign reversals", "HYPOTHESIS", "E10", "", "Shapley causality"),
    ("C-H6s", "Static graph scores predict contextual sign / total sensitivity", "HYPOTHESIS", "E11", "", "novelty of Laplacian decomposition"),
    ("C-H6m", "Policy changes modal mixing at fixed L; mixing correlates with contextuality", "HYPOTHESIS", "E12", "", "visual correlation"),
    ("C-H6r", "A few graph/dynamic modes preserve decisions (spectral closure)", "HYPOTHESIS", "E13", "N08", ""),
    ("C-H7", "A genuinely reduced model preserves decisions with speedup and certificate/abstention", "HYPOTHESIS", "E13/E14 GOLD-C", "", "rigor of sampled certificates"),
    ("C-H8", "Dominant observed modal energy is often not the limiting mode (phasor)", "HYPOTHESIS", "E16", "SQ1 observation", "EMT result"),
    ("C-H9", "Findings transfer to an alternative converter model", "HYPOTHESIS", "E23", "V28", "universal transfer"),
    ("C-AF", "Africano replication and PV hosting (GOLD-E/F)", "NOT_YET_TESTED", "E17 gate", "", "any replication without exact material"),
    ("C-NL", "Plan-level tuning improving alpha may reduce nonlinear recovery margin", "NOT_YET_TESTED", "E22 optional", "V30", "alpha improvement implies larger region of attraction"),
]


def main():
    with OUT.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(COLS)
        w.writerows(ROWS)
    print(len(ROWS), "claims")


if __name__ == "__main__":
    main()
