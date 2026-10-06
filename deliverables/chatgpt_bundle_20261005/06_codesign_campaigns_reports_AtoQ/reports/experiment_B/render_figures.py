"""Render the reproducible Experiment B-pre figures from generated CSV tables."""

from __future__ import annotations

import csv
import math
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(sys.argv[1]).resolve()
TABLES = ROOT / "reports" / "experiment_B" / "tables"
FIGURES = ROOT / "reports" / "experiment_B" / "figures"
FIGURES.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.titlesize": 13,
    "axes.labelsize": 11,
    "figure.dpi": 120,
    "savefig.dpi": 180,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.23,
})
COLORS = ["#176B87", "#D97706", "#5B8E7D", "#9B5DE5", "#D1495B", "#4C566A"]


def read_table(name: str):
    with (TABLES / name).open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def floats(rows, key):
    out = []
    for r in rows:
        try:
            v = float(r[key])
        except (TypeError, ValueError, KeyError):
            v = math.nan
        out.append(v)
    return np.asarray(out, dtype=float)


def save(fig, stem):
    fig.tight_layout()
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(FIGURES / f"{stem}.{suffix}", bbox_inches="tight")
    plt.close(fig)


def line_plot(stem, title, xlabel, ylabel, series, xlog=False, ylog=False,
              scatter=False, zero_line=False):
    fig, ax = plt.subplots(figsize=(8.4, 5.1))
    for i, (label, x, y, style) in enumerate(series):
        color = COLORS[i % len(COLORS)]
        good = np.isfinite(x) & np.isfinite(y)
        if scatter:
            ax.scatter(x[good], y[good], label=label, color=color, s=35, alpha=0.82)
        else:
            ax.plot(x[good], y[good], label=label, color=color, linewidth=2.0,
                    marker="o" if len(x) < 20 else None, markersize=4,
                    linestyle=style)
    if xlog:
        ax.set_xscale("log")
    if ylog:
        ax.set_yscale("log")
    if zero_line:
        ax.axhline(0.0, color="#333333", linewidth=0.8)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend(frameon=False, ncol=2 if len(series) > 3 else 1)
    save(fig, stem)


def heatmap(stem, title, field, cmap, label):
    rows = read_table("TABLE_B11_sigmahat_representative.csv")
    n = max(int(r["row"]) for r in rows)
    mat = np.full((n, n), np.nan)
    for r in rows:
        mat[int(r["row"]) - 1, int(r["col"]) - 1] = float(r[field])
    fig, ax = plt.subplots(figsize=(6.3, 5.5))
    im = ax.imshow(mat, cmap=cmap, interpolation="nearest", aspect="equal")
    ax.set_xticks(range(n), [str(i) for i in range(1, n + 1)])
    ax.set_yticks(range(n), [str(i) for i in range(1, n + 1)])
    ax.set_xlabel("Graph mode ℓ")
    ax.set_ylabel("Graph mode k")
    ax.set_title(title)
    fig.colorbar(im, ax=ax, label=label, shrink=0.82)
    save(fig, stem)


def main():
    spectrum = read_table("TABLE_B12_graph_spectrum.csv")
    s2_spec = [r for r in spectrum if r["case"] == "S2"]
    k = np.asarray([int(r["mode"]) for r in s2_spec])
    nu = floats(s2_spec, "nu")
    zero = np.asarray([r["is_zero"].lower() == "true" for r in s2_spec])
    fig, ax = plt.subplots(figsize=(8.4, 5.1))
    ax.scatter(k[~zero], nu[~zero], color=COLORS[0], s=54, label="Oscillatory graph modes")
    if zero.any():
        ax.scatter(k[zero], nu[zero], color=COLORS[4], marker="x", s=75, label="Common-angle zero mode")
    ax.set(title="S2 graph eigenvalue spectrum", xlabel="Mode index k", ylabel="Graph eigenvalue νₖ [s⁻²]")
    ax.legend(frameon=False)
    save(fig, "FIG_B01_graph_eigenvalue_spectrum")

    heatmap("FIG_B02_sigmahat_magnitude_heatmap",
            "S2 |Σ̂(jω)| at f = 2.1 Hz", "sigma_abs", "magma", "|Σ̂ₖℓ|")
    heatmap("FIG_B03_sigmahat_real_heatmap",
            "S2 Re{Σ̂(jω)} at f = 2.1 Hz", "sigma_real", "RdBu_r", "Re{Σ̂ₖℓ}")

    damping = read_table("TABLE_B13_graph_damping.csv")
    fig, ax = plt.subplots(figsize=(8.4, 5.1))
    for i, f in enumerate((0.2, 1.0, 5.0, 20.0)):
        subset = [r for r in damping if math.isclose(float(r["frequency_hz"]), f)]
        ax.plot(floats(subset, "nu"), floats(subset, "d_k"), marker="o",
                linewidth=1.8, color=COLORS[i], label=f"f = {f:g} Hz")
    ax.set(title="S2 graph dissipative diagonal", xlabel="Graph eigenvalue νₖ [s⁻²]",
           ylabel="dₖ(f) [self-energy units]")
    ax.legend(frameon=False)
    save(fig, "FIG_B04_graph_damping_spectrum")

    freq_metrics = read_table("TABLE_B03_sigma_frequency_metrics.csv")
    s2_metrics = [r for r in freq_metrics if r["case"] == "S2"]
    line_plot("FIG_B05_hermitian_minmax_vs_frequency", "S2 Hermitian dissipative spectrum",
        "Frequency f [Hz]", "Eigenvalue of D_G [self-energy units]",
        [("λmin(D_G)", floats(s2_metrics,"frequency_hz"), floats(s2_metrics,"DG_lambda_min"), "-"),
         ("λmax(D_G)", floats(s2_metrics,"frequency_hz"), floats(s2_metrics,"DG_lambda_max"), "-")], xlog=True, zero_line=True)
    cases = ["S1", "S2", "S3", "S4"]
    line_plot("FIG_B06_offdiag_ratio_vs_frequency", "Off-diagonal graph self-energy ratio",
        "Frequency f [Hz]", "r_OD(f)",
        [(c, floats([r for r in freq_metrics if r["case"]==c],"frequency_hz"),
          floats([r for r in freq_metrics if r["case"]==c],"offdiag_ratio"), "-") for c in cases], xlog=True)
    line_plot("FIG_B07_commutator_vs_frequency", "Graph stiffness/self-energy mismatch",
        "Frequency f [Hz]", "χ_comm(f)",
        [(c, floats([r for r in freq_metrics if r["case"]==c],"frequency_hz"),
          floats([r for r in freq_metrics if r["case"]==c],"chi_comm"), "-") for c in ("S0","S1","S2","S3")], xlog=True)

    gamma = read_table("TABLE_B14_gamma_frequency.csv")
    gamma = [r for r in gamma if r["case"] == "S2"]
    line_plot("FIG_B08_gamma_real_imag_vs_frequency", "S2 collective self-energy, mode k = 3",
        "Frequency f [Hz]", "Γ₃(jω) [operator units, s⁻² for normalized M]",
        [("Re Γ₃", floats(gamma,"frequency_hz"), floats(gamma,"gamma_real"), "-"),
         ("Im Γ₃", floats(gamma,"frequency_hz"), floats(gamma,"gamma_imag"), "--")], xlog=True, zero_line=True)

    scaling = read_table("TABLE_B06_gamma_scaling.csv")
    sweep = read_table("TABLE_B05_coupling_sweep.csv")
    series = []
    for i, mode in enumerate(sorted({int(r["mode_k"]) for r in sweep})):
        sr = sorted([r for r in sweep if int(r["mode_k"]) == mode], key=lambda r: float(r["epsilon"]))
        fit = next((r for r in scaling if int(r["mode_k"]) == mode), None)
        label = f"mode {mode}" + (f", slope {float(fit['fitted_loglog_slope']):.2f}" if fit else "")
        series.append((label, floats(sr,"epsilon"), floats(sr,"gamma_abs"), "-"))
    line_plot("FIG_B09_gamma_scaling_epsilon", "S3 quadratic coupling check",
        "Coupling scale ε [-]", "|Γₖ(sₖ⁰)|", series, xlog=True, ylog=True)

    mode = min(int(r["mode_k"]) for r in sweep)
    sr = sorted([r for r in sweep if int(r["mode_k"]) == mode], key=lambda r: float(r["epsilon"]))
    fig, ax = plt.subplots(figsize=(7.4, 6.0))
    ax.plot(floats(sr,"uncoupled_pole_real"), floats(sr,"uncoupled_pole_imag"), "k--", label="Uncoupled pole")
    ax.plot(floats(sr,"predicted_pole_real"), floats(sr,"predicted_pole_imag"), "o-", color=COLORS[1], label="Perturbative prediction")
    ax.plot(floats(sr,"exact_pole_real"), floats(sr,"exact_pole_imag"), "s-", color=COLORS[0], label="Augmented-system pole")
    ax.set(title=f"S3 pole locus, graph mode k = {mode}", xlabel="Re{s} [s⁻¹]", ylabel="Im{s} [rad/s]")
    ax.legend(frameon=False)
    save(fig, "FIG_B10_exact_vs_predicted_pole_locus")

    errors = [("mode " + str(m), floats(sorted([r for r in sweep if int(r["mode_k"])==m], key=lambda r:float(r["epsilon"])),"epsilon"),
               floats(sorted([r for r in sweep if int(r["mode_k"])==m], key=lambda r:float(r["epsilon"])),"predictor_rel_error"), "-")
              for m in sorted({int(r["mode_k"]) for r in sweep})]
    line_plot("FIG_B11_predictor_error_vs_epsilon", "S3 perturbative pole prediction error",
        "Coupling scale ε [-]", "Relative pole error [-]", errors, xlog=True, ylog=True)
    line_plot("FIG_B12_predictor_error_vs_offdiag_ratio", "Predictor error versus r_OD",
        "r_OD(sₖ⁰) [-]", "Relative pole error [-]",
        [(f"mode {m}", floats([r for r in sweep if int(r["mode_k"])==m],"offdiag_ratio"),
          floats([r for r in sweep if int(r["mode_k"])==m],"predictor_rel_error"), "o")
         for m in sorted({int(r["mode_k"]) for r in sweep})], scatter=True)
    line_plot("FIG_B13_predictor_error_vs_commutator", "Predictor error versus χ_comm",
        "χ_comm(sₖ⁰) [-]", "Relative pole error [-]",
        [(f"mode {m}", floats([r for r in sweep if int(r["mode_k"])==m],"chi_comm"),
          floats([r for r in sweep if int(r["mode_k"])==m],"predictor_rel_error"), "o")
         for m in sorted({int(r["mode_k"]) for r in sweep})], scatter=True)

    resonance = read_table("TABLE_B08_resonance_sweep.csv")
    resonance.sort(key=lambda r: float(r["detuning_parameter"]))
    line_plot("FIG_B14_near_resonance_amplification", "S4 complementary-mode resonance",
        "Frequency detuning Δf [Hz]", "|Γ₂| [operator units, s⁻² for normalized M]",
        [("|Γ₂|", floats(resonance,"detuning_parameter"), floats(resonance,"gamma_abs"), "-")],
        xlog=True, ylog=True)

    pair = read_table("TABLE_B07_pairwise_contributions.csv")
    pair.sort(key=lambda r: int(r["rank"]))
    fig, ax = plt.subplots(figsize=(8.4, 5.1))
    ax.bar([f"ℓ={r['mode_l']}" for r in pair], floats(pair,"gamma_pair_abs"), color=COLORS[0])
    ax.set(title="S3 leading pairwise Γ contributions", xlabel="Complementary graph mode",
           ylabel="|γₖ←ℓ⁽²⁾| ε²")
    save(fig, "FIG_B15_pairwise_gamma_contributions")

    chi_freq = read_table("TABLE_B15_chiG_frequency.csv")
    line_plot("FIG_B16_candidate_small_gain_metric", "S2 finite-grid candidate χ_G diagnostic",
        "Frequency f [Hz]", "χ_G(f) [-]",
        [("S2", floats(chi_freq,"frequency_hz"), floats(chi_freq,"chiG"), "-")], xlog=True)

    chi_stab = read_table("TABLE_B09_chiG_stability.csv")
    fig, ax = plt.subplots(figsize=(8.4, 5.1))
    for i, case in enumerate(sorted({r["case"] for r in chi_stab})):
        subset = [r for r in chi_stab if r["case"] == case]
        ax.scatter(floats(subset,"chiG_sup"), floats(subset,"spectral_abscissa"),
                   label=case, color=COLORS[i % len(COLORS)], s=42, alpha=0.85)
    ax.axhline(0.0, color="#333333", linewidth=0.8)
    ax.set(title="Synthetic spectral abscissa versus sampled χ_G",
           xlabel="Sampled χ_G supremum [-]", ylabel="Spectral abscissa [s⁻¹]")
    ax.legend(frameon=False, ncol=3)
    save(fig, "FIG_B17_spectral_abscissa_vs_chiG")

    if any(r["case"] == "S5" for r in spectrum):
        s5 = read_table("TABLE_B17_optional_s5.csv")
        fig, ax = plt.subplots(figsize=(6.6, 4.8))
        labels = [f"dₖ, k={r['mode']}" for r in s5] + ["λmin(D_G)"]
        values = [float(r["diagonal_dissipative"]) for r in s5] + [float(s5[0]["lambda_min"])]
        colors = [COLORS[2]] * len(s5) + [COLORS[4]]
        ax.bar(labels, values, color=colors)
        ax.axhline(0, color="#333333", linewidth=0.8)
        ax.set(title="S5 positive diagonal entries with indefinite D_G", ylabel="Dissipative value [self-energy units]")
        ax.tick_params(axis="x", rotation=35)
        save(fig, "FIG_B18_positive_diagonal_indefinite_DG")


if __name__ == "__main__":
    main()
