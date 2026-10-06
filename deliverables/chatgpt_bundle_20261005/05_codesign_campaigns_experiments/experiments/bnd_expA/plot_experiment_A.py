#!/usr/bin/env python
"""Generate Experiment A figures from the saved CSV tables."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def read_csv(path: Path):
    with path.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def floats(rows, key):
    out = []
    for row in rows:
        try:
            out.append(float(row[key]))
        except (ValueError, KeyError, TypeError):
            out.append(np.nan)
    return np.asarray(out, dtype=float)


def save(fig, outdir: Path, stem: str, title: str):
    fig.suptitle(f"{title}\nOperating case: IEEE-39, one SimpleGFLDC at bus 33", y=1.02)
    fig.tight_layout()
    description = "Operating case: IEEE-39 bus 33 SG-to-GFL; see REPORT_EXP_A.md."
    fig.savefig(outdir / f"{stem}.png", dpi=220, bbox_inches="tight", metadata={"Title": title, "Description": description})
    fig.savefig(outdir / f"{stem}.pdf", bbox_inches="tight", metadata={"Title": title, "Subject": description})
    fig.savefig(outdir / f"{stem}.svg", bbox_inches="tight", metadata={"Title": title, "Description": description})
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--report-dir", required=True, type=Path)
    args = parser.parse_args()
    out = args.report_dir / "figures"
    tables = args.report_dir / "tables"
    out.mkdir(parents=True, exist_ok=True)

    # A01: full finite spectrum with the 20 least-damped poles highlighted.
    eig = read_csv(tables / "full_eigenvalues.csv")
    least = {int(r["mode"]) for r in read_csv(tables / "TABLE_A04_least_damped_poles.csv")}
    x, y = floats(eig, "lambda_real"), floats(eig, "lambda_imag") / (2 * np.pi)
    keep = np.asarray([int(r["mode"]) in least for r in eig])
    fig, ax = plt.subplots(figsize=(7.4, 5.2))
    ax.scatter(x[~keep], y[~keep], s=15, alpha=.55, label="Other finite poles")
    ax.scatter(x[keep], y[keep], s=34, marker="x", label="20 least-damped poles")
    ax.axvline(0, color="black", lw=.8)
    ax.set(xlabel=r"Real part of pole [s$^{-1}$]", ylabel="Imaginary part / 2π [Hz]", title="Full finite eigenvalue map")
    ax.grid(True, alpha=.25); ax.legend()
    save(fig, out, "FIG_A01_full_eigenvalue_map", "Full finite eigenvalue map")

    # A02: honest pole certificate: full poles classified by the Schur singular-value test.
    poles = read_csv(tables / "TABLE_A05_schur_pole_validation.csv")
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(10.4, 4.6))
    classes = sorted({r["classification"] for r in poles})
    for cls in classes:
        rr = [r for r in poles if r["classification"] == cls]
        ax.scatter(floats(rr, "full_pole_real"), floats(rr, "full_pole_imag")/(2*np.pi), s=24, label=cls)
    ax.axvline(0, color="black", lw=.8)
    ax.set(xlabel=r"Real part of full pole [s$^{-1}$]", ylabel="Imaginary part / 2π [Hz]", title="Full poles and Schur classification")
    ax.grid(True, alpha=.25); ax.legend(fontsize=7)
    visible = [r for r in poles if r["classification"] == "retained-visible pole"]
    bx.semilogy(np.arange(1, len(visible)+1), np.maximum(floats(visible, "normalized_schur_pole_residual"), 1e-18), "o")
    bx.axhline(1e-8, color="black", ls="--", lw=.9, label="Pole residual threshold")
    bx.set(xlabel="Retained-visible pole rank", ylabel=r"σmin(Sr(λ)) / ||Sr(λ)||", title="Schur zero residual at full poles")
    bx.grid(True, which="both", alpha=.25); bx.legend(fontsize=8)
    save(fig, out, "FIG_A02_full_vs_schur_poles", "Full poles and exact Schur pole certificates")

    # A03: block identity residual on the imaginary axis.
    residual = read_csv(tables / "TABLE_A06_frequency_residuals.csv")
    imag = [r for r in residual if r["point_kind"] == "imaginary_axis"]
    f, e = floats(imag, "frequency_hz"), floats(imag, "schur_residual")
    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    ax.loglog(f, np.maximum(e, 1e-18), label="Normalized block reconstruction")
    ax.axvspan(.1, 5, alpha=.12, label="Primary synchronization band")
    ax.set(xlabel="Frequency [Hz]", ylabel="Relative Frobenius residual", title="Exact Schur identity residual")
    ax.grid(True, which="both", alpha=.25); ax.legend()
    save(fig, out, "FIG_A03_schur_identity_residual", "Schur identity residual")

    # A04: second-order bridge residual.
    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    f, e = floats(imag, "frequency_hz"), floats(imag, "bnd_bridge_residual")
    ax.loglog(f, np.maximum(e, 1e-18), label="Independent Schur/angle bridge residual")
    ax.axvspan(.1, 5, alpha=.12, label="Primary synchronization band")
    ax.axhline(1e-8, color="black", ls="--", lw=.9, label="PASS-EXACT threshold")
    ax.axhline(1e-2, color="gray", ls=":", lw=.9, label="PASS-APPROX threshold")
    ax.set(xlabel="Frequency [Hz]", ylabel=r"εBND", title="Second-order bridge residual")
    ax.grid(True, which="both", alpha=.25); ax.legend()
    save(fig, out, "FIG_A04_BND_bridge_residual", "Second-order bridge residual")

    # A05: frequency-varying self-energy norms (Πq where a global Σ is unavailable).
    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    f = floats(imag, "frequency_hz")
    ax.loglog(f, np.maximum(floats(imag, "sigma_norm2"), 1e-18), label=r"2-norm")
    ax.loglog(f, np.maximum(floats(imag, "sigma_normF"), 1e-18), label="Frobenius norm")
    ax.set(xlabel="Frequency [Hz]", ylabel="Self-energy norm [operator units]", title="Frequency-dependent self-energy")
    ax.grid(True, which="both", alpha=.25); ax.legend()
    save(fig, out, "FIG_A05_self_energy_norms", "Self-energy norms")

    # A06: first five singular values from the tidy table.
    sing = read_csv(tables / "sigma_singular_spectrum.csv")
    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    for rank in sorted({int(r["rank"]) for r in sing}):
        rr = [r for r in sing if int(r["rank"]) == rank]
        ax.loglog(floats(rr, "frequency_hz"), np.maximum(floats(rr, "singular_value"), 1e-18), label=f"σ{rank}")
    ax.set(xlabel="Frequency [Hz]", ylabel="Singular value [operator units]", title="Leading self-energy singular values")
    ax.grid(True, which="both", alpha=.25); ax.legend()
    save(fig, out, "FIG_A06_self_energy_singular_values", "Self-energy singular values")

    # A07: Hermitian spectrum under the declared q/v port convention.
    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    kind = imag[0].get("self_energy_kind", "self-energy") if imag else "self-energy"
    ax.semilogx(floats(imag, "frequency_hz"), floats(imag, "hermitian_min"), label="Minimum eigenvalue")
    ax.semilogx(floats(imag, "frequency_hz"), floats(imag, "hermitian_max"), label="Maximum eigenvalue")
    ax.axvspan(.1, 5, alpha=.12, label="Primary synchronization band")
    ax.axhline(0, color="black", lw=.8)
    ax.set(xlabel="Frequency [Hz]", ylabel="Hermitian-part eigenvalue [operator units]", title=f"Hermitian part of {kind} under declared port convention")
    ax.grid(True, which="both", alpha=.25); ax.legend()
    save(fig, out, "FIG_A07_hermitian_spectrum", "Hermitian-part spectrum")

    # A08: critical-frequency real part of Σ or Π.
    mat = read_csv(tables / "TABLE_A07_sigma_critical_frequency.csv")
    kind = mat[0].get("self_energy_kind", "self-energy") if mat else "self-energy"
    names = sorted({r["row_mode_or_state"] for r in mat})
    values = np.full((len(names), len(names)), np.nan)
    lookup = {name: i for i, name in enumerate(names)}
    for r in mat:
        values[lookup[r["row_mode_or_state"]], lookup[r["col_mode_or_state"]]] = float(r["real"])
    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    im = ax.imshow(values, aspect="auto", interpolation="nearest", cmap="coolwarm")
    ax.set_xticks(range(len(names)), [n.replace("VIndex(", "") for n in names], rotation=90, fontsize=6)
    ax.set_yticks(range(len(names)), [n.replace("VIndex(", "") for n in names], fontsize=6)
    ax.set(xlabel="Synchronization coordinate", ylabel="Synchronization coordinate", title=f"Real part of {kind} at critical-mode frequency")
    fig.colorbar(im, ax=ax, label="Real part [operator units]")
    save(fig, out, "FIG_A08_sigma_heatmap_critical_frequency", f"Critical-frequency {kind} matrix")

    # A09: conditional on a verified generalized graph-modal transform.
    gpath = tables / "graph_modal_diagnostics.csv"
    shpath = args.report_dir / "matrices" / "bus33_graph_modal_Sigmahat_critical.csv"
    if gpath.exists() and shpath.exists():
        with shpath.open(newline="", encoding="utf-8-sig") as fobj:
            raw = list(csv.reader(fobj))
        arr = np.asarray([[float(x) for x in row[1:]] for row in raw[1:]], dtype=float)
        fig, ax = plt.subplots(figsize=(6.5, 5.5))
        im = ax.imshow(np.abs(arr), aspect="auto", interpolation="nearest")
        ax.set(xlabel="Graph mode", ylabel="Graph mode", title=r"$|\hat{\Sigma}(j\omega_c)|$")
        fig.colorbar(im, ax=ax, label="Magnitude [operator units]")
        save(fig, out, "FIG_A09_graph_modal_sigma_heatmap", "Graph-modal self-energy at critical frequency")

    # A10: nonlinear and full linearized response comparison.
    trace = read_csv(tables / "TDS_trace.csv")
    t = floats(trace, "time_s")
    t_plot = t.copy()
    t_plot[(np.abs(t-1.0)<=1e-9) | (np.abs(t-1.1)<=1e-9)] = np.nan
    fig, axes = plt.subplots(2, 2, figsize=(10, 6.2), sharex=True)
    angle_signal = next(k[:-10] for k in trace[0] if k.startswith("bus") and k.endswith("_angle_rad_nonlinear"))
    signals = ["bus30_speed_hz", angle_signal, "gfl33_pll_angle_rad", "gfl33_terminal_power_pu"]
    for ax, signal in zip(axes.flat, signals):
        ax.plot(t_plot, floats(trace, signal+"_nonlinear"), label="Nonlinear PowerDynamics")
        ax.plot(t_plot, floats(trace, signal+"_linear"), ls="--", label="Full linearized descriptor")
        ax.axvspan(1.0, 1.1, alpha=.1)
        ax.set(xlabel="Time [s]", ylabel=signal.replace("_", " "), title=signal.replace("_", " "))
        ax.grid(True, alpha=.25)
    axes[0,0].legend(fontsize=8)
    save(fig, out, "FIG_A10_linear_vs_nonlinear_TDS", "Small-signal nonlinear versus linear TDS")

    # A11: required four-case summary.
    cross = read_csv(tables / "TABLE_A09_cross_bus_validation.csv")
    if cross:
        buses = [str(r["bus"]) for r in cross]
        xx = np.arange(len(buses))
        fig, axes = plt.subplots(2, 2, figsize=(10, 6.2))
        metrics = [("spectral_abscissa", "Spectral abscissa [s⁻¹]"),
                   ("critical_frequency_hz", "Critical frequency [Hz]"),
                   ("bnd_p95_error", "BND p95 residual"),
                   ("schur_p95_error", "Schur p95 residual")]
        for ax, (key, label) in zip(axes.flat, metrics):
            vals = floats(cross, key)
            if key.endswith("error"):
                vals = np.maximum(vals, 1e-18)
                ax.set_yscale("log")
            ax.plot(xx, vals, "o-")
            ax.set_xticks(xx, buses)
            ax.set(xlabel="GFL replacement bus", ylabel=label, title=label)
            ax.grid(True, which="both", alpha=.25)
        save(fig, out, "FIG_A11_cross_bus_summary", "Cross-bus extraction summary")


if __name__ == "__main__":
    main()
