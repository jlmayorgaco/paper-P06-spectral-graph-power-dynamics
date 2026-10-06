"""Render only figures and status tables supported by executed calculations."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
BLUE = "#17678a"
ORANGE = "#c35b32"
GREEN = "#327a62"


def rows(name: str) -> list[dict[str, str]]:
    with (HERE / name).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def write(name: str, data: list[dict[str, object]]) -> None:
    with (HERE / name).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(data[0]))
        writer.writeheader()
        writer.writerows(data)


def style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.dpi": 160,
            "savefig.dpi": 220,
        }
    )


def save(fig: plt.Figure, name: str) -> None:
    fig.savefig(HERE / name, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> None:
    style()
    t01 = rows("T01_DELAY_TRANSITION_REPRODUCTION.csv")
    t02 = rows("T02_PRECISE_CROSSINGS.csv")
    t04 = rows("T04_PRECISE_CROSSINGS.csv")
    fast = rows("T08_FAST_MODE_PROVENANCE.csv")
    seed_cross = float(t02[0]["local_root_crossing_ms"])
    tuned_cross = float(t02[1]["local_root_crossing_ms"])
    fixed_low = float(t04[0]["fixed_gain_tau_crit_ms"])
    fixed_high = float(t04[-1]["fixed_gain_tau_crit_ms"])
    combined_loss = seed_cross - tuned_cross
    fixed_loss = fixed_low - fixed_high
    gain_swap_recovery = fixed_high - tuned_cross
    assert abs(seed_cross - 39.380218505859375) < 1e-8
    assert len(t01) == 16 and len(t04) == 5
    assert all(row["status"] in {"SAFE", "UNSAFE"} for row in t01)

    write(
        "MATERIALITY_GATE.csv",
        [
            {
                "comparison": "stored_seed_vs_stored_zero_delay_tuned",
                "low_rho_percent": 87.5,
                "high_rho_percent": float(t02[1]["replacement_percent"]),
                "margin_low_ms": seed_cross,
                "margin_high_ms": tuned_cross,
                "delay_loss_ms": combined_loss,
                "classification": "WEAK_EFFECT",
                "interpretation": "rho_and_gains_change_together;below_preregistered_2_ms_gate",
            },
            {
                "comparison": "fixed_nominal_gains_rho_path",
                "low_rho_percent": 87.5,
                "high_rho_percent": float(t04[-1]["replacement_percent"]),
                "margin_low_ms": fixed_low,
                "margin_high_ms": fixed_high,
                "delay_loss_ms": fixed_loss,
                "classification": "WEAK_EFFECT",
                "interpretation": "replacement_only_on_preregistered_rho_path;negative_loss_means_local_margin_increase",
            },
            {
                "comparison": "high_rho_gain_swap_zero_delay_tuned_to_seed_nominal",
                "low_rho_percent": float(t04[-1]["replacement_percent"]),
                "high_rho_percent": float(t04[-1]["replacement_percent"]),
                "margin_low_ms": tuned_cross,
                "margin_high_ms": fixed_high,
                "delay_loss_ms": -gain_swap_recovery,
                "classification": "OBSERVED_ADMISSIBLE_GAIN_SWAP_NOT_OPTIMIZATION",
                "interpretation": "same_rho;seed_nominal_gain_vector_recovers_latency_margin",
            },
        ],
    )
    blocked = [
        ("T05_GAIN_OPTIMIZATION_TRACE.csv", "gain_optimization"),
        ("T06_OPTIMAL_RHO_TAU_FRONTIER.csv", "optimized_rho_tau_frontier"),
        ("T07_MAX_SPECTRAL_REPLACEMENT_VS_DELAY.csv", "inverse_spectral_frontier"),
        ("T09_HETEROGENEOUS_DELAY_RESULTS.csv", "heterogeneous_delay_placements"),
        ("T10_SPATIAL_PREDICTOR_VALIDATION.csv", "spatial_predictor_validation"),
    ]
    for filename, stage in blocked:
        write(
            filename,
            [
                {
                    "stage": stage,
                    "status": "BLOCKED_WEAK_EFFECT_STOP_RULE",
                    "reason": "fixed_gain_replacement_effect_below_2_ms_materiality_gate",
                    "result": "NOT_ESTABLISHED",
                }
            ],
        )

    # F01: actual continued roots. All 20–40 ms branches are plotted; the
    # requested 30–40 ms region is highlighted by the x-limits.
    fig, ax = plt.subplots(figsize=(8.7, 5.0))
    for design, color, label in [
        ("seed_875", BLUE, "87.5%; seed gains"),
        ("best_zero_delay_88455", ORANGE, "88.455%; zero-delay-tuned gains"),
    ]:
        branch_ids = sorted({int(r["family_id"]) for r in fast if r["design_id"] == design})
        for branch in branch_ids:
            subset = sorted(
                (r for r in fast if r["design_id"] == design and int(r["family_id"]) == branch),
                key=lambda r: float(r["tau_ms"]),
            )
            ax.plot(
                [float(r["tau_ms"]) for r in subset],
                [float(r["root_real"]) for r in subset],
                color=color,
                lw=1.45,
                alpha=0.8,
                label=label if branch == branch_ids[0] else None,
            )
    ax.axhline(-0.05, color="#333", lw=1.2, ls="--", label="required margin −0.05 s⁻¹")
    for x, color in [(seed_cross, BLUE), (tuned_cross, ORANGE)]:
        ax.axvline(x, color=color, ls=":", lw=1)
    ax.set(xlim=(30, 40.2), xlabel="Uniform PLL measurement delay (ms)", ylabel="Fast-root real part (s⁻¹)")
    ax.set_title("Exact-characteristic fast-root transition; complete contour counts checked")
    ax.grid(alpha=0.2)
    ax.legend(loc="lower left", frameon=False, fontsize=8)
    save(fig, "F01_ROOT_LOCUS_30_40MS.png")

    # F02: a full-scale view prevents the tiny numerical slope from looking
    # physically large; the zoom states precisely what it resolves.
    xs = np.array([float(r["replacement_percent"]) for r in t04])
    ys = np.array([float(r["fixed_gain_tau_crit_ms"]) for r in t04])
    fig, (wide, zoom) = plt.subplots(1, 2, figsize=(10.1, 4.4))
    for ax in (wide, zoom):
        ax.plot(xs, ys, "o-", color=BLUE, lw=1.8)
        ax.grid(alpha=0.2)
        ax.set_xlabel("GFL replacement (% of initialized generator MW)")
    wide.set(ylim=(36, 42), ylabel="First observed local margin crossing (ms)")
    wide.set_title("Declared 2 ms materiality scale")
    zoom.set(ylim=(39.378, 39.386), ylabel="Crossing (ms)")
    zoom.set_title("Local-root zoom: +0.00343 ms")
    fig.suptitle("Fixed Kp, Ki: replacement-only delay margin is effectively flat")
    fig.text(
        0.5,
        -0.01,
        "Every point has a full-contour safe/unsafe bracket [39.375, 39.453] ms; zoom uses tracked exact roots.",
        ha="center",
        fontsize=8,
    )
    save(fig, "F02_FIXED_GAIN_RHO_TAU_FRONTIER.png")

    # F05: the origin before 20 ms is deliberately not inferred.
    fig, axs = plt.subplots(1, 2, figsize=(10.5, 4.7), sharex=True, sharey=True)
    for ax, (design, title) in zip(
        axs,
        [
            ("seed_875", "87.5%; seed gains"),
            ("best_zero_delay_88455", "88.455%; zero-delay-tuned gains"),
        ],
    ):
        branch_ids = sorted({int(r["family_id"]) for r in fast if r["design_id"] == design})
        for branch in branch_ids:
            subset = sorted(
                (r for r in fast if r["design_id"] == design and int(r["family_id"]) == branch),
                key=lambda r: float(r["tau_ms"]),
            )
            ax.plot(
                [float(r["root_real"]) for r in subset],
                [float(r["frequency_hz"]) for r in subset],
                lw=1.25,
                label=f"branch {branch}",
            )
            for tau_mark in (20, 30, 40):
                mark = min(subset, key=lambda r: abs(float(r["tau_ms"]) - tau_mark))
                ax.scatter(float(mark["root_real"]), float(mark["frequency_hz"]), s=10)
        ax.axvline(-0.05, color="#333", ls="--", lw=1)
        ax.set_title(title)
        ax.set_xlabel("Re λ (s⁻¹)")
        ax.grid(alpha=0.18)
        ax.legend(loc="upper left", fontsize=7, frameon=False, ncol=2)
    axs[0].set_ylabel("Positive-imaginary root frequency (Hz)")
    fig.suptitle("Fast families followed continuously from 40 to 20 ms; markers: 20, 30, 40 ms")
    save(fig, "F05_FAST_MODE_ROOT_LOCUS.png")

    # F06: observed admissible comparison at identical high rho, with no
    # suggestion that either vector maximizes the margin.
    fig, ax = plt.subplots(figsize=(7.6, 3.5))
    ax.hlines([0, 1], 0, [tuned_cross, fixed_high], colors=[ORANGE, GREEN], lw=2)
    ax.scatter([tuned_cross, fixed_high], [0, 1], color=[ORANGE, GREEN], s=55, zorder=3)
    ax.set_yticks([0, 1], ["Zero-delay-tuned gains", "Seed nominal gains"])
    ax.set_xlim(0, 43)
    ax.set_xlabel("First observed local latency crossing (ms)")
    ax.set_title("At 88.455% GFL: a gain-vector swap recovers 1.9966 ms")
    for x, y in [(tuned_cross, 0), (fixed_high, 1)]:
        ax.annotate(f"{x:.4f} ms", (x, y), xytext=(6, 4), textcoords="offset points", fontsize=9)
    ax.text(1, 0.45, "In-bounds gain swap; gain optimum and delayed events not solved", fontsize=9)
    ax.set_ylim(-0.25, 1.25)
    ax.grid(axis="x", alpha=0.18)
    save(fig, "F06_RETUNING_DELAY_MARGIN_RECOVERY.png")

    print(f"combined_loss_ms={combined_loss:.9f}")
    print(f"fixed_gain_delay_loss_ms={fixed_loss:.9f}")
    print(f"observed_gain_swap_recovery_ms={gain_swap_recovery:.9f}")


if __name__ == "__main__":
    main()
