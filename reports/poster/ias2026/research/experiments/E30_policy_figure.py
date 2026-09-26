"""Figure and comparison note for E30, from the saved lattice only."""

from __future__ import annotations

import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

from _overnight import run_directory  # noqa: E402

RUN = run_directory()
DIRECTORY = RUN / "E30_rebalanced_policy_lattice"
LABELS = {
    "matched": "P1 Q-matched (mechanism isolation)",
    "unity_pf": "P2 unity power factor",
    "vreg_pf095": "P3 voltage regulation, 0.95 pf capability",
    "vreg_scap": "P4 voltage regulation, apparent-power capability",
}
COLOURS = {
    "matched": "#c1121f",
    "unity_pf": "#e07a00",
    "vreg_pf095": "#0353a4",
    "vreg_scap": "#3f8efc",
}


def main() -> int:
    lattice = pd.read_csv(DIRECTORY / "E30_rebalanced_policy_lattice.csv")
    summary = pd.read_csv(DIRECTORY / "E30_policy_summary.csv")
    accepted = lattice[lattice.status == "ACCEPTED"]

    figure, axes = plt.subplots(1, 2, figsize=(13.0, 5.2))

    ax = axes[0]
    for policy, colour in COLOURS.items():
        block = accepted[accepted.policy == policy].sort_values(["size", "members"])
        jitter = 0.10 * (list(COLOURS).index(policy) - 1.5)
        ax.scatter(
            block["size"] + jitter,
            block.alpha_IA,
            s=42,
            color=colour,
            label=LABELS[policy],
            zorder=3,
            edgecolor="white",
            linewidth=0.6,
        )
        by_size = block.groupby("size").alpha_IA.max()
        ax.plot(by_size.index + jitter, by_size.values, color=colour, lw=1.2, alpha=0.5)
    ax.axhline(0.0, color="black", lw=1.2)
    ax.set_xlabel("number of machines replaced")
    ax.set_ylabel(r"$\alpha_{IA}$   (inter-area family envelope)")
    ax.set_title("Every subset of the flagship, by reactive policy")
    ax.set_xticks(range(5))
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8, loc="upper left")
    ax.text(
        0.02,
        0.02,
        "above the line = unstable",
        transform=ax.transAxes,
        fontsize=8,
        style="italic",
    )

    ax = axes[1]
    orders = range(5)
    for policy, colour in COLOURS.items():
        row = summary[summary.policy == policy]
        if row.empty:
            continue
        row = row.iloc[0]
        values = [row[f"alpha_le_{k}"] if k else row["alpha_base"] for k in orders]
        ax.plot(
            list(orders),
            values,
            marker="o",
            color=colour,
            label=f"{LABELS[policy]}  [{row.classification}]",
        )
    ax.axhline(0.0, color="black", lw=1.2)
    ax.set_xlabel("interaction order retained in the reconstruction")
    ax.set_ylabel(r"reconstructed $\alpha_{IA}$")
    ax.set_title("Where the reconstruction first predicts instability")
    ax.set_xticks(list(orders))
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8, loc="upper left")

    figure.suptitle(
        "Track A under four reactive-power policies: the instability is not "
        "irreducibly fourth order under every dispatch",
        fontsize=11,
    )
    figure.tight_layout()
    target = DIRECTORY / "E30_family_boundary_by_Qpolicy.png"
    figure.savefig(target, dpi=200)
    plt.close(figure)

    # exact figure source data
    accepted[
        ["policy", "members", "size", "alpha_IA", "freq_IA_worst_hz", "band_rhp_count"]
    ].to_csv(DIRECTORY / "E30_figure_source_panel_a.csv", index=False)
    summary[
        [
            "policy",
            "classification",
            "alpha_base",
            "alpha_le_1",
            "alpha_le_2",
            "alpha_le_3",
            "alpha_le_4",
            "mu4",
            "effective_order",
        ]
    ].to_csv(DIRECTORY / "E30_figure_source_panel_b.csv", index=False)
    print(f"figure -> {target}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
