"""
Create an IEEE39 topology result map for the manuscript.

The figure overlays the IEEE39 topology with:
- Mix60 SG/GFL/GFM replacement classes,
- top damping weak nodes from the IEEE39-derived surrogate,
- the best damping-aware weak link,
- the strongest Braess-like line reversals.

All data are read from versioned outputs; this script does not rerun ANDES.
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs" / "ieee39_result_map"
FIG_OUT = ROOT / "paper_ieee_transactions" / "figures"


def read_json(path: Path):
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def gen_to_bus_map(generator_buses: list[int]) -> dict[int, int]:
    return {idx + 1: int(bus) for idx, bus in enumerate(generator_buses)}


def line_key(a: int, b: int) -> tuple[int, int]:
    return tuple(sorted((int(a), int(b))))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG_OUT.mkdir(parents=True, exist_ok=True)

    meta = read_json(ROOT / "outputs" / "ieee39_rational_filter" / "network_metadata.json")
    cases = read_json(ROOT / "outputs" / "ieee39_ibr_cases" / "case_generation_summary.json")
    surrogate = read_json(ROOT / "outputs" / "ieee39_rational_filter" / "surrogate_summary.json")
    weak_nodes = pd.read_csv(ROOT / "outputs" / "ieee39_rational_filter" / "weak_node_table.csv")
    weak_links = pd.read_csv(ROOT / "outputs" / "ieee39_rational_filter" / "weak_link_table.csv")

    bus_ids = [int(b) for b in meta["bus_ids"]]
    generator_buses = [int(b) for b in meta["generator_buses"]]
    gmap = gen_to_bus_map(generator_buses)
    mix60 = next(c for c in cases if c["name"] == "ieee39_ibr_mix60")
    gfl_buses = {gmap[int(g)] for g in mix60["gfl"]}
    gfm_buses = {gmap[int(g)] for g in mix60["gfm"]}
    sg_buses = set(generator_buses) - gfl_buses - gfm_buses

    graph = nx.Graph()
    graph.add_nodes_from(bus_ids)
    for line in meta["lines"]:
        a = int(line.get("bus1", line.get("from_bus")))
        b = int(line.get("bus2", line.get("to_bus")))
        graph.add_edge(a, b, weight=float(line["weight"]))
    pos = nx.spring_layout(graph, seed=39, weight="weight", k=0.42, iterations=800)

    top_nodes = weak_nodes.sort_values("delta_zeta", ascending=False).head(3)
    top_node_buses = [int(x) for x in top_nodes["bus"].tolist()]
    best_link = weak_links.sort_values("delta_zeta", ascending=False).iloc[0]
    best_link_key = line_key(best_link["from_bus"], best_link["to_bus"])
    worst_reversals = (
        weak_links[weak_links["braess_like_reversal"] == True]
        .sort_values("delta_zeta", ascending=True)
        .head(5)
    )
    reversal_keys = {line_key(r["from_bus"], r["to_bus"]) for _, r in worst_reversals.iterrows()}

    fig, (ax, ax_text) = plt.subplots(
        1, 2, figsize=(13.5, 6.8), gridspec_kw={"width_ratios": [2.25, 1.0]}
    )
    fig.subplots_adjust(left=0.02, right=0.985, top=0.88, bottom=0.06, wspace=0.03)
    fig.suptitle("IEEE 39-bus simulation evidence map", fontsize=14, y=0.965)

    regular_edges = []
    best_edges = []
    reversal_edges = []
    for a, b in graph.edges():
        key = line_key(a, b)
        if key == best_link_key:
            best_edges.append((a, b))
        elif key in reversal_keys:
            reversal_edges.append((a, b))
        else:
            regular_edges.append((a, b))

    nx.draw_networkx_edges(graph, pos, ax=ax, edgelist=regular_edges, edge_color="#b8bec7", width=0.8, alpha=0.6)
    nx.draw_networkx_edges(graph, pos, ax=ax, edgelist=reversal_edges, edge_color="#d95f02", width=2.0, alpha=0.9)
    nx.draw_networkx_edges(graph, pos, ax=ax, edgelist=best_edges, edge_color="#1a9850", width=3.4, alpha=0.95)

    other_buses = sorted(set(bus_ids) - set(generator_buses))
    nx.draw_networkx_nodes(graph, pos, nodelist=other_buses, node_color="#d7dbe2", node_size=105,
                           edgecolors="#666666", linewidths=0.35, ax=ax)
    nx.draw_networkx_nodes(graph, pos, nodelist=sorted(sg_buses), node_color="#303030", node_shape="s",
                           node_size=260, edgecolors="white", linewidths=0.8, ax=ax)
    nx.draw_networkx_nodes(graph, pos, nodelist=sorted(gfl_buses), node_color="#1f77b4", node_shape="^",
                           node_size=320, edgecolors="white", linewidths=0.9, ax=ax)
    nx.draw_networkx_nodes(graph, pos, nodelist=sorted(gfm_buses), node_color="#d62728", node_shape="D",
                           node_size=320, edgecolors="white", linewidths=0.9, ax=ax)

    # Highlight top weak nodes with rings without hiding their device class.
    nx.draw_networkx_nodes(graph, pos, nodelist=top_node_buses, node_color="none", node_size=620,
                           edgecolors="#fdae61", linewidths=2.8, ax=ax)

    labels = {b: str(b) for b in bus_ids}
    nx.draw_networkx_labels(graph, pos, labels=labels, font_size=6, font_color="#111111", ax=ax)

    for gen, bus in gmap.items():
        x, y = pos[bus]
        ax.text(
            x, y - 0.075, f"G{gen}", ha="center", va="top", fontsize=6.3,
            color="#111111", bbox=dict(facecolor="white", alpha=0.72, edgecolor="none", pad=0.6)
        )

    ax.axis("off")
    ax.set_title("Topology overlay", fontsize=11)

    ax_text.axis("off")
    text_lines = [
        "Simulation evidence",
        "",
        "C2 estimator:",
        "diag median 0.468%",
        "2nd-order median 0.032%",
        "rQEP p95 5.44%",
        "",
        "Weak nodes:",
        f"best bus {top_node_buses[0]}",
        "top ring = damping placement",
        "",
        "Weak links:",
        f"best line {int(best_link['from_bus'])}-{int(best_link['to_bus'])}",
        f"Delta zeta = {best_link['delta_zeta']:.4f}",
        "",
        "Braess-like reversals:",
        f"{surrogate['braess_like_count']}/{surrogate['n_lines_tested']} lines",
        "orange = strongest negative",
        "",
        "ANDES Phase-1:",
        "50 case/profile EIG audits",
        "full estimator gate still blocked",
    ]
    ax_text.text(
        0.0, 0.98, "\n".join(text_lines), va="top", ha="left",
        fontsize=8.5, transform=ax_text.transAxes
    )

    # Manual legend.
    legend_items = [
        ("#303030", "s", "SG retained"),
        ("#1f77b4", "^", "GFL replacement"),
        ("#d62728", "D", "GFM replacement"),
        ("#fdae61", "o", "top weak node ring"),
        ("#1a9850", "-", "best damping link"),
        ("#d95f02", "-", "Braess-like reversal"),
    ]
    y0 = 0.08
    for i, (color, marker, label) in enumerate(legend_items):
        y = y0 + i * 0.045
        if marker == "-":
            ax_text.plot([0.02, 0.12], [y, y], color=color, linewidth=3, transform=ax_text.transAxes)
        else:
            ax_text.scatter([0.07], [y], s=90, c=color if label != "top weak node ring" else "none",
                            edgecolors=color, marker=marker, linewidths=2, transform=ax_text.transAxes)
        ax_text.text(0.16, y, label, va="center", fontsize=8.2, transform=ax_text.transAxes)

    pdf = OUT / "fig_ieee39_result_map.pdf"
    png = OUT / "fig_ieee39_result_map.png"
    fig.savefig(pdf)
    fig.savefig(png, dpi=220)
    plt.close(fig)

    paper_pdf = FIG_OUT / "fig6_ieee39_result_map.pdf"
    paper_png = FIG_OUT / "fig6_ieee39_result_map.png"
    shutil.copyfile(pdf, paper_pdf)
    shutil.copyfile(png, paper_png)

    summary = {
        "top_weak_nodes": top_nodes.to_dict(orient="records"),
        "best_weak_link": best_link.to_dict(),
        "worst_braess_like_reversals": worst_reversals.to_dict(orient="records"),
        "figure_pdf": str(pdf.relative_to(ROOT)),
        "paper_figure_pdf": str(paper_pdf.relative_to(ROOT)),
    }
    with (OUT / "ieee39_result_map_summary.json").open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    print(f"Wrote {pdf}")
    print(f"Wrote {OUT / 'ieee39_result_map_summary.json'}")


if __name__ == "__main__":
    main()
