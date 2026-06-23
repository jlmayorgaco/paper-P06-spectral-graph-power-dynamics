"""
Phase-0 IEEE 39-bus IBR benchmark artifacts.

This script does not claim final estimator validation against full ANDES IBR
poles. It consolidates the already generated ANDES baseline, IBR case
construction, positive-mode diagnostic, and IEEE39-derived surrogate artifacts
into a reproducible gate table and a topology diagram showing the SG-to-IBR
replacement map used by the Mix60 pilot case.
"""

from __future__ import annotations

import csv
import json
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs" / "ieee39_phase0_bridge"
FIG_OUT = ROOT / "paper_ieee_transactions" / "figures"


def read_json(path: Path):
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def write_json(path: Path, data) -> None:
    with path.open("w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)


def gen_to_bus_map(generator_buses: list[int]) -> dict[int, int]:
    return {idx + 1: int(bus) for idx, bus in enumerate(generator_buses)}


def build_phase0_gates(
    baseline: dict,
    diagnostics: list[dict],
    cases: list[dict],
    physical_bridge: dict | None = None,
) -> list[dict]:
    stable_cases = []
    ok_cases = []
    for case in cases:
        test = case.get("test", {})
        if test.get("ok"):
            ok_cases.append(case["name"])
        poles = test.get("poles", {})
        if poles.get("small_signal_stable"):
            stable_cases.append(case["name"])

    no_pss = next((r for r in diagnostics if r["variant"] == "no_pss"), None)
    positive_diag_pass = bool(no_pss and no_pss["n_positive"] == 0)
    physical_status = (physical_bridge or {}).get("status")
    physical_results = (physical_bridge or {}).get("results", [])
    physical_crit = physical_results[0].get("critical_match", {}) if physical_results else {}
    if physical_status:
        freq_error = physical_crit.get("relative_freq_error")
        zeta_error = physical_crit.get("relative_zeta_error")
        physical_evidence = (
            f"Strict bridge audit {physical_status}: "
            f"critical frequency error {100 * float(freq_error):.1f}% and damping-ratio error {100 * float(zeta_error):.1f}%"
            if freq_error is not None and zeta_error is not None
            else f"Strict bridge audit {physical_status}"
        )
    else:
        physical_evidence = "No calibrated mapping from full IBR controller states to scalar M,D has been validated."

    return [
        {
            "gate": "G0.1 baseline ANDES load/PFlow/EIG",
            "status": "PASS",
            "evidence": (
                f"{baseline['inventory']['Bus']} buses, "
                f"{baseline['inventory']['Line']} lines, "
                f"{baseline['eigs']['n_eigs']} eigenvalues"
            ),
            "next_action": "Use as synchronous extraction benchmark, not as IBR validation.",
        },
        {
            "gate": "G0.2 packaged-case positive-mode diagnostic",
            "status": "PASS-DIAGNOSTIC" if positive_diag_pass else "OPEN",
            "evidence": (
                "Removing IEEEST/PSS removes all positive-real modes"
                if positive_diag_pass
                else "Positive-real mode source is not isolated"
            ),
            "next_action": "Report positive modes as a model-quality warning.",
        },
        {
            "gate": "G0.3 SG-to-GFL/GFM case generation",
            "status": "PASS" if len(ok_cases) == len(cases) else "PARTIAL",
            "evidence": f"{len(ok_cases)}/{len(cases)} derived cases solved PFlow and EIG",
            "next_action": "Use generated cases as reproducible benchmark artifacts.",
        },
        {
            "gate": "G0.4 generic IBR small-signal stability",
            "status": "PARTIAL",
            "evidence": f"{len(stable_cases)}/{len(cases)} cases stable under generic gains: {', '.join(stable_cases) or 'none'}",
            "next_action": "Calibrate REGCP1/REECA1/REPCA1/PLL1 and REGF1 gains before final validation.",
        },
        {
            "gate": "G0.5 physical ANDES-to-(L,M,D) IBR bridge",
            "status": "BLOCKED",
            "evidence": physical_evidence,
            "next_action": "Extract a matched reduced model from the full ANDES Jacobian and compare poles first.",
        },
        {
            "gate": "G0.6 final estimator validation against full ANDES IBR poles",
            "status": "BLOCKED",
            "evidence": "Requires G0.4 and G0.5 before estimator errors can be claimed.",
            "next_action": "Compare diagonal, second-order, reduced QEP, and adaptive estimates against full eigensolve.",
        },
    ]


def draw_ieee39_mix60(meta: dict, cases: list[dict]) -> dict:
    bus_ids = [int(b) for b in meta["bus_ids"]]
    generator_buses = [int(b) for b in meta["generator_buses"]]
    gmap = gen_to_bus_map(generator_buses)
    mix60 = next(c for c in cases if c["name"] == "ieee39_ibr_mix60")
    gfl_gens = [int(g) for g in mix60["gfl"]]
    gfm_gens = [int(g) for g in mix60["gfm"]]
    gfl_buses = {gmap[g] for g in gfl_gens}
    gfm_buses = {gmap[g] for g in gfm_gens}
    sg_buses = set(generator_buses) - gfl_buses - gfm_buses

    graph = nx.Graph()
    graph.add_nodes_from(bus_ids)
    for line in meta["lines"]:
        bus1 = line.get("bus1", line.get("from_bus"))
        bus2 = line.get("bus2", line.get("to_bus"))
        graph.add_edge(int(bus1), int(bus2), weight=float(line["weight"]))

    pos = nx.spring_layout(graph, seed=39, weight="weight", k=0.42, iterations=800)

    fig, ax = plt.subplots(figsize=(9.0, 6.8))
    ax.set_title("IEEE 39-bus Mix60 SG-to-IBR replacement map", fontsize=13, pad=12)
    nx.draw_networkx_edges(graph, pos, ax=ax, edge_color="#9aa0a6", width=0.85, alpha=0.65)

    other_buses = sorted(set(bus_ids) - set(generator_buses))
    nx.draw_networkx_nodes(
        graph, pos, nodelist=other_buses, node_color="#d7dbe2",
        node_size=115, edgecolors="#666666", linewidths=0.35, ax=ax, label="network/load bus"
    )
    nx.draw_networkx_nodes(
        graph, pos, nodelist=sorted(sg_buses), node_color="#303030",
        node_shape="s", node_size=275, edgecolors="white", linewidths=0.8, ax=ax, label="SG retained"
    )
    nx.draw_networkx_nodes(
        graph, pos, nodelist=sorted(gfl_buses), node_color="#1f77b4",
        node_shape="^", node_size=330, edgecolors="white", linewidths=0.9, ax=ax, label="GFL replacement"
    )
    nx.draw_networkx_nodes(
        graph, pos, nodelist=sorted(gfm_buses), node_color="#d62728",
        node_shape="D", node_size=330, edgecolors="white", linewidths=0.9, ax=ax, label="GFM replacement"
    )

    labels = {b: str(b) for b in bus_ids}
    nx.draw_networkx_labels(graph, pos, labels=labels, font_size=6, font_color="#111111", ax=ax)

    for gen, bus in gmap.items():
        x, y = pos[bus]
        ax.text(
            x, y - 0.075, f"G{gen}", ha="center", va="top", fontsize=6.5,
            color="#111111", bbox=dict(facecolor="white", alpha=0.72, edgecolor="none", pad=0.7)
        )

    ax.legend(loc="lower left", frameon=True, fontsize=8)
    ax.axis("off")
    fig.tight_layout()

    pdf = OUT / "fig_ieee39_mix60_ibr_map.pdf"
    png = OUT / "fig_ieee39_mix60_ibr_map.png"
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(png, dpi=220, bbox_inches="tight")
    plt.close(fig)

    FIG_OUT.mkdir(parents=True, exist_ok=True)
    paper_pdf = FIG_OUT / "fig5_ieee39_ibr_map.pdf"
    paper_png = FIG_OUT / "fig5_ieee39_ibr_map.png"
    shutil.copyfile(pdf, paper_pdf)
    shutil.copyfile(png, paper_png)

    return {
        "case": mix60["name"],
        "generator_to_bus": gmap,
        "gfl_gens": gfl_gens,
        "gfl_buses": sorted(gfl_buses),
        "gfm_gens": gfm_gens,
        "gfm_buses": sorted(gfm_buses),
        "sg_retained_buses": sorted(sg_buses),
        "figure_pdf": str(pdf.relative_to(ROOT)),
        "figure_png": str(png.relative_to(ROOT)),
        "paper_figure_pdf": str(paper_pdf.relative_to(ROOT)),
    }


def write_gate_csv(path: Path, gates: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["gate", "status", "evidence", "next_action"])
        writer.writeheader()
        writer.writerows(gates)


def write_report(path: Path, gates: list[dict], diagram: dict, surrogate: dict) -> None:
    lines = [
        "# Phase-0 IEEE39 IBR Bridge Artifacts",
        "",
        "This artifact records what is currently validated before any final full-ANDES IBR estimator claim.",
        "",
        "## Mix60 replacement map",
        f"- GFL generators {diagram['gfl_gens']} at buses {diagram['gfl_buses']}.",
        f"- GFM generators {diagram['gfm_gens']} at buses {diagram['gfm_buses']}.",
        f"- Retained synchronous generator buses {diagram['sg_retained_buses']}.",
        f"- Figure: `{diagram['figure_pdf']}`.",
        "",
        "## Phase-0 gates",
    ]
    for gate in gates:
        lines.append(f"- {gate['status']}: {gate['gate']} -- {gate['evidence']}")
    lines.extend(
        [
            "",
            "## Current reduced-model evidence",
            (
                "- IEEE39-derived surrogate estimator errors against the surrogate full QEP: "
                f"diagonal median {surrogate['errors']['diag']['median']:.6f}, "
                f"second-order median {surrogate['errors']['second']['median']:.6f}, "
                f"reduced-QEP median {surrogate['errors']['rqep']['median']:.6f}, "
                f"adaptive median {surrogate['errors']['adaptive']['median']:.6f}."
            ),
            (
                "- Braess-like surrogate line reversals: "
                f"{surrogate['braess_like_count']}/{surrogate['n_lines_tested']} tested lines."
            ),
            "",
            "## Blocking item",
            "The physical bridge from the full calibrated ANDES IBR Jacobian to a scalar `(L,M,D)` model is not yet validated. "
            "The manuscript must therefore present these as benchmark-construction and surrogate-validation artifacts, not as final full-ANDES IBR estimator validation.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    baseline = read_json(ROOT / "outputs" / "ieee39_rational_filter" / "andes_baseline_summary.json")
    diagnostics = read_json(ROOT / "outputs" / "ieee39_mode_diagnostics" / "baseline_positive_mode_diagnostics.json")
    cases = read_json(ROOT / "outputs" / "ieee39_ibr_cases" / "case_generation_summary.json")
    surrogate = read_json(ROOT / "outputs" / "ieee39_rational_filter" / "surrogate_summary.json")
    meta = read_json(ROOT / "outputs" / "ieee39_rational_filter" / "network_metadata.json")
    physical_path = ROOT / "outputs" / "phase0_physical_bridge" / "phase0_physical_bridge_status.json"
    physical_bridge = read_json(physical_path) if physical_path.exists() else None

    gates = build_phase0_gates(baseline, diagnostics, cases, physical_bridge)
    diagram = draw_ieee39_mix60(meta, cases)

    phase0 = {
        "diagram": diagram,
        "gates": gates,
        "surrogate_snapshot": {
            "n_trials": surrogate["n_trials"],
            "errors": surrogate["errors"],
            "braess_like_count": surrogate["braess_like_count"],
            "n_lines_tested": surrogate["n_lines_tested"],
        },
    }

    write_json(OUT / "phase0_bridge_status.json", phase0)
    write_gate_csv(OUT / "phase0_gate_table.csv", gates)
    write_report(OUT / "phase0_bridge_report.md", gates, diagram, surrogate)

    print(f"Wrote {OUT / 'phase0_bridge_status.json'}")
    print(f"Wrote {OUT / 'phase0_gate_table.csv'}")
    print(f"Wrote {OUT / 'phase0_bridge_report.md'}")
    print(f"Wrote {OUT / 'fig_ieee39_mix60_ibr_map.pdf'}")


if __name__ == "__main__":
    main()
