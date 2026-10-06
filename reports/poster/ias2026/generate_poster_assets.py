#!/usr/bin/env python3
"""Build editable poster values and the IEEE-39 local TikZ diagrams from audits."""
from __future__ import annotations

import csv
import json
import math
import re
import shutil
import sys
from pathlib import Path

import networkx as nx
from PIL import Image, ImageDraw


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RUN = HERE / "research" / "bnd_h4_mechanism" / "results" / "20260927T144629Z_d0fecb32_ias26_060_operating_v1"
OUT = HERE / "generated"
TABLES = OUT / "tables"
FIGURES = OUT / "figures"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def tex_number(value: float, digits: int = 3, signed: bool = False) -> str:
    text = f"{value:+.{digits}f}" if signed else f"{value:.{digits}f}"
    return text.rstrip("0").rstrip(".") if "." in text else text


def tex_sci(value: str | float, sig: int = 3) -> str:
    number = float(value)
    mantissa, exponent = f"{number:.{sig - 1}e}".split("e")
    return rf"{mantissa}\times10^{{{int(exponent)}}}"


def write_tikz_network(bus_ids: list[int]) -> None:
    data = json.loads((HERE / "research" / "configs" / "ias2026" / "ieee39_network.json").read_text(encoding="utf-8"))
    nodes = [int(row["idx"]) for row in data["buses"]]
    edges = [(int(row["bus1"]), int(row["bus2"])) for row in data["lines"] if int(row.get("u", 1)) == 1]
    graph = nx.Graph()
    graph.add_nodes_from(nodes)
    graph.add_edges_from(edges)
    positions = nx.spring_layout(graph, seed=20260927, k=0.50, iterations=300)
    xs = [float(positions[n][0]) for n in nodes]
    ys = [float(positions[n][1]) for n in nodes]
    xmin, xmax, ymin, ymax = min(xs), max(xs), min(ys), max(ys)
    # Portrait spring layout to match the IEEE-39 one-line figure in the approved mockup.
    # The topology and bus positions still come from the frozen network/configuration.
    coords = {
        n: (5.0 + (float(positions[n][1]) - ymin) / (ymax - ymin) * 110.0,
            5.0 + (float(positions[n][0]) - xmin) / (xmax - xmin) * 140.0)
        for n in nodes
    }
    lines = ["% Generated from ieee39_network.json; spring layout seed fixed in generator."]
    for n in nodes:
        x, y = coords[n]
        lines.append(f"\\coordinate (b{n}) at ({x:.2f}mm,{y:.2f}mm);")
    for a, b in edges:
        lines.append(f"\\draw[IASDarkText!32,line width=0.8pt] (b{a})--(b{b});")
    for n in nodes:
        x, y = coords[n]
        if n in bus_ids:
            style = "circle,draw=IASVermillion,fill=IASVermillion,line width=1.8pt,text=white,minimum size=11.5mm,font=\\bfseries\\fontsize{30pt}{30pt}\\selectfont"
        elif n in {int(row["bus"]) for row in data["machines"]}:
            style = "circle,draw=IASMidGreen!60,fill=IASIvory,line width=1.0pt,text=IASGray!90,minimum size=7.2mm,font=\\fontsize{24.2pt}{24.2pt}\\selectfont"
        else:
            style = "circle,draw=IASGray!60,fill=IASIvory,line width=0.8pt,text=IASGray!90,minimum size=7.2mm,font=\\fontsize{24.2pt}{24.2pt}\\selectfont"
        lines.append(
            f"\\node[{style},inner sep=0pt] "
            f"at (b{n}) {{{n}}};"
        )
    (FIGURES / "network_ieee39.tikz").write_text("\n".join(lines) + "\n", encoding="utf-8")


def f10_baseline_macros() -> dict[str, str]:
    """Read the re-audited F10 exact-H counts from the frozen baseline table."""
    rows = read_csv(HERE / "research" / "results" / "20260911_BASELINE_COMPARISON.csv")
    methods = {row["method_id"]: row for row in rows}
    expected = {"BaselineModalExactCount": "B3", "BaselineAdditiveExactCount": "B4", "BaselinePairwiseExactCount": "B5"}
    parsed: dict[str, tuple[int, int]] = {}
    for macro, method_id in expected.items():
        token = methods[method_id]["taskB_F10_52_H_exact"].strip()
        match = re.fullmatch(r"(\d+)/(\d+)", token)
        if not match:
            raise ValueError(f"Unexpected F10 exact-H value for {method_id}: {token!r}")
        parsed[macro] = (int(match.group(1)), int(match.group(2)))
    denominators = {total for _, total in parsed.values()}
    if len(denominators) != 1:
        raise ValueError("The F10 baseline methods do not share a common benchmark size.")
    result = {macro: str(count) for macro, (count, _) in parsed.items()}
    result["BaselineTotal"] = str(denominators.pop())
    return result


def update_f10_baseline_macros() -> None:
    """Refresh only the F10 macros in results.tex without rebuilding other assets."""
    path = OUT / "results.tex"
    source = path.read_text(encoding="utf-8")
    macros = f10_baseline_macros()
    legacy_names = ("F10ModalExactCount", "F10AdditiveExactCount", "F10PairwiseExactCount", "F10BaselineTotal")
    for name in (*macros, *legacy_names):
        source = re.sub(rf"^\\newcommand\{{\\{name}\}}\{{[^\n]*\}}\n?", "", source, flags=re.MULTILINE)
    lines = [f"\\newcommand{{\\{name}}}{{{value}}}" for name, value in macros.items()]
    first_line, separator, rest = source.partition("\n")
    path.write_text(first_line + "\n" + "\n".join(lines) + (separator + rest if separator else "\n"), encoding="utf-8")


def write_hasse(bus_ids: list[int]) -> None:
    from itertools import combinations

    ordered = tuple(bus_ids)
    subsets: list[tuple[int, ...]] = []
    for size in range(len(ordered) + 1):
        subsets.extend(combinations(ordered, size))
    coords: dict[tuple[int, ...], tuple[float, float]] = {}
    for size in range(len(ordered) + 1):
        rank = [s for s in subsets if len(s) == size]
        step = 18.0
        for i, subset in enumerate(rank):
            coords[subset] = ((i - (len(rank) - 1) / 2) * step, size * 24.0)
    lines = ["% 16 audited H4 subsets; proper subsets stable; full H4 unstable."]
    for subset in subsets:
        if len(subset) == len(ordered):
            continue
        for item in ordered:
            if item not in subset:
                superset = tuple(sorted((*subset, item)))
                x1, y1 = coords[subset]
                x2, y2 = coords[superset]
                lines.append(f"\\draw[IASRule,line width=1.0pt] ({x1:.1f}mm,{y1:.1f}mm)--({x2:.1f}mm,{y2:.1f}mm);")
    for subset in subsets:
        x, y = coords[subset]
        if len(subset) == len(ordered):
            lines.append(f"\\node[circle,draw=IASVermillion,fill=IASVermillion,line width=1pt,minimum size=7.0mm,inner sep=0pt] at ({x:.1f}mm,{y:.1f}mm) {{}};")
        else:
            lines.append(f"\\node[circle,draw=IASBlue,fill=IASBlue,minimum size=6.0mm,inner sep=0pt] at ({x:.1f}mm,{y:.1f}mm) {{}};")
    blocker_label = r"\{" + ",".join(str(item) for item in ordered) + r"\}"
    lines.append(
        f"\\node[font=\\bfseries\\fontsize{{24.2pt}}{{26.2pt}}\\selectfont,text=IASVermillion,"
        f"align=center,anchor=south] at (0mm,101mm) {{${blocker_label}$\\\\(unstable)}};"
    )
    (FIGURES / "hasse_h4.tikz").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_policy_atlas() -> None:
    source = HERE / "research" / "results" / "FINAL_CLOSURE" / "figures" / "FIG1_policy_map_source.csv"
    rows = read_csv(source)
    if not rows or set(rows[0]) < {"g", "k", "kappa"}:
        raise ValueError("The audited F7A policy-map data are missing required columns.")
    g_min = min(float(row["g"]) for row in rows)
    g_max = max(float(row["g"]) for row in rows)
    k_min = min(float(row["k"]) for row in rows)
    k_max = max(float(row["k"]) for row in rows)
    # The sweep is sampled on a regular log grid. Paint one tile per sampled
    # coordinate so the raster stays continuous at poster scale and missing
    # combinations remain visibly unresolved.
    width, height = 3000, 2600
    image = Image.new("RGB", (width, height), "#D8DEDA")
    draw = ImageDraw.Draw(image)
    colors = {
        -1: "#EEF5F0",  # no hyperedge / transversely composable
         1: "#1768AC",
         2: "#006B4A",
         3: "#C99A20",
         4: "#C93624",
    }
    g_values = sorted({float(row["g"]) for row in rows})
    k_values = sorted({float(row["k"]) for row in rows}, reverse=True)
    g_index = {value: index for index, value in enumerate(g_values)}
    k_index = {value: index for index, value in enumerate(k_values)}
    def cell_bounds(index: int, count: int, extent: int) -> tuple[int, int]:
        start = round(index * extent / count)
        end = round((index + 1) * extent / count) - 1
        return start, end
    for row in rows:
        g, k = float(row["g"]), float(row["k"])
        order = int(float(row["kappa"]))
        if order not in colors:
            raise ValueError(f"Unexpected F7A blocker order: {order}")
        x0, x1 = cell_bounds(g_index[g], len(g_values), width)
        y0, y1 = cell_bounds(k_index[k], len(k_values), height)
        draw.rectangle((x0, y0, x1, y1), fill=colors[order])
    image.save(FIGURES / "policy_atlas_heatmap.png", optimize=True)
    shutil.copy2(source, TABLES / "policy_atlas_source.csv")


def write_header_logos() -> None:
    """Color the supplied official marks for direct placement on the green header."""
    project_root = HERE.parents[2]
    ias_source = Image.open(project_root / "temp" / "poster" / "artifacts" / "logos" / "ias_annual_2026_logo.png").convert("RGBA")
    ias_pixels = ias_source.load()
    for y in range(ias_source.height):
        for x in range(ias_source.width):
            r, g, b, a = ias_pixels[x, y]
            if not a:
                continue
            if r < 12 and 100 <= g <= 150 and b < 90:
                ias_pixels[x, y] = (0, 59, 45, a)
            elif max(r, g, b) < 52:
                ias_pixels[x, y] = (255, 255, 255, a)
    ias_source.save(FIGURES / "ias_annual_2026_header.png")

    uniandes_source = Image.open(project_root / "tmp" / "pdfs" / "uniandes_logo_color-1.png").convert("RGB")
    uniandes_out = Image.new("RGBA", uniandes_source.size)
    src, dst = uniandes_source.load(), uniandes_out.load()
    for y in range(uniandes_source.height):
        for x in range(uniandes_source.width):
            r, g, b = src[x, y]
            alpha = 255 - min(r, g, b)
            if alpha < 12:
                dst[x, y] = (255, 255, 255, 0)
            elif r > 150 and g > 125 and b < 130 and r - g < 115:
                dst[x, y] = (229, 190, 69, alpha)
            else:
                dst[x, y] = (255, 255, 255, alpha)
    uniandes_out.resize((uniandes_out.width * 3, uniandes_out.height * 3), Image.Resampling.LANCZOS).save(
        FIGURES / "uniandes_logo_white_header.png"
    )


def write_editable_poster_plot_data() -> None:
    """Export the exact scenario and phasor-DAE rows used by the two poster plots."""
    intervention_rows = read_csv(RUN / "tables" / "IAS26-FINAL_P4_MC_DATA.csv")
    valid = [r for r in intervention_rows if r["scenario_valid"].strip().upper() == "TRUE"]
    categories = {
        "improved": [r for r in valid if r["intervention_improves_h4"].strip().upper() == "TRUE"],
        "deteriorated": [r for r in valid if r["intervention_deteriorates_h4"].strip().upper() == "TRUE"],
    }
    for category, rows in categories.items():
        path = TABLES / f"intervention_{category}.csv"
        with path.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=("alpha_perp_original", "alpha_perp_intervention"))
            writer.writeheader()
            for row in rows:
                writer.writerow({key: row[key] for key in writer.fieldnames})

    trace_rows = read_csv(RUN / "tables" / "IAS26-FINAL_P4_TDS_TRACES.csv")
    selected = {
        "h4_original": "CANONICAL_CANONICAL_H4_ORIGINAL_D2_A1",
        "h4_retuned": "CANONICAL_CANONICAL_H4_RETUNED_D2_A1",
        "proper_triple": "CANONICAL_CANONICAL_PROPER_30_33_35_D2_A1",
    }
    for key, trajectory in selected.items():
        rows = [r for r in trace_rows if r["trajectory_id"] == trajectory]
        if not rows:
            raise ValueError(f"The required phasor-DAE trace is missing: {trajectory}")
        rows.sort(key=lambda r: float(r["time_s"]))
        with (TABLES / f"tds_{key}.csv").open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=("time_s", "omega_sg34_minus_omega_sg39_pu"))
            writer.writeheader()
            writer.writerows({key: row[key] for key in writer.fieldnames} for row in rows)


def write_f2_plot_data() -> None:
    """Export the aligned, audited local and collective singular-value sweep."""
    source = HERE / "research" / "bnd_h4_mechanism" / "results" / "20260926T162252Z_d0fecb32_ias26_030_f2_closure_v1" / "tables" / "F2_ALIGNED_LOCAL_COLLECTIVE_ALPHA.csv"
    rows = read_csv(source)
    by_gain: dict[float, list[float]] = {}
    collective_by_gain: dict[float, float] = {}
    for row in rows:
        gain = float(row["g"])
        by_gain.setdefault(gain, []).append(float(row["physical_local_sigma_min_I_plus_Mii"]))
        collective_by_gain[gain] = float(row["collective_sigma_min_I_plus_QH"])
    if len(by_gain) < 10:
        raise ValueError("The audited F2 local/collective sweep is incomplete.")
    with (TABLES / "f2_local_collective_closure.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["g", "local_sigma_min", "collective_sigma_min"])
        for gain in sorted(by_gain):
            writer.writerow([f"{gain:.12g}", f"{min(by_gain[gain]):.12g}", f"{collective_by_gain[gain]:.12g}"])
    shutil.copy2(source, TABLES / "F2_ALIGNED_LOCAL_COLLECTIVE_ALPHA.csv")


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    numbers = {r["metric"]: r["value"] for r in read_csv(RUN / "tables" / "POSTER_NUMBERS_FINAL.csv")}
    readiness = json.loads((RUN / "reports" / "IAS2026_POSTER_READINESS.json").read_text(encoding="utf-8"))
    events = {r["event"]: r for r in read_csv(RUN / "tables" / "MC_EVENT_RATES.csv")}
    crosscode = read_csv(RUN / "tables" / "IAS26-FINAL_EVIDENCE_STRIP.csv")[0]
    status = json.loads((HERE / "research" / "bnd_h4_mechanism" / "tickets" / "STATUS.json").read_text(encoding="utf-8"))
    holdout_path = HERE.parents[2] / "research" / "ias2026_final_closure" / "reports" / "P6_GENUINE_HOLDOUT_STATUS.md"
    holdout_text = holdout_path.read_text(encoding="utf-8")
    holdout_match = re.search(r"accuracy:\s*1\.0\s*\((\d+)/(\d+)\)", holdout_text)
    if not holdout_match or "SINGLE_CLASS_HOLDOUT_REVEALED" not in holdout_text:
        raise ValueError("The audited alternative-model holdout is absent or has changed status.")
    holdout_correct, holdout_total = map(int, holdout_match.groups())
    b02 = next(r for r in read_csv(HERE / "research" / "results" / "20260911_CLAIM_MATRIX.csv") if r["id"] == "B02")
    b01 = next(r for r in read_csv(HERE / "research" / "results" / "20260911_CLAIM_MATRIX.csv") if r["id"] == "B01")
    b05 = next(r for r in read_csv(HERE / "research" / "results" / "20260911_CLAIM_MATRIX.csv") if r["id"] == "B05")
    m1_line = next(x["reason"] for x in status["blocked"] if x["ticket"] == "IAS26-010")
    residual, gate = re.search(r"residual ([0-9.eE+-]+) exceeds the frozen ([0-9.eE+-]+) gate", m1_line).groups()
    blocker = re.search(r"Flagship ([0-9+]+) ", b02["claim"]).group(1).split("+")
    policy_counts = [int(x) for x in re.search(r"\((\d+)/(\d+)/(\d+) distinct H", b01["claim"]).groups()]
    policy_path = re.search(r"witness ([^)]*)", b01["claim"]).group(1)
    mw, mva = re.search(r"\(([0-9.]+) MW / ([0-9.]+) MVA\)", b02["claim"]).groups()
    condenser_mva = re.search(r"\(([0-9.]+) MVA at P4\)", b05["claim"]).group(1)
    config = json.loads((HERE / "research" / "bnd_h4_mechanism" / "configs" / "IAS26-050_INTERVENTION_V1.json").read_text(encoding="utf-8"))
    mc_valid = readiness["MC"]["n_valid"]
    mc_total = readiness["MC"]["n_total"]
    mc_invalid = readiness["MC"]["n_infeasible"]
    tds = readiness["TDS"]

    scenario_rows = read_csv(RUN / "derived" / "SCENARIO_METRICS.csv")
    valid_rows = [r for r in scenario_rows if r["scenario_valid"].strip().upper() == "TRUE"]
    h4_stable = sum(r["h4_class"].strip().upper() == "STABLE" for r in valid_rows)
    h4_minimal = sum(r["h4_minimal_blocker_persists"].strip().upper() == "TRUE" for r in valid_rows)
    alternative_witness = sum(
        r["h4_class"].strip().upper() == "UNSTABLE"
        and r["h4_minimal_blocker_persists"].strip().upper() == "FALSE"
        for r in valid_rows
    )
    if len(valid_rows) != mc_valid or len(scenario_rows) != mc_total:
        raise ValueError("The frozen scenario table no longer matches the poster-readiness totals.")
    if h4_stable + h4_minimal + alternative_witness != mc_valid:
        raise ValueError("The frozen H4 scenario classes do not form the audited valid-scenario partition.")
    if h4_minimal != int(events["E1_H4_MINIMAL_BLOCKER_PERSISTS"]["numerator"]):
        raise ValueError("Scenario table and event-rate table disagree on H4 minimality.")

    boundary_path = (
        HERE / "research" / "bnd_h4_mechanism" / "results"
        / "20260926T162252Z_d0fecb32_ias26_030_f2_closure_v1" / "tables" / "F2_BOUNDARY_PORT_AUDIT.csv"
    )
    boundary_rows = read_csv(boundary_path)
    boundary_frequencies = [float(r["frequency_hz"]) for r in boundary_rows]
    boundary_roots = [float(r["g_root"]) for r in boundary_rows]
    if not boundary_rows or max(boundary_frequencies) - min(boundary_frequencies) > 1e-6:
        raise ValueError("The four-port boundary table no longer has a common audited frequency.")
    if max(boundary_roots) - min(boundary_roots) > 1e-6:
        raise ValueError("The four-port boundary table no longer has a common audited gain.")
    boundary_frequency = sum(boundary_frequencies) / len(boundary_frequencies)
    boundary_alpha = sum(float(r["s_real"]) for r in boundary_rows) / len(boundary_rows)
    local_sigma_by_bus = {
        int(r["device_bus"]): float(r["physical_local_sigma_min_I_plus_Mii"])
        for r in boundary_rows
    }

    alpha = float(numbers["canonical original H4 α⊥"])
    # F1's nominal target-band frequency, retained with its scope label in the poster.
    freq = float(config["metrics_and_decisions"]["nominal_pre_run_prediction"]["frequency_original_hz"])

    stable_count = int(re.search(r"(\d+) proper subsets stable", readiness["safe_claims"][0]).group(1))
    m = config["system"]["unchanged_parameters"]
    macro_rows = {
        "FlagshipBlocker": r"\{ " + ", ".join(blocker) + r" \}",
        "FlagshipAlpha": tex_number(alpha, 4, signed=True),
        "FlagshipFrequency": tex_number(freq, 4),
        "NumBlockerUnits": str(len(blocker)),
        "NumProperStable": str(stable_count),
        "NumProperTotal": str(stable_count),
        "BaselineGain": tex_number(float(config["treatment"]["control"]["g"]), 5),
        "RetunedGain": tex_number(float(config["treatment"]["primary_intervention"]["g"]), 2),
        "ControllerK": tex_number(float(m["k"]), 3),
        "ControllerT": tex_number(float(m["t"]), 1),
        "ControllerH": tex_number(float(m["h"]), 1),
        "BoundaryGain": tex_number(sum(boundary_roots) / len(boundary_roots), 6),
        "BoundaryAlpha": tex_number(boundary_alpha, 4),
        "BoundaryFrequency": tex_number(boundary_frequency, 4),
        **{
            f"LocalSigma{letter}": tex_number(value, 6)
            for letter, value in zip("ABCD", local_sigma_by_bus.values())
        },
        **{
            f"LocalSigmaLog{letter}": tex_number(math.log10(value), 6)
            for letter, value in zip("ABCD", local_sigma_by_bus.values())
        },
        "LocalSigmaBoundary": tex_number(float(numbers["minimum physical local σmin at boundary"]), 6),
        "CollectiveSigmaBoundary": tex_sci(numbers["collective σmin(I+QH) at boundary"]),
        "CollectiveSigmaPlotValue": f"{float(numbers['collective σmin(I+QH) at boundary']):.8e}",
        "CollectiveSigmaLog": tex_number(math.log10(float(numbers["collective σmin(I+QH) at boundary"])), 6),
        "DisplacedMW": mw,
        "ConverterMVA": mva,
        "CondenserMVA": condenser_mva,
        "NumValidMC": str(mc_valid),
        "NumScenarios": str(mc_total),
        "NumInfeasibleMC": str(mc_invalid),
        "NumStableHFourScenarios": str(h4_stable),
        "NumAlternativeWitnessScenarios": str(alternative_witness),
        "NumCollectiveBlockerScenarios": str(int(events["E2_COLLECTIVE_PHENOMENON_PERSISTS"]["numerator"])),
        "ValidScenarioPercent": tex_number(100 * mc_valid / mc_total, 1),
        "InfeasibleScenarioPercent": tex_number(100 * mc_invalid / mc_total, 1),
        "ValidBarWidth": tex_number(426 * mc_valid / mc_total, 4),
        "InvalidBarWidth": tex_number(426 * mc_invalid / mc_total, 4),
        "ValidBarCenter": tex_number(213 * mc_valid / mc_total, 4),
        "InvalidBarCenter": tex_number(426 - 213 * mc_invalid / mc_total, 4),
        "HFourMinimalPercent": tex_number(100 * h4_minimal / mc_valid, 1),
        "CollectiveBlockerPercent": tex_number(100 * int(events["E2_COLLECTIVE_PHENOMENON_PERSISTS"]["numerator"]) / mc_valid, 1),
        "AlternativeWitnessPercent": tex_number(100 * alternative_witness / mc_valid, 1),
        "BlockerPersistenceCount": str(int(events["E1_H4_MINIMAL_BLOCKER_PERSISTS"]["numerator"])),
        "BlockerPersistenceRate": tex_number(100*float(events["E1_H4_MINIMAL_BLOCKER_PERSISTS"]["estimate"]), 1),
        "CollectivePersistenceCount": str(int(events["E2_COLLECTIVE_PHENOMENON_PERSISTS"]["numerator"])),
        "InterventionImproves": str(int(events["E3_INTERVENTION_IMPROVES"]["numerator"])),
        "InterventionRescues": str(int(events["E4_INTERVENTION_RESCUES_H4"]["numerator"])),
        "InterventionDeteriorates": str(int(events["E5_INTERVENTION_DETERIORATES"]["numerator"])),
        "InterventionImprovesPercent": tex_number(100 * float(events["E3_INTERVENTION_IMPROVES"]["estimate"]), 1),
        "InterventionRescuesPercent": tex_number(100 * float(events["E4_INTERVENTION_RESCUES_H4"]["estimate"]), 1),
        "InterventionDeterioratesPercent": tex_number(100 * float(events["E5_INTERVENTION_DETERIORATES"]["estimate"]), 1),
        "TDSAgreement": f"{tds['n_sign_agreement']}/{tds['n_completed_finite_primary_fits']}",
        "TDSAlphaMAE": tex_number(tds["alpha_mae_completed_s-1"], 4),
        "BlindHoldoutCorrect": str(holdout_correct),
        "BlindHoldoutTotal": str(holdout_total),
        "TDSCompleted": str(tds["n_completed_finite_primary_fits"]),
        "TDSTotal": str(tds["n_planned"]),
        "TDSSignDisagreements": str(tds["n_sign_disagreement"]),
        "TDSZeroHzDisagreements": str(readiness["falsification_N8"]["mismatches_zero_frequency_n"]),
        "TDSSolverFailures": str(tds["run_status_counts"]["SOLVER_FAILED"]),
        "CrossCodeAgreement": f"{crosscode['verdict_agreement']}/{crosscode['N_case_comparisons']}",
        "NumDistinctPolicyA": str(policy_counts[0]),
        "NumDistinctPolicyB": str(policy_counts[1]),
        "NumDistinctPolicyC": str(policy_counts[2]),
        "PolicyPath": policy_path.replace("->", r"\to").replace("EMPTY", r"\varnothing"),
        "MOneResidual": tex_sci(residual),
        "MOneGate": tex_sci(gate),
        # The current release does not expose a port-ratio residual; do not alias the M1 gate residual.
        "PortResidual": r"\PosterPending",
    }
    macro_rows.update(f10_baseline_macros())
    results_lines = ["% GENERATED from current IAS26 poster audit and cited frozen results. Do not edit by hand."]
    results_lines += [f"\\newcommand{{\\{name}}}{{{value}}}" for name, value in macro_rows.items()]
    (OUT / "results.tex").write_text("\n".join(results_lines) + "\n", encoding="utf-8")

    claims = [
        "% Generated claim-state macros; claim provenance is recorded in claim_audit.md.",
        r"\newcommand{\FlagshipEvidenceState}{VERIFIED}",
        r"\newcommand{\MechanismEvidenceState}{CONDITIONAL}",
        r"\newcommand{\MOneEvidenceState}{BLOCKED}",
        r"\newcommand{\PolicyAtlasEvidenceState}{CONDITIONAL}",
        r"\newcommand{\MonteCarloEvidenceState}{DESCRIPTIVE\ SYNTHETIC\ ENSEMBLE}",
        r"\newcommand{\MonteCarloReviewStatus}{PENDING\ USER\ REVIEW}",
        r"\newcommand{\CrossCodeEvidenceState}{SAME-MODEL\ REPRODUCTION}",
        r"\newcommand{\TdsEvidenceState}{SAME-MODEL\ PHASOR\ DAE}",
        r"\newcommand{\AlternativeConverterState}{NOT\ RUN}",
        r"\newcommand{\GovernorState}{CONDITIONAL}",
        r"\newcommand{\KundurState}{CONDITIONAL}",
        r"\newcommand{\IEEEsixtyEightState}{NOT\ RUN}",
        r"\newcommand{\BlindHoldoutState}{CONDITIONAL}",
        r"\newcommand{\TopologyState}{PENDING}",
        r"\newcommand{\PlanningState}{PENDING}",
        r"\newcommand{\SupportState}{CONDITIONAL}",
    ]
    (OUT / "claims.tex").write_text("\n".join(claims) + "\n", encoding="utf-8")

    with (TABLES / "policy_atlas.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["policy", "distinct_hypergraphs"])
        for policy, count in zip(("F7A", "F7B", "F7C"), policy_counts):
            w.writerow([policy, count])
    shutil.copy2(RUN / "tables" / "MC_EVENT_RATES.csv", TABLES / "MC_EVENT_RATES.csv")
    write_tikz_network([int(x) for x in blocker])
    write_hasse([int(x) for x in blocker])
    write_header_logos()
    write_policy_atlas()
    write_editable_poster_plot_data()
    write_f2_plot_data()
    order_tikz = [
        "% Reserved order-truncation plot geometry. No M1 order-by-order result is released.",
        r"\draw[IASRule,line width=1.5pt] (28,17)--(316,17);",
        r"\draw[IASRule,line width=1.5pt] (28,17)--(28,84);",
    ]
    for x, order in zip((64, 132, 200, 268), (1, 2, 3, 4)):
        order_tikz.append(f"\\draw[IASRule] ({x},14)--({x},20);")
        order_tikz.append(f"\\node[font=\\fontsize{{28pt}}{{30pt}}\\selectfont,text=IASDarkText] at ({x},7) {{{order}}};")
    order_tikz += [
        r"\node[font=\fontsize{26.2pt}{29pt}\selectfont,text=IASDarkText] at (174,1) {retained interaction order};",
        r"\node[font=\fontsize{33pt}{36pt}\selectfont,text=IASGray,align=center] at (174,63) {\UnresolvedMark\\ PENDING};",
        r"\node[font=\fontsize{26.2pt}{29pt}\selectfont,text=IASGray,align=center] at (174,42) {M1 order-by-order outputs\\ not released};",
    ]
    (FIGURES / "order_cliff.tikz").write_text("\n".join(order_tikz) + "\n", encoding="utf-8")
    print(f"Generated {OUT / 'results.tex'} and local IEEE-39/H4 TikZ figures from audited sources.")


if __name__ == "__main__":
    if sys.argv[1:] == ["--f10-baselines-only"]:
        update_f10_baseline_macros()
    else:
        main()
