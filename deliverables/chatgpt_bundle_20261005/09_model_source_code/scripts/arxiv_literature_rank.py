"""Rank harvested arXiv records for literature-review screening.

This script is deliberately simple and transparent: it applies keyword-based
screening to the harvested metadata and emits ranked candidate lists by topic.
The output is not a final citation list; it is a reproducible first pass for
human screening.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter
from pathlib import Path
from typing import Any


CORE_WEIGHTS = {
    "power system": 5,
    "power grid": 5,
    "power network": 5,
    "inverter": 5,
    "converter": 5,
    "grid-forming": 6,
    "grid following": 6,
    "grid-following": 6,
    "inverter-based": 6,
    "low inertia": 5,
    "pll": 5,
    "phase locked loop": 5,
    "small-signal stability": 5,
    "damping": 4,
    "eigenvalue": 4,
    "nonlinear eigenvalue": 5,
    "schur complement": 5,
    "matrix polynomial": 5,
    "self-energy": 6,
    "self energy": 6,
    "resolvent": 5,
    "pseudospectrum": 4,
    "nonnormal": 4,
    "non-normal": 4,
    "passivity": 4,
    "small gain": 4,
    "small phase": 4,
    "spectral graph": 5,
    "laplacian": 3,
    "algebraic connectivity": 4,
    "effective resistance": 4,
    "kirchhoff": 4,
    "hosting capacity": 6,
    "dynamic hosting": 6,
    "virtual inertia": 6,
    "inertia allocation": 5,
    "inertia placement": 5,
    "braess": 6,
    "line reinforcement": 5,
    "transmission expansion": 4,
    "krein": 5,
    "sign characteristic": 6,
    "dissipation-induced": 5,
    "indefinite damping": 5,
    "impedance": 4,
    "loewner": 4,
    "modal identification": 5,
    "frequency response": 3,
}

TOPIC_BOOSTS = {
    "power_system_ibr_stability": [
        "grid-forming",
        "grid-following",
        "inverter-based",
        "pll",
        "weak grid",
        "low inertia",
    ],
    "schur_nep_self_energy": [
        "nonlinear eigenvalue",
        "rational eigenvalue",
        "schur complement",
        "matrix polynomial",
        "resolvent",
    ],
    "spectral_graph_power_networks": [
        "spectral graph",
        "laplacian",
        "algebraic connectivity",
        "effective resistance",
        "power grid",
    ],
    "dynamic_hosting_capacity": [
        "hosting capacity",
        "dynamic hosting",
        "distributed energy resources",
        "inverter hosting",
    ],
    "virtual_inertia_damping_placement": [
        "virtual inertia",
        "inertia allocation",
        "inertia placement",
        "damping placement",
        "rocof",
    ],
    "braess_reinforcement_power_networks": [
        "braess",
        "line reinforcement",
        "transmission expansion",
        "topology control",
        "power grid",
    ],
    "networked_control_resolvent_passivity": [
        "networked control",
        "resolvent",
        "nonnormal",
        "pseudospectrum",
        "passivity",
        "small gain",
        "small phase",
    ],
    "krein_sign_dissipation_matrix_polynomials": [
        "krein",
        "sign characteristic",
        "matrix polynomial",
        "dissipation-induced",
        "indefinite damping",
    ],
    "blackbox_converter_modal_identification": [
        "impedance",
        "black-box",
        "converter",
        "loewner",
        "modal identification",
        "frequency response",
    ],
}


def normalize(value: str | None) -> str:
    if not value:
        return ""
    return re.sub(r"\s+", " ", value).strip()


def load_jsonl(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if line:
                rows.append(json.loads(line))
    return rows


def count_phrase(text: str, phrase: str) -> int:
    return text.count(phrase.lower())


def score_record(record: dict[str, Any], topic: str) -> tuple[int, list[str]]:
    text = " ".join(
        [
            normalize(record.get("title")),
            normalize(record.get("summary")),
            " ".join(record.get("categories", [])),
            normalize(record.get("primary_category")),
        ]
    ).lower()
    score = 0
    hits: list[str] = []
    for phrase, weight in CORE_WEIGHTS.items():
        n_hits = count_phrase(text, phrase)
        if n_hits:
            score += weight * min(n_hits, 3)
            hits.append(phrase)
    for phrase in TOPIC_BOOSTS.get(topic, []):
        if phrase.lower() in text:
            score += 5
            if phrase not in hits:
                hits.append(phrase)

    categories = set(record.get("categories", []))
    primary = record.get("primary_category", "")
    if "eess.SY" in categories or primary == "eess.SY":
        score += 8
    if {"math.OC", "cs.SY", "math.NA", "math.SP"} & categories:
        score += 3
    if record.get("published_year") and record["published_year"] >= 2020:
        score += 2
    if record.get("doi"):
        score += 1
    return score, hits[:12]


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            flat = dict(row)
            for key, value in list(flat.items()):
                if isinstance(value, list):
                    flat[key] = "; ".join(str(v) for v in value)
            writer.writerow(flat)


def markdown_escape(value: str) -> str:
    return value.replace("|", "\\|")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("data/arxiv_literature/index/all_metadata.jsonl"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/arxiv_literature/analysis"),
    )
    parser.add_argument("--top-n", type=int, default=40)
    args = parser.parse_args()

    rows = load_jsonl(args.input)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    ranked_rows: list[dict[str, Any]] = []
    topic_counters: Counter[str] = Counter()
    year_counters: Counter[int] = Counter()

    for record in rows:
        year = record.get("published_year")
        if year:
            year_counters[int(year)] += 1
        for topic in record.get("matched_topics", []):
            topic_counters[topic] += 1
            score, hits = score_record(record, topic)
            ranked_rows.append(
                {
                    "topic": topic,
                    "score": score,
                    "hits": hits,
                    "arxiv_id": record.get("arxiv_id"),
                    "published_year": record.get("published_year"),
                    "title": record.get("title"),
                    "authors": record.get("authors", []),
                    "primary_category": record.get("primary_category"),
                    "categories": record.get("categories", []),
                    "doi": record.get("doi"),
                    "abs_url": record.get("abs_url"),
                    "summary": record.get("summary"),
                }
            )

    ranked_rows.sort(
        key=lambda row: (
            row["topic"],
            row["score"],
            row["published_year"] or 0,
            row["arxiv_id"] or "",
        ),
        reverse=True,
    )

    top_rows: list[dict[str, Any]] = []
    for topic in sorted(topic_counters):
        topic_rows = [row for row in ranked_rows if row["topic"] == topic]
        top_rows.extend(topic_rows[: args.top_n])

    write_csv(
        args.output_dir / "top_candidates_by_topic.csv",
        top_rows,
        [
            "topic",
            "score",
            "hits",
            "arxiv_id",
            "published_year",
            "title",
            "authors",
            "primary_category",
            "categories",
            "doi",
            "abs_url",
            "summary",
        ],
    )

    write_csv(
        args.output_dir / "topic_counts.csv",
        [{"topic": key, "records": value} for key, value in topic_counters.items()],
        ["topic", "records"],
    )
    write_csv(
        args.output_dir / "year_counts.csv",
        [{"year": key, "records": value} for key, value in sorted(year_counters.items())],
        ["year", "records"],
    )

    lines = [
        "# arXiv screening report",
        "",
        f"Input records: {len(rows)}",
        f"Top candidates per topic: {args.top_n}",
        "",
        "## Topic coverage",
        "",
        "| Topic | Records |",
        "|---|---:|",
    ]
    for topic, count in topic_counters.most_common():
        lines.append(f"| `{topic}` | {count} |")

    lines.extend(["", "## Top candidates", ""])
    for topic in sorted(topic_counters):
        lines.extend([f"### {topic}", "", "| Score | Year | Title | arXiv | Hits |", "|---:|---:|---|---|---|"])
        topic_rows = [row for row in top_rows if row["topic"] == topic]
        for row in topic_rows[: min(args.top_n, 15)]:
            title = markdown_escape(normalize(row["title"]))
            hits = markdown_escape(", ".join(row["hits"]))
            lines.append(
                f"| {row['score']} | {row['published_year']} | {title} | "
                f"[{row['arxiv_id']}]({row['abs_url']}) | {hits} |"
            )
        lines.append("")

    (args.output_dir / "screening_report.md").write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )
    print(f"wrote {args.output_dir / 'screening_report.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
