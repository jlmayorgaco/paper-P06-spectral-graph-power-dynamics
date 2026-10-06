"""Download PDFs for ranked arXiv literature candidates.

Use this after `scripts/arxiv_literature_rank.py`. It deduplicates arXiv IDs,
respects a delay between downloads, and stores PDFs by topic/year.
"""

from __future__ import annotations

import argparse
import csv
import re
import time
import urllib.error
import urllib.request
from collections import defaultdict
from pathlib import Path


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")


def fetch_bytes(url: str, user_agent: str, timeout: float, retries: int) -> bytes:
    last_error: Exception | None = None
    for attempt in range(retries + 1):
        request = urllib.request.Request(url, headers={"User-Agent": user_agent})
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                return response.read()
        except (urllib.error.HTTPError, urllib.error.URLError) as exc:
            last_error = exc
            if attempt < retries:
                time.sleep(10 * (attempt + 1))
                continue
            raise
    raise RuntimeError(f"failed to fetch {url}") from last_error


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("data/arxiv_literature/analysis/top_candidates_by_topic.csv"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/arxiv_literature/pdfs/top_ranked"),
    )
    parser.add_argument("--max-per-topic", type=int, default=10)
    parser.add_argument("--delay-seconds", type=float, default=3.2)
    parser.add_argument("--timeout", type=float, default=90.0)
    parser.add_argument("--retries", type=int, default=3)
    parser.add_argument(
        "--user-agent",
        default="spectral-ibr-literature-review/0.1 (mailto:research@example.com)",
    )
    args = parser.parse_args()

    if args.delay_seconds < 3.0:
        parser.error("Use --delay-seconds >= 3.0 to be polite to arXiv")

    by_topic: dict[str, list[dict[str, str]]] = defaultdict(list)
    with args.input.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            by_topic[row["topic"]].append(row)

    selected: list[dict[str, str]] = []
    seen_ids: set[str] = set()
    for topic, rows in sorted(by_topic.items()):
        rows.sort(key=lambda row: int(row.get("score") or 0), reverse=True)
        topic_count = 0
        for row in rows:
            arxiv_id = row["arxiv_id"]
            if arxiv_id in seen_ids:
                continue
            seen_ids.add(arxiv_id)
            selected.append(row)
            topic_count += 1
            if topic_count >= args.max_per_topic:
                break

    args.output_dir.mkdir(parents=True, exist_ok=True)
    downloaded = 0
    skipped = 0
    errors: list[dict[str, str]] = []

    for index, row in enumerate(selected, start=1):
        arxiv_id = row["arxiv_id"]
        topic = safe_name(row["topic"])
        year = row.get("published_year") or "unknown"
        target = args.output_dir / topic / str(year) / f"{safe_name(arxiv_id)}.pdf"
        if target.exists() and target.stat().st_size > 0:
            skipped += 1
            continue
        target.parent.mkdir(parents=True, exist_ok=True)
        pdf_url = f"https://arxiv.org/pdf/{arxiv_id}"
        print(f"[{index}/{len(selected)}] {arxiv_id} -> {target}", flush=True)
        try:
            payload = fetch_bytes(pdf_url, args.user_agent, args.timeout, args.retries)
            tmp = target.with_suffix(".pdf.tmp")
            tmp.write_bytes(payload)
            tmp.replace(target)
            downloaded += 1
        except Exception as exc:  # noqa: BLE001
            errors.append({"arxiv_id": arxiv_id, "topic": topic, "error": str(exc)})
        time.sleep(args.delay_seconds)

    if errors:
        error_path = args.output_dir / "download_errors.csv"
        with error_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["arxiv_id", "topic", "error"])
            writer.writeheader()
            writer.writerows(errors)

    print(
        f"selected={len(selected)} downloaded={downloaded} "
        f"skipped={skipped} errors={len(errors)}"
    )
    return 0 if not errors else 2


if __name__ == "__main__":
    raise SystemExit(main())
