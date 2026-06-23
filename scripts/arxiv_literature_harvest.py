"""Harvest arXiv metadata for the thesis literature review.

The script queries arXiv's legacy API by topic and year, stores topic/year
metadata, builds a deduplicated global index, and optionally downloads PDFs.

It intentionally defaults to metadata-only harvesting. Full PDF downloads can
be large and should be run overnight with the arXiv rate limit enabled.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any


ARXIV_API_URL = "https://export.arxiv.org/api/query"
ATOM_NS = "{http://www.w3.org/2005/Atom}"
ARXIV_NS = "{http://arxiv.org/schemas/atom}"
OPENSEARCH_NS = "{http://a9.com/-/spec/opensearch/1.1/}"
MAX_API_RESULTS = 30000
MAX_PAGE_SIZE = 2000


def normalize_ws(value: str | None) -> str:
    if not value:
        return ""
    return re.sub(r"\s+", " ", value).strip()


def arxiv_id_parts(entry_id: str) -> tuple[str, str]:
    raw = entry_id.rstrip("/").split("/")[-1]
    base = re.sub(r"v\d+$", "", raw)
    return raw, base


def mkdir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def write_json(path: Path, data: Any) -> None:
    mkdir(path.parent)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")
    tmp.replace(path)


def read_jsonl(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if line:
                rows.append(json.loads(line))
    return rows


def write_jsonl(path: Path, rows: list[dict[str, Any]]) -> None:
    mkdir(path.parent)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8", newline="\n") as handle:
        for row in rows:
            handle.write(json.dumps(row, ensure_ascii=False, sort_keys=True) + "\n")
    tmp.replace(path)


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    mkdir(path.parent)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            flat = dict(row)
            for key, value in list(flat.items()):
                if isinstance(value, list):
                    flat[key] = "; ".join(str(v) for v in value)
                elif isinstance(value, dict):
                    flat[key] = json.dumps(value, sort_keys=True)
            writer.writerow(flat)
    tmp.replace(path)


def parse_year(value: str) -> int | None:
    if not value:
        return None
    try:
        return datetime.fromisoformat(value.replace("Z", "+00:00")).year
    except ValueError:
        return None


def build_query(topic_query: str, year: int) -> str:
    start = f"{year}01010000"
    end = f"{year}12312359"
    return f"({topic_query}) AND submittedDate:[{start} TO {end}]"


def build_api_url(
    query: str,
    start: int,
    max_results: int,
    sort_by: str,
    sort_order: str,
) -> str:
    params = {
        "search_query": query,
        "start": str(start),
        "max_results": str(max_results),
        "sortBy": sort_by,
        "sortOrder": sort_order,
    }
    return ARXIV_API_URL + "?" + urllib.parse.urlencode(params)


def fetch_bytes(url: str, user_agent: str, timeout: float, retries: int) -> bytes:
    last_error: Exception | None = None
    for attempt in range(retries + 1):
        request = urllib.request.Request(url, headers={"User-Agent": user_agent})
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                return response.read()
        except urllib.error.HTTPError as exc:
            last_error = exc
            if exc.code in {429, 500, 502, 503, 504} and attempt < retries:
                time.sleep(10 * (attempt + 1))
                continue
            raise
        except urllib.error.URLError as exc:
            last_error = exc
            if attempt < retries:
                time.sleep(10 * (attempt + 1))
                continue
            raise
    raise RuntimeError(f"failed to fetch {url}") from last_error


def parse_feed(
    payload: bytes,
    topic_slug: str,
    topic_label: str,
    query_year: int,
    query: str,
) -> tuple[int | None, list[dict[str, Any]]]:
    root = ET.fromstring(payload)
    total_node = root.find(f"{OPENSEARCH_NS}totalResults")
    total = int(total_node.text) if total_node is not None and total_node.text else None
    records: list[dict[str, Any]] = []

    for entry in root.findall(f"{ATOM_NS}entry"):
        entry_id = normalize_ws(entry.findtext(f"{ATOM_NS}id"))
        id_version, id_base = arxiv_id_parts(entry_id)
        authors = [
            normalize_ws(author.findtext(f"{ATOM_NS}name"))
            for author in entry.findall(f"{ATOM_NS}author")
        ]
        authors = [author for author in authors if author]

        categories = [
            category.attrib.get("term", "")
            for category in entry.findall(f"{ATOM_NS}category")
            if category.attrib.get("term")
        ]
        primary_category_node = entry.find(f"{ARXIV_NS}primary_category")
        primary_category = ""
        if primary_category_node is not None:
            primary_category = primary_category_node.attrib.get("term", "")

        links: list[dict[str, str]] = []
        pdf_url = ""
        for link in entry.findall(f"{ATOM_NS}link"):
            link_data = {k: v for k, v in link.attrib.items()}
            links.append(link_data)
            if link.attrib.get("title") == "pdf" or link.attrib.get("type") == "application/pdf":
                pdf_url = link.attrib.get("href", "")

        published = normalize_ws(entry.findtext(f"{ATOM_NS}published"))
        updated = normalize_ws(entry.findtext(f"{ATOM_NS}updated"))
        records.append(
            {
                "arxiv_id": id_base,
                "arxiv_id_version": id_version,
                "entry_id": entry_id,
                "abs_url": f"https://arxiv.org/abs/{id_version}",
                "pdf_url": pdf_url or f"https://arxiv.org/pdf/{id_version}",
                "title": normalize_ws(entry.findtext(f"{ATOM_NS}title")),
                "authors": authors,
                "summary": normalize_ws(entry.findtext(f"{ATOM_NS}summary")),
                "comment": normalize_ws(entry.findtext(f"{ARXIV_NS}comment")),
                "journal_ref": normalize_ws(entry.findtext(f"{ARXIV_NS}journal_ref")),
                "doi": normalize_ws(entry.findtext(f"{ARXIV_NS}doi")),
                "primary_category": primary_category,
                "categories": categories,
                "published": published,
                "published_year": parse_year(published),
                "updated": updated,
                "updated_year": parse_year(updated),
                "links": links,
                "query_topic": topic_slug,
                "query_topic_label": topic_label,
                "query_year": query_year,
                "query": query,
                "harvested_at_utc": datetime.now(UTC).isoformat(),
            }
        )
    return total, records


def download_pdf(
    record: dict[str, Any],
    target_path: Path,
    user_agent: str,
    timeout: float,
    retries: int,
) -> bool:
    if target_path.exists() and target_path.stat().st_size > 0:
        return False
    mkdir(target_path.parent)
    payload = fetch_bytes(record["pdf_url"], user_agent, timeout, retries)
    tmp = target_path.with_suffix(".pdf.tmp")
    tmp.write_bytes(payload)
    tmp.replace(target_path)
    return True


@dataclass(frozen=True)
class HarvestArgs:
    config: Path
    output_dir: Path
    from_year: int
    to_year: int
    max_results_per_topic_year: int
    page_size: int
    delay_seconds: float
    sort_by: str
    sort_order: str
    user_agent: str
    timeout: float
    retries: int
    resume: bool
    download_pdfs: bool
    pdf_limit_per_topic_year: int
    confirm_mass_download: bool
    write_abstracts: bool


def parse_args() -> HarvestArgs:
    parser = argparse.ArgumentParser(
        description="Harvest arXiv metadata by thesis topic and year."
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("configs/research/arxiv_topics.json"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/arxiv_literature"),
    )
    parser.add_argument("--from-year", type=int, default=2010)
    parser.add_argument(
        "--to-year",
        type=int,
        default=datetime.now(UTC).year,
    )
    parser.add_argument(
        "--max-results-per-topic-year",
        type=int,
        default=1000,
        help="Use 30000 for the maximum arXiv API cap per topic-year query.",
    )
    parser.add_argument("--page-size", type=int, default=1000)
    parser.add_argument(
        "--delay-seconds",
        type=float,
        default=3.2,
        help="arXiv requires no more than one legacy API request every 3 seconds.",
    )
    parser.add_argument(
        "--sort-by",
        choices=["relevance", "lastUpdatedDate", "submittedDate"],
        default="submittedDate",
    )
    parser.add_argument(
        "--sort-order",
        choices=["ascending", "descending"],
        default="descending",
    )
    parser.add_argument(
        "--user-agent",
        default="spectral-ibr-literature-review/0.1 (mailto:research@example.com)",
    )
    parser.add_argument("--timeout", type=float, default=60.0)
    parser.add_argument("--retries", type=int, default=3)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--download-pdfs", action="store_true")
    parser.add_argument(
        "--pdf-limit-per-topic-year",
        type=int,
        default=0,
        help="0 means no PDF cap when --download-pdfs is used.",
    )
    parser.add_argument(
        "--confirm-mass-download",
        action="store_true",
        help="Required for uncapped PDF downloads.",
    )
    parser.add_argument("--write-abstracts", action="store_true")
    ns = parser.parse_args()

    if ns.max_results_per_topic_year > MAX_API_RESULTS:
        parser.error(f"--max-results-per-topic-year cannot exceed {MAX_API_RESULTS}")
    if ns.page_size > MAX_PAGE_SIZE:
        parser.error(f"--page-size cannot exceed {MAX_PAGE_SIZE}")
    if ns.delay_seconds < 3.0:
        parser.error("Use --delay-seconds >= 3.0 to respect arXiv API policy")
    if (
        ns.download_pdfs
        and ns.pdf_limit_per_topic_year <= 0
        and not ns.confirm_mass_download
    ):
        parser.error(
            "Uncapped PDF downloads require --confirm-mass-download. "
            "Use --pdf-limit-per-topic-year for a bounded run."
        )
    return HarvestArgs(**vars(ns))


def flatten_for_index(record: dict[str, Any]) -> dict[str, Any]:
    return {
        "arxiv_id": record["arxiv_id"],
        "arxiv_id_version": record["arxiv_id_version"],
        "published": record["published"],
        "published_year": record["published_year"],
        "updated": record["updated"],
        "title": record["title"],
        "authors": record["authors"],
        "primary_category": record["primary_category"],
        "categories": record["categories"],
        "doi": record["doi"],
        "journal_ref": record["journal_ref"],
        "abs_url": record["abs_url"],
        "pdf_url": record["pdf_url"],
        "matched_topics": record.get("matched_topics", []),
        "matched_topic_years": record.get("matched_topic_years", []),
        "summary": record["summary"],
    }


def main() -> int:
    args = parse_args()
    config = json.loads(args.config.read_text(encoding="utf-8"))
    topics = config["topics"]
    run_id = datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")

    mkdir(args.output_dir)
    mkdir(args.output_dir / "index")
    mkdir(args.output_dir / "runs" / run_id)

    all_topic_records: list[dict[str, Any]] = []
    topic_year_summary: list[dict[str, Any]] = []
    request_count = 0
    pdf_download_count = 0

    for topic in topics:
        topic_slug = topic["slug"]
        topic_label = topic["label"]
        for year in range(args.from_year, args.to_year + 1):
            topic_year_dir = args.output_dir / "topics" / topic_slug / str(year)
            metadata_path = topic_year_dir / "metadata.jsonl"
            query_path = topic_year_dir / "query.json"

            if args.resume and metadata_path.exists():
                records = read_jsonl(metadata_path)
                total = len(records)
                source = "resume"
                print(f"[resume] {topic_slug} {year}: {len(records)} records")
            else:
                query = build_query(topic["query"], year)
                write_json(
                    query_path,
                    {
                        "topic": topic,
                        "year": year,
                        "query": query,
                        "sort_by": args.sort_by,
                        "sort_order": args.sort_order,
                    },
                )
                records = []
                total = None
                max_results = args.max_results_per_topic_year
                start = 0
                source = "api"
                while start < max_results:
                    page_size = min(args.page_size, max_results - start)
                    url = build_api_url(
                        query=query,
                        start=start,
                        max_results=page_size,
                        sort_by=args.sort_by,
                        sort_order=args.sort_order,
                    )
                    print(
                        f"[api] {topic_slug} {year}: "
                        f"start={start} page_size={page_size}",
                        flush=True,
                    )
                    payload = fetch_bytes(
                        url=url,
                        user_agent=args.user_agent,
                        timeout=args.timeout,
                        retries=args.retries,
                    )
                    request_count += 1
                    page_total, page_records = parse_feed(
                        payload=payload,
                        topic_slug=topic_slug,
                        topic_label=topic_label,
                        query_year=year,
                        query=query,
                    )
                    if total is None:
                        total = page_total
                    if not page_records:
                        break
                    records.extend(page_records)
                    start += len(page_records)
                    if total is not None and start >= total:
                        break
                    if start >= max_results:
                        break
                    time.sleep(args.delay_seconds)

                write_jsonl(metadata_path, records)

            if args.write_abstracts:
                abstracts_dir = topic_year_dir / "abstracts"
                mkdir(abstracts_dir)
                for record in records:
                    abstract_path = abstracts_dir / f"{record['arxiv_id']}.txt"
                    if not abstract_path.exists():
                        abstract_path.write_text(
                            (
                                f"{record['title']}\n\n"
                                f"Authors: {', '.join(record['authors'])}\n"
                                f"arXiv: {record['arxiv_id_version']}\n\n"
                                f"{record['summary']}\n"
                            ),
                            encoding="utf-8",
                        )

            if args.download_pdfs:
                limit = args.pdf_limit_per_topic_year
                pdf_records = records if limit <= 0 else records[:limit]
                for record in pdf_records:
                    pdf_path = topic_year_dir / "pdfs" / f"{record['arxiv_id']}.pdf"
                    changed = download_pdf(
                        record=record,
                        target_path=pdf_path,
                        user_agent=args.user_agent,
                        timeout=args.timeout,
                        retries=args.retries,
                    )
                    if changed:
                        pdf_download_count += 1
                    time.sleep(args.delay_seconds)

            all_topic_records.extend(records)
            topic_year_summary.append(
                {
                    "topic_slug": topic_slug,
                    "topic_label": topic_label,
                    "query_year": year,
                    "source": source,
                    "api_total_results": total,
                    "stored_records": len(records),
                    "metadata_path": str(metadata_path),
                }
            )
            time.sleep(args.delay_seconds)

    unique: dict[str, dict[str, Any]] = {}
    for record in all_topic_records:
        key = record["arxiv_id"]
        if key not in unique:
            merged = dict(record)
            merged["matched_topics"] = []
            merged["matched_topic_years"] = []
            unique[key] = merged
        topic_hit = record["query_topic"]
        topic_year_hit = f"{record['query_topic']}:{record['query_year']}"
        if topic_hit not in unique[key]["matched_topics"]:
            unique[key]["matched_topics"].append(topic_hit)
        if topic_year_hit not in unique[key]["matched_topic_years"]:
            unique[key]["matched_topic_years"].append(topic_year_hit)

    index_rows = [flatten_for_index(row) for row in unique.values()]
    index_rows.sort(
        key=lambda row: (
            row["published_year"] or 0,
            row["published"],
            row["arxiv_id"],
        ),
        reverse=True,
    )

    by_year = Counter(row.get("published_year") for row in index_rows)
    by_primary_category = Counter(row.get("primary_category", "") for row in index_rows)
    by_topic = Counter()
    for row in index_rows:
        by_topic.update(row.get("matched_topics", []))

    write_jsonl(args.output_dir / "index" / "all_metadata.jsonl", index_rows)
    write_csv(
        args.output_dir / "index" / "all_metadata.csv",
        index_rows,
        [
            "arxiv_id",
            "arxiv_id_version",
            "published",
            "published_year",
            "updated",
            "title",
            "authors",
            "primary_category",
            "categories",
            "doi",
            "journal_ref",
            "abs_url",
            "pdf_url",
            "matched_topics",
            "matched_topic_years",
            "summary",
        ],
    )
    write_csv(
        args.output_dir / "index" / "summary_by_topic_year.csv",
        topic_year_summary,
        [
            "topic_slug",
            "topic_label",
            "query_year",
            "source",
            "api_total_results",
            "stored_records",
            "metadata_path",
        ],
    )

    manifest = {
        "run_id": run_id,
        "created_at_utc": datetime.now(UTC).isoformat(),
        "config": str(args.config),
        "output_dir": str(args.output_dir),
        "from_year": args.from_year,
        "to_year": args.to_year,
        "topics": [topic["slug"] for topic in topics],
        "request_count": request_count,
        "topic_year_records": len(all_topic_records),
        "unique_records": len(index_rows),
        "pdf_download_count": pdf_download_count,
        "by_year": {str(k): v for k, v in sorted(by_year.items())},
        "by_primary_category": dict(by_primary_category.most_common()),
        "by_topic": dict(by_topic.most_common()),
        "arxiv_policy": {
            "delay_seconds": args.delay_seconds,
            "page_size": args.page_size,
            "max_results_per_topic_year": args.max_results_per_topic_year,
            "note": "arXiv legacy API policy: one request every three seconds, single connection.",
        },
    }
    write_json(args.output_dir / "manifest.json", manifest)
    write_json(args.output_dir / "runs" / run_id / "manifest.json", manifest)

    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
