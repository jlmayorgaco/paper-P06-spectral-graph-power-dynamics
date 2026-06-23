# arXiv literature harvest

This folder documents the reproducible arXiv harvest used to build the
literature-review corpus for the spectral graph, IBR stability, rational
self-energy, and dynamic hosting thesis.

## Policy guardrail

The script follows the arXiv legacy API policy: one connection and at least one
request every three seconds. The API can return up to 30,000 results per query
in pages of at most 2,000 records, but broad queries should be narrowed when
they return very large result sets.

Official policy references:

- https://info.arxiv.org/help/api/user-manual.html
- https://info.arxiv.org/help/api/tou.html

## Files

- `configs/research/arxiv_topics.json`: topic families and arXiv queries.
- `scripts/arxiv_literature_harvest.py`: metadata/PDF harvester.
- `data/arxiv_literature/topics/<topic>/<year>/metadata.jsonl`: raw topic-year
  metadata.
- `data/arxiv_literature/index/all_metadata.jsonl`: deduplicated global index.
- `data/arxiv_literature/index/all_metadata.csv`: spreadsheet-friendly index.
- `data/arxiv_literature/index/summary_by_topic_year.csv`: per-query coverage.
- `data/arxiv_literature/manifest.json`: run summary.

## Smoke test

```powershell
python scripts/arxiv_literature_harvest.py `
  --from-year 2025 `
  --to-year 2026 `
  --max-results-per-topic-year 20 `
  --page-size 20 `
  --write-abstracts
```

## Deep metadata harvest

This is the recommended first full pass. It gathers a large corpus without
downloading PDFs.

```powershell
python scripts/arxiv_literature_harvest.py `
  --from-year 2010 `
  --to-year 2026 `
  --max-results-per-topic-year 1000 `
  --page-size 1000 `
  --write-abstracts `
  --resume
```

## Maximum topic-year metadata harvest

This respects the arXiv API cap per query, but it can take a long time.

```powershell
python scripts/arxiv_literature_harvest.py `
  --from-year 2010 `
  --to-year 2026 `
  --max-results-per-topic-year 30000 `
  --page-size 2000 `
  --write-abstracts `
  --resume
```

## Bounded PDF download

Use this after inspecting `all_metadata.csv`. This downloads a bounded number
of PDFs per topic-year and avoids duplicate-heavy uncontrolled downloads.

```powershell
python scripts/arxiv_literature_harvest.py `
  --from-year 2020 `
  --to-year 2026 `
  --max-results-per-topic-year 500 `
  --page-size 500 `
  --download-pdfs `
  --pdf-limit-per-topic-year 10 `
  --resume
```

Uncapped PDF downloads require `--confirm-mass-download` by design.

## Ranked PDF download

After running the screening ranker, download the highest-scoring PDFs by topic:

```powershell
python scripts/arxiv_literature_rank.py --top-n 50

python scripts/arxiv_download_ranked_pdfs.py `
  --max-per-topic 10
```

The PDFs are stored under:

```text
data/arxiv_literature/pdfs/top_ranked/<topic>/<year>/<arxiv_id>.pdf
```

The first verified run used `--max-per-topic 5`, yielding 45 PDFs with zero
download errors.
