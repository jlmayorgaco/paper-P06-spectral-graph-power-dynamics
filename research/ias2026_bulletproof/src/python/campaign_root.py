"""Locate the campaign root in the source tree or a freshly extracted bundle."""

from __future__ import annotations

import sys
from pathlib import Path


def resolve_campaign_root(explicit: str | Path | None = None) -> Path:
    if explicit is not None:
        root = Path(explicit).expanduser().resolve()
        if not (root / "raw").is_dir() or not (root / "reports").is_dir():
            raise FileNotFoundError(f"not an IAS2026 campaign root: {root}")
        return root
    start = Path(__file__).resolve()
    for parent in (start, *start.parents):
        if (parent / "raw").is_dir() and (parent / "reports").is_dir():
            if (parent / "src").is_dir() or (parent / "code").is_dir():
                return parent
    raise FileNotFoundError("could not locate campaign root from script path")


def campaign_root_from_argv() -> Path:
    args = sys.argv
    if "--campaign-root" in args:
        index = args.index("--campaign-root")
        if index + 1 >= len(args):
            raise ValueError("--campaign-root requires a path")
        return resolve_campaign_root(args[index + 1])
    return resolve_campaign_root()
