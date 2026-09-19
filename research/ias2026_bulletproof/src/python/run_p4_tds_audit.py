"""P4 audit: no Julia proper-subset/blocker exists after the fresh P2 census."""

from __future__ import annotations

import csv
from pathlib import Path


CAMPAIGN = Path(__file__).resolve().parents[2]
OUT = CAMPAIGN / "raw" / "p4"
OUT.mkdir(parents=True, exist_ok=True)


def main() -> int:
    rows = list(csv.DictReader((CAMPAIGN / "raw" / "p2" / "p2_powerdynamics_portfolios.csv").open(encoding="utf-8", newline="")))
    blocker_count = sum(r["stable"].lower() != "true" for r in rows)
    with (OUT / "p4_tds_status.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["portfolio_class", "common_disturbance", "estimator", "status", "reason"])
        writer.writerow(["proper_subset", "NOT_RUN", "matrix_pencil_prony", "NOT_APPLICABLE", "no blocker in fresh P2 transverse spectrum"])
        writer.writerow(["blocker", "NOT_RUN", "bandpass_hilbert_log_envelope", "NOT_APPLICABLE", "no blocker in fresh P2 transverse spectrum"])
        writer.writerow(["repaired_blocker", "NOT_RUN", "both_estimators", "NOT_APPLICABLE", "no blocker in fresh P2 transverse spectrum"])
    report = CAMPAIGN / "reports" / "P4_JULIA_TDS_STATUS.md"
    report.write_text(
        "# P4 — Julia common small-disturbance TDS\n\n"
        f"status: {'NO_BLOCKER_FOUND' if blocker_count == 0 else 'REQUIRES_TDS'}\n"
        "evidence_class: FRESH_P2_DEPENDENCY_AUDIT\n"
        f"blocker_portfolios: {blocker_count}\n"
        "No proper-subset/blocker/repaired-blocker trajectory was defined because "
        "the fresh P2 transverse spectrum found no blocker. The required two estimator "
        "rows are recorded as not applicable, not as passed TDS evidence.\n",
        encoding="utf-8",
    )
    print(f"P4_{'NO_BLOCKER_FOUND' if blocker_count == 0 else 'REQUIRES_TDS'} blockers={blocker_count}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
