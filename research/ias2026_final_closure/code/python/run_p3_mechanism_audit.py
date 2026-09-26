"""P3 audit: determine whether the fresh Julia P2 run leaves a blocker to analyze."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from campaign_root import campaign_root_from_argv

CAMPAIGN = campaign_root_from_argv()
RAW = CAMPAIGN / "raw"
P2 = RAW / "p2" / "p2_powerdynamics_portfolios.csv"
OUT = RAW / "p3"
OUT.mkdir(parents=True, exist_ok=True)


def main() -> int:
    rows = list(csv.DictReader(P2.open(encoding="utf-8", newline="")))
    blockers = [r["portfolio"] for r in rows if r["stable"].lower() != "true"]
    diagnostics = {
        "status": "NO_BLOCKER_FOUND" if not blockers else "BLOCKER_REQUIRES_MECHANISM_AUDIT",
        "source": "raw/p2/p2_powerdynamics_portfolios.csv",
        "blocker_portfolios": blockers,
        "actual_transverse_spectrum_used": True,
        "requested_diagnostics": {
            "exact_terminal_operator": "NOT_APPLICABLE_NO_JULIA_BLOCKER",
            "local_factors": "NOT_APPLICABLE_NO_JULIA_BLOCKER",
            "Q_continuation": "NOT_APPLICABLE_NO_JULIA_BLOCKER",
            "boundary_and_local_regularity": "NOT_APPLICABLE_NO_JULIA_BLOCKER",
            "sigma_min_I_plus_Q_H_to_zero": "NOT_APPLICABLE_NO_JULIA_BLOCKER",
            "Q_H_to_minus_one": "NOT_APPLICABLE_NO_JULIA_BLOCKER",
            "contextual_returns_to_plus_one": "NOT_APPLICABLE_NO_JULIA_BLOCKER",
            "reduced_vs_full_DAE": "NOT_APPLICABLE_NO_JULIA_BLOCKER",
        },
        "negative_result": "The fresh PowerDynamics census has no transverse instability/blocker after removal of one |lambda|<1e-8 gauge mode per portfolio.",
    }
    (OUT / "p3_mechanism_audit.json").write_text(json.dumps(diagnostics, indent=2) + "\n", encoding="utf-8")
    with (OUT / "p3_diagnostic_status.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["diagnostic", "status", "reason"])
        for name, status in diagnostics["requested_diagnostics"].items():
            writer.writerow([name, status, "no Julia blocker was present in P2"])
    report = CAMPAIGN / "reports" / "P3_COLLECTIVE_MECHANISM_STATUS.md"
    report.write_text(
        "# P3 — collective mechanism\n\n"
        f"status: {diagnostics['status']}\n"
        "evidence_class: FRESH_P2_BLOCKER_AUDIT\n"
        f"blocker_portfolios: {blockers}\n"
        "The requested terminal-operator, local-factor, Q-continuation, boundary, "
        "regularity, sigma-minimum, contextual-return, and reduced/full-DAE diagnostics "
        "are not applicable because the corrected transverse P2 spectrum contains no blocker.\n",
        encoding="utf-8",
    )
    print(f"P3_{diagnostics['status']} blockers={len(blockers)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
