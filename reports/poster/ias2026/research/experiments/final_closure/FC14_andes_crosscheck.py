"""FC14 (step 17): ANDES cross-checks, only where the comparison is meaningful.

Run with .venv/tx3-andes (ANDES 2.0.0; the frozen analysis venv is not touched).

What ANDES can independently check:
  1. the AC power flow of the documented case: bus voltages and every generator's
     P and Q, hence the MW / Mvar numbers of the unit audit;
  2. the base-case electromechanical modes of the documented dynamic models
     (GENROU 6th order, full IEEEX1, IEEEST, TGOV1N). This is NOT equation-
     equivalent to L0 (two-axis, harmonized first-order AVR, no governor), so the
     comparison is qualitative;
  3. the direction of the governor effect: ANDES eigenvalues with TGOV1N in
     service vs switched off (u = 0), an implementation-independent test of the
     governed-model finding (FC03).

What ANDES cannot check: anything with the grid-following converter of L0. ANDES
does not implement those equations, so no portfolio result is claimed as
ANDES-validated.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[6]
SRC = REPO / "data" / "raw" / "ieee39_full.xlsx"
RUN = Path(
    (REPO / "reports/poster/ias2026/research/outputs/ias2026/FINAL_CLOSURE_CURRENT_RUN")
    .read_text(encoding="utf-8")
    .strip()
)
OUT = RUN / "FC14_andes_crosscheck"
OUT.mkdir(parents=True, exist_ok=True)


def run(governors: bool):
    import andes

    andes.config_logger(stream_level=40)
    ss = andes.load(str(SRC), setup=False, no_output=True, default_config=True)
    if not governors:
        # u = 0 makes the TGOV1N initialization fail (LL_y residual ~1), so the
        # linearization would not be at an equilibrium. A droop R = 1e6 keeps
        # every governor state consistent and makes Pm constant to 1e-6.
        for idx in list(ss.TGOV1N.idx.v):
            ss.TGOV1N.alter("R", idx, 1e6)
    ss.setup()
    ss.PFlow.run()
    ss.TDS.init()
    ok = (
        bool(getattr(ss.TDS, "initialized", False))
        and float(np.abs(ss.dae.g).max()) < 1e-6
    )
    ss.EIG.run()
    mu = np.asarray(ss.EIG.mu)
    pf_ok = bool(ss.PFlow.converged)
    pf = {
        "bus": [int(b) for b in ss.Bus.idx.v],
        "v": [float(v) for v in ss.Bus.v.v],
        "a": [float(a) for a in ss.Bus.a.v],
        "gen_bus": [int(b) for b in list(ss.PV.bus.v) + list(ss.Slack.bus.v)],
        "p": [float(p) for p in list(ss.PV.p.v) + list(ss.Slack.p.v)],
        "q": [float(q) for q in list(ss.PV.q.v) + list(ss.Slack.q.v)],
        "pf_converged": pf_ok,
        "tds_init_ok": ok,
    }
    return mu, pf


def band(mu, lo=0.3, hi=1.5):
    f = np.abs(mu.imag) / (2 * np.pi)
    sel = mu[(f >= lo) & (f <= hi) & (mu.imag > 0)]
    return sorted(
        ((float(m.real), float(abs(m.imag) / (2 * np.pi))) for m in sel),
        key=lambda t: -t[0],
    )


def main() -> int:
    mu_g, pf = run(True)
    mu_n, _ = run(False)
    out = {
        "andes_power_flow": pf,
        "eig_with_tgov1n": {
            "rightmost_band": band(mu_g)[:6],
            "n_positive_real": int((mu_g.real > 1e-6).sum()),
            "near_zero": [complex(m).__repr__() for m in mu_g[np.abs(mu_g) < 1e-3]],
        },
        "eig_tgov1n_off": {
            "rightmost_band": band(mu_n)[:6],
            "n_positive_real": int((mu_n.real > 1e-6).sum()),
            "near_zero": [complex(m).__repr__() for m in mu_n[np.abs(mu_n) < 1e-3]],
        },
    }
    (OUT / "FC14_andes.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(
        json.dumps({k: v for k, v in out.items() if k != "andes_power_flow"}, indent=1)
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
