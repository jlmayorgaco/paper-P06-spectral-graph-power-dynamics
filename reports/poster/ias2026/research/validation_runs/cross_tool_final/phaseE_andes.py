"""Phase E1, ANDES side: d lambda / d gamma_e on the frozen holdout, same cases.

Runs under the tx3-andes interpreter and imports NONE of the project model code. The
equation-equivalent R3 case is built by the unchanged F1 worker functions
(experiments/F1_andes_equivalent_worker.py: build_case, run_case): SEXS with
TA/TB = 1 as the first-order AVR, GENROU degraded to fourth order (tau 1e-4),
IEEEST/TGOV1N/Toggler/ACEc removed, the flagship machines replaced by constant-power
PQ injections, PQ loads constant power.

Line scaling. gamma_e multiplies the complete branch two-port admittance of the
canonical line e (series and charging, tap unchanged). In the ANDES Line sheet (same
row order as the canonical JSON; verified by bus1/bus2/tap in Phase C) this is
    r -> r / gamma,  x -> x / gamma,  b, g, b1, b2, g1, g2 -> gamma * (value).
ANDES evaluates u/((r+1e-8) + 1j(x+1e-8)) in its line equations (line.py:200); the
regularization perturbs the scaled admittance by < 1e-8/|z| relative, far below the
finite-difference resolution.

Same-case semantics. Internally ("static_power", q_policy "matched") the replaced
unit takes the (P, Q) that the machine produces in the all-PV power flow of the
SCALED network. To solve exactly the same case, the injections here are re-read at
each gamma from the ANDES all-PV power flow of the scaled full case (F1 read them
once at gamma = 1; at gamma = 1 the two are identical).

Tracking. At each gamma the eigenvalue nearest the gamma = 1 critical band mode is
used; d lambda / d gamma = (lambda(1+eps) - lambda(1-eps)) / (2 eps), eps = 0.002.

Writes phaseE/andes_E1.csv and phaseE/andes_E1.json.
Usage (tx3-andes venv): python phaseE_andes.py <research_root> <holdout_lines.json>
"""

from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
from pathlib import Path

import andes
import numpy as np
import pandas as pd

repo = Path(sys.argv[1]).resolve()
holdout = json.loads(Path(sys.argv[2]).read_text())
SEL = [int(x) for x in holdout["lines"]]
assert holdout["seed"] == 20260911 and SEL == [
    3,
    9,
    15,
    17,
    22,
    26,
    31,
    32,
    35,
    39,
    42,
    44,
]
spec = importlib.util.spec_from_file_location(
    "f1worker", repo / "experiments/F1_andes_equivalent_worker.py"
)
w = importlib.util.module_from_spec(spec)
spec.loader.exec_module(w)

HERE = Path(__file__).resolve().parent
OUT = HERE / "phaseE"
OUT.mkdir(exist_ok=True)
EPS = 0.002
CORE = w.CORE
payload = json.loads((repo / "configs/ias2026/ieee39_network.json").read_text())
sheets = pd.read_excel(w.case_path(), sheet_name=None)
line_sheet = sheets["Line"]
canon = [(int(x["bus1"]), int(x["bus2"])) for x in payload["lines"]]
assert (
    list(zip(line_sheet.bus1.astype(int), line_sheet.bus2.astype(int), strict=False))
    == canon
)


def scaled_sheets(li: int | None, gamma: float) -> dict:
    out = {k: v.copy() for k, v in sheets.items()}
    if li is None:
        return out
    ln = out["Line"]
    row = ln.index[li]
    for col in ("r", "x"):
        ln.loc[row, col] = float(ln.loc[row, col]) / gamma
    for col in ("b", "g", "b1", "b2", "g1", "g2"):
        if col in ln.columns and pd.notna(ln.loc[row, col]):
            ln.loc[row, col] = float(ln.loc[row, col]) * gamma
    return out


def generation(sh: dict, workdir: str, tag: str):
    """All-PV power flow of the (scaled) full case -> machine P, Q per bus."""
    path = str(Path(workdir) / f"full_{tag}.xlsx")
    with pd.ExcelWriter(path) as writer:
        for name, frame in sh.items():
            if len(frame):
                frame.to_excel(writer, sheet_name=name, index=False)
    andes.config_logger(stream_level=40)
    ss = andes.load(path, setup=True, no_output=True)
    ss.PFlow.run()
    assert ss.PFlow.converged
    return dict(zip(ss.PV.bus.v, ss.PV.p.v, strict=False)), dict(
        zip(ss.PV.bus.v, ss.PV.q.v, strict=False)
    )


def spectrum(li, gamma, workdir, tag):
    sh = scaled_sheets(li, gamma)
    pgen, qgen = generation(sh, workdir, tag)
    path = w.build_case(sh, CORE, pgen, qgen, workdir, tag, avr=True)
    r = w.run_case(path, fourth_order=True, tau=1e-4)
    assert r["status"] == "OK", r
    return r


rows = []
with tempfile.TemporaryDirectory(prefix="phaseE_andes_") as work:
    r0 = spectrum(None, 1.0, work, "g1")
    band = w.band_modes(r0["all_eigenvalues"])
    lam0 = complex(band[0])
    info = {
        "andes_version": andes.__version__,
        "case_file": w.case_path(),
        "lambda0_real": lam0.real,
        "lambda0_imag": lam0.imag,
        "lambda0_freq_hz": lam0.imag / (2 * np.pi),
        "alpha_band_gamma1": r0["alpha_band"],
        "holdout_seed": holdout["seed"],
        "holdout_lines": SEL,
        "eps": EPS,
    }
    print(json.dumps(info, indent=2), flush=True)
    for li in SEL:
        lams = []
        for k, gam in enumerate((1 - EPS, 1 + EPS)):
            ev = spectrum(li, gam, work, f"L{li}_{k}")["all_eigenvalues"]
            lams.append(complex(ev[np.argmin(np.abs(ev - lam0))]))
        dl = (lams[1] - lams[0]) / (2 * EPS)
        rows.append(
            {
                "line_index": li,
                "line": f"L{li:02d}:{canon[li][0]}-{canon[li][1]}",
                "andes_dlambda_real": dl.real,
                "andes_dlambda_imag": dl.imag,
                "andes_lambda_minus": str(lams[0]),
                "andes_lambda_plus": str(lams[1]),
            }
        )
        print(rows[-1], flush=True)

pd.DataFrame(rows).to_csv(OUT / "andes_E1.csv", index=False)
(OUT / "andes_E1.json").write_text(json.dumps(info, indent=2))
