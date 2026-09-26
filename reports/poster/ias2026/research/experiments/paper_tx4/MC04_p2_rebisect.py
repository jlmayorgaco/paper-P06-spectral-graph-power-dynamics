"""MC04: POST-HOC re-localization of the P2 events (not a preregistered pass).

The preregistered P2 check (MC02) failed as implemented: its bisection used the
three-state label (STABLE / UNSTABLE / BOUNDARY_OR_UNRESOLVED), so it stopped at
the edge of the classifier's unresolved band, where the crossing eigenvalue
still has abs(Re) of 4e-6 to 7e-5 > 1e-6. The frozen threshold is unchanged here;
only the bisection variable changes: the sign of alpha_perp (the continuous
rightmost real part of the transverse spectrum, no classifier), 50 steps, on
the same segments, grid brackets and subsets recorded by MC02. The result is
reported beside the preregistered failure, never instead of it.
"""

from __future__ import annotations

import json
import sys
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd
from _mc import WORKERS, MCExperiment, out_dir, write_json
from MC02_ieee39 import _init, solve  # noqa: E402

from ibr_cycles.certification.symmetry import (  # noqa: E402
    frequency_partner,
    rotation_generator,
)
from ibr_cycles.certification.transverse import transverse_operator  # noqa: E402

OUT = out_dir("MC04_p2_rebisect")
SRC = out_dir("MC02_ieee39")


def spectrum(members, v):
    case = solve(members, v)
    r_x, _ = rotation_generator(case.dae, case.equilibrium.z)
    tr = transverse_operator(case.system.A, r_x, frequency_partner(case.dae).w)
    return np.linalg.eigvals(tr.a_perp)


def task(item):
    seg, ev = item
    va, vb = np.array(seg["theta_a"]), np.array(seg["theta_b"])
    npts = len(seg["H_sequence"])
    grid = np.linspace(0.0, 1.0, npts)
    members = (
        () if ev["subset"] == "BASE" else tuple(int(b) for b in ev["subset"].split("+"))
    )
    lo, hi = grid[ev["j"]], grid[ev["j"] + 1]
    a_lo = spectrum(members, va + lo * (vb - va)).real.max()
    a_hi = spectrum(members, va + hi * (vb - va)).real.max()
    if np.sign(a_lo) == np.sign(a_hi):
        return {**ev, "sign_change": False, "re": np.nan, "im": np.nan}
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        a = spectrum(members, va + mid * (vb - va)).real.max()
        if np.sign(a) == np.sign(a_lo):
            lo, a_lo = mid, a
        else:
            hi = mid
    lam = spectrum(members, va + 0.5 * (lo + hi) * (vb - va))
    k = int(np.argmax(lam.real))
    return {
        **ev,
        "sign_change": True,
        "re": float(lam[k].real),
        "im": float(abs(lam[k].imag)),
        "abs_lambda": float(abs(lam[k])),
    }


def main(argv) -> int:
    exp = MCExperiment(
        name="MC04_p2_rebisect",
        question="Post hoc: are the P2 hypergraph changes bracketed by axis crossings?",
        config={"post_hoc": True, "bisection": "sign of alpha_perp, 50 steps"},
        workers=WORKERS,
    )
    started = time.time()
    segs = json.loads((SRC / "P2_segments.json").read_text(encoding="utf-8"))
    items = [(s, e) for s in segs for e in s["events"] if e.get("subset")]
    empty = [e for s in segs for e in s["events"] if not e.get("subset")]
    with Pool(WORKERS, initializer=_init) as pool:
        rows = pool.map(task, items, chunksize=1)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "MC04_events.csv", index=False)
    ok = df.sign_change & (df.re.abs() <= 1e-6)
    res = {
        "events_with_subset": len(df),
        "events_without_changed_subset": len(empty),
        "sign_change_confirmed": int(df.sign_change.sum()),
        "located_abs_re_le_1e-6": int(ok.sum()),
        "abs_re_max": float(df.re.abs().max()),
        "freq_hz_range": [
            float(df.im.min() / (2 * np.pi)),
            float(df.im.max() / (2 * np.pi)),
        ],
        "real_crossings": int((df.im < 1e-6).sum()),
        "elapsed_s": round(time.time() - started, 1),
    }
    write_json(OUT / "MC04_summary.json", res)
    exp.finish("COMPUTED", **res)
    print(json.dumps(res, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
