"""F7 check - do the hypergraph labels depend on the frozen regulator leak?

The leak w = 0.05 rad/s is what makes the reactive-policy coordinate regular at
g = 0 (Safeguard A). It is a modelling constant, so the labels it produces must
not hinge on it. A stratified sample of already-evaluated F7 points (half drawn
from points adjacent to a boundary, half from region interiors, g > 0 only since
g = 0 is leak-independent by construction) is re-labelled at w = 0.02 and 0.20.
"""

from __future__ import annotations

import json
from multiprocessing import Pool

import numpy as np
import pandas as pd

import _f7_common
from _bootstrap import RESULTS  # noqa: I001  (sets the import path)
from _f7_common import Theta
from _overnight import choose_workers, pin_blas_threads
from F7_report import points_file

F7 = RESULTS / "F7"
LEAKS = (0.02, 0.20)
PER_MAP = 120
SEED = 20260910
BAD = {"BASE_UNSTABLE", "INFEASIBLE"}


def _init(leak):
    pin_blas_threads()
    _f7_common.LEAK = leak


def _label(theta):
    row = _f7_common.evaluate(theta, fast=False)
    return row.get("label", row.get("status"))


def sample(name: str, rng) -> pd.DataFrame:
    points = pd.read_csv(points_file(name), low_memory=False)
    events = pd.read_csv(F7 / f"{name}_events.csv", low_memory=False)
    ok = points[~points.label.isin(BAD) & (points.g > 0)]
    near_keys = set(zip(events.I_a, events.J_a, strict=True)) | set(
        zip(events.I_b, events.J_b, strict=True)
    )
    near = ok[[k in near_keys for k in zip(ok.I, ok.J, strict=True)]]
    far = ok.drop(near.index)
    take = PER_MAP // 2
    picked = pd.concat(
        [
            near.sample(
                min(take, len(near)), random_state=rng.integers(1 << 31)
            ).assign(stratum="boundary"),
            far.sample(min(take, len(far)), random_state=rng.integers(1 << 31)).assign(
                stratum="interior"
            ),
        ]
    )
    picked["map"] = name
    return picked


def main() -> int:
    rng = np.random.default_rng(SEED)
    frames = [
        sample(n, rng)
        for n in ("F7A", "F7B", "F7C")
        if (F7 / f"{n}_events.csv").exists()
    ]
    picked = pd.concat(frames, ignore_index=True)
    thetas = [Theta(r.g, r.k, r.t, r.h) for r in picked.itertuples()]
    for leak in LEAKS:
        with Pool(choose_workers(), initializer=_init, initargs=(leak,)) as pool:
            picked[f"label_w{leak:g}"] = pool.map(_label, thetas, chunksize=2)
    report = {
        "leaks": LEAKS,
        "reference_leak": _f7_common.LEAK,
        "sampled": int(len(picked)),
    }
    for leak in LEAKS:
        agree = picked.label == picked[f"label_w{leak:g}"]
        report[f"w={leak:g}"] = {
            "agreement": float(agree.mean()),
            "by_stratum": picked.assign(a=agree).groupby("stratum").a.mean().to_dict(),
            "by_map": picked.assign(a=agree).groupby("map").a.mean().to_dict(),
            "disagreements": picked.loc[
                ~agree, ["map", "g", "k", "t", "h", "label", f"label_w{leak:g}"]
            ].to_dict("records")[:40],
        }
    picked.to_csv(F7 / "F7_leak_sensitivity.csv", index=False)
    (F7 / "F7_leak_sensitivity.json").write_text(
        json.dumps(report, indent=2, default=str), encoding="utf-8"
    )
    print(
        json.dumps(
            {
                k: (
                    v
                    if not isinstance(v, dict)
                    else {a: b for a, b in v.items() if a != "disagreements"}
                )
                for k, v in report.items()
            },
            indent=1,
            default=str,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
