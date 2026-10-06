# ruff: noqa: E501
"""Prereg §3 determinism: rerun the H1 census at HARDENING_H01 and the H6 link condition at
HARDENING_H01/H4 into separate stores (raw/H_DET_*) and compare with the originals on every
summary field except wall-clock. Writes results/hardening/CDWH_DETERMINISM.json."""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
for _p in (str(HERE), str(HERE.parent / "cdw")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import _infra as I  # noqa: E402
import _hinfra as HI  # noqa: E402
import H01_census as H1  # noqa: E402
import H06_links as H6  # noqa: E402


def strip(o):
    if isinstance(o, dict):
        return {k: strip(v) for k, v in o.items() if k not in ("wall_s",)}
    if isinstance(o, list):
        return [strip(v) for v in o]
    return o


def main():
    t1 = [t for t in H1.tasks() if t["pid"] == "HARDENING_H01"]
    t6 = [t for t in H6.tasks() if t["pid"] == "HARDENING_H01" and t["target"] == "H4" and t["family"] == "link"]
    I.run_tasks("H_DET_H01", "H01_census", "run_task", t1, 4)
    I.run_tasks("H_DET_H06", "H06_links", "run_task", t6, 4)
    out = {}
    for store_o, store_d, ts in (("H_H01", "H_DET_H01", t1), ("H_H06", "H_DET_H06", t6)):
        so, sd = I.Store(store_o), I.Store(store_d)
        n = same = 0
        diffs = []
        for t in ts:
            k = I.task_key(t)
            a, b = strip(so.get(k)), strip(sd.get(k))
            n += 1
            if a == b:
                same += 1
            else:
                diffs.append(k)
        out[store_o] = {"n_tasks": n, "identical": same, "differing_keys": diffs[:10]}
    out["all_identical"] = all(v["identical"] == v["n_tasks"] for k, v in out.items() if isinstance(v, dict))
    HI.write_json("CDWH_DETERMINISM.json", out)
    print(out)


if __name__ == "__main__":
    main()
