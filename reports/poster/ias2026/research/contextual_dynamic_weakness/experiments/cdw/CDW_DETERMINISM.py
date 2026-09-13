# ruff: noqa: E501
"""Phase 28: rerun headline deterministic experiments into separate stores and compare outputs.

Prereg: E1 census at D01, D11, H01; E4 at D01/H4 (SPR links); E9 T1 (all families/methods).
Comparison excludes wall-clock fields.
"""

from __future__ import annotations

import json

import _infra as I
import E01_census as E1
import E09_design as E9
import E34_sens as E34


def strip(o):
    if isinstance(o, dict):
        return {k: strip(v) for k, v in o.items() if k not in ("wall_s", "wall_design_s", "t_full", "t_eig_full", "t_eig", "t_contour")}
    if isinstance(o, list):
        return [strip(v) for v in o]
    return o


def compare(phase, tasks):
    orig = I.Store(phase)
    rep = I.Store("DET_" + phase)
    same = diff = missing = 0
    examples = []
    for t in tasks:
        k = I.task_key(t)
        if not orig.done(k) or not rep.done(k):
            missing += 1
            continue
        a, b = strip(orig.get(k)), strip(rep.get(k))
        if json.dumps(a, sort_keys=True) == json.dumps(b, sort_keys=True):
            same += 1
        else:
            diff += 1
            if len(examples) < 3:
                examples.append(k)
    return {"phase": phase, "n": len(tasks), "identical": same, "different": diff, "missing": missing, "examples": examples}


def main(workers=None):
    res = I.resources(workers or 12)
    t1 = [t for t in E1.tasks() if t["pid"] in ("D01", "D11", "H01")]
    t34 = [t for t in E34.tasks() if t["pid"] == "D01" and t["target"] == "H4" and t["family"] == "link" and t["semantics"] == "SPR" and t["env"] is None]
    t9 = [t for t in E9.tasks() if t["target"] == "T1"]
    I.run_tasks("DET_E01", "E01_census", "run_task", t1, res["workers"])
    I.run_tasks("DET_E34", "E34_sens", "run_task", t34, res["workers"])
    I.run_tasks("DET_E09", "E09_design", "run_task", t9, res["workers"])
    out = [compare("E01", t1), compare("E34", t34), compare("E09", t9)]
    summary = {"environment": I.environment(), "workers": res["workers"], "comparisons": out,
               "all_identical": all(c["different"] == 0 and c["missing"] == 0 for c in out)}
    I.atomic_write_json(I.RESULTS / "CDW_DETERMINISM.json", summary)
    print(json.dumps(out, indent=1))
    return summary


if __name__ == "__main__":
    main()
