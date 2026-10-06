"""Preserve the first-boundary catalogue and add later modal crossings."""
from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

HERE=Path(__file__).resolve().parent

def rows(path:Path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))

def main():
    cid=sys.argv[1]
    out=HERE/"evaluations"/cid
    c=rows(out/"ROOT_COVERAGE_HIGH.csv")[0]
    assert c["status"]=="COMPLETE_NUMERICAL_ROOT_COVERAGE_AT_UNSAFE_CONTOUR"
    rootrows=[]
    for source in (out/"DISCOVERED_ROOTS.csv",out/"DISCOVERED_ROOTS_HIGH.csv"):
        for r in rows(source):
            tau=float(r["local_crossing_ms"])
            freq=float(r["frequency_at_crossing_hz"])
            if not math.isfinite(tau):continue
            if any(abs(tau-x[0])<1e-5 and abs(freq-x[1])<1e-3 for x in rootrows):continue
            rootrows.append((tau,freq,float(r["crossing_root_residual"]),source.name))
    rootrows.sort()
    formatted=[dict(candidate_id=cid,rank=j,seed_id=x[3],local_crossing_ms=x[0],
                    critical_real=-0.05,critical_imag=2*math.pi*x[1],
                    frequency_hz=x[1],small_factor_residual=x[2],
                    status="FULL_CHARACTERISTIC_HIGH_CONTOUR_COVERAGE")
               for j,x in enumerate(rootrows,1)]
    assert len(formatted)>=int(c["distinct_positive_imag_roots"])
    with (out/"ROOTS_EXPANDED.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(formatted[0]));w.writeheader();w.writerows(formatted)
    print(cid,"expanded roots",len(formatted),"first",formatted[0]["local_crossing_ms"])

if __name__=="__main__":main()
