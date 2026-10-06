"""Replace a seed-only spectral estimate after full contour/root coverage.

The seed estimate is copied to SPECTRAL_INITIAL.csv; correction is traceable.
"""

from __future__ import annotations

import csv
import shutil
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent


def read(path: Path) -> list[dict[str,str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def main() -> None:
    assert len(sys.argv)==2,"usage: correct_spectral_from_coverage.py candidate_id"
    cid=sys.argv[1]
    folder=HERE/"evaluations"/cid
    coverage=read(folder/"ROOT_COVERAGE.csv")[0]
    assert coverage["status"]=="COMPLETE_NUMERICAL_ROOT_COVERAGE_AT_UNSAFE_CONTOUR"
    roots=read(folder/"ROOTS.csv")
    first=min(roots,key=lambda r:float(r["local_crossing_ms"]))
    spectral=folder/"SPECTRAL.csv"
    original=folder/"SPECTRAL_INITIAL.csv"
    if not original.exists():
        shutil.copy2(spectral,original)
    old=read(original)[0]
    new=dict(old)
    new["local_tau_crit_ms"]=first["local_crossing_ms"]
    new["critical_frequency_hz"]=first["frequency_hz"]
    new["critical_root_real"]=first["critical_real"]
    new["critical_root_imag"]=first["critical_imag"]
    later=sorted((float(r["local_crossing_ms"]) for r in roots if r is not first))
    new["second_local_crossing_ms"]=later[0] if later else "nan"
    new["third_local_crossing_ms"]=later[1] if len(later)>1 else "nan"
    new["fast_roots_discovered"]=len(roots)
    new["status"]="FULL_CONTOUR_ROOT_COVERAGE_CORRECTED_INITIAL_SEED"
    with spectral.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(new));w.writeheader();w.writerow(new)
    print("CORRECTED",cid,"seed",old["local_tau_crit_ms"],"complete",new["local_tau_crit_ms"])


if __name__=="__main__":
    main()
