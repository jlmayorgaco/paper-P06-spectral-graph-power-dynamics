"""Append generated evidence to the existing open document, preserving its input."""
from pathlib import Path
import hashlib
import json
import shutil

OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[1]
THEORY=ROOT/"experiments/theory_collective_damping_20261003/THEORY.tex"
def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

summary=json.loads((OUT/"SUMMARY.json").read_text())
assert summary["status"]=="TARGETED_VALIDATION_PASSED"
manifest=json.loads((OUT/"INPUT_MANIFEST.json").read_text())
key=next(k for k in manifest["inputs"] if k.endswith("THEORY.tex"))
snapshot=OUT/"source_snapshot/THEORY_before_validation.tex"
snapshot.parent.mkdir(exist_ok=True)
if not snapshot.exists():
    assert sha(THEORY)==manifest["inputs"][key],"Unexpected theory modification; inspect before merging."
    shutil.copyfile(THEORY,snapshot)
assert sha(snapshot)==manifest["inputs"][key]
previous=OUT/"THEORY_UPDATE_PROVENANCE.json"
if previous.exists():
    known=json.loads(previous.read_text())
    assert sha(THEORY)==known["after_sha256"],"The editor changed the document; merge manually."
else:
    assert sha(THEORY)==sha(snapshot)
source=snapshot.read_text(encoding="utf-8")
replacements=[
    (r"\date{3 October 2026 --- theory first; no new dynamic simulations}",
     r"\date{3 October 2026 --- theory and one preregistered validation}"),
    ("No new grid experiment is executed.",
     "A subsequent preregistered IEEE--39 trial validates the all-PLL gain\n"
     "construction for one finite replacement step, including the numerical\n"
     "full spectral margin and five nonlinear delayed events. This does not\n"
     "demonstrate a replacement maximum or superiority over fixed gains."),
    ("parity, all-root checks and event validation remain unexecuted here.",
     "parity, all-root checks and event validation were the subsequent gates.\n"
     "Section~\\ref{sec:mapvalidation} reports their targeted execution."),
    ("A first experiment should test whether the formula reconstructs known",
     "The preregistered first experiment tests whether the formula reconstructs known"),
    ("& SUPPORTED\\_LOCAL formulation; no grid run or convergence evidence\\\\",
     "& NUMERICALLY\\_VALIDATED for one prescribed step; no convergence claim\\\\"),
    ("they are not IEEE--39 dynamic validation. No simulation, optimizer\n"
     "or DDE spectral sweep is run. New files stay in\n"
     "\\path{experiments/theory_collective_damping_20261003/}.",
     "they are not themselves IEEE--39 dynamic validation. The subsequent\n"
     "one-step trial in Section~\\ref{sec:mapvalidation} is stored separately in\n"
     "\\path{experiments/all_pll_gain_map_validation_20261003/}.\n"
     "It runs no replacement optimizer or delay sweep. Original theory files\n"
     "stay in \\path{experiments/theory_collective_damping_20261003/}."),
]
for before,after in replacements:
    assert source.count(before)==1,f"Expected one occurrence: {before}"
    source=source.replace(before,after)
marker=r"\section{Claim ledger, novelty and missing gates}"
assert source.count(marker)==1
source=source.replace(marker,(OUT/"GENERATED_VALIDATION_SECTION.tex").read_text(encoding="utf-8")+marker)
THEORY.write_text(source,encoding="utf-8")
previous.write_text(json.dumps({
    "before_sha256":sha(snapshot),"after_sha256":sha(THEORY),
    "before_snapshot":str(snapshot),"source":str(THEORY),
    "summary_sha256":sha(OUT/"SUMMARY.json"),
    "generated_section_sha256":sha(OUT/"GENERATED_VALIDATION_SECTION.tex"),
    "scope":"append executed validation and reconcile stale theory-only statements; no theorem changes"
},indent=2)+"\n",encoding="utf-8")
print("Updated the existing THEORY.tex in place; input snapshot and prior provenance preserved.")
