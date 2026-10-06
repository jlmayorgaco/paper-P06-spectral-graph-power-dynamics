"""Record environment + sha256 of every frozen model input. Read-only on history."""
import hashlib, json, platform, subprocess, pathlib, sys
import numpy, scipy, pandas
ROOT = pathlib.Path(__file__).resolve().parents[3]
CAMP = pathlib.Path(__file__).resolve().parents[1]
SRC = ROOT / 'experiments' / 'graph_gsp_codesign_20261003'
files = sorted((SRC / 'model').glob('*.csv')) + [SRC / 'baseline.toml', SRC / 'TABLE_02_BASELINE_ROOTS.csv', SRC / 'DELAY_MODEL_CONTRACT.md',
         ROOT / 'experiments/nonlinear_codesign_20261001/ReducedDAE.jl', SRC / 'DelayedEvents.jl', ROOT / 'Project.toml', ROOT / 'Manifest.toml',
         ROOT / 'experiments/interaction_decision_20261004/model.py']
def sh(*a):
    r = subprocess.run(a, cwd=ROOT, text=True, capture_output=True); return r.stdout.strip()
rec = dict(head=sh('git', 'rev-parse', 'HEAD'), branch=sh('git', 'branch', '--show-current'),
           dirty_tracked=sh('git', 'status', '--porcelain', '--untracked-files=no').splitlines(),
           julia=sh('julia', '--version'), python=platform.python_version(), numpy=numpy.__version__, scipy=scipy.__version__, pandas=pandas.__version__,
           inputs={p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in files})
(CAMP / 'derived' / 'ENV_RECORD.json').write_text(json.dumps(rec, indent=2))
print(json.dumps({k: v for k, v in rec.items() if k != 'inputs'}, indent=1)); print(len(rec['inputs']), 'inputs hashed')
