"""Restore frozen dependencies on a fresh clone; never overwrite existing files."""
from pathlib import Path
import hashlib,json,zipfile
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[1]
REG=ROOT/'experiments/regional_paper_closure_20261003'
sha=lambda b:hashlib.sha256(b).hexdigest()
restored=[];identical=[];skipped=[]
def put(rel,data):
    target=(ROOT/rel).resolve()
    if not target.is_relative_to(ROOT.resolve()):raise ValueError('Unsafe archive member')
    if target.exists():
        if target.read_bytes()!=data:raise ValueError('Existing input differs; preserved: '+rel)
        identical.append(rel);return
    target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(data);restored.append(rel)
for archive in sorted((REG/'evidence').glob('*.zip')):
    with zipfile.ZipFile(archive) as z:
        for n in z.namelist():
            if n.endswith('/'):continue
            if n in ['Project.toml','Manifest.toml','experiments/theory_collective_damping_20261003/THEORY.tex']:
                skipped.append(n);continue
            put(n,z.read(n))
for p in (REG/'inputs').rglob('*'):
    if p.is_file():put(p.relative_to(REG/'inputs').as_posix(),p.read_bytes())
print(json.dumps({'restored':len(restored),'identical':len(identical),
                  'skipped_live_manuscript_and_root_environment':skipped},indent=2))
