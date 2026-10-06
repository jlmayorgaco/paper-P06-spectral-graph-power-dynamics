"""Copy immutable inputs and archive prior evidence without changing originals."""
from pathlib import Path
import hashlib,json,shutil,zipfile
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[1]
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
rels=[p.relative_to(ROOT) for p in (ROOT/'experiments/graph_gsp_codesign_20261003/model').glob('*.csv')]
rels += [Path(x) for x in ['experiments/graph_gsp_codesign_20261003/baseline.toml',
 'experiments/theory_collective_damping_20261003/regional_certificate_20261003/FROZEN_ALGEBRA_POINT.json',
 'experiments/all_pll_gain_map_validation_20261003/SUMMARY.json']]
records={}
for rel in rels:
    source=ROOT/rel;dest=OUT/'inputs'/rel;dest.parent.mkdir(parents=True,exist_ok=True)
    if dest.exists(): assert sha(dest)==sha(source)
    else: shutil.copyfile(source,dest)
    records[str(rel).replace('\\','/')]=sha(source)
evidence=OUT/'evidence';evidence.mkdir(exist_ok=True)
for dirname in ['all_pll_gain_map_validation_20261003','theory_collective_damping_20261003/regional_certificate_20261003']:
    folder=ROOT/'experiments'/dirname
    name=dirname.split('/')[-1]+'_frozen.zip';target=evidence/name
    if not target.exists():
        with zipfile.ZipFile(target,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
            for p in sorted(folder.rglob('*')):
                if p.is_file() and '__pycache__' not in p.parts:
                    z.write(p,p.relative_to(ROOT).as_posix())
    records[str(target.relative_to(OUT))]=sha(target)
(OUT/'BUNDLED_INPUTS.json').write_text(json.dumps(records,indent=2)+'\n')
print('Bundled',len(rels),'direct inputs and two frozen evidence archives.')

# Transitive device equations, network input data and Julia environment.
# These are archived copies, never edits to the frozen sources.
snapshot=ROOT/'experiments/physical_collective_damping_20261003/source_snapshot'
entries={p.relative_to(snapshot).as_posix():p for p in snapshot.rglob('*') if p.is_file()}
for name in ['Project.toml','Manifest.toml']:
    entries['experiments/physical_collective_damping_20261003/source_snapshot/'+name]=snapshot/name
prior=json.loads((ROOT/'experiments/all_pll_gain_map_validation_20261003/INPUT_MANIFEST.json').read_text())
for rel,expected in prior['inputs'].items():
    rel=rel.replace('\\','/')
    source=ROOT/rel
    if rel.endswith('theory_collective_damping_20261003/THEORY.tex'):
        source=ROOT/'experiments/all_pll_gain_map_validation_20261003/source_snapshot/THEORY_before_validation.tex'
    assert sha(source)==expected,(rel,sha(source),expected)
    entries[rel]=source
target=evidence/'HISTORICAL_REPRO_SOURCE.zip'
if not target.exists():
    with zipfile.ZipFile(target,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for rel,p in sorted(entries.items()): z.write(p,rel)
source_manifest={rel:sha(p) for rel,p in entries.items()}
(OUT/'HISTORICAL_SOURCE_MANIFEST.json').write_text(json.dumps(source_manifest,indent=2)+'\n')
records[str(target.relative_to(OUT))]=sha(target)
(OUT/'BUNDLED_INPUTS.json').write_text(json.dumps(records,indent=2)+'\n')
print('Archived',len(entries),'historical source/environment inputs.')
