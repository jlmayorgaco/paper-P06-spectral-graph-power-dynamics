"""Portable source/evidence archive, with exact-byte inventory and no historical mutation."""
from moments import *
import zipfile

def main():
    archive=OUT/'NETWORK_MOMENT_CODESIGN_20261004_COMPLETE.zip'
    skip={archive.name,archive.name+'.sha256','BUNDLE_MANIFEST.json'}
    def allowed(p):return p.is_file() and p.name not in skip and p.suffix not in {'.aux','.out','.pyc'} and not any(x in p.parts for x in ['__pycache__','qa'])
    paths={p for p in OUT.rglob('*') if allowed(p)}
    sibling=OUT.parent/'graph_gain_geometry_20261004'
    paths.update(p for p in sibling.rglob('*') if allowed(p))
    # Original dependency archive entries are exact repository paths.
    dep=json.loads((OLD/'DEPENDENCY_MANIFEST.json').read_text())['files']
    paths.update(ROOT/name for name in dep)
    extra=['model.py','designs/corrected.toml','TABLE_08_MULTIBRANCH_CONTINUATION.csv','TABLE_12_NONLINEAR_EVENTS.csv',
        'FIG_01_CERTIFIED_DECISION.pdf','FIG_02_COMPENSATION_RULE.pdf','POSTER_CLAIMS.md','REPORT_ES.md']
    paths.update(OLD/name for name in extra)
    paths.add(OLD/'INTERACTION_DECISION_20261004_COMPLETE.zip')
    paths.update(OLD/f'CERT_{name}.json' for name in ['anchor','single30','single37','joint','complex_pair','corrected'])
    causal=tomllib.loads((sibling/'causal_audit/MANIFEST.toml').read_text())
    paths.update(ROOT/name for name in causal['source_sha256'])
    paths.add(ROOT/'experiments/theory_collective_damping_20261003/THEORY.tex')
    manifest={'files':{str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
        'environment':'Python packages and Julia package artifacts installed separately from pinned manifests',
        'scope':'new moment campaign plus prior negative geometry campaign and required read-only inputs'}
    dump('BUNDLE_MANIFEST.json',manifest);paths.add(OUT/'BUNDLE_MANIFEST.json')
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
        for p in sorted(paths):z.write(p,str(p.relative_to(ROOT)).replace('\\','/'))
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None
        for name,digest in manifest['files'].items():assert hashlib.sha256(z.read(name)).hexdigest()==digest,name
    (OUT/(archive.name+'.sha256')).write_text(hashlib.sha256(archive.read_bytes()).hexdigest()+'  '+archive.name+'\n')
    print('Verified',len(paths),'files; archive bytes',archive.stat().st_size)

if __name__=='__main__':main()
