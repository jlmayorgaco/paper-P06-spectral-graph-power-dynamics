"""Verify CRC, all input hashes, then reproduce model export from extracted zip."""
import json, zipfile, hashlib, subprocess
from pathlib import Path
OUT=Path(__file__).resolve().parent
bundle=OUT/'physical_collective_damping_20261003_bundle.zip'
root=OUT/'repro'
with zipfile.ZipFile(bundle) as z:
    assert z.testzip() is None
    for name in z.namelist():
        assert not Path(name).is_absolute() and '..' not in Path(name).parts
    z.extractall(root)
manifest=json.loads((OUT/'EXPERIMENT_MANIFEST.json').read_text())
for rel,sha in manifest['source_hashes'].items():
    assert hashlib.sha256((root/rel).read_bytes()).hexdigest()==sha,rel
exp='experiments/physical_collective_damping_20261003'
result=subprocess.run(['julia','--project=.',exp+'/export_model.jl'],cwd=root,text=True,capture_output=True)
(OUT/'BUNDLE_SMOKE_TEST.log').write_text(result.stdout+'\n'+result.stderr,encoding='utf-8')
assert result.returncode==0,(result.returncode,result.stderr)
(OUT/'BUNDLE_VERIFICATION.json').write_text(json.dumps(dict(crc='PASS',source_hashes='PASS',
    independent_extracted_model_export='PASS',tested_bundle_sha256=hashlib.sha256(bundle.read_bytes()).hexdigest(),
    scope='The final package may additionally include this verification log; no numerical outputs or source inputs change.'),indent=2))
print('PASS: CRC, source hashes, independent Julia model export from extracted package',flush=True)
