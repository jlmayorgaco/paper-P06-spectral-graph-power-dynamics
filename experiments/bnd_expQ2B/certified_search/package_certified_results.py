"""New archive; preserves the original Q2B bundle and its historical status."""
from pathlib import Path
import hashlib,zipfile,json
ROOT=Path(__file__).resolve().parents[3]
archive=ROOT/'reports/experiment_Q2B_certified_local_bundle.zip'
roots=['reports/experiment_Q2B','experiments/bnd_expQ2B','src/bnd_expQ2B','test/bnd_expQ2B',
       'src/bnd_model_expN','src/bnd_design','src/bnd_design_e','src/bnd_expQ',
       'src/pd39','reports/experiment_D/inputs']
files=set()
for name in roots:
    files.update(p for p in (ROOT/name).rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix not in ('.pyc','.pyo'))
mf=json.loads((ROOT/'reports/experiment_N/MODEL_FREEZE.json').read_text())
files.update(ROOT/p for p in mf['inputs'])
for name in ['Project.toml','Manifest.toml','reports/experiment_N/MODEL_FREEZE.json',
             'reports/experiment_N/TABLE_N01_original_operating_point.csv.sha256',
             'reports/experiment_A/matrices/bus33_baseline_equilibrium.csv']:
    files.add(ROOT/name)
entries=sorted(files,key=lambda p:p.relative_to(ROOT).as_posix())
manifest={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in entries}
readme='''# Cierre local Q2B

Leer primero reports/experiment_Q2B/CERTIFIED_SEARCH/RESULTADO_LOCAL_Y_BRECHA_GLOBAL.md.
Los informes anteriores conservan sus estados históricos FAIL_OPTIMIZATION; no se reescriben.
Este paquete añade la búsqueda corregida, certificado local numérico, robustez de banda completa,
validación PD independiente, figuras, CSV, sensores, código e insumos congelados.
La globalidad sigue abierta. No hubo commit/push. Dependencias Julia externas según Manifest.toml.
'''
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p in entries:z.write(p,p.relative_to(ROOT).as_posix())
    z.writestr('START_HERE.md',readme)
    z.writestr('SHA256_MANIFEST.json',json.dumps(manifest,indent=2))
with zipfile.ZipFile(archive) as z:
    assert z.testzip() is None
    for name,h in manifest.items():assert hashlib.sha256(z.read(name)).hexdigest()==h
h=hashlib.sha256(archive.read_bytes()).hexdigest()
archive.with_suffix('.zip.sha256').write_text(h+'  '+archive.name+'\n')
print(json.dumps(dict(path=str(archive),files=len(entries),bytes=archive.stat().st_size,sha256=h,CRC_and_manifest_pass=True),indent=2))
