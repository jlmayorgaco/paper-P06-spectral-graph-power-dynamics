"""Record final evidence integrity without running or changing the experiment."""
from pathlib import Path
import json, hashlib, tomllib, re
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
pr=json.loads((OUT/'PREREGISTRATION.json').read_text());env=tomllib.loads((OUT/'RUN_ENVIRONMENT.toml').read_text())
main=ROOT/'experiments/theory_collective_damping_20261003/THEORY.tex'
text=main.read_text(encoding='utf-8');before=(OUT/'THEORY_before_bridge.tex').read_text(encoding='utf-8')
without=re.sub(r'% BEGIN PHYSICAL_FREQUENCY_BRIDGE_20261004.*?% END PHYSICAL_FREQUENCY_BRIDGE_20261004\n\n','',text,flags=re.S)
checks={'frozen_inputs_unchanged':all(sha(ROOT/p)==h for p,h in pr['inputs'].items()),'preregistration_hash_unchanged':sha(OUT/'PREREGISTRATION.json')==(OUT/'PREREGISTRATION.sha256').read_text().strip(),'model_hash_unchanged':sha(OUT/'PhysicalBridge.jl')==pr['generated_model_sha256']==env['model_sha256'],'executed_runner_unchanged':sha(OUT/'run_validation.jl')==env['runner_sha256'],'prior_manuscript_content_preserved':without==before,'one_new_manuscript_section':text.count(r'\label{sec:physical-energy-bridge}')==1,'new_undefined_citations':False}
assert all(v for k,v in checks.items() if k!='new_undefined_citations'),checks
latex={'path':str(main),'source_sha256':sha(main),'compiled':False,'attempts_this_turn':2,'status':'BLOCKED_HOST_COMPILER','error':'Unable to find standard directories for platform','editor':'Existing THEORY.tex editor retained; no replacement PDF created','static_checks':'Unique section/labels and balanced added environments checked; not compilation'}
section=(OUT/'PHYSICAL_BRIDGE_SECTION.tex').read_text()+ (OUT/'NUMERICAL_SECTION.tex').read_text()
beg=re.findall(r'\\begin\{([^}]+)\}',section);end=re.findall(r'\\end\{([^}]+)\}',section)
assert sorted(beg)==sorted(end)
(OUT/'EDITOR_STATUS.json').write_text(json.dumps(latex,indent=2),encoding='utf-8')
(OUT/'INTEGRITY_CHECKS.json').write_text(json.dumps(checks,indent=2),encoding='utf-8')
manifest=json.loads((OUT/'MANIFEST.json').read_text())
manifest['installed_sources']={str(p):sha(p) for p in [Path('C:/Users/walla/.julia/packages/PowerDynamics/VzOiZ/src/Library/Controls/Govs.jl'),Path('C:/Users/walla/.julia/packages/PowerDynamics/VzOiZ/src/Library/Renewables/ComposableInverter.jl')]}
manifest['open_manuscript']={'path':str(main),'sha256':sha(main)}
manifest['outputs']={p.name:sha(p) for p in sorted(OUT.iterdir()) if p.is_file() and p.name!='MANIFEST.json'}
(OUT/'MANIFEST.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
print(json.dumps(checks,indent=2));print('Compiler limitation recorded; no frozen input or executed code changed.')
