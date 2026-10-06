import sys,json,hashlib
from pathlib import Path
source=Path(sys.argv[1]);target=Path(sys.argv[2]);d=json.loads(source.read_text())
lines=['dc_convention = "physical_supply"',
       'status = "CANDIDATE_REQUIRING_INDEPENDENT_VALIDATION"',
       f'source_json_sha256 = "{hashlib.sha256(source.read_bytes()).hexdigest()}"',
       f'retained_SG_MW = {d["retained_MW"]!r}',f'replacement_percent = {d["replacement_percent"]!r}']
for key in ('rho','Kp','Ki'):lines.append(key+' = '+repr(d[key]))
target.write_text('\n'.join(lines)+'\n')
print(target)
