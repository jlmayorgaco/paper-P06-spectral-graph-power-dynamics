"""Freeze an isolated compiled-PowerDynamics reference with explicit DC supply signs."""
from pathlib import Path
import hashlib,json
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
reference=ROOT/'src/bnd_model_expN/PDReferenceN.jl'
model=ROOT/'src/pd39/model.jl'
s=reference.read_text(encoding='utf-8')
s=s.replace('module PDReferenceN','module PDPhysicalReference',1)
s=s.replace('include(joinpath(@__DIR__, "..", "pd39", "PD39.jl"))',
            'include(joinpath(@__DIR__, "..", "..", "src", "pd39", "PD39.jl"))',1)
source=model.read_text(encoding='utf-8')
start=source.index('@component function WeightedSimpleGFLDC')
end=source.index('\nend',start)+len('\nend')
component=source[start:end].replace('WeightedSimpleGFLDC','PhysicalGFLDC',1)
component=component.replace('ScaledPortLFilter(;','PD39.PD39Model.ScaledPortLFilter(;')
assert component.count('(V_dc-v_dc_state)')==2
assert component.count('(p_ac-P_dc)')==1
component=component.replace('(V_dc-v_dc_state)','(v_dc_state-V_dc)').replace('(p_ac-P_dc)','(P_dc-p_ac)')
imports='using ModelingToolkitBase\nusing PowerDynamics.Library\nusing ModelingToolkitBase: @component, @variables, @parameters, t_nounits, D_nounits, System\n'
s=s.replace('using .PD39','using .PD39\n'+imports+'\n'+component,1)
s=s.replace('PD39.PD39Model.WeightedSimpleGFLDC(','PhysicalGFLDC(',1)
target=HERE/'PDPhysicalReference.jl';target.write_text(s,encoding='utf-8')
out=ROOT/'reports/nonlinear_codesign_20261001';out.mkdir(exist_ok=True)
manifest={'source_reference_sha256':hashlib.sha256(reference.read_bytes()).hexdigest(),
          'source_model_sha256':hashlib.sha256(model.read_bytes()).hexdigest(),
          'physical_reference_sha256':hashlib.sha256(target.read_bytes()).hexdigest(),
          'changes':['DC energy balance uses positive Pdc supply minus converter-to-AC power.',
                     'Both DC PI errors use measured voltage minus reference.',
                     'No shared core source or installed library was modified.']}
(out/'physical_pd_reference_manifest.json').write_text(json.dumps(manifest,indent=2))
print(target)
