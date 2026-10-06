"""Freeze the protocol and create an isolated positive-DC-energy model."""
from pathlib import Path
import hashlib, json, subprocess, datetime

OUT = Path(__file__).resolve().parent
ROOT = OUT.parent.parent
if (OUT/'PREREGISTRATION.json').exists():
    raise RuntimeError('Protocol exists: use the saved model and run_validation.jl; do not regenerate it')
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
source = ROOT / 'src/pd39/model.jl'
reference = ROOT / 'src/bnd_model_expN/PDReferenceN.jl'
candidate = ROOT / 'reports/experiment_Q2B/CERTIFIED_SEARCH/Z_LOCAL_SECURE_FINAL.toml'
assert sha(candidate) == '66513a8d3a1b4cc37c9d6f2b1a5824371fed2d0c5621d2ccdb7540490fe59124'
text = source.read_text(encoding='utf-8')
component = text[text.index('@component function WeightedSimpleGFLDC'):text.index('"Return copies of the official CSV tables')]
for old, new in [
 ('WeightedSimpleGFLDC', 'PhysicalGFLDC'),
 ('(V_dc-v_dc_state)*kp_v_dc', '(v_dc_state-V_dc)*kp_v_dc'),
 ('(p_ac-P_dc)/v_dc_state', '(P_dc-p_ac)/v_dc_state'),
 ('(V_dc-v_dc_state)*ki_v_dc', '(v_dc_state-V_dc)*ki_v_dc')]:
    assert component.count(old) == 1
    component = component.replace(old, new)
ref = reference.read_text(encoding='utf-8')
builder = ref[ref.index('function _gfl'):ref.index('"""Closed-form trim')]
builder = builder.replace('PD39.PD39Model.WeightedSimpleGFLDC', 'PhysicalGFLDC').replace('base::FrozenBaseline','base::P.FrozenBaseline').replace('mdl=PD39.PD39Model','mdl=P.PD39.PD39Model')
header = '''# Generated only in this new experiment; frozen sources remain unchanged.
module PhysicalBridge
using Graphs, LinearAlgebra, NetworkDynamics, PowerDynamics, PowerDynamics.Library
using ModelingToolkitBase
using ModelingToolkitBase: @component, @variables, @parameters, t_nounits, D_nounits, System
const P = Main.PDReferenceN
const ScaledPortLFilter = P.PD39.PD39Model.ScaledPortLFilter
const NOMINAL_KP = P.NOMINAL_KP
const NOMINAL_KI = P.NOMINAL_KI
'''
(OUT/'PhysicalBridge.jl').write_text(header+component+builder+'\nend\n',encoding='utf-8')
inputs = [source, reference, candidate, ROOT/'Project.toml', ROOT/'Manifest.toml', ROOT/'experiments/bnd_expQ2B/certified_search/ConsistentEventObservation.jl', ROOT/'reports/experiment_D/inputs/machine.csv', ROOT/'reports/experiment_D/inputs/gov.csv']
protocol = {
 'registered_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
 'question':'Does the positive-DC-energy variant preserve the frozen trim and linear AC dynamics, and close its physical energy balance in selected nonlinear IEEE39 events?',
 'status':'PREREGISTERED_BEFORE_NEW_SIMULATIONS',
 'candidate':str(candidate.relative_to(ROOT)), 'candidate_sha256':sha(candidate),
 'tuning':'NONE: fixed rho, Kp, Ki, ratings, DC and current gains, PQ sharing',
 'pure_PLL_delay_s':0, 'measurement_window_s':0.5,
 'events':[{'name':'bus16_1MW','bus':16,'MW':1.0},{'name':'bus16_100MW','bus':16,'MW':100.0},{'name':'bus8_100MW','bus':8,'MW':100.0},{'name':'bus29_100MW','bus':29,'MW':100.0}],
 'new_runs':'all four physical cases plus legacy bus16_1MW independent reproduction; compare physical to archived legacy traces for all four',
 'event_contract':'Sustained ZIP Pset decrement MW/100 at t=1s, Qset unchanged; realized consumption remains voltage dependent',
 'horizon_s':61.0,'saveat_s':0.005,'solver':'Rodas5P','abstol':1e-10,'reltol':1e-10,'maxiters':2000000,'wall_limit_each_s':600,
 'gates':{'equilibrium_max_residual':1e-10,'PQ_max_pu':1e-10,'Jacobian_similarity_relative':1e-10,'finite_pole_matching_absolute':1e-5,'pointwise_energy_balance_MW':1e-6,'frequency_limit_Hz':0.5,'RoCoF_limit_Hz_s':0.5},
 'energy_verification':'State-gradient times model RHS independently compared with terminal powers, measured DC source, load and network/filter losses. Integrated trapezoid errors at 5/10/20ms reported, not certified.',
 'nonclaims':['No new optimization','No pure-delay validation','No current-limiter guarantee','No asymptotic endpoint assumed from 61s','No automatic transfer of the lossless constant-voltage theorem to the full IEEE39 model','No global replacement bound'],
 'git_HEAD':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
 'git_dirty':subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True).splitlines(),
 'inputs':{str(p.relative_to(ROOT)):sha(p) for p in inputs},
 'generated_model_sha256':sha(OUT/'PhysicalBridge.jl')
}
path=OUT/'PREREGISTRATION.json'
path.write_text(json.dumps(protocol,indent=2),encoding='utf-8')
(OUT/'PREREGISTRATION.sha256').write_text(sha(path)+'\n',encoding='ascii')
print('Frozen',path,'sha256',sha(path))
