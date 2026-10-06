"""Audit saved evidence and write the final machine-readable closure ledger."""
from pathlib import Path
import sys,json,hashlib,re,zipfile,subprocess
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[1]
sys.path.insert(0,str(OUT/'vendor'))
from flint import arb,acb,ctx
ctx.prec=128
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def save(n,o): (OUT/n).write_text(json.dumps(o,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
def ball(d): return acb(arb(d['real']),arb(d['imag']))
result=json.loads((OUT/'CONTOUR_RESULT.json').read_text())
panels=json.loads((OUT/'PANELS.json').read_text())
assert len(panels)==485==result['panel_count']
assert all(panels[i]['b']==panels[(i+1)%len(panels)]['a'] for i in range(len(panels)))
assert all(p['a']!=p['b'] for p in panels)
integral=sum((ball(p['count']) for p in panels),acb(0))/(2*arb.pi()*acb(0,1))
assert integral.real>arb('.5') and integral.real<arb('1.5') and integral.contains(acb(1))
assert result['count_one_proved']
archived=0
for name,prefix in [
 ('all_pll_gain_map_validation_20261003_frozen.zip','experiments/all_pll_gain_map_validation_20261003/'),
 ('regional_certificate_20261003_frozen.zip','experiments/theory_collective_damping_20261003/regional_certificate_20261003/')]:
    with zipfile.ZipFile(OUT/'evidence'/name) as z:
        manifest=json.loads(z.read(prefix+'FINAL_MANIFEST.json'))
        for rel,expected in manifest['artifacts'].items():
            actual=hashlib.sha256(z.read(prefix+rel.replace('\\','/'))).hexdigest()
            assert actual==expected,rel
            archived+=1
bundled=json.loads((OUT/'BUNDLED_INPUTS.json').read_text())
for rel,expected in bundled.items():
    rel=rel.replace('\\','/')
    p=OUT/rel if rel.startswith('evidence/') else OUT/'inputs'/rel
    assert sha(p)==expected,rel
source=ROOT/'experiments/theory_collective_damping_20261003/THEORY.tex'
tex=source.read_text(encoding='utf-8')
stack=[]
for action,name in re.findall(r'\\(begin|end)\{([^}]+)\}',tex):
    if action=='begin': stack.append(name)
    else: assert stack and stack.pop()==name,(name,stack)
assert not stack
labels=re.findall(r'\\label\{([^}]+)\}',tex)
assert len(labels)==len(set(labels))
assert set(re.findall(r'\\(?:eqref|ref)\{([^}]+)\}',tex))<=set(labels)
save('AUDIT.json',{'archived_artifact_hashes_verified':archived,
 'bundled_input_hashes_verified':len(bundled),'panels_closed_chain':True,
 'count_reintegrated_from_saved_balls':{'real':integral.real.str(35),'imag':integral.imag.str(35),
 'real_lower':integral.real.lower().str(35),'real_upper':integral.real.upper().str(35)},
 'unique_integer_count':1,'latex_environment_and_reference_structure':'PASS',
 'native_pdf_compile':'BLOCKED_ENVIRONMENT','full_spectrum_verified':False})
save('STATUS.json',{'manuscript':'COMPLETE_RESEARCH_DRAFT_WITH_EXPLICIT_SUBMISSION_GATES',
 'claim_status':'NUMERICALLY_VALIDATED','mathematical_derivations':'INDEPENDENTLY_REVIEWED',
 'one_cluster_continuous_certificate':True,'whole_parameter_box_for_this_cluster':True,
 'full_DDE_spectral_certificate':False,'new_nonlinear_simulations':0,
 'new_IEEE39_replacement_optimization':False,
 'informative_physical_replacement_upper_bound':None,
 'verified_coarse_remainder_upper_per_second':4.458179,
 'historical_feasible_GFL_MW':4732.81871482147,'historical_retained_SG_MW':669.942375157377,
 'historical_GFL_percent':87.6,'historical_claim_status':'NUMERICALLY_VALIDATED',
 'optimal_rho_Kp_Ki':False,'global_or_near_global_gap':None,
 'poster_ready_for_maximum_replacement_claim':False,'submission_ready':False,
 'remaining_gates':['Informative verified bound or strict controller-exclusion witness',
 'Fair same-contract co-design comparison and paired nonlinear comparator',
 'Material full-model decision explained by collective terms',
 'Native LaTeX compilation and visual review'],
 'date':'2026-10-03'})
head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
paths=[p for p in OUT.rglob('*') if p.is_file() and not any(x in p.parts for x in ['vendor','__pycache__'])
       and p.name!='FINAL_MANIFEST.json']
save('FINAL_MANIFEST.json',{'git_base_HEAD':head,'manuscript_path':str(source.relative_to(ROOT)),
 'manuscript_sha256':sha(source),'artifacts':{p.relative_to(OUT).as_posix():sha(p) for p in sorted(paths)}})
print((OUT/'AUDIT.json').read_text())
