"""Final evidence gates and source integrity; no repeated simulation."""
from moments import *
from pypdf import PdfReader
import subprocess
def main():
    baseline=json.loads((OUT/'BASELINE_MANIFEST.json').read_text())
    for name,expected in baseline['inputs'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==expected,name
    validation=json.loads((OUT/'VALIDATION_LOCK.json').read_text())
    for name,expected in validation['files'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==expected,name
    for name,expected in json.loads((OUT/'MOMENT_INPUT_LOCK_V1.json').read_text()).items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==expected,name
    lock=json.loads((OUT/'globalized/SOURCE_LOCK.json').read_text())
    assert hashlib.sha256((OUT/'codesign_globalized.py').read_bytes()).hexdigest()==lock['code_sha256']
    spectra=pd.read_csv(OUT/'TABLE_11_FULL_SPECTRUM.csv')
    assert spectra['count'].tolist()==[0,6,6,0]
    ev=pd.read_csv(OUT/'TABLE_12_NONLINEAR_EVENTS.csv');assert len(ev)==5 and ev['pass'].all()
    assert ev.F.max()<=.5 and ev.R.max()<=.5 and ev.Vmin.min()>=.9 and ev.Vmax.max()<=1.1 and ev.slack.min()>=.002
    co=pd.read_csv(OUT/'globalized/TABLE_13_GLOBALIZED.csv');cc=co[co.method=='collective']
    assert len(cc)==14 and cc.catalog_margin_pass.all() and cc.gain_bounds_pass.all()
    assert not json.loads((OUT/'globalized/FAILURES.json').read_text())
    mom=pd.read_csv(OUT/'TABLE_17_DESIGN_MOMENT_ENCLOSURES.csv');cm=mom[mom.design.str.endswith('_collective')]
    assert len(cm)==14 and cm.neutral_within_1e_12.all()
    assert json.loads((OUT/'MOMENT_ASSUMPTION_CERTIFICATE.json').read_text())['all_weight_signs_enclosed']
    pdf=OUT/'Beyond_Nodal_Damping_Network_Moment_Laws_20261004.pdf'
    theory=ROOT/'experiments/theory_collective_damping_20261003/THEORY.tex'
    log=(OUT/'latex_build_second.log').read_text(errors='replace')
    assert 'Output written' in log and 'Warning' not in log and 'Overfull' not in log and 'undefined' not in log
    tex={'native_compiler':'UNAVAILABLE: Unable to find standard directories for platform',
        'fallback':'existing MiKTeX pdflatex, two successful passes, no warnings',
        'editor_source':str(theory),'source_sha256':hashlib.sha256(theory.read_bytes()).hexdigest(),
        'pdf_sha256':hashlib.sha256(pdf.read_bytes()).hexdigest(),'pages':len(PdfReader(pdf).pages)}
    dump('LATEX_COMPILATION.json',tex)
    status={'frozen_input_integrity':True,'headline_contour_counts':spectra['count'].tolist(),
        'nonlinear_events_passed':len(ev),'catalog_only_collective_patterns_passed':len(cc),
        'all14_rounded_moment_residuals_enclosed_below_1e_12':True,'warnings_in_latex':False,
        'journal_ready':False,'optimality':'UNPROVED','communication_controller':'THEOREM_ONLY',
        'headline_optimizer':'MAX_ITERATIONS','modal_only_optimizer':'LINE_SEARCH_STALLED'}
    dump('FINAL_AUDIT.json',status)
    paths=[p for p in OUT.rglob('*') if p.is_file() and p.suffix in {'.py','.jl','.toml','.csv','.tex','.md'}
        and not any(x in p.parts for x in ['__pycache__','qa'])]
    dump('EXPERIMENT_MANIFEST.json',{'git_parent':baseline['git_HEAD'],'claim_status':'NUMERICALLY_VALIDATED',
        'files':{str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
        'main_source_sha256':tex['source_sha256'],'publication_limit':'model-specific theory and targeted example; not maximum replacement'})
    print(json.dumps(status,indent=2))
if __name__=='__main__':main()
