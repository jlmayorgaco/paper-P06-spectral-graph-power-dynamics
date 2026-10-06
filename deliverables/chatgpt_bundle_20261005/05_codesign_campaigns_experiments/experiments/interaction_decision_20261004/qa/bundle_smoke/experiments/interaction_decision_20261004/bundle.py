"""Portable dependency archive and complete source/evidence bundle."""
from model import *
import zipfile
def main():
    rels=['experiments/graph_gsp_codesign_20261003/baseline.toml','experiments/graph_gsp_codesign_20261003/TABLE_02_BASELINE_ROOTS.csv',
        'experiments/graph_gsp_codesign_20261003/DelayedEvents.jl',
        'experiments/nonlinear_codesign_20261001/ReducedDAE.jl',
        'experiments/analytical_delay_codesign_mega_20261002/m3_a_trace_integral.jl',
        'experiments/delay_dressed_replacement_frontier_20261002/DelayCharacteristic.jl',
        'src/bnd_model_expN/PDExactDesignN.jl','src/bnd_design_e/CollectiveModel.jl',
        'src/bnd_design/AnalyticSG.jl','src/bnd_design/AnalyticGFLPLL.jl',
        'reports/experiment_N/TABLE_N01_original_operating_point.csv','reports/experiment_N/TABLE_N01_original_operating_point.csv.sha256',
        'reports/experiment_A/matrices/bus33_baseline_equilibrium.csv',
        'experiments/physical_collective_damping_20261003/source_snapshot/Project.toml',
        'experiments/physical_collective_damping_20261003/source_snapshot/Manifest.toml',
        'reports/poster/ias2026/collective_interaction_20261003/generated/figures/rank_two.pdf',
        'reports/poster/ias2026/collective_interaction_20261003/generated/figures/events.pdf']
    inputs=[ROOT/x for x in rels]+list((SOURCE/'model').glob('*.csv'))+list((ROOT/'reports/experiment_D/inputs').glob('*.csv'))
    save('DEPENDENCY_MANIFEST.json',{'files':{str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}})
    with zipfile.ZipFile(OUT/'DEPENDENCIES.zip','w',zipfile.ZIP_DEFLATED) as z:
        for p in inputs:z.write(p,str(p.relative_to(ROOT)).replace('\\','/'))
    # All scientific intermediates and logs, excluding TeX scratch and visual QA.
    skip={'.aux','.out','.pyc'};archive=OUT/'INTERACTION_DECISION_20261004_COMPLETE.zip'
    files=[p for p in OUT.rglob('*') if p.is_file() and p.suffix not in skip and '__pycache__' not in p.parts and 'qa' not in p.parts
        and p.name not in [archive.name,archive.name+'.sha256','BUNDLE_MANIFEST.json']]
    files.append(ROOT/'experiments/theory_collective_damping_20261003/THEORY.tex')
    save('BUNDLE_MANIFEST.json',{'files':{str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}})
    files.append(OUT/'BUNDLE_MANIFEST.json')
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
        for p in files:z.write(p,str(p.relative_to(ROOT)).replace('\\','/'))
    (OUT/(archive.name+'.sha256')).write_text(hashlib.sha256(archive.read_bytes()).hexdigest()+'  '+archive.name+'\n')
    print('Bundle',archive.name,'files',len(files),'bytes',archive.stat().st_size)
if __name__=='__main__':main()
