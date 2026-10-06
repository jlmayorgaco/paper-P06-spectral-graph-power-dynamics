"""Snapshot exact scientific inputs without mutating historical source/results."""
from moments import *
import subprocess,platform,scipy,cvxpy
def main():
    import zipfile
    paths={ROOT/x for x in json.loads((OLD/'DEPENDENCY_MANIFEST.json').read_text())['files']}
    paths.update([OLD/'model.py',OLD/'designs/corrected.toml',OLD/'TABLE_08_MULTIBRANCH_CONTINUATION.csv',OLD/'TABLE_12_NONLINEAR_EVENTS.csv'])
    paths.update((OUT/'model').glob('*.csv'))
    hashes={str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}
    data={'inputs':hashes,'git_HEAD':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
          'git_dirty_at_snapshot':subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True),
          'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'cvxpy':cvxpy.__version__,
          'model':'204-state full physical-supply IEEE39, ten filtered type-II PLLs, pure detector delay',
          'fixed_GFL_MW':float(PORTS.P0@P[:10]),'total_dispatch_MW':float(PORTS.P0.sum()),
          'input_event':'constant-impedance load parameter, not an exactly constant-MW trajectory',
          'no_current_limiter_claim':True}
    dest=OUT/'BASELINE_MANIFEST.json'
    if dest.exists():assert json.loads(dest.read_text())['inputs']==hashes
    else:dump(dest.name,data)
    val=[OUT/'validate_julia.jl',OUT/'globalized/designs/uniform_1ms_collective.toml']
    val+=[OUT/'designs'/f'{x}.toml' for x in ['base','uniform_1ms_unchanged','uniform_1ms_sitewise']]
    lock={'files':{str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in val},
          'spectrum_gamma':-.05,'event_horizon_s':60,'window_s':.5,'events':[[8,-100],[16,100],[16,-100],[29,100],[29,-100]],
          'limits':{'F_Hz':.5,'R_Hz_s':.5,'Vmin':.9,'Vmax':1.1,'SG_normalized_actuator_slack':.002},
          'solver':'exact exponential contour; nonlinear method of steps with adaptive Rosenbrock, tol1e-9,dtmax.01'}
    v=OUT/'VALIDATION_LOCK.json'
    if v.exists():assert json.loads(v.read_text())==lock
    else:dump(v.name,lock)
    print('INPUTS',len(hashes),'HEAD',data['git_HEAD'])
if __name__=='__main__':main()
