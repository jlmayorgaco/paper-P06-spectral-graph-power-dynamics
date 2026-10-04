"""Enclose moment-neutrality residual of the actual rounded gain vectors."""
from moments import *
sys.path.insert(0,str(ROOT/'experiments/regional_paper_closure_20261003/vendor'))
from flint import arb,ctx
ctx.prec=192
def main():
    frame=pd.read_csv(OUT/'TABLE_16_INTERVAL_WEIGHTS.csv')
    weights=[arb(s) for s in frame.interval]
    rows=[]
    for path in sorted((OUT/'globalized/designs').glob('*.toml')):
        d=tomllib.loads(path.read_text());kappa=arb(0)
        for i,w in enumerate(weights):
            ki=arb(float(P[20+i]));dp=arb(float(d['Kp'][i]))-arb(float(P[10+i]))
            dt=arb(float(d['tau_vector'][i]))-arb(float(.04))
            kappa+=w*(dt/ki-dp/ki**2)
        rows.append(dict(design=path.stem,kappa_interval=kappa.str(30),
            abs_upper=str(abs(kappa).upper()),neutral_within_1e_12=bool(abs(kappa)<arb('1e-12')),
            sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    pd.DataFrame(rows).to_csv(OUT/'TABLE_17_DESIGN_MOMENT_ENCLOSURES.csv',index=False)
    collective=[r for r in rows if r['design'].endswith('_collective')]
    assert len(collective)==14 and all(r['neutral_within_1e_12'] for r in collective)
    print('All14 collective finite designs satisfy enclosed |kappa|<1e-12; exact zero is not claimed.')
if __name__=='__main__':main()
