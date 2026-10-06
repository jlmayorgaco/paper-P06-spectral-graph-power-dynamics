"""Run E3 + E4 cases in parallel (<=12 workers); each case -> raw/runs/<label>.json (+ designs toml if success). Resumable."""
import sys, pathlib, json, os
sys.path.insert(0, str(pathlib.Path(__file__).parent))
from repair_nhop import *
from multiprocessing import Pool
RUNS = CAMP / 'raw' / 'runs'; RUNS.mkdir(parents=True, exist_ok=True); (CAMP / 'raw' / 'designs').mkdir(parents=True, exist_ok=True)


def cases():
    out = []
    for g, br in (('frozen', False), ('bridge', True)):
        diam = comm_graph(br)['diam']
        for n in range(diam + 1):
            for tau in (44, 48, 52): out.append((g, n, tau, 0.875))
            for rho in (0.90, 0.925, 0.95): out.append((g, n, 44, rho))   # (44 ms, 0.875) shared by E3 and E4
    return out


def job(c):
    g, n, tau, rho = c; lab = f'{g}_n{n}_tau{tau}_rho{int(round(rho * 1000))}'
    f = RUNS / f'{lab}.json'
    if f.exists(): return lab
    dist = comm_graph(g == 'bridge')['dist']
    Kp, KI, hist, fin = repair_nhop(dist, n, tau / 1000, rho, label=lab)
    fin.update(label=lab, graph=g); f.write_text(json.dumps(dict(final=fin, hist=hist)))
    if fin['success']: write_toml(CAMP / 'raw' / 'designs' / f'{lab}.toml', rho, tau / 1000, Kp, KI)
    print('FINAL', fin, flush=True); return lab


if __name__ == '__main__':
    cs = cases()
    if len(sys.argv) > 1 and sys.argv[1] == 'test': cs = [c for c in cs if c == ('frozen', 1, 44, 0.875)]
    cs = sorted(cs, key=lambda c: -c[1])
    with Pool(min(12, len(cs))) as P: list(P.imap_unordered(job, cs, chunksize=1))
