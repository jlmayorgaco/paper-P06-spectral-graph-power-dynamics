"""Pilot: run the first task of each phase in-process (not checkpointed) to catch bugs."""
import sys, time, json, importlib
import _infra  # noqa
only = sys.argv[1].split(",") if len(sys.argv) > 1 else None
from CDW_MASTER_RUN import PHASES
for name, module, tfn, rfn, afn, deps, pre in PHASES:
    if not tfn or (only and name not in only):
        continue
    if name in ("E08b",):
        continue
    mod = importlib.import_module(module)
    t0 = time.time(); ts = getattr(mod, tfn)(); t1 = time.time()
    task = ts[0] if name != "E09" else ts[0]
    try:
        r = getattr(mod, rfn)(task); ok = True
    except Exception as e:
        import traceback; traceback.print_exc(); r = str(e); ok = False
    print(f"{name}: n_tasks={len(ts)} gen={t1-t0:.1f}s run1={time.time()-t1:.1f}s ok={ok} out={json.dumps(r, default=str)[:300]}", flush=True)
