"""BC02: which representation of the binary family actually exists?

Statements A, B and C are in src/ibr_cycles/certification/binary.py.
Families (configs/binary_certification_v1/BC00_config.yaml, bc02):
  IEEE-39  core {30,33,35,37}, 16 vertices, at the six F8 points
  Kundur   candidates {2,3,4}, 8 vertices, at K12A (g, k) = (0, 0.75),
           (0.08, 1.25), (0.11, 1.25), the G2 cases, declared here before the run
  IEEE-68  candidates {9,6,3,4}, 16 vertices at k = 1, g in {0, 0.25, 1};
           twelve plants G1-G12: affinity on the empty vertex, the 12
           singletons and 64 random vertices (seed 20260911), and vertex
           equivalence on 12 of those random vertices

For every vertex:
  A  control affinity for fixed S: A(theta) from anchors vs the direct path at
     a random theta (Kundur; IEEE-39 and IEEE-68 were verified in F7 and G3)
  B  port locality: T_S(s) against T_0(s) + sum_{i in S} [T_{i}(s) - T_0(s)] at
     three values of s, the support of each increment, and
     det P = h det T
  C  common realization: J(delta) = J0 + sum delta_i J_i (C1 ghost and C2
     filler), f = g = 0 at every vertex, vertex equivalence
     (spectrum = actual portfolio + ghost or filler), ghost Hurwitz, and the
     second difference of the REDUCED A (non-affinity)
  ledger: RHP count of every device block (h_S), to decide whether dY_i may be
     assumed stable (a BC03 prerequisite)
"""

from __future__ import annotations

import json
import sys
import time
from multiprocessing import Pool

import numpy as np
from _bc import WORKERS, BCExperiment, out_dir, write_json
from scipy.optimize import linear_sum_assignment

from _f7_common import LEAK, Theta  # noqa: E402
from _overnight import pin_blas_threads  # noqa: E402
from F8_service_attribution import r_configs, solve_config  # noqa: E402
from F12_kundur import net as k_net  # noqa: E402
from F12_kundur import solve as k_solve  # noqa: E402
from G3_ieee68 import NETWORK as NET68  # noqa: E402
from G3_ieee68 import solve as s68  # noqa: E402
from ibr_cycles.certification.binary import (  # noqa: E402
    all_vertices,
    build_common_realization,
    reduced,
    stacked,
)
from ibr_cycles.dynamics.linearize import central_difference_jacobians  # noqa: E402
from ibr_cycles.models.ieee39_devices import ConverterParameters  # noqa: E402
from ibr_cycles.models.ieee39_network import load_network  # noqa: E402

OUT = out_dir("BC02")
RES = OUT.parent.parent
SEED = 20260911
S_VALUES = (0.5 + 2.0j, 1e-1 + 4.0j, 2.0 + 0.5j)


def _dist(a, b):
    if a.size != b.size:
        return float("nan")
    c = np.abs(a[:, None] - b[None, :])
    r, cc = linear_sum_assignment(c)
    return float(c[r, cc].max())


def port_t(jac, s):
    n = jac.fx.shape[0]
    return jac.gz + jac.gx @ np.linalg.solve(s * np.eye(n) - jac.fx, jac.fz)


def det_identity(jac, s):
    """|log det P(s) - log[det(sI - f_x) det T(s)]| (complex log, modulo 2 pi i)."""

    n = jac.fx.shape[0]
    p = np.block([[s * np.eye(n) - jac.fx, -jac.fz], [jac.gx, jac.gz]])
    sp, lp = np.linalg.slogdet(p)
    sh, lh = np.linalg.slogdet(s * np.eye(n) - jac.fx)
    st, lt = np.linalg.slogdet(port_t(jac, s))
    return float(abs(lp - lh - lt) + abs(sp - sh * st))


def device_ledger(case):
    jac = central_difference_jacobians(
        case.dae, case.equilibrium.x, case.equilibrium.z, {}
    )
    out = []
    for slot in case.dae.slots:
        block = jac.fx[slot.start : slot.stop, slot.start : slot.stop]
        ev = np.linalg.eigvals(block)
        out.append(
            (slot.bus, slot.kind, float(ev.real.max()), int((ev.real > 1e-9).sum()))
        )
    return jac, out


# --------------------------------------------------------------- one family --


def _spec(spec):
    kind = spec[0]
    if kind == "39":
        return _b39(spec[1])
    if kind == "K":
        return _bk(spec[1], spec[2])
    return _b68(spec[1])


def family(task):
    name, candidates, spec, sample = task
    builder, realization_kwargs = _spec(spec)
    started = time.time()
    rows = []
    cr = build_common_realization(candidates, **realization_kwargs)
    build_s = time.time() - started
    vertices = list(all_vertices(candidates)) if sample is None else sample
    jac_c = {}
    for delta, s in vertices:
        for var in ("C1", "C2"):
            jac, res, _ = cr.jacobian(delta, var)
            jac_c[(s, var)] = (stacked(jac), reduced(jac), res)
    j0 = {v: jac_c[((), v)][0] for v in ("C1", "C2")}
    a0 = jac_c[((), "C1")][1]
    t0_jac = None
    for delta, s in vertices:
        row = {
            "family": name,
            "vertex": "+".join(map(str, s)) or "EMPTY",
            "size": len(s),
        }
        for var in ("C1", "C2"):
            jmat, amat, res = jac_c[(s, var)]
            singles = [((b,), var) for b in s]
            if all(k in jac_c for k in singles):
                pred = j0[var] + sum(jac_c[k][0] - j0[var] for k in singles)
                row[f"{var}_affinity_rel"] = float(
                    np.linalg.norm(jmat - pred) / np.linalg.norm(jmat)
                )
            row[f"{var}_equilibrium_residual"] = res
        # reduced non-affinity (C1): second difference over pairs
        if len(s) >= 2 and all(((b,), "C1") in jac_c for b in s):
            pred = a0 + sum(jac_c[((b,), "C1")][1] - a0 for b in s)
            amat = jac_c[(s, "C1")][1]
            row["reduced_nonaffinity_rel"] = float(
                np.linalg.norm(amat - pred) / np.linalg.norm(amat)
            )
        if builder is not None and (
            sample is None or row["vertex"] in sample_equiv(sample)
        ):
            case = builder(s)
            jac_s, ledger = device_ledger(case)
            actual = np.linalg.eigvals(case.system.A)
            ghosts = np.concatenate([e for _, _, e in cr.absent_blocks(delta)])
            ev_c1 = np.linalg.eigvals(jac_c[(s, "C1")][1])
            ev_c2 = np.linalg.eigvals(jac_c[(s, "C2")][1])
            row["C1_vertex_equiv_dist"] = _dist(ev_c1, np.concatenate([actual, ghosts]))
            row["C2_vertex_equiv_dist"] = _dist(
                ev_c2, np.concatenate([actual, -np.ones(ev_c2.size - actual.size)])
            )
            row["ghost_max_re"] = (
                float(ghosts.real.max()) if ghosts.size else float("nan")
            )
            row["device_rhp_blocks"] = sum(1 for _, _, _, n in ledger if n)
            row["det_identity_logerr"] = det_identity(jac_s, S_VALUES[0])
            row["device_max_re"] = max(m for _, _, m, _ in ledger)
            # B: port locality against the empty vertex and singletons
            if not s:
                t0_jac = jac_s
                row["port_locality_rel"] = 0.0
            row["_jac"] = jac_s
        rows.append(row)
    # B needs T_0 and the singleton ports at the same s values
    port = {r["vertex"]: r.pop("_jac") for r in rows if "_jac" in r}
    if t0_jac is not None:
        for r in rows:
            if r["vertex"] not in port or r["size"] < 2:
                continue
            members = r["vertex"].split("+")
            if not all(m in port for m in members):
                continue
            errs = []
            for sv in S_VALUES:
                t0 = port_t(t0_jac, sv)
                ts = port_t(port[r["vertex"]], sv)
                pred = t0 + sum(port_t(port[m], sv) - t0 for m in members)
                errs.append(np.linalg.norm(ts - pred) / np.linalg.norm(ts))
            r["port_locality_rel"] = float(max(errs))
        for m, jac in port.items():
            if "+" in m or m == "EMPTY":
                continue
            d = port_t(jac, S_VALUES[0]) - port_t(t0_jac, S_VALUES[0])
            big = np.abs(d) > 1e-8 * np.abs(d).max()
            rows_idx, cols_idx = np.nonzero(big)
            for r in rows:
                if r["vertex"] == m:
                    r["port_increment_support_size"] = int(
                        len(set(rows_idx.tolist()) | set(cols_idx.tolist()))
                    )
    return {
        "family": name,
        "rows": rows,
        "build_s": build_s,
        "elapsed_s": time.time() - started,
    }


def sample_equiv(sample):
    return {"+".join(map(str, s)) or "EMPTY" for _, s in sample[: 1 + 12 + 12]}


# -------------------------------------------------------------- families --


def _b39(point):
    th = Theta(point["g"], point["k"], point["t"], point["h"])
    cfg = {c["name"]: c for c in r_configs()}["R_none"]
    conv = ConverterParameters(
        voltage_control=True, voltage_gain=th.g, voltage_leak=LEAK
    )
    kwargs = {"converter": conv, "machine_scaling": {"ka": th.k, "ta": th.t}}

    def build(s):
        return solve_config(s, th, cfg)

    return build, kwargs


def _bk(g, k):
    th = {"g": g, "k": k, "t": 1.0}
    conv = ConverterParameters(voltage_control=True, voltage_gain=g, voltage_leak=0.05)
    kwargs = {
        "network": k_net(),
        "converter": conv,
        "machine_scaling": {"ka": k, "ta": 1.0},
    }

    def build(s):
        return k_solve(s, th)

    return build, kwargs


def _b68(g, k=1.0):
    conv = ConverterParameters(voltage_control=True, voltage_gain=g, voltage_leak=0.05)
    kwargs = {
        "network": load_network(NET68),
        "converter": conv,
        "machine_scaling": {"ka": k},
    }

    def build(s):
        return s68(s, g, k)

    return build, kwargs


def tasks():
    f8 = json.loads((RES / "F8" / "F8_points.json").read_text(encoding="utf-8"))
    out = []
    for n, p in f8.items():
        out.append((f"IEEE-39 {n}", (30, 33, 35, 37), ("39", p), None))
    for g, k in ((0.0, 0.75), (0.08, 1.25), (0.11, 1.25)):
        out.append((f"Kundur K12A g={g} k={k}", (2, 3, 4), ("K", g, k), None))
    for g in (0.0, 0.25, 1.0):
        out.append((f"IEEE-68 cand g={g}", (3, 4, 6, 9), ("68", g), None))
    rng = np.random.default_rng(SEED)
    plants = tuple(range(1, 13))
    sample = [({b: 0 for b in plants}, ())]
    sample += [({b: int(b == p) for b in plants}, (p,)) for p in plants]
    seen = {()} | {(p,) for p in plants}
    while len(sample) < 1 + 12 + 64:
        bits = rng.integers(0, 2, size=12)
        s = tuple(p for p, bit in zip(plants, bits, strict=True) if bit)
        if s in seen:
            continue
        seen.add(s)
        sample.append(({b: int(b in s) for b in plants}, s))
    out.insert(0, ("IEEE-68 twelve g=0", plants, ("68", 0.0), sample))
    return out


def main(argv) -> int:
    exp = BCExperiment(
        name="BC02_binary_representation",
        question="Is there an exact common realization affine in the indicators?",
        config={"seed": SEED, "s_values": [str(s) for s in S_VALUES]},
    )
    started = time.time()
    with Pool(min(WORKERS, 13), initializer=pin_blas_threads) as pool:
        results = pool.map(family, tasks(), chunksize=1)
    import pandas as pd

    frame = pd.DataFrame([r for res in results for r in res["rows"]])
    frame.to_csv(OUT / "BC02_vertices.csv", index=False)
    summary = {}
    for res in results:
        f = frame[frame.family == res["family"]]
        summary[res["family"]] = {
            "vertices": int(len(f)),
            "build_s": round(res["build_s"], 1),
            "elapsed_s": round(res["elapsed_s"], 1),
            "C1_affinity_max": float(
                f.get("C1_affinity_rel", pd.Series(dtype=float)).max()
            ),
            "C2_affinity_max": float(
                f.get("C2_affinity_rel", pd.Series(dtype=float)).max()
            ),
            "equilibrium_residual_max": float(
                f[["C1_equilibrium_residual", "C2_equilibrium_residual"]]
                .to_numpy()
                .max()
            ),
            "C1_vertex_equiv_max": float(
                f.get("C1_vertex_equiv_dist", pd.Series(dtype=float)).max()
            ),
            "C2_vertex_equiv_max": float(
                f.get("C2_vertex_equiv_dist", pd.Series(dtype=float)).max()
            ),
            "ghost_max_re": float(f.get("ghost_max_re", pd.Series(dtype=float)).max()),
            "device_rhp_blocks_max": float(
                f.get("device_rhp_blocks", pd.Series(dtype=float)).max()
            ),
            "reduced_nonaffinity_max": float(
                f.get("reduced_nonaffinity_rel", pd.Series(dtype=float)).max()
            ),
            "port_locality_max": float(
                f.get("port_locality_rel", pd.Series(dtype=float)).max()
            ),
            "port_increment_support": sorted(
                set(
                    f.get("port_increment_support_size", pd.Series(dtype=float))
                    .dropna()
                    .astype(int)
                )
            ),
        }
    write_json(OUT / "BC02_summary.json", summary)
    exp.finish("COMPUTED", elapsed_s=time.time() - started, **summary)
    print(json.dumps(summary, indent=1, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
