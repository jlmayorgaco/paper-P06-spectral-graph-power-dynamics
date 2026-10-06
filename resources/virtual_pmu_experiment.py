#!/usr/bin/env python3
"""
virtual_pmu_experiment.py
=========================
End-to-end test of the sparse-PMU "virtual PMU" estimator on the IEEE-39 grid.

  1. Simulates an event with Julia/PowerDynamics v5.0.0 (native IEEE-39: Sauer-Pai
     machines + AVR + governors, constant-Z loads). The event runs from T0 to TF and is
     then removed, so the grid returns to normal.
  2. Keeps ONLY what the 8 PMUs would measure: bus voltage + one branch current each
     (branches identified in the SGSMA data), plus PMU-like noise.
  3. Estimates the 39 bus voltages, all branch currents/flows and the bus frequencies
     from those 8 PMUs with the injection-coordinate estimator:
        V = V0 + M_G g (+ z_j a)   with   M_G = A0^-1 C_G,   A0 = Y_N + diag(y_L0)
        - H0 : 10 generator-current deviations g (nominal loads), plus the physical
               pseudo-measurement for the G4/G5 pocket (K34 g33 - K33 g34 = 0)
        - H1 : H0 + an unknown current injection a at the event bus (event location known)
        Solved frame-by-frame (GLS/MMSE) and with a random-walk Kalman filter + RTS smoother.
  4. Compares against the simulated truth at the 31 non-PMU buses: plots, RMSE, angle
     error, 95% coverage, branch flows, frequency, and two baselines (nominal power flow,
     nearest PMU). A pre-registered rule then says whether it qualifies as an estimator
     of the non-PMU nodes.

Requirements
  python >= 3.9 : numpy scipy pandas matplotlib
  julia  >= 1.10 on PATH (or set JULIA=/path/to/julia). The first run installs
  PowerDynamics 5.0.0 and friends in ./julia_env (this can take 20-40 min, then it is cached).

Run
  python virtual_pmu_experiment.py
  (REUSE_SIM=1 python virtual_pmu_experiment.py  re-runs only the estimator on the last simulation)
"""
import json, os, subprocess, sys, textwrap, shutil
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ============================ CONFIGURATION ===================================
EVENT_TYPE = "load"      # "load" (constant-Z load scaled at EVENT_BUS) | "gen" (governor p_ref at a generator bus 30..38)
EVENT_BUS = 7            # a NON-PMU bus is the hard case
EVENT_SCALE = 0.20       # +20 % load admittance (or +20 % p_ref for "gen")
T0, TF, TEND = 5.0, 10.0, 20.0   # event on at T0, removed at TF, simulation until TEND [s]
FPS = 30                 # PMU reporting rate
ADD_NOISE = True         # PMU noise (values measured on the SGSMA data)
SIGMA_VMAG = 0.0027      # pu, voltage magnitude noise
SIGMA_VANG_DEG = 0.004   # deg, phasor angle noise (V and I)
SIGMA_IMAG_REL = 0.003   # relative current magnitude noise
SIGMA_POCKET = 0.02      # pu current, uncertainty of the G4/G5 droop-sharing prior
SIGMA_X_PRIOR = 10.0     # pu current, (weak) prior on generator-current deviations
SIGMA_A_PRIOR = 10.0     # pu current, (weak) prior on the H1 event injection
KF_Q_STD = 0.02          # pu current per frame, random-walk std of the Kalman variant
SEED = 1
OUT = os.path.abspath("vpmu_out")
JULIA = os.environ.get("JULIA", "julia")
JULIA_ENV = os.environ.get("JULIA_PROJECT_DIR", os.path.abspath("julia_env"))
SIM_JL = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sim_ieee39_event.jl")

PMU = [2, 5, 6, 10, 19, 22, 29, 39]
# branch current measured by each PMU (current leaving the PMU bus into this branch)
PMU_BRANCH = {2: 3, 5: 6, 6: 7, 10: 11, 19: 16, 22: 21, 29: 28, 39: 1}
POCKET = [20, 33, 34]    # {20,33,34}: only seen through line 16-19 -> needs the droop prior
NB = 39
# ==============================================================================


def setup_julia():
    if os.path.exists(os.path.join(JULIA_ENV, "Manifest.toml")):
        return
    print(f"[setup] creating Julia environment in {JULIA_ENV} (first run only)...")
    os.makedirs(JULIA_ENV, exist_ok=True)
    code = ('using Pkg; Pkg.activate(ARGS[1]); Pkg.add([PackageSpec(name="PowerDynamics", version="5.0.0"),'
            'PackageSpec(name="NetworkDynamics"), PackageSpec(name="ModelingToolkitBase"), PackageSpec(name="CSV"),'
            'PackageSpec(name="DataFrames"), PackageSpec(name="OrdinaryDiffEqRosenbrock"),'
            'PackageSpec(name="OrdinaryDiffEqNonlinearSolve"), PackageSpec(name="SciMLBase"), PackageSpec(name="JSON")]);'
            'Pkg.precompile()')
    subprocess.run([JULIA, "-e", code, JULIA_ENV], check=True)


def simulate():
    if not os.path.exists(SIM_JL):
        sys.exit(f"missing {SIM_JL} (put sim_ieee39_event.jl next to this script)")
    os.makedirs(OUT, exist_ok=True)
    cmd = [JULIA, f"--project={JULIA_ENV}", SIM_JL, OUT, EVENT_TYPE, str(EVENT_BUS), str(EVENT_SCALE),
           str(T0), str(TF), str(TEND), str(FPS)]
    print("[julia]", " ".join(cmd))
    subprocess.run(cmd, check=True)


# ----------------------------- network model -----------------------------------
def network(datadir):
    bus = pd.read_csv(os.path.join(datadir, "bus.csv"))
    br = pd.read_csv(os.path.join(datadir, "branch.csv"))
    ld = pd.read_csv(os.path.join(datadir, "loads.csv"))          # initialised Pset, Qset, Vset (from Julia)
    mach = pd.read_csv(os.path.join(datadir, "machine.csv"))
    gov = pd.read_csv(os.path.join(datadir, "gov.csv"))
    blocks = []                                                      # PowerDynamics PiLine convention
    for _, r in br.iterrows():
        rs, rd = r.r_src, (r.r_dst if "r_dst" in br.columns else 1.0)
        y = 1 / complex(r.R, r.X); ys = complex(r.G_src, r.B_src); yd = complex(r.G_dst, r.B_dst)
        Yl = np.array([[rs * rs * (y + ys), -rs * rd * y], [-rs * rd * y, rd * rd * (y + yd)]])
        blocks.append((int(r.src_bus) - 1, int(r.dst_bus) - 1, Yl))
    YN = np.zeros((NB, NB), complex)
    for i, j, Yl in blocks:
        YN[np.ix_([i, j], [i, j])] += Yl
    # power flow (bus-level P,Q,V from bus.csv), Newton with a numerical Jacobian
    typ = bus.bus_type.values; P = bus.P.values.astype(float); Q = bus.Q.values.astype(float)
    Vm = np.where(np.isnan(bus.V.values.astype(float)), 1.0, bus.V.values.astype(float)); Va = np.zeros(NB)
    ai = np.where(typ != "Slack")[0]; mi = np.where(typ == "PQ")[0]
    def mis(x):
        a = Va.copy(); m = Vm.copy(); a[ai] = x[:len(ai)]; m[mi] = x[len(ai):]
        V = m * np.exp(1j * a); S = V * np.conj(YN @ V)
        return np.r_[(P - S.real)[ai], (Q - S.imag)[mi]]
    x = np.r_[Va[ai], Vm[mi]]
    for _ in range(30):
        F = mis(x)
        if np.max(abs(F)) < 1e-12: break
        J = np.array([(mis(x + 1e-7 * e) - mis(x - 1e-7 * e)) / 2e-7 for e in np.eye(len(x))]).T
        x = x - np.linalg.solve(J, F)
    Va[ai] = x[:len(ai)]; Vm[mi] = x[len(ai):]; V0 = Vm * np.exp(1j * Va)
    yL = np.zeros(NB, complex)
    for _, r in ld.iterrows():                                       # consumption admittance
        yL[int(r.bus) - 1] = -np.conj(complex(r.Pset, r.Qset)) / r.Vset ** 2
    A0 = YN + np.diag(yL)
    gb = [int(b) - 1 for b in mach.bus]
    CG = np.zeros((NB, len(gb))); CG[gb, range(len(gb))] = 1
    MG = np.linalg.solve(A0, CG)
    IG0 = (A0 @ V0)[gb]
    nong = [i for i in range(NB) if i not in gb]
    gate = dict(pf_kcl_nongen=float(np.max(abs((A0 @ V0)[nong]))), V0_vs_MG_IG0=float(np.max(abs(V0 - MG @ IG0))))
    # governor gains for the G4/G5 pocket prior  K = Sn / R
    K = {int(r.bus): float(mach[mach.bus == r.bus].Sn.iloc[0]) / float(r.R) for _, r in gov.iterrows()}
    return dict(YN=YN, blocks=blocks, V0=V0, A0=A0, MG=MG, IG0=IG0, gb=gb, K=K, gate=gate, bus=bus)


def branch_row(blocks, frm, to):
    """row r such that r @ V = current leaving bus `frm` into branch frm-to."""
    for i, j, Yl in blocks:
        if {i, j} == {frm - 1, to - 1}:
            k = 0 if i == frm - 1 else 1
            r = np.zeros(NB, complex); r[i] = Yl[k, 0]; r[j] = Yl[k, 1]; return r
    raise KeyError((frm, to))


def realify(M):
    return np.block([[M.real, -M.imag], [M.imag, M.real]])


def nearest_pmu(blocks):
    adj = {i: set() for i in range(NB)}
    for i, j, _ in blocks: adj[i].add(j); adj[j].add(i)
    near = {}
    for s in range(NB):
        seen, front = {s}, [s]
        while front:
            hit = [b for b in front if b + 1 in PMU]
            if hit: near[s] = hit[0] + 1; break
            nxt = []
            for b in front:
                for c in adj[b]:
                    if c not in seen: seen.add(c); nxt.append(c)
            front = nxt
    return near


# ----------------------------- estimator ----------------------------------------
def estimator_matrices(net, with_event_bus):
    MG, A0, V0 = net["MG"], net["A0"], net["V0"]
    Hrows = [np.eye(NB)[b - 1] for b in PMU] + [branch_row(net["blocks"], b, PMU_BRANCH[b]) for b in PMU]
    H = np.array(Hrows)                                               # 16 x 39 complex
    M = MG.copy()
    if with_event_bus:                                                # H1: + current injection at the event bus
        M = np.hstack([MG, np.linalg.solve(A0, np.eye(NB)[:, [EVENT_BUS - 1]])])
    G = realify(H @ M)                                                # 32 x nx
    nc = M.shape[1]
    # pocket pseudo-measurement:  K34*g33 - K33*g34 = 0  (complex -> 2 real rows)
    gi = [b + 1 for b in net["gb"]]
    c = np.zeros(nc, complex); c[gi.index(33)] = net["K"][34]; c[gi.index(34)] = -net["K"][33]
    c = c / np.abs(c).max()
    D = realify(c[None, :])
    prior = np.full(2 * nc, SIGMA_X_PRIOR ** 2)
    if with_event_bus:
        prior[[nc - 1, 2 * nc - 1]] = SIGMA_A_PRIOR ** 2
    return H, M, G, D, prior


def noise_cov(H, V0):
    """real 32x32 covariance of [Re z; Im z] for polar (magnitude, angle) PMU noise."""
    z0 = H @ V0; n = len(z0); R = np.zeros((2 * n, 2 * n)); sa = np.deg2rad(SIGMA_VANG_DEG)
    for c in range(n):
        sm = SIGMA_VMAG if c < 8 else SIGMA_IMAG_REL * abs(z0[c])
        th = np.angle(z0[c]); Rot = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
        C2 = Rot @ np.diag([sm ** 2, (abs(z0[c]) * sa) ** 2]) @ Rot.T if ADD_NOISE else np.eye(2) * 1e-10
        idx = [c, c + n]; R[np.ix_(idx, idx)] = C2 + np.eye(2) * 1e-12
    return R


def gls(G, D, prior, R, Y):
    """frame-by-frame MMSE: x = argmin |y-Gx|^2_R + |Dx|^2/sp^2 + |x|^2_prior."""
    Ri = np.linalg.inv(R)
    Lam = G.T @ Ri @ G + D.T @ D / SIGMA_POCKET ** 2 + np.diag(1 / prior)
    P = np.linalg.inv(Lam)
    return (P @ G.T @ Ri @ Y), P


def kalman_rts(G, D, prior, R, Y):
    Ck = np.vstack([G, D]); Rk = np.block([[R, np.zeros((R.shape[0], D.shape[0]))],
                                           [np.zeros((D.shape[0], R.shape[0])), np.eye(D.shape[0]) * SIGMA_POCKET ** 2]])
    Yk = np.vstack([Y, np.zeros((D.shape[0], Y.shape[1]))])
    n, N = Ck.shape[1], Y.shape[1]; Q = np.eye(n) * KF_Q_STD ** 2
    m, P = np.zeros(n), np.diag(prior)
    mf, Pf, mp, Pp = np.zeros((N, n)), np.zeros((N, n, n)), np.zeros((N, n)), np.zeros((N, n, n))
    nis = np.zeros(N); Rki = np.linalg.inv(Rk)
    for k in range(N):
        Pm = P + Q if k else P
        mp[k], Pp[k] = m, Pm
        S = Ck @ Pm @ Ck.T + Rk; nu = Yk[:, k] - Ck @ m; nis[k] = nu @ np.linalg.solve(S, nu)
        Lam = np.linalg.inv(Pm) + Ck.T @ Rki @ Ck
        P = np.linalg.inv(Lam); m = P @ (np.linalg.solve(Pm, m) + Ck.T @ Rki @ Yk[:, k])
        mf[k], Pf[k] = m, P
    ms, Ps = mf.copy(), Pf.copy()
    for k in range(N - 2, -1, -1):
        J = Pf[k] @ np.linalg.inv(Pp[k + 1])
        ms[k] = mf[k] + J @ (ms[k + 1] - mp[k + 1]); Ps[k] = Pf[k] + J @ (Ps[k + 1] - Pp[k + 1]) @ J.T
    return ms.T, Ps, nis


# ----------------------------- metrics -------------------------------------------
def bus_voltage_posterior(M, x, P):
    """voltage deviations (39 x N complex) and per-bus 2x2 real covariances."""
    nc = M.shape[1]; g = x[:nc] + 1j * x[nc:]
    Mr = realify(M)
    covs = []
    for Pk in (P if P.ndim == 3 else [P]):
        C = Mr @ Pk @ Mr.T
        covs.append(np.stack([np.array([[C[i, i], C[i, i + NB]], [C[i + NB, i], C[i + NB, i + NB]]]) for i in range(NB)]))
    return M @ g, np.array(covs)


def freq_from_phase(V, t):
    ph = np.unwrap(np.angle(V), axis=1)
    return 60 + np.gradient(ph, t, axis=1) / (2 * np.pi)


def main():
    rng = np.random.default_rng(SEED)
    if os.environ.get("REUSE_SIM") and os.path.exists(os.path.join(OUT, "truth_V.csv")):
        print("[sim] reusing", OUT)
    else:
        setup_julia(); simulate()
    meta = json.load(open(os.path.join(OUT, "meta.json")))
    tv = pd.read_csv(os.path.join(OUT, "truth_V.csv")); t = tv.t.values
    Vt = (tv[[f"Vr{b}" for b in range(1, 40)]].values + 1j * tv[[f"Vi{b}" for b in range(1, 40)]].values).T   # 39 x N
    net = network(OUT); V0 = net["V0"]
    print("[gates]", net["gate"], " |V_sim(t=0) - V0_pf| =", float(np.max(abs(Vt[:, 0] - V0))))

    # ------------- what the 8 PMUs measure -------------
    Hrows = [np.eye(NB)[b - 1] for b in PMU] + [branch_row(net["blocks"], b, PMU_BRANCH[b]) for b in PMU]
    H = np.array(Hrows); Zt = H @ Vt                                   # 16 x N true PMU phasors
    if ADD_NOISE:
        sm = np.r_[np.full(8, SIGMA_VMAG), SIGMA_IMAG_REL * abs(Zt[8:, 0])][:, None]
        Zm = (abs(Zt) + sm * rng.standard_normal(Zt.shape)) * np.exp(1j * (np.angle(Zt) + np.deg2rad(SIGMA_VANG_DEG) * rng.standard_normal(Zt.shape)))
    else:
        Zm = Zt.copy()
    # per-frame angle reference (bus-2 PMU), as in the SGSMA data pipeline
    psi = np.angle(Zm[0] / V0[1])
    Zd = Zm * np.exp(-1j * psi)
    Y = np.vstack([(Zd - (H @ V0)[:, None]).real, (Zd - (H @ V0)[:, None]).imag])        # 32 x N

    hidden = [b for b in range(1, 40) if b not in PMU]
    ident = [b for b in hidden if b not in POCKET]
    rot = np.exp(1j * psi)
    results, estimates = {}, {}
    variants = [("H0", False)] + ([("H1(event bus known)", True)] if EVENT_TYPE == "load" else [])
    for name, h1 in variants:
        _, M, G, D, prior = estimator_matrices(net, h1)
        R = noise_cov(H, V0)
        s = np.linalg.svd(G, compute_uv=False); rankG = int(np.sum(s > 1e-10 * s[0]))
        s2 = np.linalg.svd(np.vstack([G, D]), compute_uv=False); rankGD = int(np.sum(s2 > 1e-10 * s2[0]))
        for solver in ("frame-by-frame", "Kalman+RTS"):
            if solver == "frame-by-frame":
                x, P = gls(G, D, prior, R, Y); dV, cov = bus_voltage_posterior(M, x, P); cov = np.repeat(cov, len(t), 0)
            else:
                x, Ps, nis = kalman_rts(G, D, prior, R, Y); dV, cov = bus_voltage_posterior(M, x, Ps)
            Vh = (V0[:, None] + dV) * rot                               # back to the absolute (GPS) frame
            key = f"{name} | {solver}"; estimates[key] = (Vh, cov)
            results[key] = dict(rank_G=rankG, rank_with_pocket_prior=rankGD, n_states=G.shape[1])
    # baselines
    near = nearest_pmu(net["blocks"])
    Vnom = V0[:, None] * rot
    Vnear = np.array([Zm[PMU.index(near[i])] for i in range(NB)])

    ev = (t >= T0) & (t <= TF + 2.0); norm = ~ev
    def metrics(Vh, cov=None):
        out = {}
        for b in hidden:
            i = b - 1; e = Vh[i] - Vt[i]
            m = dict(rmse_mag=float(np.sqrt(np.mean((abs(Vh[i]) - abs(Vt[i])) ** 2))),
                     rmse_ang_deg=float(np.rad2deg(np.sqrt(np.mean(np.angle(Vh[i] / Vt[i]) ** 2)))),
                     rmse_complex=float(np.sqrt(np.mean(abs(e) ** 2))),
                     rmse_complex_event=float(np.sqrt(np.mean(abs(e[ev]) ** 2))),
                     rmse_complex_normal=float(np.sqrt(np.mean(abs(e[norm]) ** 2))))
            if cov is not None:
                ed = e * np.conj(rot)                                    # error in the estimator frame
                E = np.stack([ed.real, ed.imag], 1); C = cov[:, i] + np.eye(2) * 1e-14
                d2 = np.einsum("ki,kij,kj->k", E, np.linalg.inv(C), E)
                m["coverage95"] = float(np.mean(d2 < 5.991)); m["std_pred"] = float(np.sqrt(np.mean(np.trace(C, axis1=1, axis2=2))))
            out[b] = m
        return out
    allm = {k: metrics(Vh, cov) for k, (Vh, cov) in estimates.items()}
    allm["baseline: nominal PF"] = metrics(Vnom); allm["baseline: nearest PMU"] = metrics(Vnear)

    # branch active-power flows (all 46 branches, from-end) and bus frequencies
    def flows(V): return np.array([V[i] * np.conj(Yl[0, 0] * V[i] + Yl[0, 1] * V[j]) * 100 for i, j, Yl in net["blocks"]]).real
    Pt = flows(Vt)
    ft = freq_from_phase(Vt, t)
    mask = np.ones(len(t), bool)
    for ts in (T0, TF): mask &= abs(t - ts) > 3.0 / FPS
    mask[:2] = mask[-2:] = False
    summary = {}
    for k in allm:
        Vh = estimates[k][0] if k in estimates else (Vnom if "nominal" in k else Vnear)
        rc = np.array([allm[k][b]["rmse_complex"] for b in ident])
        summary[k] = dict(
            median_rmse_complex_identifiable=float(np.median(rc)), max_rmse_complex_identifiable=float(rc.max()),
            median_rmse_mag=float(np.median([allm[k][b]["rmse_mag"] for b in ident])),
            median_rmse_ang_deg=float(np.median([allm[k][b]["rmse_ang_deg"] for b in ident])),
            median_rmse_event_window=float(np.median([allm[k][b]["rmse_complex_event"] for b in ident])),
            pocket_rmse_complex={b: allm[k][b]["rmse_complex"] for b in POCKET},
            branch_P_rmse_MW=float(np.sqrt(np.mean((flows(Vh) - Pt) ** 2))),
            freq_rmse_mHz_hidden=float(1e3 * np.sqrt(np.mean((freq_from_phase(Vh, t)[np.array(hidden) - 1][:, mask] - ft[np.array(hidden) - 1][:, mask]) ** 2))))
        if "coverage95" in allm[k][ident[0]]:
            summary[k]["median_coverage95_identifiable"] = float(np.median([allm[k][b]["coverage95"] for b in ident]))

    # ------------- pre-registered decision rule -------------
    noise_ref = np.sqrt(SIGMA_VMAG ** 2 + (np.deg2rad(SIGMA_VANG_DEG)) ** 2) if ADD_NOISE else 1e-3
    decisions = {}
    for k in estimates:
        s = summary[k]
        beats = np.mean([min(allm[k][b]["rmse_complex"] < allm["baseline: nominal PF"][b]["rmse_complex"],
                             allm[k][b]["rmse_complex"] < allm["baseline: nearest PMU"][b]["rmse_complex"]) for b in ident])
        crit = {"C1 median RMSE (identifiable buses) < 1 x single-PMU noise": s["median_rmse_complex_identifiable"] < 1.0 * noise_ref,
                "C2 beats BOTH baselines on >= 90 % of identifiable buses": beats >= 0.90,
                "C3 95 % coverage in [0.80, 1.00]": 0.80 <= s.get("median_coverage95_identifiable", 0) <= 1.0}
        decisions[k] = dict(criteria=crit, fraction_beating_baselines=float(beats), verdict=all(crit.values()))

    json.dump(dict(meta=meta, gates=net["gate"], estimator=results, summary=summary, decisions=decisions,
                   per_bus=allm), open(os.path.join(OUT, "metrics.json"), "w"), indent=1, default=str)

    # ------------- plots -------------
    fig_dir = os.path.join(OUT, "figs"); os.makedirs(fig_dir, exist_ok=True)
    col = {"truth": "#111111", "H0 | frame-by-frame": "#2a78d6", "H0 | Kalman+RTS": "#1baf7a",
           "H1(event bus known) | frame-by-frame": "#eb6834", "H1(event bus known) | Kalman+RTS": "#e87ba4",
           "baseline: nominal PF": "#9a9993", "baseline: nearest PMU": "#eda100"}
    show = [EVENT_BUS if EVENT_BUS not in PMU else 7, 18, 27, 20]
    for q, lab, fn in [("mag", "|V| [pu]", abs), ("ang", "angle vs bus 2 [deg]", None)]:
        fig, axs = plt.subplots(len(show), 1, figsize=(9, 2.3 * len(show)), sharex=True)
        for ax, b in zip(axs, show):
            i = b - 1
            ref = np.angle(Vt[1]) if q == "ang" else None
            val = (lambda V: abs(V[i])) if q == "mag" else (lambda V: np.rad2deg(np.angle(V[i]) - ref))
            ax.plot(t, val(Vt), color=col["truth"], lw=2.4, label="truth (simulation)")
            for k in ["H0 | Kalman+RTS"] + [k for k in estimates if k.startswith("H1") and "Kalman" in k]:
                Vh, cov = estimates[k]; ax.plot(t, val(Vh), color=col[k], lw=1.2, label=k)
                if q == "mag":
                    sd = np.sqrt(cov[:, i, 0, 0] + cov[:, i, 1, 1]); ax.fill_between(t, abs(Vh[i]) - 2 * sd, abs(Vh[i]) + 2 * sd, color=col[k], alpha=.15, lw=0)
            ax.plot(t, val(Vnom), color=col["baseline: nominal PF"], lw=1, ls="--", label="baseline: nominal PF")
            ax.axvspan(T0, TF, color="#eb6834", alpha=.06)
            ax.set_ylabel(f"bus {b}\n{lab}", fontsize=8); ax.grid(alpha=.3)
        axs[0].legend(fontsize=7, ncol=2, loc="lower left", bbox_to_anchor=(0, 1.02), frameon=False)
        axs[-1].set_xlabel("t [s]"); fig.tight_layout(); fig.savefig(os.path.join(fig_dir, f"1_{q}_hidden_buses.png"), dpi=150); plt.close(fig)
    fig, ax = plt.subplots(figsize=(11, 3.6)); x = np.arange(len(hidden)); w = 0.2
    for n, k in enumerate(["H0 | Kalman+RTS"] + [k for k in estimates if k.startswith("H1") and "Kalman" in k] + ["baseline: nominal PF", "baseline: nearest PMU"]):
        ax.bar(x + (n - 1.5) * w, [allm[k][b]["rmse_complex"] for b in hidden], w, color=col[k], label=k)
    ax.set_yscale("log"); ax.set_xticks(x); ax.set_xticklabels(hidden, fontsize=7); ax.set_xlabel("non-PMU bus")
    ax.set_ylabel("RMSE complex V [pu]"); ax.legend(fontsize=7, frameon=False); ax.grid(alpha=.3, axis="y")
    for b in POCKET: ax.axvspan(hidden.index(b) - .5, hidden.index(b) + .5, color="#9a9993", alpha=.12)
    ax.set_title("Reconstruction error per non-PMU bus (grey = pocket 20/33/34, prior-dependent)", fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(fig_dir, "2_rmse_per_bus.png"), dpi=150); plt.close(fig)
    fig, ax = plt.subplots(figsize=(11, 3))
    for k in [k for k in estimates if "Kalman" in k]:
        ax.plot(hidden, [allm[k][b]["coverage95"] for b in hidden], "o-", color=col[k], label=k, ms=4)
    ax.axhline(.95, color="k", lw=.8, ls=":"); ax.set_ylim(0, 1.02); ax.set_xticks(hidden); ax.tick_params(labelsize=7)
    ax.set_ylabel("95 % band coverage"); ax.set_xlabel("non-PMU bus"); ax.legend(fontsize=7, frameon=False); ax.grid(alpha=.3)
    fig.tight_layout(); fig.savefig(os.path.join(fig_dir, "3_coverage.png"), dpi=150); plt.close(fig)
    kbest = [k for k in estimates if "Kalman" in k][-1]
    fig, axs = plt.subplots(1, 2, figsize=(11, 3.6))
    for b in show:
        axs[0].plot(t[mask], ft[b - 1][mask], lw=2, label=f"bus {b} truth")
        axs[0].plot(t[mask], freq_from_phase(estimates[kbest][0], t)[b - 1][mask], lw=1, ls="--", label=f"bus {b} est.")
    axs[0].set_xlabel("t [s]"); axs[0].set_ylabel("f [Hz]"); axs[0].legend(fontsize=6, ncol=2, frameon=False); axs[0].grid(alpha=.3)
    Ph = flows(estimates[kbest][0])
    axs[1].plot(Pt[:, ev].ravel(), Ph[:, ev].ravel(), ".", ms=1, color=col[kbest])
    lim = [Pt.min(), Pt.max()]; axs[1].plot(lim, lim, "k", lw=.8)
    axs[1].set_xlabel("true branch P [MW]"); axs[1].set_ylabel("estimated branch P [MW]"); axs[1].grid(alpha=.3)
    axs[1].set_title(f"all 46 branches during the event ({kbest})", fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(fig_dir, "4_frequency_and_flows.png"), dpi=150); plt.close(fig)

    # ------------- report -------------
    print("\n" + "=" * 100)
    print(f"EVENT: {EVENT_TYPE} at bus {EVENT_BUS}, scale {EVENT_SCALE:+.2f}, on {T0}-{TF} s, sim to {TEND} s, noise={ADD_NOISE}")
    for k, v in results.items():
        print(f"  {k:36s} rank(G)={v['rank_G']}/{v['n_states']}  rank with pocket prior={v['rank_with_pocket_prior']}")
    print("-" * 100)
    hdr = f"{'method':38s}{'RMSE|V|':>9s}{'RMSEang°':>9s}{'RMSEcplx':>10s}{'max':>8s}{'event win':>10s}{'cov95':>7s}{'P [MW]':>8s}{'f [mHz]':>9s}"
    print("identifiable non-PMU buses (28), medians over buses:"); print(hdr)
    for k, s in summary.items():
        print(f"{k:38s}{s['median_rmse_mag']:9.4f}{s['median_rmse_ang_deg']:9.3f}{s['median_rmse_complex_identifiable']:10.4f}"
              f"{s['max_rmse_complex_identifiable']:8.4f}{s['median_rmse_event_window']:10.4f}{s.get('median_coverage95_identifiable', float('nan')):7.2f}"
              f"{s['branch_P_rmse_MW']:8.1f}{s['freq_rmse_mHz_hidden']:9.2f}")
    print("\npocket buses 20/33/34 (depend on the droop-sharing prior):")
    for k in estimates:
        print(f"  {k:38s}", {b: round(v, 4) for b, v in summary[k]['pocket_rmse_complex'].items()})
    print("\nDECISION (pre-registered rule):")
    for k, d in decisions.items():
        print(f"  {k:38s} -> {'ESTIMATOR: YES' if d['verdict'] else 'ESTIMATOR: NO'}  (beats both baselines on {100*d['fraction_beating_baselines']:.0f} % of buses)")
        for c, ok in d["criteria"].items(): print(f"       [{'x' if ok else ' '}] {c}")
    print(f"\nfigures: {fig_dir}\nmetrics: {os.path.join(OUT, 'metrics.json')}")


if __name__ == "__main__":
    main()
