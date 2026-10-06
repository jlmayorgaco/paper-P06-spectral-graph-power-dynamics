# ruff: noqa: E501
"""H12 modal mixing at fixed frequencies (prereg H12/H13): mu(theta, omega) for H4 at every
policy, and mu(theta_ref, omega_c(theta)) at theta_ref = P4 (D01). Old E12 construction."""

from __future__ import annotations

import _hinfra as HI
import numpy as np

import _cdw as C
import _infra as I
import E11_spectral as E11

PHASE = "H_H12"
PHASE_B = "H_H12b"
FIXED_HZ = (0.60, 0.70, 0.80)
ALL_POL = HI.DISC + HI.OLD_HOLD + HI.HPOL


def _T(jac, s):
    n = jac.fx.shape[0]
    return jac.gz + jac.gx @ np.linalg.solve(s * np.eye(n) - jac.fx, jac.fz)


def mu_at(case, hz_list):
    _, _, jac = C.physical_matrices(case)
    _, U, _ = E11.graph_basis()
    return [E11.mixing(_T(jac, 1j * 2 * np.pi * hz), U) for hz in hz_list]


def tasks():
    return [{"phase": PHASE, "pid": p} for p in ALL_POL]


def run_task(task):
    theta = HI.theta_of(task["pid"])
    case = C.solve(C.V4, theta)
    ev = C.transverse_eval(case, modes=False)
    wc = ev["lam_hz"] if ev["lam_hz"] > 1e-3 else 0.7  # old E12 rule
    hz = [HI.OMEGA_REF_HZ, *FIXED_HZ, wc]
    mus = mu_at(case, hz)
    return {"status": ev["status"], "alpha": ev["alpha"], "lam_hz": ev["lam_hz"], "wc_hz": wc,
            "mu_ref": mus[0], "mu_060": mus[1], "mu_070": mus[2], "mu_080": mus[3], "mu_wc": mus[4]}


def tasks_b():
    wcs = {}
    for r in I.Store(PHASE).all():
        if r.get("ok"):
            wcs[r["task"]["pid"]] = r["wc_hz"]
    return [{"phase": PHASE_B, "pids": sorted(wcs), "wc": [wcs[p] for p in sorted(wcs)]}]


def run_task_b(task):
    case = C.solve(C.V4, C.policy("D01"))
    mus = mu_at(case, task["wc"])
    return {"mu_ref_at_wc": dict(zip(task["pids"], mus, strict=True)), "mu00": mu_at(case, [HI.OMEGA_REF_HZ])[0]}
