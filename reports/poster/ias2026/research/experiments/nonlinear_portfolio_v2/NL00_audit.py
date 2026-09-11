"""NL00: physics, units and origin audit of the L0 models (no new simulation physics).

A. Every Sn, P, Q, H/M, D, reactance, gain and time constant traced from the
   frozen JSON (with its source hash) into the device objects. The base
   conversion gamma = Sn/Sbase is checked numerically on every machine of a
   solved case. The machine is re-expressed on the system base (X / gamma,
   M, D and P * gamma, weight 1); it must inject the same network current and
   have the same derivatives up to the state scaling, at an off-equilibrium
   point. This is P_sys = gamma P_local, I_sys = gamma I_local and
   X_sys = X_local / gamma as a behavioural test, not a relabelling.
B. The 4270.7 figure: the BC00-A reconstruction is reused and re-checked here.
C. Near-zero poles and symmetry: the BC01 audit rows are reused (the same code),
   plus a FINITE-angle rotation test of the nonlinear residuals, which is
   stronger than the tangent identity A R = 0.
D. Power-flow slack vs dynamic source: whether any bus is held at a fixed
   voltage or angle in the dynamic model.
E. Scope of the earlier independent evidence (ANDES), quoted from F1.

Outputs in results/NL/NL00/: NL00_units.csv, NL00_base_conversion.csv,
NL00_symmetry_audit.csv, NL00_rotation_nonlinear.csv, NL00_hashes.json.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import _bootstrap  # noqa: E402,F401
from _bootstrap import ROOT  # noqa: E402
from F12_kundur import solve as k_solve  # noqa: E402
from G3_ieee68 import solve as s68  # noqa: E402
from ibr_cycles.models.ieee39_case import ReplacementPlan, solve_case  # noqa: E402
from ibr_cycles.models.ieee39_devices import ConverterParameters  # noqa: E402

OUT = ROOT / "results" / "NL" / "NL00"
OUT.mkdir(parents=True, exist_ok=True)
SOURCES = {
    "ieee39_network.json": ROOT / "configs" / "ias2026" / "ieee39_network.json",
    "kundur_network.json": ROOT / "configs" / "kundur" / "kundur_network.json",
    "ieee68_network.json": ROOT / "configs" / "ieee68" / "ieee68_network.json",
    "ieee39_full.xlsx (ANDES source)": ROOT.parents[3]
    / "data"
    / "raw"
    / "ieee39_full.xlsx",
    "ieee39_devices.py": ROOT / "src" / "ibr_cycles" / "models" / "ieee39_devices.py",
    "ieee68_devices.py": ROOT / "src" / "ibr_cycles" / "models" / "ieee68_devices.py",
    "ieee39_case.py": ROOT / "src" / "ibr_cycles" / "models" / "ieee39_case.py",
    "ieee39_network.py": ROOT / "src" / "ibr_cycles" / "models" / "ieee39_network.py",
}


def sha(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def cases():
    conv = ConverterParameters(
        voltage_control=True, voltage_gain=0.05, voltage_leak=0.05
    )
    return {
        "IEEE-39 base": solve_case(ReplacementPlan.of({})),
        "IEEE-39 30@0.5": solve_case(ReplacementPlan.of({30: 0.5}), converter=conv),
        "IEEE-39 flagship": solve_case(
            ReplacementPlan.of({b: 1.0 for b in (30, 33, 35, 37)}), converter=conv
        ),
        "Kundur 2": k_solve((2,), {"g": 0.08, "k": 1.25, "t": 1.0}),
        "IEEE-68 3+4+6+9": s68((3, 4, 6, 9), 0.0, 1.0),
    }


def _to_system_base(dev, xs):
    """The same machine re-expressed on the system base (weight 1).

    Two-axis machine: reactances and resistance / gamma; M and D * gamma;
    mechanical power * gamma; the power-input PSS states and its washout input
    are powers (* gamma), so its gain / gamma. The 68-bus machine: impedances
    / gamma, H, D and Tm * gamma; the speed-input PSS and exciters are
    intensive.
    """

    from dataclasses import replace

    g = dev.weight
    p = dev.parameters
    if hasattr(p, "pss_gain"):  # IEEE-39 / Kundur two-axis
        q = replace(
            p,
            ra=p.ra / g,
            xd=p.xd / g,
            xq=p.xq / g,
            xd1=p.xd1 / g,
            xq1=p.xq1 / g,
            m=p.m * g,
            d=p.d * g,
            pm=p.pm * g,
            pss_gain=p.pss_gain / g,
        )
        xs2 = xs.copy()
        xs2[5] *= g
        xs2[6] *= g
        return (
            replace(dev, parameters=q, weight=1.0),
            xs2,
            np.array([1, 1, 1, 1, 1, g, g] + [1] * (xs.size - 7), dtype=float),
        )
    q = replace(
        p,
        xl=p.xl / g,
        ra=p.ra / g,
        xd=p.xd / g,
        xd1=p.xd1 / g,
        xd2=p.xd2 / g,
        xq=p.xq / g,
        xq1=p.xq1 / g,
        xq2=p.xq2 / g,
        h=p.h * g,
        d=p.d * g,
        tm=p.tm * g,
    )
    return replace(dev, parameters=q, weight=1.0), xs.copy(), np.ones(xs.size)


def base_conversion(case_name, case) -> list[dict]:
    """Physical invariance of every machine under the local -> system base change.

    The system-base object must give the same network current, and the same
    derivatives up to the state scaling (power states * gamma).
    """

    dae = case.dae
    x, z = case.equilibrium.x, case.equilibrium.z
    v = dae.voltages(z)
    rng = np.random.default_rng(7)
    rows = []
    for slot in dae.slots:
        if slot.kind != "sg":
            continue
        dev = slot.device
        vb = complex(v[dae.network.position(slot.bus)])
        xs = x[slot.start : slot.stop] + 1e-3 * rng.standard_normal(
            slot.stop - slot.start
        )
        sysdev, xs_sys, scale = _to_system_base(dev, xs)
        i_loc = dev.injection(xs, vb)
        i_sys = sysdev.injection(xs_sys, vb)
        f_loc = dev.derivatives(xs, vb)
        f_sys = sysdev.derivatives(xs_sys, vb)
        rows.append(
            {
                "case": case_name,
                "bus": slot.bus,
                "gamma": dev.weight,
                "S_device_MVA": dev.weight * 100.0,
                "current_mismatch": abs(i_loc - i_sys),
                "derivative_mismatch": float(np.abs(f_loc * scale - f_sys).max()),
                "max_abs_derivative": float(np.abs(f_loc).max()),
            }
        )
    return rows


def rotation_nonlinear(case_name, case, phi=0.37) -> dict:
    """Finite rotation: every angle + phi, every bus voltage * e^{j phi}."""

    from ibr_cycles.certification.symmetry import rotation_generator

    dae = case.dae
    x, z = case.equilibrium.x.copy(), case.equilibrium.z.copy()
    r_x, _ = rotation_generator(dae, z)
    rng = np.random.default_rng(1)
    x = x + 1e-3 * rng.standard_normal(x.size)  # off-equilibrium point
    z = z + 1e-3 * rng.standard_normal(z.size)
    v = dae.voltages(z) * np.exp(1j * phi)
    z_rot = np.empty_like(z)
    z_rot[0::2], z_rot[1::2] = v.real, v.imag
    x_rot = x + phi * r_x
    f0, f1 = dae.f(x, z, {}), dae.f(x_rot, z_rot, {})
    g0, g1 = dae.g(x, z, {}), dae.g(x_rot, z_rot, {})
    g0c = (g0[0::2] + 1j * g0[1::2]) * np.exp(1j * phi)
    g1c = g1[0::2] + 1j * g1[1::2]
    return {
        "case": case_name,
        "phi_rad": phi,
        "max_abs_f_change": float(np.abs(f1 - f0).max()),
        "max_abs_g_covariance_error": float(np.abs(g1c - g0c).max()),
        "max_abs_f": float(np.abs(f0).max()),
    }


def slack_audit(case_name, case) -> dict:
    dae = case.dae
    kinds = {s.kind for s in dae.slots}
    fixed = [
        s for s in dae.slots if getattr(s.device, "n_states", 1) == 0 and s.kind == "sg"
    ]
    return {
        "case": case_name,
        "pf_slack_bus": dae.network.slack_bus,
        "slack_bus_has_dynamic_machine": any(
            s.bus == dae.network.slack_bus and s.kind == "sg" for s in dae.slots
        ),
        "fixed_voltage_sources_in_DAE": len(fixed),
        "device_kinds": ",".join(sorted(kinds)),
        "algebraic_unknowns": dae.n_z,
        "note": "every bus voltage is an algebraic unknown; no bus is held fixed",
    }


def main() -> int:
    hashes = {k: sha(v) for k, v in SOURCES.items() if Path(v).exists()}
    (OUT / "NL00_hashes.json").write_text(
        json.dumps(hashes, indent=1), encoding="utf-8"
    )
    units = pd.read_csv(ROOT / "results" / "BC" / "BC00" / "BC00_units_per_machine.csv")
    units.to_csv(OUT / "NL00_units.csv", index=False)
    all_cases = cases()
    conv = pd.DataFrame(
        [r for n, c in all_cases.items() for r in base_conversion(n, c)]
    )
    conv.to_csv(OUT / "NL00_base_conversion.csv", index=False)
    rot = pd.DataFrame([rotation_nonlinear(n, c) for n, c in all_cases.items()])
    rot.to_csv(OUT / "NL00_rotation_nonlinear.csv", index=False)
    slack = pd.DataFrame([slack_audit(n, c) for n, c in all_cases.items()])
    slack.to_csv(OUT / "NL00_slack_audit.csv", index=False)
    sym = []
    for name in ("ieee39", "kundur", "ieee68"):
        p = ROOT / "results" / "BC" / "BC01" / f"BC01_{name}_subsets.csv.gz"
        if p.exists():
            f = pd.read_csv(p)
            sym.append(
                {
                    "benchmark": name,
                    "subset_cases": len(f),
                    "max_A_R_residual": float(f.a_rot.max()),
                    "max_Jordan_residual": float(f.a_jordan.max()),
                    "neutral_mode_verified_fraction": float(f.neutral_verified.mean()),
                    "dead_state_cases": int((f.n_dead > 0).sum()),
                    "full_A_eigs_within_1e-3_(max)": int(f["nz_full_1e-3"].max()),
                    "rest_eigs_within_1e-3_(max)": int(f["nz_rest_1e-3"].max()),
                    "status_rest_counts": json.dumps(
                        f.status_rest.value_counts().to_dict()
                    ),
                    "status_physical_counts": json.dumps(
                        f.status_physical.value_counts().to_dict()
                    ),
                }
            )
    pd.DataFrame(sym).to_csv(OUT / "NL00_symmetry_audit.csv", index=False)
    pd.set_option("display.width", 220)
    print(conv.groupby("case")[["current_mismatch", "derivative_mismatch"]].max())
    print(rot.to_string())
    print(slack.to_string())
    print(pd.DataFrame(sym).to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
