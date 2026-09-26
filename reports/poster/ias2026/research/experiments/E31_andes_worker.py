"""ANDES-side worker for E31. Runs under the tx3-andes interpreter only.

ANDES is used as an INDEPENDENT implementation: a different machine model
(GENROU, six states with damper windings, against the four-state model used
here), a different exciter, a different stabilizer, a different network assembly
and a different eigenvalue path. Nothing in this file imports project code.

Removing a machine is done by REBUILDING the case without it, not by zeroing a
switch: zeroing leaves ANDES trying to initialize a device that no longer has a
consistent operating point, and it fails with NaNs. The bus keeps its injection
as a constant-power PQ term at the base-case generation, which is the same
negative control the project calls ``static_power``: the synchronous dynamics are
gone and the power is still delivered.

What ANDES cannot represent here is the project's grid-following converter. ANDES
ships REGCA1 and REECA1, whose equations are structurally different, so using
them would compare two different converters rather than validate one. Those cases
are reported NOT COMPARABLE instead of compared badly.

Usage
    python E31_andes_worker.py <output.csv> <workdir>
"""

from __future__ import annotations

import os
import sys
import warnings
from itertools import combinations

warnings.filterwarnings("ignore")

import andes  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

CORE = (30, 33, 35, 37)
BAND = (0.3, 1.5)
ZERO = 1e-3
SHEETS_TO_STRIP = ("GENROU", "IEEEX1", "IEEEST", "TGOV1N")


def case_path() -> str:
    return os.path.join(
        os.path.dirname(andes.__file__), "cases", "ieee39", "ieee39_full.xlsx"
    )


def band_modes(values):
    frequencies = np.abs(values.imag) / (2.0 * np.pi)
    keep = (
        (frequencies >= BAND[0])
        & (frequencies <= BAND[1])
        & (values.imag >= 0.0)
        & (np.abs(values) > ZERO)
    )
    return sorted(values[keep], key=lambda z: -z.real)


def base_generation():
    andes.config_logger(stream_level=40)
    system = andes.load(case_path(), setup=True, no_output=True)
    system.PFlow.run()
    return (
        dict(zip(system.PV.bus.v, system.PV.p.v)),
        dict(zip(system.PV.bus.v, system.PV.q.v)),
    )


def build_case(sheets, members, pgen, qgen, workdir, tag):
    """Case with the machines at ``members`` gone and their power kept as PQ."""

    out = {name: frame.copy() for name, frame in sheets.items()}
    # Toggle is a time-domain disturbance that points at a specific machine. It
    # plays no part in an eigenvalue run and it blocks setup once that machine is
    # gone, so it is dropped from every case including the base, which keeps the
    # comparison like for like.
    for _event in ("Toggler", "Toggle"):
        out.pop(_event, None)
    genrou = out["GENROU"]
    doomed = genrou[genrou.bus.isin(members)]
    machine_ids = set(doomed.idx)
    out["GENROU"] = genrou[~genrou.bus.isin(members)]

    avr_ids: set = set()
    if "IEEEX1" in out:
        exciters = out["IEEEX1"]
        avr_ids = set(exciters[exciters.syn.isin(machine_ids)].idx)
        out["IEEEX1"] = exciters[~exciters.syn.isin(machine_ids)]
    if "IEEEST" in out:
        stabilizers = out["IEEEST"]
        out["IEEEST"] = stabilizers[~stabilizers.avr.isin(avr_ids)]
    if "TGOV1N" in out:
        governors = out["TGOV1N"]
        out["TGOV1N"] = governors[~governors.syn.isin(machine_ids)]

    pv = out["PV"]
    out["PV"] = pv[~pv.bus.isin(members)]
    pq = out["PQ"]
    additions = []
    for offset, bus in enumerate(sorted(members)):
        additions.append(
            {
                "uid": len(pq) + offset,
                "idx": f"PQ_GEN_{bus}",
                "u": 1,
                "name": f"PQ_GEN_{bus}",
                "bus": bus,
                "Vn": float(pq.Vn.iloc[0]),
                "p0": -float(pgen[bus]),  # a generator is a negative load
                "q0": -float(qgen[bus]),
                "vmax": 1.4,
                "vmin": 0.6,
                "owner": np.nan,
            }
        )
    out["PQ"] = pd.concat([pq, pd.DataFrame(additions)], ignore_index=True)

    target = os.path.join(workdir, f"ieee39_{tag}.xlsx")
    with pd.ExcelWriter(target) as writer:
        for name, frame in out.items():
            if len(frame):
                frame.to_excel(writer, sheet_name=name, index=False)
    return target


def run_case(path, *, governors: bool):
    andes.config_logger(stream_level=40)
    system = andes.load(path, setup=True, no_output=True)
    # constant power, not constant impedance: the injection must not soften
    for key, value in (("p2p", 1), ("q2q", 1), ("p2z", 0), ("q2z", 0)):
        if key in system.PQ.config.__dict__:
            system.PQ.config.__dict__[key] = value
    if not governors and getattr(system, "TGOV1N", None) is not None:
        if system.TGOV1N.n:
            system.TGOV1N.u.v[:] = 0.0
    system.PFlow.run()
    if not system.PFlow.converged:
        return {"status": "PF_FAILED"}
    system.EIG.run()
    values = system.EIG.mu
    modes = band_modes(values)
    if not modes:
        return {"status": "NO_BAND_MODE"}
    worst = modes[0]
    voltages = np.array(system.Bus.v.v, dtype=float)
    return {
        "status": "OK",
        "alpha_band": float(worst.real),
        "freq_band_hz": float(abs(worst.imag) / (2.0 * np.pi)),
        "band_rhp_count": int(sum(1 for m in modes if m.real > 0.0)),
        "band_mode_count": len(modes),
        "spectral_abscissa": float(values.real.max()),
        "rhp_total": int(((values.real > 1e-6) & (values.imag >= 0)).sum()),
        "min_voltage": float(voltages.min()),
        "max_voltage": float(voltages.max()),
        "n_modes": int(len(values)),
        "band_modes": ";".join(
            f"{m.real:+.5f}@{abs(m.imag) / (2 * np.pi):.4f}" for m in modes[:4]
        ),
    }


def main() -> int:
    target, workdir = sys.argv[1], sys.argv[2]
    os.makedirs(workdir, exist_ok=True)
    source = case_path()
    sheets = pd.read_excel(source, sheet_name=None)
    pgen, qgen = base_generation()

    rows = []
    for size in range(len(CORE) + 1):
        for members in combinations(CORE, size):
            label = "+".join(map(str, members)) or "BASE"
            path = build_case(sheets, members, pgen, qgen, workdir, label)
            for governors in (True, False):
                row = {"governors": governors, "members": label, "size": size}
                try:
                    row.update(run_case(path, governors=governors))
                except Exception as error:  # noqa: BLE001 - recorded per case
                    row.update({"status": f"ERROR: {str(error)[:110]}"})
                rows.append(row)
                print(
                    "  andes %-12s gov=%-5s %-10s %s"
                    % (
                        label,
                        governors,
                        row.get("status"),
                        (
                            "a=%+.5f f=%.4f rhp=%d"
                            % (
                                row.get("alpha_band", float("nan")),
                                row.get("freq_band_hz", float("nan")),
                                row.get("band_rhp_count", -1),
                            )
                            if row.get("status") == "OK"
                            else ""
                        ),
                    ),
                    flush=True,
                )
    pd.DataFrame(rows).to_csv(target, index=False)
    print(f"andes results -> {target}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
