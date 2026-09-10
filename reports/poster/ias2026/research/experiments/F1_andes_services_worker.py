"""ANDES side of the F1 service ladder. Runs under the tx3-andes interpreter.

For each service configuration and each replacement subset, rebuild the case
without the removed machines (their power kept as a constant-power injection),
optionally without the exciters and stabilizers, and report the inter-area band.

Service configurations, matched to what the internal model can also express:

    A  AVR on, PSS on, no governor
    B  AVR on, PSS off
    C  AVR off (field held at its initial value), PSS off

Removing IEEEX1 from an ANDES case leaves the machine field voltage at its
initialized constant, which is manual excitation, the same ablation the internal
model calls ``avr_manual``.

Usage
    python F1_andes_services_worker.py <output.csv> <workdir>
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
CONFIGURATIONS = {
    "A avr+pss": {"avr": True, "pss": True},
    "B avr only": {"avr": True, "pss": False},
    "C neither": {"avr": False, "pss": False},
}


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


def build_case(sheets, members, pgen, qgen, workdir, tag, *, avr, pss):
    out = {name: frame.copy() for name, frame in sheets.items()}
    for event in ("Toggler", "Toggle"):
        out.pop(event, None)
    out.pop("TGOV1N", None)  # the internal model has no governor

    genrou = out["GENROU"]
    machine_ids = set(genrou[genrou.bus.isin(members)].idx)
    out["GENROU"] = genrou[~genrou.bus.isin(members)]

    avr_ids: set = set()
    if "IEEEX1" in out:
        exciters = out["IEEEX1"]
        avr_ids = set(exciters[exciters.syn.isin(machine_ids)].idx)
        kept = exciters[~exciters.syn.isin(machine_ids)]
        out["IEEEX1"] = kept if avr else kept.iloc[0:0]
    if "IEEEST" in out:
        stabilizers = out["IEEEST"]
        kept = stabilizers[~stabilizers.avr.isin(avr_ids)]
        if not avr:
            kept = kept.iloc[0:0]  # a stabilizer with no exciter has nowhere to act
        out["IEEEST"] = kept if pss else kept.iloc[0:0]

    pv = out["PV"]
    out["PV"] = pv[~pv.bus.isin(members)]
    pq = out["PQ"]
    additions = [
        {
            "uid": len(pq) + offset,
            "idx": f"PQ_GEN_{bus}",
            "u": 1,
            "name": f"PQ_GEN_{bus}",
            "bus": bus,
            "Vn": float(pq.Vn.iloc[0]),
            "p0": -float(pgen[bus]),
            "q0": -float(qgen[bus]),
            "vmax": 1.4,
            "vmin": 0.6,
            "owner": np.nan,
        }
        for offset, bus in enumerate(sorted(members))
    ]
    if additions:
        out["PQ"] = pd.concat([pq, pd.DataFrame(additions)], ignore_index=True)

    target = os.path.join(workdir, f"ieee39_{tag}.xlsx")
    with pd.ExcelWriter(target) as writer:
        for name, frame in out.items():
            if len(frame):
                frame.to_excel(writer, sheet_name=name, index=False)
    return target


def run_case(path):
    andes.config_logger(stream_level=40)
    system = andes.load(path, setup=True, no_output=True)
    config = system.PQ.config
    config.pq2z = 0
    config.p2p, config.q2q = 1.0, 1.0
    config.p2z, config.q2z = 0.0, 0.0
    config.p2i, config.q2i = 0.0, 0.0
    system.PFlow.run()
    if not system.PFlow.converged:
        return {"status": "PF_FAILED"}
    system.EIG.run()
    values = system.EIG.mu
    modes = band_modes(values)
    if not modes:
        return {"status": "NO_BAND_MODE"}
    worst = modes[0]
    return {
        "status": "OK",
        "alpha_band": float(worst.real),
        "freq_band_hz": float(abs(worst.imag) / (2.0 * np.pi)),
        "band_rhp_count": int(sum(1 for m in modes if m.real > 0.0)),
        "band_mode_count": len(modes),
        "spectral_abscissa": float(values.real.max()),
        "rhp_total": int(((values.real > 1e-6) & (values.imag >= 0)).sum()),
        "n_modes": int(len(values)),
        "band_modes": ";".join(
            f"{m.real:+.5f}@{abs(m.imag) / (2 * np.pi):.4f}" for m in modes[:4]
        ),
    }


def main() -> int:
    target, workdir = sys.argv[1], sys.argv[2]
    os.makedirs(workdir, exist_ok=True)
    sheets = pd.read_excel(case_path(), sheet_name=None)
    pgen, qgen = base_generation()

    subsets = [()] + [(b,) for b in CORE] + [tuple(CORE)]
    rows = []
    for name, services in CONFIGURATIONS.items():
        for members in subsets:
            label = "+".join(map(str, members)) or "BASE"
            tag = f"{name.split()[0]}_{label}"
            row = {"configuration": name, "members": label, "size": len(members)}
            try:
                path = build_case(
                    sheets,
                    members,
                    pgen,
                    qgen,
                    workdir,
                    tag,
                    avr=services["avr"],
                    pss=services["pss"],
                )
                row.update(run_case(path))
            except Exception as error:  # noqa: BLE001 - recorded per case
                row["status"] = f"ERROR: {str(error)[:110]}"
            rows.append(row)
            print(
                "  andes %-11s %-12s %-10s %s"
                % (
                    name,
                    label,
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
    print(f"andes service ladder -> {target}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
