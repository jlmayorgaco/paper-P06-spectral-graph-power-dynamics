"""F1 R5: an equation-equivalent ANDES case. Runs under the tx3-andes interpreter.

The internal model is a fourth-order machine, a first-order AVR
``efd' = (K (vref - |V|) - efd) / T``, an optional washout-plus-lag stabilizer on
electrical power, no governor, constant-power loads.

ANDES can express that exactly except for the machine order:

  * **SEXS with TA/TB = 1** collapses its lead-lag to unity and leaves
    ``K / (1 + s TE)``, which is the internal AVR term for term. Limits are
    opened to +-99 because the internal AVR has none.
  * **GENROU with xd2 = xd1, xq2 = xq1 and very small subtransient constants**
    makes the damper windings fast, which is the closest ANDES gets to a
    fourth-order machine without a new model card. The residual is quantified by
    running two subtransient constants a decade apart.
  * **no exciter at all** is the manual-excitation configuration, which both
    tools express exactly and which is therefore the cleanest equivalence point.

Ladder stages, all with no governor and constant-power loads:

    R2  manual excitation, no stabilizer          <- exact on both sides
    R3  first-order AVR via SEXS, no stabilizer   <- exact on both sides
    R3b same, with the machine left at sixth order

Usage
    python F1_andes_equivalent_worker.py <output.csv> <workdir>
"""

from __future__ import annotations

import os
import sys
import warnings

warnings.filterwarnings("ignore")

import andes  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

CORE = (30, 33, 35, 37)
BAND = (0.3, 1.5)
ZERO = 1e-3
SUBSETS = [()] + [(b,) for b in CORE] + [tuple(CORE)]
STAGES = {
    # name: (use SEXS avr, degrade GENROU to fourth order, subtransient constant)
    "R2 manual excitation": (False, True, 1e-4),
    "R3 first-order AVR": (True, True, 1e-4),
}
#: Full spectra are dumped so the reconciliation is not limited to one band mode.
SPECTRA: dict = {}


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


def build_case(sheets, members, pgen, qgen, workdir, tag, *, avr):
    out = {name: frame.copy() for name, frame in sheets.items()}
    for drop in ("Toggler", "Toggle", "TGOV1N", "IEEEST", "ACEc"):
        out.pop(drop, None)

    genrou = out["GENROU"]
    machine_ids = set(genrou[genrou.bus.isin(members)].idx)
    out["GENROU"] = genrou[~genrou.bus.isin(members)]

    exciters = out.pop("IEEEX1")
    if avr:
        kept = exciters[~exciters.syn.isin(machine_ids)].reset_index(drop=True)
        out["SEXS"] = pd.DataFrame(
            {
                "uid": range(len(kept)),
                "idx": [f"SEXS_{i + 1}" for i in range(len(kept))],
                "u": 1,
                "name": [f"SEXS_{i + 1}" for i in range(len(kept))],
                "syn": kept.syn.to_numpy(),
                "TATB": 1.0,   # collapses the lead-lag to unity
                "TB": 1.0,
                "K": kept.KA.to_numpy(),   # the internal AVR gain
                "TE": kept.TE.to_numpy(),  # the internal AVR time constant
                "EMIN": -99.0,
                "EMAX": 99.0,
            }
        )

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

    target = os.path.join(workdir, f"eq_{tag}.xlsx")
    with pd.ExcelWriter(target) as writer:
        for name, frame in out.items():
            if len(frame):
                frame.to_excel(writer, sheet_name=name, index=False)
    return target


def run_case(path, *, fourth_order, tau):
    andes.config_logger(stream_level=40)
    system = andes.load(path, setup=True, no_output=True)
    config = system.PQ.config
    config.pq2z = 0
    config.p2p, config.q2q = 1.0, 1.0
    config.p2z, config.q2z = 0.0, 0.0
    config.p2i, config.q2i = 0.0, 0.0
    if fourth_order:
        system.GENROU.xd2.v[:] = system.GENROU.xd1.v[:]
        system.GENROU.xq2.v[:] = system.GENROU.xq1.v[:]
        system.GENROU.Td20.v[:] = tau
        system.GENROU.Tq20.v[:] = tau
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
        "all_eigenvalues": values.copy(),
        "alpha_band": float(worst.real),
        "freq_band_hz": float(abs(worst.imag) / (2.0 * np.pi)),
        "band_rhp_count": int(sum(1 for m in modes if m.real > 0.0)),
        "spectral_abscissa": float(values.real.max()),
        "rhp_total": int(((values.real > 1e-6) & (values.imag >= 0)).sum()),
        "n_modes": int(len(values)),
    }


def main() -> int:
    target, workdir = sys.argv[1], sys.argv[2]
    os.makedirs(workdir, exist_ok=True)
    sheets = pd.read_excel(case_path(), sheet_name=None)
    pgen, qgen = base_generation()

    rows = []
    for stage, (avr, fourth, tau) in STAGES.items():
        for members in SUBSETS:
            label = "+".join(map(str, members)) or "BASE"
            row = {"stage": stage, "members": label, "size": len(members)}
            try:
                path = build_case(
                    sheets, members, pgen, qgen, workdir,
                    f"{stage.split()[0]}_{label}", avr=avr,
                )
                row.update(run_case(path, fourth_order=fourth, tau=tau))
            except Exception as error:  # noqa: BLE001
                row["status"] = f"ERROR: {str(error)[:110]}"
            spectrum = row.pop("all_eigenvalues", None)
            if spectrum is not None:
                SPECTRA[f"{stage}|{label}"] = spectrum
            rows.append(row)
            print(
                "  andes %-28s %-12s %-10s %s"
                % (
                    stage, label, row.get("status"),
                    (
                        "a=%+.5f f=%.4f rhp=%d"
                        % (row.get("alpha_band", float("nan")),
                           row.get("freq_band_hz", float("nan")),
                           row.get("band_rhp_count", -1))
                        if row.get("status") == "OK" else ""
                    ),
                ),
                flush=True,
            )
    pd.DataFrame(rows).to_csv(target, index=False)
    np.savez_compressed(
        os.path.join(os.path.dirname(target), "F1_andes_spectra.npz"), **SPECTRA
    )
    print(f"equation-equivalent ladder -> {target}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
