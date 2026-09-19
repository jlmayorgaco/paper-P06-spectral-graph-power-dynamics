"""P1 independent oracle for the frozen GFL11 equations.

This file intentionally contains no imports from the historical tx3 model.  It
is a small, explicit transcription of the preregistered GFL11 equations used
to compare against the Julia/PowerDynamics implementation.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np


CAMPAIGN = Path(__file__).resolve().parents[2]
RAW = CAMPAIGN / "raw" / "gfl11"
RAW.mkdir(parents=True, exist_ok=True)


PARAMS = {
    "kp_pll": 53.0,
    "ki_pll": 1400.0,
    "tau_p": 0.03,
    "kp_p": 0.20,
    "ki_p": 8.0,
    "kp_q": 0.20,
    "ki_q": 8.0,
    "kp_i": 0.25,
    "ki_i": 6.0,
    "xf": 0.15,
    "rf": 0.01,
    "kp_v": 2.0,
    "ki_v": 20.0,
    "leak": 0.05,
    "omega_b": 2.0 * math.pi * 60.0,
}


def initialize(v: float, a: float, p0: float, q0: float, w: float, g: float) -> tuple[float, ...]:
    terminal = v * np.exp(1j * a)
    current = np.conj((p0 + 1j * q0) / w / terminal)
    local = current * np.exp(-1j * a)
    i_d, i_q = float(local.real), float(local.imag)
    p_ref = v * i_d
    q_ref = -v * i_q
    return (a, 0.0, p_ref, q_ref, i_d, -i_q, i_d, i_q, PARAMS["rf"] * i_d, PARAMS["rf"] * i_q, q_ref)


def evaluate(x: np.ndarray, v: float, a: float, p_ref: float, q_ref: float, v_ref: float, g: float, w: float) -> tuple[np.ndarray, float, float, float, float]:
    theta, x_pll, p_f, q_f, x_p, x_q, i_d, i_q, x_id, x_iq, x_v = x
    v_d = v * math.cos(a - theta)
    v_q = v * math.sin(a - theta)
    power = v_d * i_d + v_q * i_q
    reactive = v_q * i_d - v_d * i_q
    err = g * (v_ref - v)
    q_cmd = PARAMS["kp_v"] * err + x_v
    id_ref = PARAMS["kp_p"] * (p_ref - p_f) + x_p
    iq_ref = -(PARAMS["kp_q"] * (q_cmd - q_f) + x_q)
    e_d = v_d + PARAMS["kp_i"] * (id_ref - i_d) + x_id - PARAMS["xf"] * i_q
    e_q = v_q + PARAMS["kp_i"] * (iq_ref - i_q) + x_iq + PARAMS["xf"] * i_d
    f = np.array(
        [
            PARAMS["kp_pll"] * v_q + x_pll,
            PARAMS["ki_pll"] * v_q,
            (power - p_f) / PARAMS["tau_p"],
            (reactive - q_f) / PARAMS["tau_p"],
            PARAMS["ki_p"] * (p_ref - p_f),
            PARAMS["ki_q"] * (q_cmd - q_f),
            (PARAMS["omega_b"] / PARAMS["xf"])
            * (e_d - v_d - PARAMS["rf"] * i_d + PARAMS["xf"] * i_q),
            (PARAMS["omega_b"] / PARAMS["xf"])
            * (e_q - v_q - PARAMS["rf"] * i_q - PARAMS["xf"] * i_d),
            PARAMS["ki_i"] * (id_ref - i_d),
            PARAMS["ki_i"] * (iq_ref - i_q),
            PARAMS["ki_v"] * err - PARAMS["leak"] * (x_v - q_ref),
        ],
        dtype=float,
    )
    current = w * (i_d + 1j * i_q) * np.exp(1j * theta)
    return f, float(current.real), float(current.imag), float(w * power), float(w * reactive)


def main() -> int:
    rng = np.random.default_rng(20260919)
    rows: list[dict[str, float | int | str]] = []
    for case in range(64):
        v = float(rng.uniform(0.90, 1.10))
        a = float(rng.uniform(-math.pi, math.pi))
        p0 = float(rng.uniform(0.10, 0.80))
        q0 = float(rng.uniform(-0.25, 0.25))
        w = float(rng.uniform(0.25, 1.25))
        g = float(rng.choice([0.03625, 0.2076814, 0.25]))
        x_eq = np.array(initialize(v, a, p0, q0, w, g), dtype=float)
        x = x_eq + 0.05 * np.maximum(1.0, np.abs(x_eq)) * rng.normal(size=11)
        p_ref, q_ref, v_ref = x_eq[2], x_eq[3], v
        f, ir, ii, power, reactive = evaluate(x, v, a, p_ref, q_ref, v_ref, g, w)
        row: dict[str, float | int | str] = {
            "case": case,
            "g": g,
            "w": w,
            "v": v,
            "a": a,
            "p_ref": p_ref,
            "q_ref": q_ref,
            "v_ref": v_ref,
        }
        row.update({f"x{i}": float(value) for i, value in enumerate(x)})
        row.update({f"f{i}": float(value) for i, value in enumerate(f)})
        row.update({"i_r": ir, "i_i": ii, "P_system": power, "Q_system": reactive})
        rows.append(row)

    out = RAW / "p1_python_oracle.csv"
    fields = list(rows[0])
    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    (RAW / "p1_python_environment.txt").write_text(
        "python_oracle=run_p1_gfl11.py\n"
        f"numpy={np.__version__}\n"
        "rng=PCG64 seed 20260919\n"
        "equations=20260911_GFL_REPRODUCTION_SPEC.md section 4.2\n",
        encoding="utf-8",
    )
    print(f"P1_PYTHON_ORACLE_PASS cases={len(rows)} path={out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
