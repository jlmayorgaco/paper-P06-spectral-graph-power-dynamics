"""Post-freeze SciPy SLSQP falsification on the exact frozen Julia closure."""
from __future__ import annotations

import csv
import math
import os
import random
import shutil
import socket
import subprocess
import time
import tomllib
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "reports" / "experiment_N"


def main() -> None:
    candidate_path = OUT / "Z_N_NOMINAL_FINAL.toml"
    candidate = tomllib.loads(candidate_path.read_text(encoding="utf-8"))
    with (OUT / "TABLE_N01_original_operating_point.csv").open(newline="") as fh:
        original = list(csv.DictReader(fh))
    weight = float(original[8]["P_gen_MW"])
    incumbent = float(candidate["retained_SG_MW"])
    julia = shutil.which("julia")
    if julia is None:
        raise SystemExit("Julia executable not found")
    flags = subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0
    proc = subprocess.Popen([julia, "--project=.",
                             "experiments/bnd_expN/falsification_server.jl"],
                            cwd=ROOT, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True,
                            creationflags=flags, bufsize=1)
    try:
        assert proc.stdout is not None
        while True:
            line = proc.stdout.readline()
            if not line:
                raise RuntimeError("Julia closure server exited before ready")
            if "FALSIFICATION_SERVER_READY" in line:
                break
        with socket.create_connection(("127.0.0.1", 32761), timeout=15) as conn:
            handle = conn.makefile("rw", encoding="utf-8", newline="\n")
            cache_key: bytes | None = None
            cache_value: np.ndarray | None = None
            evaluations = 0

            def model(x: np.ndarray) -> np.ndarray:
                nonlocal cache_key, cache_value, evaluations
                x = np.asarray(x, dtype=float)
                key = x.tobytes()
                if key == cache_key and cache_value is not None:
                    return cache_value
                handle.write(",".join(format(float(v), ".17g") for v in x) + "\n")
                handle.flush()
                response = handle.readline().strip()
                if response.startswith("ERR:"):
                    raise RuntimeError(response)
                values = np.fromstring(response, sep=",")
                if len(values) != 22 or not np.all(np.isfinite(values)):
                    raise RuntimeError("invalid Julia closure response: " + response)
                cache_key, cache_value = key, values
                evaluations += 1
                return values

            k0p = 2 * math.pi * 5
            k0i = k0p * k0p / 4
            bounds = [(1e-7, 1.0)] + [(0.25*k0p, 4*k0p)]*10 + \
                     [(0.25*k0i, 4*k0i)]*10
            c0 = np.array([1-candidate["rho"][8]] + candidate["Kp"] +
                          candidate["Ki"], dtype=float)
            with (OUT / "validation_cases" / "ExpE_corrected" /
                  "case_definition.csv").open(newline="") as fh:
                e = list(csv.DictReader(fh))
            e0 = np.array([1-float(e[8]["rho"])] +
                          [float(r["Kp"]) for r in e] +
                          [float(r["Ki"]) for r in e])
            rng = random.Random(20260930)
            seeds = [("frozen_candidate", c0), ("ExpE_retrimmed", e0)]
            for j in range(1, 7):
                eps = rng.uniform(0.001, 0.015)
                kp = [rng.uniform(*bounds[k]) for k in range(1, 11)]
                ki = [rng.uniform(*bounds[k]) for k in range(11, 21)]
                seeds.append((f"deterministic_random_{j}",
                              np.array([eps] + kp + ki)))

            constraints = [{"type": "ineq",
                            "fun": lambda x: -0.05-model(x)[0],
                            "jac": lambda x: -model(x)[1:]}]
            results: list[dict] = []
            for name, seed in seeds:
                start = time.monotonic()
                try:
                    result = minimize(lambda x: weight*x[0], seed,
                                      jac=lambda x: np.r_[weight,
                                                           np.zeros(20)],
                                      bounds=bounds, constraints=constraints,
                                      method="SLSQP",
                                      options={"maxiter": 100, "ftol": 1e-11,
                                               "disp": False})
                    alpha = float(model(result.x)[0])
                    cost = float(weight*result.x[0])
                    feasible = alpha <= -0.05 and np.all(
                        [lo-1e-9 <= val <= hi+1e-9
                         for val, (lo, hi) in zip(result.x, bounds)])
                    better = feasible and cost < incumbent-1e-5
                    row = {"seed": name, "solver": "SciPy_SLSQP",
                           "success": bool(result.success), "message": str(result.message),
                           "iterations": int(result.nit), "retained_SG_MW": cost,
                           "alpha": alpha, "feasible": feasible,
                           "better_than_frozen_by_1e-5_MW": better,
                           "elapsed_s": time.monotonic()-start,
                           "evaluations_cumulative": evaluations}
                except Exception as exc:
                    row = {"seed": name, "solver": "SciPy_SLSQP",
                           "success": False, "message": str(exc),
                           "iterations": 0, "retained_SG_MW": math.nan,
                           "alpha": math.nan, "feasible": False,
                           "better_than_frozen_by_1e-5_MW": False,
                           "elapsed_s": time.monotonic()-start,
                           "evaluations_cumulative": evaluations}
                results.append(row)
                print("FALSIFICATION", row, flush=True)
            handle.write("QUIT\n")
            handle.flush()
    finally:
        if proc.poll() is None:
            proc.terminate()
        proc.wait(timeout=30)
    path = OUT / "TABLE_N10_generic_solver_falsification.csv"
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(results[0]))
        writer.writeheader()
        writer.writerows(results)
    print("GENERIC_SOLVER_BETTER_POINT",
          any(r["better_than_frozen_by_1e-5_MW"] for r in results))


if __name__ == "__main__":
    main()
