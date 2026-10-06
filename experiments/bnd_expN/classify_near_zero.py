"""Classify near-zero modes using explicit gauge quotient and state participation."""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from scipy.linalg import eig, null_space, solve, svdvals
from scipy.optimize import linear_sum_assignment

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "reports" / "experiment_N"
OLD_G = ROOT / "reports" / "experiment_M" / "matrices" / "ExpG_candidate"


def read(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def matrix(path: Path) -> np.ndarray:
    return np.atleast_2d(np.loadtxt(path, delimiter=","))


def current(original: dict, bus: int) -> complex:
    row = original[bus]
    v = complex(float(row["V_real_pu"]), float(row["V_imag_pu"]))
    s = complex(float(row["P_gen_MW"]), float(row["Q_gen_Mvar"])) / 100
    return np.conj(s / v)


def case_modes(label: str, path: Path, original: dict, source: str) -> list[dict]:
    a = matrix(path / "PD_A.csv")
    m = matrix(path / "PD_M.csv")
    state = read(path / "PD_state_map.csv")
    dyn = np.flatnonzero(np.diag(m) == 1)
    alg = np.flatnonzero(np.diag(m) == 0)
    red = a[np.ix_(dyn, dyn)] - a[np.ix_(dyn, alg)] @ solve(
        a[np.ix_(alg, alg)], a[np.ix_(alg, dyn)])
    names = [state[j]["state_name"] for j in dyn]
    values = {r["state_name"]: float(r["equilibrium_value"])
              for r in state if "equilibrium_value" in r and r["equilibrium_value"]}
    g = np.zeros(len(dyn))
    for k, name in enumerate(names):
        if name.endswith("machine₊δ)") or name.endswith("gfl₊pll₊θ)"):
            g[k] = 1
        elif name.endswith("gfl₊filter₊i_f_i)"):
            if source == "ExpM_uncorrected":
                g[k] = values[name.replace("i_f_i)", "i_f_r)")]
            else:
                bus = int(name.split(",", 1)[0].split("(")[1])
                g[k] = current(original, bus).real
        elif name.endswith("gfl₊filter₊i_f_r)"):
            if source == "ExpM_uncorrected":
                g[k] = -values[name.replace("i_f_r)", "i_f_i)")]
            else:
                bus = int(name.split(",", 1)[0].split("(")[1])
                g[k] = -current(original, bus).imag
    if not np.linalg.norm(g):
        raise RuntimeError(f"no gauge generator in {label}")
    gres = np.linalg.norm(red @ g) / (np.linalg.norm(red) * np.linalg.norm(g))
    q = null_space(g[None, :])
    quotient = q.T @ red @ q
    raw_lambda, raw_vectors = eig(red)
    qlambda, qvectors = eig(quotient)
    near = np.flatnonzero(np.abs(raw_lambda) < 1e-5)
    # Identify the removed gauge by asking which one-to-one deletion reproduces
    # the explicitly quotiented spectrum best, instead of removing by magnitude.
    omissions = []
    for j in near:
        remain = np.delete(raw_lambda, j)
        rr, cc = linear_sum_assignment(np.abs(remain[:, None] - qlambda[None, :]))
        omissions.append((float(max(np.abs(remain[rr] - qlambda[cc]))), int(j)))
    if not omissions:
        raise RuntimeError(f"no near-zero gauge candidate in {label}")
    quotient_error, gauge_idx = min(omissions)
    singular = svdvals(red)
    sv_right = np.linalg.svd(red, full_matrices=False)[2][-1]
    null_overlap = abs(np.vdot(sv_right, g)) / (np.linalg.norm(sv_right) * np.linalg.norm(g))
    out = []
    for j in near:
        vec = raw_vectors[:, j]
        overlap = abs(np.vdot(vec, g)) / (np.linalg.norm(vec) * np.linalg.norm(g))
        qdist = float(min(np.abs(qlambda - raw_lambda[j])))
        if j == gauge_idx and gres < 1e-8 and null_overlap > 0.9:
            classification = "GAUGE"
        elif qdist < 1e-6 and j != gauge_idx:
            classification = "PHYSICAL"
        elif j != gauge_idx and qdist >= 1e-6:
            classification = "NUMERICAL_ARTIFACT"
        else:
            classification = "UNRESOLVED"
        physical_vector = q @ qvectors[:, int(np.argmin(np.abs(qlambda - raw_lambda[j])))]
        generalized_coupling = float(np.linalg.norm(red @ physical_vector))
        if classification == "PHYSICAL":
            subtype = ("GENERALIZED_ZERO_CHAIN" if generalized_coupling > 1e-8
                       else "INDEPENDENT_ZERO_EIGENVECTOR")
            top_vector = physical_vector
        else:
            subtype = "GLOBAL_ROTATION" if classification == "GAUGE" else "NONE"
            top_vector = vec
        top = np.argsort(np.abs(top_vector) ** 2)[-3:][::-1]
        out.append(dict(case=label, source=source, mode_index=int(j) + 1,
            real=float(raw_lambda[j].real), imag=float(raw_lambda[j].imag),
            classification=classification, subtype=subtype,
            generalized_coupling_norm=generalized_coupling,
            gauge_residual=float(gres),
            eigenvector_gauge_overlap=float(overlap),
            smallest_singular_value=float(singular[-1]),
            second_smallest_singular_value=float(singular[-2]),
            singular_null_gauge_overlap=float(null_overlap),
            quotient_mode_distance=qdist,
            gauge_omission_spectrum_error=quotient_error,
            top_state_1=names[top[0]], top_state_2=names[top[1]],
            top_state_3=names[top[2]]))
    return out


def main() -> None:
    original = {int(r["bus"]): r for r in
                read(OUT / "TABLE_N01_original_operating_point.csv")}
    rows = case_modes("ExpG_uncorrected_ExpM", OLD_G, original, "ExpM_uncorrected")
    exported = read(OUT / "TABLE_N05_validation_case_export.csv")
    for case in exported:
        rows.extend(case_modes(case["case"], OUT / "validation_cases" / case["case"],
                               original, "ExpN_corrected"))
    path = OUT / "TABLE_N03_near_zero_mode_classification.csv"
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print("N03", len(rows), "GAUGE", sum(r["classification"] == "GAUGE" for r in rows),
          "PHYSICAL", sum(r["classification"] == "PHYSICAL" for r in rows),
          "UNRESOLVED", sum(r["classification"] == "UNRESOLVED" for r in rows))
    old = [r for r in rows if r["case"] == "ExpG_uncorrected_ExpM"]
    print("OLD_EXPG", [(r["classification"], r["subtype"], r["real"]) for r in old])


if __name__ == "__main__":
    main()
