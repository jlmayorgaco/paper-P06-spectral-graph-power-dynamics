"""Assemble the final ExpD report from the frozen candidate and PD artifacts."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "reports" / "experiment_D"
TABLES = REPORT / "tables"


def main() -> None:
    candidate = json.loads((REPORT / "ANALYTIC_CANDIDATE.json").read_text(encoding="utf-8"))
    validation = json.loads((REPORT / "PD_VALIDATION_EXP_D.json").read_text(encoding="utf-8"))
    optima = pd.DataFrame(candidate["cross_bus_candidates"])
    pd_endpoints = pd.read_csv(TABLES / "TABLE_D11_PD_validation.csv")
    grid = pd.read_csv(TABLES / "TABLE_D17_PD_rho_PLL_grid.csv")
    spectra = pd.read_csv(TABLES / "TABLE_D19_PD_grid_full_spectra.csv")
    if validation["candidate_sha256"] != (REPORT / "ANALYTIC_CANDIDATE.sha256").read_text(encoding="utf-8").strip():
        raise RuntimeError("PowerDynamics validation does not refer to the currently frozen candidate")
    if len(pd_endpoints) != 4 or len(grid) != 132:
        raise RuntimeError("Expected four endpoint validations and the complete 132-cell grid")
    if len(spectra[["bus", "rho", "beta"]].drop_duplicates()) != 132:
        raise RuntimeError("Full-spectrum file does not cover all 132 grid cells")
    for key, group in spectra.groupby(["bus", "rho", "beta"]):
        row = grid.loc[(grid["bus"] == key[0]) & (grid["rho"] == key[1]) & (grid["beta"] == key[2])].iloc[0]
        nongauge = group.loc[~group["is_gauge"].astype(bool)]
        if len(nongauge) != int(row["non_gauge_poles"]):
            raise RuntimeError(f"Pole count mismatch in full-spectrum grid cell {key}")
        alpha = float(nongauge["lambda_real"].max())
        if abs(alpha - float(row["alpha"])) > 1e-8:
            raise RuntimeError(f"Rightmost-pole mismatch in full-spectrum grid cell {key}")
        if bool(group["mode_pass"].astype(bool).all()) != bool(row["within_sigma"]):
            raise RuntimeError(f"All-mode classification mismatch in grid cell {key}")

    cross = []
    for item in candidate["cross_bus_candidates"]:
        bus = int(item["bus"])
        endpoint = pd_endpoints.loc[pd_endpoints["bus"] == bus].iloc[0]
        surface = grid.loc[grid["bus"] == bus]
        rho1 = surface.loc[surface["rho"] == 1.0]
        cross.append(
            {
                "bus": bus,
                "P0_MW": item["P0_MW"],
                "rho_star": item["rho_star"],
                "P_GFL_max_MW": item["MW_star"],
                "Kp_star_rad_per_s": item["Kp_star"],
                "Ki_star_rad_per_s2": item["Ki_star"],
                "active_constraint": item["active_constraint"],
                "rightmost_real": item["rightmost_pole_real"],
                "rightmost_imag": item["rightmost_pole_imag"],
                "analytic_alpha": item["alpha"],
                "PD_alpha": endpoint["PD_alpha"],
                "PD_abs_pole_error": endpoint["nearest_pole_abs_error"],
                "PD_endpoint_status": endpoint["status"],
                "grid_cells": len(surface),
                "grid_feasible_cells": int(surface["within_sigma"].astype(bool).sum()),
                "rho1_feasible_gain_points": int(rho1["within_sigma"].astype(bool).sum()),
                "rho1_gain_points": len(rho1),
            }
        )
    cross_df = pd.DataFrame(cross).sort_values("bus")
    cross_df.to_csv(TABLES / "TABLE_D18_cross_bus_optima.csv", index=False)

    final_gates = {
        "D0 design/validation firewall": "PASS",
        "D1 PLL parameter structure": "PASS_EXACT_AFFINE_LOW_RANK",
        "D2 Woodbury PLL representation": "PASS_LOCAL_IDENTITY",
        "D3 continuous rho physical model": "PASS_AT_LOCAL_PORT_AND_PD_GRID",
        "D4 rho update structure": "PASS_LOCAL_AFFINE_RHO",
        "D5 global low-rank update": "NOT_USED_FOR_SATURATED_BOUND_CERTIFICATE",
        "D6 replacement closure": "NOT_REQUIRED_FOR_SATURATED_ENDPOINT_PROOF",
        "D7 algebraic PLL boundary": "NOT_ENUMERATED_NONBINDING",
        "D8 algebraic rho solution": "NOT_ENUMERATED_NONBINDING",
        "D9 analytic maximum rho": "PASS_ZERO_GAP_SATURATED_BOUND",
        "D10 multimode safety": "PASS_FULL_STATE_ENDPOINT_SPECTRA",
        "D11 independent PowerDynamics validation": "PASS_ALL_FOUR_ENDPOINTS",
        "D12 validation grid": "COMPLETE_132_CELLS_WITH_9_VIOLATIONS",
        "D13 cross-bus optimum replication": "PASS_BUSES_30_33_35_37",
        "D14 self-energy mechanism": "NOT_REQUIRED_FOR_OPTIMUM_CERTIFICATE",
        "D15 ready for ExpE": "YES_WITH_BUS37_INTERIOR_HOLE_RECORDED",
    }
    pd.DataFrame([{"gate": k, "status": v} for k, v in final_gates.items()]).to_csv(
        TABLES / "TABLE_D16_gate_summary.csv", index=False
    )

    report = r"""# Experiment D — Single-bus SG-to-GFL replacement optimum

## Result

**EXP_D_STATUS: PASS for the preregistered single-bus objective and PLL domain.** The primary bus-33 optimum is

\[
\rho^\star=1,\qquad K_p^\star=31.4159265359\;\mathrm{rad/s},\qquad
K_i^\star=246.7401100272\;\mathrm{rad/s^2},\qquad
P_{\rm GFL,max}=632\;\mathrm{MW}.
\]

The full-replacement rightmost pole is \(-0.09884700665\pm j0.43026220071\;\mathrm{s^{-1}}\), giving a spectral margin of 0.09884700665 s⁻¹ against the required 0.05 s⁻¹. The detailed PowerDynamics pole agrees exactly at the stored precision. The same saturated optimum is feasible at buses 30, 35, and 37; see `TABLE_D18_cross_bus_optima.csv`.

## Certified design domain

The frozen preregistration specifies

\[
0\le\rho\le1,\qquad 0.9\le\beta\le1.1,\qquad
K_p=\beta K_{p0},\quad K_i=\beta^2K_{i0},
\]

with \(K_{p0}=5(2\pi)\), \(K_{i0}=K_{p0}^2/4\), and required decay rate \(\sigma_{\rm req}=0.05\;\mathrm{s^{-1}}\). The proportional and integral gains therefore move together along the preregistered physical PLL-bandwidth coordinate. Nominal gains are feasible and are the unique zero-effort controller tie-break.

## Global optimality certificate

For each bus, \(P_i^0>0\) and \(0\le\rho\le1\), so every admissible point satisfies

\[
P_i^0\rho\le P_i^0.
\]

The frozen candidate has the witness \((\rho,\beta)=(1,1)\). The full-state ExpC matrix has 111 finite modes at each full-GFL endpoint, including one gauge mode; every one of the 110 non-gauge modes satisfies \(\Re\lambda\le-0.05\). A fresh PowerDynamics equilibrium and full linearization independently verified all four witnesses. Thus the objective upper bound is attained and the global optimality gap is zero. No admissible combination can replace more than 100% of the dispatch. The lexicographic controller-effort tie-break selects \(\beta=1\), hence the reported \(K_p^\star,K_i^\star\). The objective-bound argument is exact; spectral feasibility is numerically verified from the full double-precision eigenvalue sets and independently cross-checked against PowerDynamics, with the engineering-margin decision applied using a 1e-10 s⁻¹ tolerance. No interval-arithmetic eigenvalue enclosure is claimed.

This is an exact saturated-bound certificate. We did not construct a system-level polynomial \(F(s,\rho,K_p,K_i)\), compute resultants, or enumerate nonbinding pole-boundary branches: none can improve the objective past the attained physical upper bound \(\rho=1\). Consequently, this result certifies the global maximum replacement and selected tie-break point, but does not claim a complete algebraic map of every pole-boundary branch.

## Frozen candidate and independent endpoint validation

The analytical candidate was frozen before the validation adapter imported PowerDynamics. Frozen candidate SHA-256:

```text
""" + validation["candidate_sha256"] + r"""
```

PowerDynamics ran after the hash gate. All four buses pass equilibrium qualification, the 0.05 s⁻¹ all-mode margin, and rightmost-pole matching; the absolute pole and alpha errors are zero in `TABLE_D11_PD_validation.csv`. The analytic candidate computation itself uses the frozen ExpC full-state matrices and makes zero PowerDynamics calls.

## Post-freeze PowerDynamics grid

The independent validation grid contains 11 replacement shares \(\rho=0,0.1,\ldots,1\) and three preregistered PLL scales \(\beta=0.9,1,1.1\) at each of buses 30, 33, 35, and 37: 132/132 cells evaluated. Every eigenvalue for every cell is retained in `TABLE_D19_PD_grid_full_spectra.csv`; the row-level grid also stores the rightmost non-gauge pole. An audit confirms that full-spectrum mode counts and all-mode pass/fail classifications reproduce every grid row. Overall, 123 cells meet the required margin. Buses 30, 33, and 35 pass all 33 points. Bus 37 passes 24/33; all nine violations occur at \(\rho=0.7,0.8,0.9\), one for each beta, with equilibria converged. At \(\rho=1\), all three gain points pass at all four buses.

The bus-37 result is nonmonotone: the interior mixed realization has a right-half-plane or weakly damped mode at the three high-share grid columns, while the zero-share SG is structurally removed at \(\rho=1\) and the endpoint is stable. The endpoint remains the global maximum for the stated static replacement objective; the interior stability holes matter for staged replacement and must be considered in any path-dependent implementation. The figure and row-level data are `figures/FIG_D01_PD_validation_grid.png` and `tables/TABLE_D17_PD_rho_PLL_grid.csv`.

## Supporting results and limitations

- The local nine-state PLL model has exact affine dependence on its implemented gains; the Kp and Ki updates have rank one each and joint rank two. The Woodbury local-resolvent factorization and same-equilibrium local current-port update remain documented in the structural-audit artifacts.
- This ExpD optimum uses exact full-state endpoint spectra, not a graph-mode truncation or a low-rank closure. It establishes replacement capacity under the stated frozen small-signal model and controller domain.
- The grid follows the preregistered one-dimensional PLL coordinate. If the intended problem is instead an independent rectangular Kp/Ki domain, its numerical bounds must be frozen separately; the saturated objective proof is unchanged whenever the nominal gain point remains admissible, but the validation surface must be recomputed on that domain.
- The current-sharing model has no current limiter. No current-limit, nonlinear transient, EMT, uncertainty-robustness, or multi-bus co-design claim is made here.

## Decision

ExpD delivers the requested numerical optimum for bus 33 and the three cross-bus repetitions, with a zero-gap global certificate and independent endpoint PowerDynamics validation. The post-freeze grid also exposes the bus-37 interior stability hole. ExpE can begin from these endpoint results, carrying that bus-37 behavior as a constraint on any staged or partial replacement policy.
"""
    (REPORT / "REPORT_EXP_D.md").write_text(report, encoding="utf-8")

    ledger = """# Claim ledger — Experiment D

| ID | Claim | Status | Evidence / scope |
|---|---|---|---|
| D-C01 | Local PLL gain dependence is affine and low rank. | EXACT | Kp and Ki each rank one; joint rank two in the independently coded nine-state model. |
| D-C02 | Local rated SG/GFL ports form an affine rho mixture at the frozen equilibrium. | LOCAL_EXACT | Device-port identity and cross-bus local-port audit; no current-limiter claim. |
| D-C03 | A full-system determinant-lemma closure was used for the optimizer. | NOT_USED | The attained physical upper bound proves the objective without pole-branch elimination; no closure is claimed. |
| D-C04 | Bus-33 replacement maximum is globally certified on the preregistered domain. | PASS | rho is bounded above by 1 and the all-mode-feasible rho=1, beta=1 witness attains that bound. |
| D-C05 | The bus-33 optimum matches detailed PowerDynamics. | PASS | Exact stored rightmost-pole match after the candidate SHA gate. |
| D-C06 | The optimum repeats at buses 30, 35, and 37. | PASS | All full-replacement endpoints pass the all-mode margin and fresh PowerDynamics validation. |
| D-C07 | The post-freeze surface confirms feasibility over the whole design box. | PARTIAL | 132/132 cells evaluated; 123 meet the margin, with nine bus-37 interior violations. The optimum endpoint is confirmed. |
| D-C08 | Every rho boundary branch or stationary point was eliminated algebraically. | NOT_CLAIMED | This is unnecessary for the saturated global maximum and was not computed. |
| D-C09 | The result supports nonlinear current-limited deployment. | NOT_ESTABLISHED | The current-sharing model has no current limiter; only equilibrium and small-signal spectra are covered. |
"""
    (REPORT / "CLAIM_LEDGER_EXP_D.md").write_text(ledger, encoding="utf-8")

    results_path = REPORT / "RESULTS_EXP_D.json"
    results = json.loads(results_path.read_text(encoding="utf-8"))
    results["validation_calls"] = validation["powerdynamics_calls"]
    results["cross_bus_scope"] = "FOUR_BUS_FULL_REPLACEMENT_ENDPOINT_AND_PD_GRID"
    results["cross_bus"] = {str(row["bus"]): row["PD_endpoint_status"] for row in cross}
    results["replacement"].update(
        {
            "rho_model": "SATURATED_ENDPOINT_OPTIMUM_CERTIFIED_BY_PHYSICAL_UPPER_BOUND",
            "local_port_rho_model": "AFFINE_IN_RHO",
            "full_operator_dimension": 111,
            "closure_dimension": None,
            "boundary_enumeration": "NOT_REQUIRED_FOR_ATTAINED_RHO_UPPER_BOUND",
            "rho_star": 1.0,
        }
    )
    results.update(
        {
            "status": "PASS",
            "status_scope": "Saturated global replacement optimum on preregistered PLL bandwidth domain; exact branch map not claimed",
            "candidate_frozen": True,
            "candidate_sha256": validation["candidate_sha256"],
            "candidate": candidate["primary_candidate"],
            "cross_bus_optima": cross_df.to_dict(orient="records"),
            "global_certificate": candidate["global_certificate"],
            "design_domain": candidate["domain"],
            "powerdynamics": validation,
            "validation_grid": {
                "status": validation["grid_status"],
                "rows": len(grid),
                "feasible_rows": int(grid["within_sigma"].astype(bool).sum()),
                "spectral_violations": int((grid["status"] == "SPECTRAL_VIOLATION").sum()),
                "full_spectrum_rows": len(spectra),
                "csv": "tables/TABLE_D17_PD_rho_PLL_grid.csv",
                "full_spectra_csv": "tables/TABLE_D19_PD_grid_full_spectra.csv",
                "figure": "figures/FIG_D01_PD_validation_grid.png",
            },
            "ready_for_expE": True,
            "gates": final_gates,
        }
    )
    results_path.write_text(json.dumps(results, indent=2, allow_nan=False) + "\n", encoding="utf-8")

    readme = """# Experiment D reproduction

## Frozen outcome

ExpD certifies the saturated replacement optimum on the existing preregistered PLL domain: `rho in [0,1]`, `beta in [0.9,1.1]`, `Kp=beta*Kp0`, and `Ki=beta^2*Ki0`. Bus 33 reaches 632 MW at `rho*=1`, `Kp*=31.4159265359 rad/s`, and `Ki*=246.7401100272 rad/s^2`. The same endpoint is feasible at buses 30, 35, and 37. The exact objective certificate is the physical upper bound `P0*rho <= P0` plus a full-spectrum-feasible witness at `rho=1`; resultants and nonbinding pole-boundary branches are not claimed.

## Reproduction commands

Run the structural local-model audit:

```powershell
julia --project=. experiments/bnd_expD/run_expD_structure.jl
```

Freeze the analytical endpoint candidate using only frozen ExpC full-state matrices:

```powershell
julia --project=. experiments/bnd_expD/freeze_endpoint_candidates.jl
```

Only after the candidate and its SHA-256 sidecar are frozen, validate the four endpoints and evaluate the 132-cell PowerDynamics grid:

```powershell
julia --project=. experiments/bnd_expD/run_expD_validation_only.jl
```

Render the grid figure and final report artifacts:

```powershell
python experiments/bnd_expD/plot_validation_grid.py
python experiments/bnd_expD/finalize_expD_results.py
```

`run_expD_validation_only.jl` checks the candidate's frozen flag, status, and SHA-256 before importing PowerDynamics. The validation surface uses 11 values of rho and three preregistered PLL bandwidth scales at each of buses 30, 33, 35, and 37. It records unstable cells rather than removing them. Bus 37 has nine intermediate-share spectral violations at rho 0.7–0.9, despite a feasible rho=1 endpoint.

## Main artifacts

- `ANALYTIC_CANDIDATE.json` and `ANALYTIC_CANDIDATE.sha256`: frozen analytical witness and hash.
- `tables/TABLE_D09_analytic_optimum.csv`: selected optima.
- `tables/TABLE_D10_multimode_check.csv`: complete analytic endpoint spectra.
- `tables/TABLE_D11_PD_validation.csv`: fresh endpoint-to-PowerDynamics match.
- `tables/TABLE_D17_PD_rho_PLL_grid.csv`: all post-freeze grid cells.
- `tables/TABLE_D19_PD_grid_full_spectra.csv`: every PowerDynamics pole at all grid cells.
- `tables/TABLE_D18_cross_bus_optima.csv`: final four-bus summary.
- `figures/FIG_D01_PD_validation_grid.png`: spectral-margin surface.
- `REPORT_EXP_D.md`, `RESULTS_EXP_D.json`, and `CLAIM_LEDGER_EXP_D.md`: final result and scope.

The local PLL and port structural-audit run writes its preliminary artifacts with `_STRUCTURE` suffixes so that it cannot overwrite the frozen candidate, full-spectrum tables, validation results, or final report.

## Model limits

The current-share model scales terminal current and retains each device's internal equations for interior rho; zero-share devices are removed at the endpoints. The model has no current limiter. The report therefore makes a frozen-operating-point small-signal claim and documents the bus-37 interior stability hole; it does not claim nonlinear or current-limited deployment performance.
"""
    (REPORT / "README_REPRODUCE.md").write_text(readme, encoding="utf-8")
    print(f"WROTE {TABLES / 'TABLE_D18_cross_bus_optima.csv'}")
    print(f"WROTE {REPORT / 'REPORT_EXP_D.md'}")


if __name__ == "__main__":
    main()
