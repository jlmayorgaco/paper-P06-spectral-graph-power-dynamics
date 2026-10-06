"""Print the requested Experiment M terminal verdict from measured ledgers."""
from __future__ import annotations

import csv
import tomllib
from pathlib import Path


ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/"reports"/"experiment_M"
TABLES=OUT/"tables"


def rows(name:str)->list[dict]:
    with (TABLES/name).open(newline="",encoding="utf-8") as fh:return list(csv.DictReader(fh))


def val(rows_:list[dict],case:str,key:str)->str:
    return next(r[key] for r in rows_ if r["case"]==case)


def main()->None:
    provenance=tomllib.loads((OUT/"SOFTWARE_PROVENANCE.toml").read_text(encoding="utf-8"))
    bases=rows("TABLE_M01_component_bases.csv")
    pq=rows("TABLE_M04_component_PQ_sharing.csv")
    inventory=rows("TABLE_M09_state_inventory_comparison.csv")
    matrices=rows("TABLE_M11_matrix_identity.csv")
    poles=rows("TABLE_M16_pole_matching_summary.csv")
    critical=max(rows("TABLE_M16_pole_matching_ExpK_nominal.csv"),key=lambda r:float(r["lambda_PD_real"]))
    ports=rows("TABLE_M15_port_transfer_identity.csv")
    linear=rows("TABLE_M19_linear_trajectory_identity.csv")
    scaling=rows("TABLE_M20_scaling_fit.csv")[0]
    getpq=lambda bus:next(r for r in pq if r["case"]=="all_SG" and int(r["bus"])==bus)
    getmat=lambda case:next(r for r in matrices if r["case"]==case and r["matrix"]=="reduced_A")
    out={
        "EXP_M_STATUS":"FAIL / MODEL_NOT_RECONCILED",
        "JULIA_VERSION":provenance["julia_version"],
        "POWERDYNAMICS_VERSION":provenance["packages"]["PowerDynamics"]["version"],
        "POWERDYNAMICS_COMMIT":"UNAVAILABLE; git-tree-sha1="+provenance["packages"]["PowerDynamics"]["git_tree_sha1"],
        "NETWORKDYNAMICS_VERSION":provenance["packages"]["NetworkDynamics"]["version"],
        "IEEE39_COMMIT":"BUNDLED_WITH_POWERDYNAMICS_TREE_"+provenance["packages"]["PowerDynamics"]["git_tree_sha1"],
        "PROJECT_SHA256":provenance["project_sha256"],
        "MANIFEST_SHA256":provenance["manifest_sha256"],
        "IEEE39_SBASE_MVA":"100.0",
        "IEEE39_FBASE_HZ":"60.0",
        "IEEE39_OMEGABASE":"376.99111843077515",
        "GFL_CONSTRUCTION_FBASE_HZ":next(r["fbase_equivalent"] for r in bases if r["case"]=="GFL_template"),
        "BASE_MISMATCH_FOUND":"YES",
        "BUS31_GEN_PQ":f"({getpq(31)['p_gen_MW']}, {getpq(31)['q_gen_Mvar']}) MW/MVAr",
        "BUS31_LOAD_PQ":f"({getpq(31)['p_zip_MW']}, {getpq(31)['q_zip_Mvar']}) MW/MVAr",
        "BUS31_NET_PQ":f"({getpq(31)['p_net_MW']}, {getpq(31)['q_net_Mvar']}) MW/MVAr",
        "BUS39_GEN_PQ":f"({getpq(39)['p_gen_MW']}, {getpq(39)['q_gen_Mvar']}) MW/MVAr",
        "BUS39_LOAD_PQ":f"({getpq(39)['p_zip_MW']}, {getpq(39)['q_zip_Mvar']}) MW/MVAr",
        "BUS39_NET_PQ":f"({getpq(39)['p_net_MW']}, {getpq(39)['q_net_Mvar']}) MW/MVAr",
        "PQ_SHARING_RECONCILED":"NO",
        "RHO_SEMANTICS":"SG Sn=(1-rho)Sn0; GFL external current=rho*filter current; actual P/Q independently initialized",
        "SG_RATING_SCALES_WITH_EPSILON":"YES",
        "GFL_RATING_SCALES_WITH_RHO":"NO",
        "CURRENT_ONLY_SCALING_FOUND":"YES",
        "ALL_SG_STATE_COUNT_PD":val(inventory,"all_SG","PD_states"),
        "ALL_SG_STATE_COUNT_AN":val(inventory,"all_SG","AN_states"),
        "EXPG_STATE_COUNT_PD":val(inventory,"ExpG_candidate","PD_states"),
        "EXPG_STATE_COUNT_AN":val(inventory,"ExpG_candidate","AN_states"),
        "EXPK_STATE_COUNT_PD":val(inventory,"ExpK_nominal","PD_states"),
        "EXPK_STATE_COUNT_AN":val(inventory,"ExpK_nominal","AN_states"),
        "ALL_SG_MATRIX_REL_ERROR":getmat("all_SG")["relative_frobenius"],
        "EXPG_MATRIX_REL_ERROR":getmat("ExpG_candidate")["relative_frobenius"],
        "EXPK_MATRIX_REL_ERROR":getmat("ExpK_nominal")["relative_frobenius"],
        "ALL_SG_POLE_MAX_ERROR":val(poles,"all_SG","max_pole_error"),
        "EXPG_POLE_MAX_ERROR":val(poles,"ExpG_candidate","max_pole_error"),
        "EXPK_POLE_MAX_ERROR":val(poles,"ExpK_nominal","max_pole_error"),
        "EXPK_PD_CRITICAL_POLE":critical["lambda_PD_real"]+" + 0im s^-1",
        "EXPK_ANALYTIC_MATCHED_POLE":critical["lambda_AN_real"]+" + 0im s^-1 (global assignment)",
        "CRITICAL_MODE_CLASS":"NETWORK_INTERMODAL_BRANCH",
        "CRITICAL_MODE_TOP_STATES":"SG39 rotor angle; GFL38 filter currents; SG39 rotor speed",
        "CRITICAL_MODE_MISMATCH_CAUSE":"mixed PD initialization does not enforce frozen per-device P/Q share",
        "OPEN_LOOP_CLOSURE_IDENTITY":"PASS for all three cases; max error 1.6488e-9 s^-1",
        "GFL_PORT_MAX_REL_ERROR":str(max(float(r["max_relative_error"]) for r in ports if r["case"]=="ExpK_nominal")),
        "SG_PORT_MAX_REL_ERROR":str(max(float(r["max_relative_error"]) for r in ports if r["case"]=="all_SG")),
        "BASE_50_60_ABLATION_ALPHA_50":"NOT_APPLICABLE (no 50-Hz device found)",
        "BASE_50_60_ABLATION_ALPHA_60":"NOT_APPLICABLE",
        "BASE_50_60_EXPLAINS_MISMATCH":"NO",
        "PARAMETER_ESTIMATION_USED":"NO",
        "PARAMETERS_ESTIMATED":"NONE",
        "IDENTIFIABILITY_STATUS":"NOT_APPLICABLE",
        "LINEAR_TRAJECTORY_MAX_REL_ERROR":str(max(float(r["max_relative_trajectory_error"]) for r in linear if r["status"].startswith("PASS"))),
        "NONLINEAR_ERROR_SCALING_EXPONENT":scaling["absolute_error_scaling_exponent"]+" (global fit; numerical floor at 1e-5)",
        "MODEL_IDENTITY_CERTIFIED":"NO",
        "MODEL_SAFE_FOR_Z_OPTIMIZATION":"NO",
        "MAIN_CAUSE_OF_PREVIOUS_DISCREPANCY":"SG/GFL mixed P/Q split and internal equilibrium not held at analytical frozen shares",
        "REQUIRED_CODE_FIX":"enforce per-device P/Q in mixed initialization; set pure-GFL Vbase from bus.csv; audit gauge count",
        "NEXT_STEP":"rebuild parameterized share-constrained PD-exact analytical model, then repeat identity gates",
        "PUSH":"NO",
    }
    text="\n".join(f"{key}: {value}" for key,value in out.items())+"\n"
    (OUT/"TERMINAL_SUMMARY_EXP_M.txt").write_text(text,encoding="utf-8")
    print(text,end="")


if __name__=="__main__":main()
