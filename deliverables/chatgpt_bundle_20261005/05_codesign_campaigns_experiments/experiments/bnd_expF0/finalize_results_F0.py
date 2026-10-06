"""Append candidate classification/projection status to the compact F0 manifest."""
from pathlib import Path
import tomllib

import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/"reports"/"experiment_F0"
result=tomllib.loads((OUT/"RESULTS_EXP_F0.toml").read_text(encoding="utf-8"))
candidates=pd.read_csv(OUT/"tables"/"TABLE_F04_candidate_classification.csv")
c0=candidates[candidates.candidate=="ExpC C0 primary synchronizing backbone L_G"].iloc[0]
c33=candidates[candidates.candidate=="ExpC C33 primary synchronizing backbone L_G"].iloc[0]
result["synchronizing_backbone_C33_lambda_min"]=float(c33.lambda_min)
result["synchronizing_backbone_C33_lambda_max"]=float(c33.lambda_max)
result["synchronizing_backbone_C33_psd"]=bool(c33.psd)
result["mode_comparison_status"]=(OUT/"F0_MODE_OVERLAP_STATUS.md").read_text(encoding="utf-8").strip()

def toml_value(value):
    if isinstance(value,bool): return "true" if value else "false"
    if isinstance(value,str): return '"'+value.replace('\\','\\\\').replace('"','\\"')+'"'
    if isinstance(value,list): return "["+", ".join(toml_value(x) for x in value)+"]"
    return repr(value)

(OUT/"RESULTS_EXP_F0.toml").write_text(
    "\n".join(f"{key} = {toml_value(value)}" for key,value in result.items())+"\n",encoding="utf-8")
report=OUT/"REPORT_EXP_F0.md"
text=report.read_text(encoding="utf-8")
if "## Synchronizing-backbone classification" not in text:
    text += ("\n## Synchronizing-backbone classification\n\n"
        f"The frozen ExpC all-SG C0 synchronizing backbone has spectrum "
        f"[{float(c0.lambda_min):.6g}, {float(c0.lambda_max):.6g}] and is PSD within tolerance. "
        f"The bus-33 GFL C33 backbone has spectrum [{float(c33.lambda_min):.6g}, {float(c33.lambda_max):.6g}] and is "
        f"{'PSD' if bool(c33.psd) else 'indefinite'}. The C33 operator is not sent through a real square root.\n")
report.write_text(text,encoding="utf-8")
print("F0 result manifest and backbone classification updated")
