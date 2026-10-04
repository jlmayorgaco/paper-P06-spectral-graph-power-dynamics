"""Use preserved input copies for portable reproduction."""
from pathlib import Path
OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[1]
def resolve(relative):
    p=OUT/'inputs'/Path(relative)
    return p if p.exists() else ROOT/relative
