"""Record installed IEEE-39 sources and immutable audit-case definitions."""
from __future__ import annotations

import csv
import hashlib
import tomllib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "reports" / "experiment_M" / "tables"
PROVENANCE = ROOT / "reports" / "experiment_M" / "SOFTWARE_PROVENANCE.toml"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(name: str, rows: list[dict]) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / name).open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def candidate(path: Path) -> dict:
    actual = digest(path)
    expected = Path(str(path) + ".sha256").read_text(encoding="utf-8").strip()
    if actual != expected:
        raise RuntimeError(f"frozen candidate hash mismatch: {path}")
    return tomllib.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    prov = tomllib.loads(PROVENANCE.read_text(encoding="utf-8"))
    source = Path(prov["ieee39_source"]["path"])
    source_dir = source.parent / "ieee39data"
    hashes = []
    for name in ("bus", "branch", "load", "machine", "avr", "gov"):
        p = source_dir / f"{name}.csv"
        expected = prov["ieee39_inputs"][name]["sha256"]
        hashes.append(dict(input_file=p.name, installed_path=str(p), sha256=digest(p),
                           provenance_sha256=expected, matches_provenance=digest(p) == expected))
    write("TABLE_M03_input_hashes.csv", hashes)
    with (source_dir / "bus.csv").open(newline="", encoding="utf-8") as fh:
        buses = list(csv.DictReader(fh))
    inv = []
    for row in buses:
        inv.append(dict(bus=row["bus"], category=row["category"], bus_type=row["bus_type"],
            has_generator=row["has_gen"], machine_model="SauerPaiMachine" if row["has_gen"] == "true" else "",
            avr_model="AVRTypeI" if row["has_avr"] == "true" else "",
            governor_model="TGOV1" if row["has_gov"] == "true" else "",
            has_zip_load=row["has_load"], load_model="ZIPLoad" if row["has_load"] == "true" else "",
            line_model="PiLine_fault", base_kv=row["base_kv"],
            pf_model=row["bus_type"], source_sha256=digest(source)))
    write("TABLE_M02_model_inventory.csv", inv)
    gpath = ROOT / "reports" / "experiment_G" / "Z_G_FINAL.toml"
    kpath = ROOT / "reports" / "experiment_K" / "Z_K_NOMINAL_FINAL.toml"
    g, k = candidate(gpath), candidate(kpath)
    grow = {int(x["bus"]): x for x in g["generator"]}
    definitions = []
    for case in ("all_SG", "ExpG_candidate", "ExpK_nominal"):
        for bus in range(30, 40):
            i = bus - 30
            if case == "all_SG":
                rho, kp, ki, source_path, sha = 0., 0., 0., str(source), digest(source)
            elif case == "ExpG_candidate":
                x = grow[bus]
                rho, kp, ki = float(x["rho"]), float(x["Kp"]), float(x["Ki"])
                source_path, sha = str(gpath), digest(gpath)
            else:
                rho, kp, ki = float(k["rho"][i]), float(k["Kp"][i]), float(k["Ki"][i])
                source_path, sha = str(kpath), digest(kpath)
            definitions.append(dict(case=case, bus=bus, rho=rho, epsilon=1-rho, Kp=kp,
                Ki=ki, architecture="GFL" if rho == 1 else "SG" if rho == 0 else "SG+GFL",
                source_file=source_path, source_sha256=sha))
    write("TABLE_M05_case_definitions.csv", definitions)
    print("installed IEEE-39 files hashed; audit cases frozen")


if __name__ == "__main__":
    main()
