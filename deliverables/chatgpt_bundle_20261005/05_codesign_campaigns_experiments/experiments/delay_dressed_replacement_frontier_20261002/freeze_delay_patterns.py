"""Serialize and hash the pre-registered delay patterns before any DDE run."""
from __future__ import annotations
import csv, hashlib, json
from pathlib import Path

HERE = Path(__file__).resolve().parent
patterns_path = HERE / "FROZEN_DELAY_PATTERNS.csv"
metrics_path = HERE / "TABLE_D05_GSP_DELAY_METRICS.csv"
with patterns_path.open(newline="", encoding="utf-8") as f:
    patterns = list(csv.DictReader(f))
with metrics_path.open(newline="", encoding="utf-8") as f:
    metric = {r["delay_pattern_id"]: r for r in csv.DictReader(f)}
uniform = [0, 2, 5, 10, 20, 30, 40, 50]
payload = {
    "freeze_stage": "created before any delayed poles, optimization, or time-domain result",
    "candidate_buses": list(range(30, 40)),
    "uniform_grid_ms": uniform,
    "heterogeneous_multiset_rule_ms": "tau_j=(j-1)*50/9, j=1,...,10",
    "random_seed": 20261002,
    "random_permutation_count": 100,
    "heuristic_search": {"candidate_permutations": 100000, "seed": 20261002,
                         "selection": "best/worst commutator and low/high spectral roughness among sampled permutations"},
    "fixed_multiset_ms": [float(x) for x in patterns[0]["tau_ms"].split(";")],
    "multiset_sha256": patterns[0]["multiset_sha256"],
    "graph_definition": "Kron-reduced lossless physical synchronizing graph on buses 30:39; original SG inertia reference mass; see GRAPH_OPERATOR_AUDIT.toml and saved matrices",
    "patterns": []
}
for r in patterns:
    t = [float(x) for x in r["tau_ms"].split(";")]
    m = metric[r["delay_pattern_id"]]
    payload["patterns"].append({
        "id": r["delay_pattern_id"], "tau_ms_by_bus_30_to_39": t,
        "mean_tau_ms": float(r["mean_tau_ms"]), "std_tau_ms": float(r["std_tau_ms"]),
        "max_tau_ms": float(r["max_tau_ms"]), "vector_sha256": r["vector_sha256"],
        "multiset_sha256": r["multiset_sha256"], "chi_tau": float(m["chi_tau"]),
        "graph_roughness": float(m["graph_roughness"]),
        "low_frequency_score": float(m["low_frequency_score"]),
        "high_frequency_score": float(m["high_frequency_score"]),
        "operator_name": m["operator_name"],
    })
serialized = json.dumps(payload, sort_keys=True, indent=2) + "\n"
out = HERE / "FROZEN_DELAY_PATTERNS.json"
out.write_text(serialized, encoding="utf-8")
print("FROZEN", hashlib.sha256(serialized.encode()).hexdigest(), "patterns", len(payload["patterns"]))
