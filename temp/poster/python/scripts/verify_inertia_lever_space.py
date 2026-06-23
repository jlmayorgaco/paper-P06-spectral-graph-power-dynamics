from __future__ import annotations

import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
POSTER = ROOT / "poster_ieee_session"
DATA = POSTER / "data"


def connected_components(adj: np.ndarray) -> int:
    n = adj.shape[0]
    seen = np.zeros(n, dtype=bool)
    count = 0
    for start in range(n):
        if seen[start]:
            continue
        count += 1
        stack = [start]
        seen[start] = True
        while stack:
            node = stack.pop()
            for nbr in np.where(adj[node])[0]:
                if not seen[nbr]:
                    seen[nbr] = True
                    stack.append(int(nbr))
    return count


def random_connected_graph(rng: np.random.Generator, n: int, density: float) -> list[tuple[int, int, float]]:
    while True:
        edges: list[tuple[int, int, float]] = []
        adj = np.zeros((n, n), dtype=bool)
        for i in range(n):
            for j in range(i + 1, n):
                if rng.random() < density:
                    weight = float(rng.uniform(0.25, 2.5))
                    edges.append((i, j, weight))
                    adj[i, j] = True
                    adj[j, i] = True
        if edges and connected_components(adj) == 1:
            return edges


def laplacian(n: int, edges: list[tuple[int, int, float]]) -> np.ndarray:
    lap = np.zeros((n, n), dtype=float)
    for i, j, weight in edges:
        lap[i, j] -= weight
        lap[j, i] -= weight
        lap[i, i] += weight
        lap[j, j] += weight
    return lap


def edge_laplacian(n: int, i: int, j: int) -> np.ndarray:
    b = np.zeros(n, dtype=float)
    b[i] = 1.0
    b[j] = -1.0
    return np.outer(b, b)


def diagonal_basis(n: int, i: int) -> np.ndarray:
    mat = np.zeros((n, n), dtype=float)
    mat[i, i] = 1.0
    return mat


def relative_residual(delta_l: np.ndarray, basis: list[np.ndarray]) -> tuple[float, float]:
    target = delta_l.reshape(-1)
    design = np.column_stack([mat.reshape(-1) for mat in basis])
    coeffs, *_ = np.linalg.lstsq(design, target, rcond=None)
    residual = target - design @ coeffs
    abs_res = float(np.linalg.norm(residual))
    rel_res = abs_res / max(float(np.linalg.norm(target)), 1.0)
    return abs_res, rel_res


def main() -> None:
    rng = np.random.default_rng(20260609)
    cases = []
    max_abs = 0.0
    max_rel = 0.0
    graph_count = 40

    for graph_index in range(graph_count):
        n = int(rng.integers(4, 12))
        density = float(rng.uniform(0.33, 1.0))
        edges = random_connected_graph(rng, n, density)
        lap = laplacian(n, edges)
        edge_basis = [edge_laplacian(n, i, j) for i, j, _ in edges]
        diag_basis = [diagonal_basis(n, i) for i in range(n)]
        basis = edge_basis + diag_basis

        graph_max_abs = 0.0
        graph_max_rel = 0.0
        for node in range(n):
            eps = np.zeros((n, n), dtype=float)
            eps[node, node] = 1.0
            delta_l = 0.5 * (eps @ lap + lap @ eps)
            abs_res, rel_res = relative_residual(delta_l, basis)
            graph_max_abs = max(graph_max_abs, abs_res)
            graph_max_rel = max(graph_max_rel, rel_res)

        max_abs = max(max_abs, graph_max_abs)
        max_rel = max(max_rel, graph_max_rel)
        cases.append(
            {
                "graph_index": graph_index,
                "n": n,
                "edge_count": len(edges),
                "density": density,
                "max_abs_residual": graph_max_abs,
                "max_relative_residual": graph_max_rel,
            }
        )

    DATA.mkdir(parents=True, exist_ok=True)
    summary = {
        "seed": 20260609,
        "graph_count": graph_count,
        "n_range": [4, 11],
        "density_range": [0.33, 1.0],
        "max_abs_residual": max_abs,
        "max_relative_residual": max_rel,
        "interpretation": "The inertia-induced perturbation lies in existing-edge Laplacians plus diagonal shunts; no new-edge component is needed.",
        "cases": cases,
    }
    out = DATA / "inertia_lever_space_verification.json"
    with out.open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print(
        json.dumps(
            {
                "graph_count": graph_count,
                "max_abs_residual": max_abs,
                "max_relative_residual": max_rel,
                "output": str(out),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
