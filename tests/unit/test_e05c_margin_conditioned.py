from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
MODULE = ROOT / "experiments" / "tx3" / "E05C_margin_conditioned" / "margin_conditioned.py"
SPEC = importlib.util.spec_from_file_location("margin_conditioned", MODULE)
assert SPEC is not None and SPEC.loader is not None
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)


def test_pair_and_triple_lower_order_closure() -> None:
    pair = {(): -1 + 2j, (1,): -0.9 + 2.1j, (2,): -0.8 + 1.9j, (1, 2): -0.5 + 2.3j}
    pair_delta = M.complex_externality(pair, (1, 2))
    assert abs(M.stable_lower_order_prediction(pair, (1, 2)) + pair_delta - pair[(1, 2)]) < 1e-14
    triple = {
        (): -1 + 2j,
        (1,): -0.9 + 2.1j,
        (2,): -0.8 + 1.9j,
        (3,): -1.1 + 2.2j,
        (1, 2): -0.5 + 2.3j,
        (1, 3): -0.7 + 2.0j,
        (2, 3): -0.6 + 2.4j,
        (1, 2, 3): -0.2 + 2.8j,
    }
    triple_delta = M.complex_externality(triple, (1, 2, 3))
    assert abs(M.stable_lower_order_prediction(triple, (1, 2, 3)) + triple_delta - triple[(1, 2, 3)]) < 1e-14


def test_candidate_union_reuses_shared_vertices() -> None:
    coalitions = [(7, 8), (2, 7, 8), (2, 5, 8)]
    vertices = M.union_vertices(coalitions)
    assert () in vertices and (7, 8) in vertices and (2, 7, 8) in vertices
    assert len(vertices) < sum(2 ** len(value) for value in coalitions)


def test_frozen_participation_is_normalized_and_positive() -> None:
    assert abs(float(M.STRESS_PARTICIPATION.sum()) - 1.0) < 1e-15
    assert (M.STRESS_PARTICIPATION > 0.0).all()
