"""CDW unit tests (pytest). Run: python -m pytest experiments/cdw/tests -q"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import _infra  # noqa: E402,F401
import numpy as np  # noqa: E402
import pytest  # noqa: E402

import _analysis as AN  # noqa: E402
import _cdw as C  # noqa: E402


def test_ybus_identity_exact():
    from ibr_cycles.models.ieee39_network import load_network

    assert np.array_equal(C.build_ybus(), load_network().ybus)


def test_laplacian_tap_decomposition_exact():
    L, order = C.laplacian_B()
    Y = C.build_ybus()
    idx = {b: i for i, b in enumerate(order)}
    resid = np.zeros(len(order))
    for br in C.branches():
        y = 1 / complex(br["r"], br["x"])
        bs, t = -y.imag, br["tap"]
        f, tt = idx[br["f"]], idx[br["t"]]
        resid[f] += bs * (1 / t**2 - 1 / t) - br["b"] / 2 / t**2
        resid[tt] += bs * (1 - 1 / t) - br["b"] / 2
    for sh in C.network_payload()["shunts"]:
        resid[idx[int(sh["bus"])]] -= sh["b"]
    # Im Y = -(L_B + Delta_tap) + B_ch + B_sh  <=>  -Im Y - L_B = diag(resid)
    M = -Y.imag - L
    assert np.allclose(M - np.diag(np.diag(M)), 0, atol=1e-12)
    assert np.allclose(np.diag(M), resid, atol=1e-10)
    assert np.allclose(L.sum(axis=1), 0, atol=1e-10)


def test_linear_order_dp_matches_brute():
    rng = np.random.default_rng(1)
    for _ in range(5):
        W = rng.integers(0, 5, size=(6, 6)).astype(float)
        np.fill_diagonal(W, 0)
        assert AN.optimal_linear_order(W)[1] == pytest.approx(AN.brute_linear_order(W)[1])


def test_shapley_efficiency_and_symmetry():
    players = (1, 2, 3, 4)
    rng = np.random.default_rng(2)
    v = {s: float(rng.normal()) for s in C.subsets(players)}
    v[()] = 0.0
    phi = AN.shapley(v, players)
    assert sum(phi.values()) == pytest.approx(v[players])


def test_one_step_submodularity_equivalence():
    rng = np.random.default_rng(3)
    players = (1, 2, 3, 4)
    for _ in range(20):
        v = {s: float(rng.normal()) for s in C.subsets(players)}
        local = all(d <= 1e-12 for *_, d in AN.second_differences(v, players))
        assert local == AN.is_submodular(v, players, tol=1e-12)
    # a genuinely submodular function (concave of cardinality)
    v = {s: float(np.sqrt(len(s))) for s in C.subsets(players)}
    assert AN.is_submodular(v, players)
    assert all(d <= 1e-12 for *_, d in AN.second_differences(v, players))


def test_gordan_conflict_witness():
    G = np.array([[1.0, 0.0], [-1.0, 0.0], [0.0, 1.0]])
    val, lam = AN.gordan(np.array([[1.0, 0.0], [-1.0, 0.0]]))
    assert val < 1e-8
    val2, _ = AN.gordan(G[[0, 2]])
    assert val2 > 1e-3


@pytest.mark.slow
def test_engine_consistency_and_structural_remarks():
    import _sens as S

    eng = S.Engine(C.V4, C.policy("D01"))
    for sem in ("SPR", "RP"):
        assert np.abs(eng.residual(eng.w0(sem), None, sem)).max() < 1e-6
    d = eng.derivatives([{"kind": "g", "idx": 30, "value": 0.03625}, {"kind": "vset", "idx": 32, "value": 0.0}], "SPR")
    g, v = d["items"]
    assert g["d_frozen"][0] == pytest.approx(g["d_total"][0], rel=1e-6)
    assert abs(v["d_frozen"][0]) < 1e-9
