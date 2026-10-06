"""Audit algebra in THEORY_GRAPH_ROBUST_CONTROL.txt; not a grid certificate.

Compares independently assembled port interconnections, graph-coordinate blocks,
ellipsoid support values, Schur signs, and permutation invariance. No SG/GFL
trajectory or robust-feasibility claim follows from this floating-point audit.
"""
from pathlib import Path
import hashlib
import json
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "reports/nonlinear_codesign_20261001"
THEORY = REPORT / "THEORY_GRAPH_ROBUST_CONTROL.txt"
rng = np.random.default_rng(20261001)
errors = {k: 0.0 for k in (
    "port_assembly", "spectral_blocks", "spectral_eigenvalues",
    "ellipsoid_attainment", "ellipsoid_boundary", "permutation",
    "spectral_endpoint_excess",
)}
schur_checks = 0
positive_tests = 0
negative_tests = 0


def track(name, value):
    errors[name] = max(errors[name], float(abs(value)))


for case in range(120):
    n = int(rng.integers(2, 10))
    edges = [(i - 1, i) for i in range(1, n)]
    edges += [(i, j) for i in range(n) for j in range(i + 2, n)
              if rng.random() < 0.3]
    incidence = np.zeros((n, len(edges)))
    for e, (i, j) in enumerate(edges):
        incidence[i, e], incidence[j, e] = 1.0, -1.0
    L = (incidence * rng.uniform(0.1, 1.0, len(edges))) @ incidence.T
    eye, zero = np.eye(n), np.zeros((n, n))
    T = np.block([[-L, eye], [eye, zero]])

    # Heterogeneous supply: compare substitution to the explicit formula (29).
    Q = np.diag(rng.uniform(-0.5, 0.5, n))
    S = np.diag(rng.uniform(-1.0, 1.0, n))
    R = np.diag(rng.uniform(-2.0, 0.0, n))
    rsum = np.diag(rng.uniform(0.0, 0.5, n))
    mu = float(rng.uniform(0.5, 3.0))
    rw = np.diag(rng.uniform(0.5, 2.0, n))
    supply = np.block([[Q, S], [S.T, R]])
    direct = T.T @ supply @ T + np.block([[zero, zero], [zero, rsum-mu*rw]])
    explicit = np.block([
        [L.T@Q@L-L.T@S-S.T@L+R, S.T-L.T@Q],
        [S-Q@L, Q+rsum-mu*rw],
    ])
    track("port_assembly", np.max(np.abs(direct-explicit)))

    # A simultaneous relabeling must preserve this matrix up to congruence.
    perm = eye[rng.permutation(n)]
    fullperm = np.block([[perm, zero], [zero, perm]])
    Lp, Qp, Sp, Rp = [perm@a@perm.T for a in (L, Q, S, R)]
    Tp = np.block([[-Lp, eye], [eye, zero]])
    Pp = np.block([[Qp, Sp], [Sp.T, Rp]])
    relabeled = Tp.T@Pp@Tp + np.block([
        [zero, zero], [zero, perm@(rsum-mu*rw)@perm.T]])
    track("permutation", np.max(np.abs(relabeled-fullperm@direct@fullperm.T)))

    # Homogeneous contracts: the full matrix is exactly an orthogonal sum
    # of the 2x2 graph blocks, even though the device dynamics need not be linear.
    q, s, sigma, r = 0.02, 0.1, 0.05, 1.3
    d = 2.0 if case % 2 == 0 else 0.001
    delta = mu*r-q-sigma
    lam, U = np.linalg.eigh(L)
    homogeneous = T.T @ np.block([[q*eye, s*eye], [s*eye, -d*eye]]) @ T
    homogeneous += np.block([[zero, zero], [zero, (sigma-mu*r)*eye]])
    transform = np.block([[U, zero], [zero, U]])
    expected = np.block([
        [np.diag(q*lam**2-2*s*lam-d), np.diag(s-q*lam)],
        [np.diag(s-q*lam), -delta*eye],
    ])
    track("spectral_blocks", np.max(np.abs(transform.T@homogeneous@transform-expected)))
    block_eigenvalues = []
    for ll in lam:
        aa, cc = q*ll**2-2*s*ll-d, s-q*ll
        eig = np.linalg.eigvalsh(np.array([[aa, cc], [cc, -delta]]))
        block_eigenvalues.extend(eig)
        schur = aa+cc**2/delta
        if min(abs(eig[-1]), abs(schur)) > 1e-9:
            assert (eig[-1] < 0) == (schur < 0)
            schur_checks += 1
            positive_tests += int(schur < 0)
            negative_tests += int(schur > 0)
    track("spectral_eigenvalues", np.max(np.abs(
        np.linalg.eigvalsh(homogeneous)-np.sort(block_eigenvalues))))
    grid = np.linspace(0.0, max(0.0, lam[-1]), 301)
    vals = q*grid**2-2*s*grid-d+(s-q*grid)**2/delta
    track("spectral_endpoint_excess", max(0.0, np.max(vals)-max(vals[0], vals[-1])))

    # Arbitrary ellipsoidal input directions; not confined to one node.
    A = rng.normal(size=(n, n))
    W = A.T@A+0.5*eye
    E = rng.normal(size=(n+2, n))
    a = rng.normal(size=n+2)
    b = E.T@a
    support = np.sqrt(b@np.linalg.solve(W, b))
    wstar = np.linalg.solve(W, b)/support
    track("ellipsoid_attainment", a@E@wstar-support)
    track("ellipsoid_boundary", wstar@W@wstar-1.0)

assert max(errors.values()) < 1e-10, errors
assert positive_tests > 0 and negative_tests > 0
result = {
    "status": "ALGEBRAIC_IDENTITIES_AUDITED_NOT_A_POWER_SYSTEM_CERTIFICATE",
    "seed": 20261001,
    "random_graphs": 120,
    "schur_sign_checks": schur_checks,
    "negative_semidefinite_blocks": positive_tests,
    "blocks_failing_condition": negative_tests,
    "max_absolute_errors": errors,
    "theory_sha256": hashlib.sha256(THEORY.read_bytes()).hexdigest(),
    "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "numpy_version": np.__version__,
    "scope": "Floating-point algebra checks; the conditional proof is in the theory document. This audit does not certify SG/GFL dynamics; the later coupled result has its own verifier and restricted scope.",
}
out = REPORT / "graph_control_identity_audit.json"
out.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
print(json.dumps(result, indent=2))
