"""BC00-B / BC01: physical spectrum of a solved case, with the full zero ledger.

Pipeline for one solved ``ReplacementCase``:

1. Jacobians at the equilibrium at two finite-difference scales. Their
   difference estimates the matrix error that the classifier uses.
2. Dead states: rows of A that are identically zero, i.e. states frozen by a
   service switch at 0 (flux_blend = 0, avr_blend = 0). A zero row makes A
   block triangular, so det(sI - A) = s^k det(sI - A'). They are removed
   exactly and ledgered as DEAD_STATE, not deleted by size.
3. Rotation generator (R_x, R_z), derived from the device frames. The DAE and
   reduced identities are checked.
4. Frequency partner w (uniform frequency shift) when every machine has zero
   damping, with the Jordan identity A w = omega_B R_x checked.
5. A_q, the quotient by R (T1), is the physical model, and its classification
   is the physical status. The neutral-frequency mode, if present, is an exact
   zero of A_q and makes that status BOUNDARY_OR_UNRESOLVED.
6. A_qq, which is A_q deflated by the exact neutral eigenvector Z^T w, is the
   rest of the spectrum. It is used to count right-half-plane eigenvalues
   without a disc around zero. The neutral mode stays in the ledger.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np

from ..dynamics.linearize import central_difference_jacobians
from .classify import SAFETY, classify_spectrum, near_zero_table
from .symmetry import (
    dae_identities,
    deflate_eigenvector,
    frequency_partner,
    quotient,
    rotation_generator,
)


def _reduced(jac):
    return jac.fx - jac.fz @ np.linalg.solve(jac.gz, jac.gx)


@dataclass
class PhysicalReport:
    n_states: int
    n_dead: int
    dead_labels: list
    matrix_error: float
    a_norm: float
    identities: dict
    partner: str
    neutral_verified: bool
    neutral_value: float | None
    status_physical: str  # classification of A_q (the physical model)
    status_rest: str  # classification of A_qq (A_q minus the declared neutral mode)
    n_positive_rest: int
    n_positive_real_rest: int
    n_rhp_raw_rest: int
    abscissa_rest: float
    d_axis_rest: float
    threshold: float
    near_zero_full: dict
    near_zero_q: dict
    near_zero_rest: dict
    near_axis_rest: list
    positive_rest: list

    def as_row(self) -> dict:
        row = asdict(self)
        for key in ("identities", "near_zero_full", "near_zero_q", "near_zero_rest"):
            for k, v in row.pop(key).items():
                row[f"{key}:{k}"] = v
        row["dead_labels"] = ";".join(self.dead_labels)
        row["near_axis_rest"] = ";".join(f"{v:.3e}" for v in self.near_axis_rest)
        row["positive_rest"] = ";".join(f"{v:.3e}" for v in self.positive_rest)
        return row


def physical_matrices(case):
    """(A, dA_norm, jac) at the case equilibrium; A from the unit-scale Jacobian."""

    dae = case.dae
    x, z = case.equilibrium.x, case.equilibrium.z
    jac1 = central_difference_jacobians(dae, x, z, {})
    jac2 = central_difference_jacobians(dae, x, z, {}, scale_x=2.0, scale_z=2.0)
    a1 = _reduced(jac1)
    return a1, a1 - _reduced(jac2), jac1


def physical_report(case, *, safety: float = SAFETY) -> PhysicalReport:
    dae = case.dae
    a_full, d_mat, jac = physical_matrices(case)
    r_x, r_z = rotation_generator(dae, case.equilibrium.z)
    partner = frequency_partner(dae)
    ident = dae_identities(jac, a_full, r_x, r_z, partner.w)

    dead = np.flatnonzero(~np.any(a_full != 0.0, axis=1))
    keep = np.setdiff1d(np.arange(a_full.shape[0]), dead)
    a = a_full[np.ix_(keep, keep)]
    d_a = float(np.linalg.norm(d_mat, 2))
    rx = r_x[keep]
    w = None if partner.w is None else partner.w[keep]

    q = quotient(a, rx)
    d_k = d_mat[np.ix_(keep, keep)]
    norm = float(np.linalg.norm(a_full, 2))
    tol = safety * (d_a / norm + 1e-12)
    verified = (
        w is not None
        and ident.a_jordan is not None
        and ident.a_jordan <= tol
        and ident.a_rot <= tol
    )
    rep_q = classify_spectrum(q.a_q, q.z.T @ d_k @ q.z, safety, known_zero=verified)
    neutral_value = None
    rest, d_rest = q.a_q, q.z.T @ d_k @ q.z
    if verified:
        v_n = q.z.T @ w
        v_n = v_n / np.linalg.norm(v_n)
        neutral_value = float(v_n @ q.a_q @ v_n)
        qq = deflate_eigenvector(q.a_q, v_n)
        rest, d_rest = qq.a_q, qq.z.T @ d_rest @ qq.z
    rep_rest = classify_spectrum(rest, d_rest, safety)
    return PhysicalReport(
        n_states=int(a_full.shape[0]),
        n_dead=int(dead.size),
        dead_labels=[dae.labels[k] for k in dead],
        matrix_error=d_a,
        a_norm=norm,
        identities=asdict(ident),
        partner=partner.reason,
        neutral_verified=bool(verified),
        neutral_value=neutral_value,
        status_physical=rep_q.status,
        status_rest=rep_rest.status,
        n_positive_rest=rep_rest.n_positive,
        n_positive_real_rest=rep_rest.n_positive_real,
        n_rhp_raw_rest=rep_rest.n_rhp_raw,
        abscissa_rest=rep_rest.abscissa,
        d_axis_rest=rep_rest.d_axis,
        threshold=rep_rest.threshold,
        near_zero_full=near_zero_table(np.linalg.eigvals(a_full)),
        near_zero_q=near_zero_table(np.array([r.value for r in rep_q.records])),
        near_zero_rest=near_zero_table(np.array([r.value for r in rep_rest.records])),
        near_axis_rest=[r.value for r in rep_rest.near_axis()],
        positive_rest=[r.value for r in rep_rest.records if r.sign == "POSITIVE"],
    )
