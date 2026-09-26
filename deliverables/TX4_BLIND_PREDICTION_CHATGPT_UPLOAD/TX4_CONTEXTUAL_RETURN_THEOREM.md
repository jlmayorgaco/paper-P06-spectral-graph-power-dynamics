# Contextual return identity used in TX4

## Scope

This note records the exact block-determinant identity used to interpret the
four-bus TX4 witness. It is standard Schur-complement algebra applied to the
frozen port-closure matrix; it is not claimed as a new theorem.

## Definitions

For a selected portfolio `H`, let

```text
C_H(s) = I + Q_H(s)
```

and choose one device `i` in `H`. Let `R = H \ {i}` and partition `Q_H` in
the corresponding block coordinates. When `I + Q_RR(s)` is nonsingular,
define

```text
G_R(s)       = (I + Q_RR(s))^{-1}
R_i|R(s)     = Q_iR(s) G_R(s) Q_Ri(s)
```

The second quantity is the contextual return from device `i` through the
remaining devices and network closure.

## Identity

The block determinant formula gives

```text
det(I + Q_H)
  = det(I + Q_RR) det(I - R_i|R).
```

Therefore, provided the proper-subset factor is regular, a zero of the full
collective closure is equivalent to a unit eigenvalue of the contextual
return (with the sign convention used by the frozen TX4 port model, the
reported eigenvalue approaches `-1`).

## TX4 numerical check

At the reduced/full H4 controller boundary, the nearest return eigenvalue in
the archived convention was `-0.9999999649 - 5.1e-9 j`. The maximum recorded
Schur residual was `3.44e-16`. The local factor remained regular: the reduced
boundary run reported local minimum singular value `0.4892355`, while the
collective factor approached singularity at `1.38e-8`.

These values support a collective closure interpretation for this frozen
model and operating policy. They do not establish universal local regularity,
a global stability-radius theorem, or a probabilistic claim.
