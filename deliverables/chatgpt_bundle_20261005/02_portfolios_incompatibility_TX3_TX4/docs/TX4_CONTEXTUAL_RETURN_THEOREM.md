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
return. The signs are different because the closure uses `I + Q_H`, whereas
the contextual return uses `I - R_i|R`: the closure eigenvalue approaches
`-1` and the return eigenvalue approaches `+1`.

## TX4 numerical check

At the reduced/full H4 controller boundary, the fresh sign audit found a
nearest `Q_H` eigenvalue at `-1` within `3.51e-8` and contextual-return
eigenvalues at `+1` within at most `2.28e-7`. The maximum Schur residual was
`3.63e-16`. The physical local factors `I+M_ii` remained regular, with a
minimum singular value of `0.2973` over the frozen sweep, while the
collective factor approached singularity at `1.26e-8`.

These values support a collective closure interpretation for this frozen
model and operating policy. They do not establish universal local regularity,
a global stability-radius theorem, or a probabilistic claim.
