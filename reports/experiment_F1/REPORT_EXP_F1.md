# Experiment F1 — Ideal second-order graph damping

**F1_STATUS: PASS**

For `M>0`, `K>=0`, mass-normalization gives `K̃=M^{-1/2}KM^{-1/2}`. In the commuting modal basis, each scalar equation is `s²+dₖs+νₖ=0`. For underdamping, the decay rate is `dₖ/2`; after critical damping, it is `(dₖ-√(dₖ²-4νₖ))/2`, which decreases as `dₖ` increases. Therefore the isolated-mode maximum occurs at `dₖ*=2√νₖ`, yielding `D*=2M^{1/2}K̃^{1/2}M^{1/2}`. This is not a claim of global optimality over arbitrary noncommuting damping matrices.

A deterministic synthetic 5-mode SPD case compared underdamped, critical, and overdamped modal damping. Maximum analytic-to-numeric pole error: `2.787913920690064e-7`; modal transformation residual: `1.2351231148954867e-14`. The acceptance tolerance is `1e-6 s^-1`; the largest discrepancy occurs at a repeated critical pole, whose direct eigenvalue is ill-conditioned.

No nonproportional optimum or Weyl–Horn matrix is proposed: this experiment limits itself to the proved modal result and explicitly marks the broader reference unimplemented. No GFL claim is made.
