# Weak-link application report

Case: Mix60/no-PSS
Line: Line_42 (20-34)
Perturbation: 1% line strengthening

Observed effect:
- Global margin zeta: 0.123566404 -> 0.123558419
- Tracked network-mode zeta: 0.133884200 -> 0.133815595
- Distance to matched control pole: 3.332268 -> 3.327898

Why static analysis fails:
- Retained stiffness term A: -1.091e-09
- Retained rotation term B: 7.974e-09
- Controller self-energy term C: -6.874e-03
- Control denominator fraction: 0.959

Engineering interpretation:
The line is weak for this action because strengthening it pulls a network-family
mode toward a PLL-dominated condensed-control pole. A graph-only or stiffness-only
ranking does not contain that pole distance and therefore misses the hidden margin.

Recommended action:
Do not apply line-only reinforcement as the first fix. Retune the implicated PLL/control
pole, add damping, or choose a reinforcement that increases damping without reducing
the network-control pole distance.
