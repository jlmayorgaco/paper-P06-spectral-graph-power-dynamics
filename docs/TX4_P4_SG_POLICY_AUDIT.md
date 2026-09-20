# TX4 P4 Surviving-SG Policy Audit

## Frozen policy

The exact P4 tuple is `(g,k,t,h)=(0.03625, 1.425, 1.5, 1.0)`.

- `g=0.03625` is the common GFL Q/V control gain.
- `k=1.425` multiplies the AVR gain `KA` of every surviving SG.
- `t=1.5` multiplies the AVR time constant `TE` of every surviving SG.
- `h=1.0` is the frozen heterogeneity exponent; with `h=1`, no additional bus-dependent perturbation is applied.

## Python implementation

The exact Python run constructs `ConverterParameters(voltage_control=True, voltage_gain=0.03625, voltage_leak=0.05)` and passes `machine_scaling={"ka": 1.425, "ta": 1.5}` into `solve_case`. The `machine_scaling` path multiplies the surviving SG `KA` and `TE` fields before initialization and linearization. The replacement policy is `q_policy="matched"` at all 16 subsets.

## Julia implementation

The exact Julia run declares the same constants and applies `ka=avr[2]*1.425` and `ta=avr[3]*1.5` in `machine_from_row`. The GFL implementation uses `P4_G=0.03625`, `kp_v=2`, `ki_v=20`, and `leak=0.05`. The Julia census uses the same bus order and the same matched AC operating point.

## Audit result

The P4 constants are explicit in both reproduction sources, and the resulting all-SG parity row and 16-case cross-code census are reported in:

- `results/TX4_EXACT_P4_ALLSG_PARITY.csv`;
- `results/TX4_EXACT_P4_V4_CROSSCODE.csv`;
- `results/TX4_EXACT_P4_CLAIMS.csv`.

No parameter was fitted to the Julia output, and no alternate P4 value was substituted.
