# Experiment A → Experiment B operator interface

## Bundle contract

Experiment B accepts a `BNDOperatorBundle` with:

- `M`: real symmetric positive-definite retained inertia matrix.
- `L`: real symmetric positive-semidefinite synchronizing stiffness on the retained coordinates.
- `sigma(s)`: callable returning a dense `Matrix{ComplexF64}` of size `m × m` for complex `s` in the declared valid domain.
- `sigma_derivative(s)`: callable returning the analytic derivative with the same shape. If unavailable, mark it unavailable; do not silently use a finite-difference derivative in the pole predictor.
- `metadata`: `source`, `retained_state_names`, `units`, `sign_convention`, `operating_point`, `model_version`, `exact_or_approximate`, and `valid_frequency_band`.

The sign convention must define which side of the retained equation the self-energy force occupies. Units must identify the interpretation of `s`, angular frequency `ω` in rad/s, and plotted frequency `f=ω/(2π)` in Hz. `retained_state_names` must follow the row and column ordering in `M`, `L`, and `Σ`.

## Version 1 portable format

`load_expA_bundle(path)` reads a TOML file with `schema_version = "1.0"`, `M`, `L`, `metadata`, and `sigma_available`. When `sigma_available = true`, `[realization]` contains real matrices `D0`, `Ac`, `B`, and `C` describing `Σ(s)=D0+C(sI−Ac)⁻¹B`. The serialized realization is an interchange option, not a requirement on how ExpA computes its operator.

When ExpA exports only `T_exact(s)` or another object from which it has not derived a defensible `Σ(s)`, set `sigma_available = false` and omit `[realization]`. The loader then returns a bundle with `sigma === nothing`. B-final must stop before graph self-energy calculations; it must not infer or fabricate `Σ(s)`.

## Validation before B-final

Check matrix dimensions, positive definiteness of `M`, Hermitian/PSD properties and declared construction of `L`, finite response values, analytic derivative consistency, state ordering, units, sign convention, operating point, model version, exactness flag, and validity band. The generalized graph mathematics itself accepts the validated bundle without changes.
