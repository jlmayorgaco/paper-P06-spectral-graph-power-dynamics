# Common Uncertainty Specification

All portfolio radii use one physical normalization. The planned block structure is `Delta_all = blkdiag(Delta_30,...,Delta_38)` with fixed left/right weights declared in `prereg/uncertainty_weights_manifest.json`. A portfolio selects blocks; it does not refit an uncertainty scale.

The transverse robust calculation is valid only after removal of the structural zero/Jordan chain and a well-posed nominal loop. If the available Julia package supports only diagonal complex perturbations, the result is labeled as that conservative/limited structure and never called a general full-block mu result.

The base-radius and compatibility definitions are only evaluated for `0 <= epsilon < r_base`. A nominally transverse-unstable portfolio has radius zero by definition.
