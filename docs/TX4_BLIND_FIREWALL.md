# TX4 operational blind-prediction firewall

Literal human blindness is impossible because historical TX4 answer tables
are in the repository. The operational separation is enforced by process and
file paths.

## Allowed before reveal

The prediction command may import only the frozen model/equilibrium builders,
all-SG network and operating-point data, local SG/GFL device equations,
controller parameters, candidate bus/rating metadata, and the exact audited
port-closure implementation. It may write only new files under the
`TX4_BLIND_PREDICTION_` prefix plus the preregistered V4/V9 prediction names.

The baseline kernel and local device library are serialized with source-path
and SHA-256 metadata. Every prediction records that it did not load an answer
key. The prediction files are hashed and committed before the reveal command
is allowed to read historical full-order tables.

## Forbidden before reveal

The prediction phase must not read `FC01_structure.csv`, `PCV02_core_portfolios`,
`PCV03_counterfactual`, `FC10_census_lattice`, V4/V9 answer tables, historical
H/kappa labels, or historical g-boundary files. It must not call any helper
that returns full-order portfolio alpha/verdict data. It may use no full-order
portfolio eigenanalysis; full-order spectra are reveal-stage ground truth.

## Reveal boundary

The reveal is a separate command and commit after prediction artifacts exist.
It reads frozen full-order tables and produces new comparison files. Any
prediction/full mismatch is reported rather than used to change the blind
algorithm. The firewall is a reproducibility protocol, not a claim of
cryptographic blindness.
