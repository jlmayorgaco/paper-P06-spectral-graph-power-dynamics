# Experiment C backbone sensitivity audit

Conclusion: **PARTIALLY_ROBUST**

The primary backbone was preregistered as `Mhalf * Pg * Herm(Minvhalf*L0*Minvhalf) * Pg * Mhalf`. The secondary backbone is the Euclidean gauge-preserving Hermitian part of `L0`. No pole-fitting criterion was used.

Maximum principal angle: 0.4541531103915557 degrees. Maximum observed Gamma magnitude ratio: 4.559220446739163. Results are case-specific; inspect TABLE_C14 and TABLE_C04 for spectra and modal mappings.

Classification criteria were fixed in expC.toml: ROBUST requires maximum principal angle <= 10.0 degrees and maximum Gamma ratio <= 2.0; PARTIALLY_ROBUST requires maximum angle <= 45 degrees; otherwise BASIS_SENSITIVE.

