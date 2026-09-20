# TX4 Robustness Deviations Log

1. A single non-recorded H4 timing probe was run before preregistration to
   size the campaign. It produced no retained result, was not used for
   parameter selection, and is excluded from every table and claim.
2. The prior audit contains no calibrated physical uncertainty weights. The MC
   campaign therefore uses explicitly assumed independent uniforms on the
   bounded engineering intervals and is labeled a bounded-box engineering
   probability, not a physical population probability.
3. The exact custom DAE exposes a reduced spectral closure directly, but the
   prior PowerDynamics network gate stopped before same-model network parity.
   The Julia part is consequently a declared cross-code spot-check and not a
   fresh PowerDynamics certificate.
4. `q_H4` is operationalized as a minimum H4 non-gauge distance to the
   imaginary axis. It is not the archived collective-return `|1+q|` quantity
   and is not interpreted as a physical robust radius.
5. The required report names use CSV/Parquet-compatible files. If the local
   Parquet runtime is unavailable, the complete CSV remains the authoritative
   result and the limitation is recorded in the execution log rather than
   silently omitting rows.
