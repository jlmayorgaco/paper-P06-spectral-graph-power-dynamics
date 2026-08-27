# E05B connected spectral shift and contour moments

`Xi_S(s)` was evaluated directly as the signed trace sum of finite-pole resolvents. Exact residue moments
were compared with the preregistered 128-node contour integral. Candidate decisions are:

| Coalition | Valid resolved holdout | Median Re mu1 | Median Im mu1 | Decision |
|---|---:|---:|---:|---|
| A7-A8 | 16/16 | 0.00465632 | -0.000899434 | SUPPORTED |
| A3-A6 | 0/16 | -0.00120829 | 0.0130478 | REJECTED |
| A2-A7-A8 | 16/16 | -0.000837924 | -0.000335143 | SUPPORTED |
| A2-A5-A8 | 16/16 | 0.00158268 | 0.00202859 | SUPPORTED |

A3-A6 is deliberately retained as a negative numerical localization result: its frozen circle is too close
to action-vertex poles for the 128-node tolerance and has 0/16 valid resolved numerical contours. A7-A8 and
both frozen triples pass. The existential D3 gate (at least one pair and one triple) is
**SUPPORTED**. This establishes localized pole motion, not material damping margin.
