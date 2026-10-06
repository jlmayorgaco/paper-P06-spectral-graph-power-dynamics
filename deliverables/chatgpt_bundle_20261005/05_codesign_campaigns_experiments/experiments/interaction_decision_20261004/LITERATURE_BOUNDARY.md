# Closest verified literature and fair scope

1. Grosdidier and Morari, Interaction measures for systems under decentralized
control, Automatica22(3),309–319,1986, DOI10.1016/0005-1098(86)90029-4.
Primary exposition: https://thesis.caltech.edu/1014/1/Grosdidier-P_1986.pdf,
ChapterIII. Normalized off-diagonal error, interaction determinants and signed
decentralized-controller obstructions are established. Our algebra is not a
new general result in decentralized control.

2. Chen, Tan, Javaid, Angeli and Chaudhuri, Decentralized Stability of
IBR-dominated Power Grids Using Block Diagonal Dominance,2026 preprint:
https://arxiv.org/html/2606.28023v1 . Full block rows include prospective
multi-IBR interaction; sufficient certificates include decay rate. Their BDD
criterion is not the singleton-additive modal surrogate used here. Failure
of BDD would be inconclusive. Its rational-LTI formulation requires justification
before adaptation to pure DDE or approximation-error control for fitted models.

3. Huang, Xin, Dong and Dörfler, Impacts of Grid Structure on PLL-Synchronization
Stability of Converter-Integrated Power Systems:
https://arxiv.org/html/1903.05489v2 . Kron reduction, network-modal decomposition,
PLL retuning and IEEE39 are antecedents. Detailed homogeneous-dynamics/network
assumptions differ from our fixed exported heterogeneous full-device model.

4. Michiels and Gumussoy, fixed-order stabilizing controllers for interconnected
systems with time delays: https://arxiv.org/abs/2003.05496 . Spectral optimization
of prescribed controller structures in DDAEs is established. No comparison
against its implementation has been executed in this campaign.

5. Borgioli et al., real stability radius for time-delay systems:
https://www.unige.ch/math/vandereycken/papers/published_Borgioli_MLV.pdf .
Structured robust-distance problems have stronger global guarantees in their
declared domains. They are not equivalent to SG-replacement maximization.

6. Bindel and Hood, nonlinear eigenvalue localization:
https://www.cs.cornell.edu/~bindel/papers/2013-simax.pdf . Contour counting and
localization are antecedents, not discoveries claimed here.

Capability added in this execution: a finite exact-delay example in which two
actual SG-replacement/local-gain actions individually retain a decay margin,
but their joint model provably violates it; a scalar interaction correction
recovers the margin and passes the fixed nonlinear tests. Both the certificate
and its application scope matter. A world-first or journal-ready claim remains
unsupported; maximum replacement is still open.
