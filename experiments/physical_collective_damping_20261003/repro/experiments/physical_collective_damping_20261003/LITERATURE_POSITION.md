# Novelty boundary, checked 2026-10-03

This experiment does not establish novelty by itself. The following primary
sources already address important neighboring claims:

- Huang, Xin, Dong and Doerfler, *Impacts of Grid Structure on PLL-Synchronization
  Stability of Converter-Integrated Power Systems*, [arXiv:1903.05489](https://arxiv.org/abs/1903.05489).
  Their abstract reports Kron reduction, network modal decomposition, grounded
  Laplacian stability margins, sensitivities, PLL retuning and IEEE-39 validation.
  Therefore graph coupling plus PLL retuning on IEEE-39 is not a new claim here.
- Zhao et al., *Corrected complex torque coefficient method for synchronization
  stability analysis of grid-connected VSC*, IET Renewable Power Generation 18(3),
  412--425 (2024), online September 2023,
  [doi:10.1049/rpg2.12857](https://ietresearch.onlinelibrary.wiley.com/doi/10.1049/rpg2.12857).
  The paper identifies a limitation of classical torque-based damping criteria
  due to the PLL proportional term and derives a corrected criterion for its
  single-converter model. This directly cautions against claiming that a scalar
  damping diagnostic universally determines PLL stability. Our port calculation
  is not a reproduction of its corrected coefficient, nor a demonstrated
  generalization of its theorem.

The executed new repository contribution is a physical-port diagnostic and
nonlinear pure-delay validation for this specific heterogeneous IEEE-39 model.
A publishable novelty claim still needs a useful, validated prediction/design
law beyond these established results, and evidence that it improves constrained
replacement rather than only retuning at fixed rho.
