# F2 preregistration

The graph basis is fixed by the F0 primary `Lc = Re(Yport)` before inspecting
GFL closed-loop spectra. Exact modal separation requires both real and
imaginary parts of the passive port admittance to be diagonal in this basis.
We classify residual off-diagonal modal coupling at or below 5% of the full
operator Frobenius norm as `APPROX_MODAL_WITH_SMALL_COUPLING`; larger coupling
is `NON_MODAL`. This tolerance is a structural gate, not a controller-fit
criterion. The heterogeneous comparison uses the frozen analytic IEEE-39
equilibrium and nominal GFL gains; it imports no PowerDynamics module.
