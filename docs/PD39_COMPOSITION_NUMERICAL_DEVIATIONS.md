# PD39 follow-up deviations and implementation notes

This note records deviations from the initially drafted short follow-up
protocol.  It was updated before the successful numerical audit completed.

1. A blanket rerun of the 189 selected portfolio/scenario rows was started
   and stopped after compilation overhead became clear.  It produced no new
   result file.  The successful campaign therefore joins the selected set to
   the frozen discovery rows for the nine existing scenarios and reruns only
   the new numerical checks.  This preserves the discovery record and avoids
   presenting repeated calculations as independent evidence.

2. The installed Julia/LinearAlgebra LAPACK interface has no `eigen!` method
   for the requested `BigFloat` matrix keyword combination.  The BigFloat
   arithmetic check is consequently recorded as unavailable.  Float32,
   Float64, complex Float64, and descriptor generalized-eigenvalue checks
   remain separately reported.

3. Perturbations around the high-PLL corner lie on the boundary of the
   discovery controller box.  The prescribed +/-0.5% finite perturbations
   were evaluated in the preregistered secondary physical box so both signs
   remain admissible; no perturbation magnitude was expanded.

4. The numerical solver sweep changes `tol`, `abstol`, and `reltol` together.
   The component initializer's own residual is retained in every row.  A
   loose 1e-8 setting is therefore not silently treated as equivalent to the
   converged 1e-10/1e-12 settings.

No holdout, radius, weak-component, co-design, topology, or TDS campaign was
added to this short follow-up.

