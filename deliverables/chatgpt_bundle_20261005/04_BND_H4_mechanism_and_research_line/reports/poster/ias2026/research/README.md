# IAS AM 2026 — controller-mediated interaction between grid actions

Research framework behind the IEEE IAS Annual Meeting 2026 poster in the parent
directory. It asks one question:

> Is the dynamic security of a portfolio of interventions inferable from the
> dynamic security of its individual components or lower-order subsets?

and it is built so that the answer can be *no* without that being an artefact.

## One-command reproduction

```powershell
# from this directory
..\..\..\..\.venv\tx3-analysis\Scripts\python.exe -m pytest tests -q
..\..\..\..\.venv\tx3-analysis\Scripts\python.exe experiments\E00_regression_toy.py
```

Anything that is an algebraic identity is checked to machine precision; anything
empirical carries an uncertainty budget. `experiments/E00_regression_toy.py`
exits nonzero if GATE 0 regresses.

## Layout

```
src/ibr_cycles/
  dynamics/     DAE, equilibrium certification, Jacobians, index-1 reduction,
                eigen-analysis, MAC-based modal tracking
  reduction/    Sigma(s), Schur determinant identity, pole/residue expansion
  actions/      low-rank factorization, Action Green operator K(s),
                individual/collective split Q(s), portfolio bookkeeping
  cycles/       action graph and SCCs, holonomies H_gamma, connected cumulants,
                finite-amplitude Moebius, winding with admissibility test
  diagnosis/    mode provenance, minimum-norm repair
  models/       toy5 (frozen regression case), ieee9_gfl, ieee39
  io/           run manifests (git commit, seeds, versions, config hash)
experiments/    E00 ... E20, each runnable from the CLI
tests/          the regression contract; every identity is a hard gate
results/        raw/ tables/ figures/ manifests/ reports/
docs/           METHODS.md, CLAIMS.md, EXPERIMENTS.md, RESULTS.md,
                FAILED_EXPERIMENTS.md, POSTER_RESULTS.md
```

## Gates

| Gate | Condition | Status |
|---|---|---|
| 0 | five-state regression envelope reproduces | **PASSED** (max deviation `4.2e-07`) |
| 1 | IEEE-9 GFL regression | **PASSED** (frozen `<4.6e-07`; outside envelope corroborated to `5.3e-4`) |
| 2 | Schur/Sigma identity `< 1e-10` | **PASSED** (`8.2e-14`) |
| 3 | action determinant identity `< 1e-8` | **PASSED** (`3.8e-12`) |
| 4 | a realistic IEEE-39 portfolio shows a robust collective effect | not started |
| 5 | provenance winding converges under contour refinement | **PASSED on toy**, admissibility enforced |
| 6 | negative controls behave correctly | not started |
| 7 | frozen flagship survives held-out operating-point uncertainty | not started |
| 8 | targeted repair beats asset removal and global retuning | toy only |

Gates 4 and beyond require IEEE-39 and are not authorized to produce poster
wording until Gate 1 and the negative controls exist.

## What GATE 0 established

At the frozen action amplitudes of the five-state case: every single action is
stable, every pair is stable, the triple is unstable at `Re(lambda) = +0.046781`,
`Im = 2.212347`, so `kappa = 3`.

The instability is **not** an accumulation of pairwise effects. Truncating the
Moebius decomposition of the spectral abscissa at order two predicts `-0.023223`
(stable); the irreducible third-order term of `+0.070005` is what crosses the
axis.

On a region isolating the created mode, the winding splits as full `+1`,
individual `0`, collective `+1` — the collective factor carries the new mode.
That statement is conditional: a larger region that also encloses the
single-action mode of A flips the split to individual `1`, collective `0`, and
the framework reports `INADMISSIBLE_REGION` instead of a provenance label.

A controller retune to `k = 1.259557` restores a `-0.02` abscissa without
removing either network reinforcement.

## What GATE 1 established

The same structure appears on a converter-dominated IEEE-9 case with a real AC
operating point: base `-0.105890`, every single and every pair stable, the triple
unstable at `+0.003217`, `1.3577 Hz`. Cross-checked by the argument principle on
`Re(s) > -0.002` with a clearance-tied step: proper subsets `0`, portfolio `+2`.
At the created eigenvalue the collective factor vanishes (`2.2e-12`) while the
individual factor does not (`0.0489`).

**And one result that cuts against the intended narrative.** The Moebius
decomposition says the IEEE-9 crossing is driven by accumulated *pairwise*
interaction, not by an irreducible three-way term: the order-2 truncation is
already unstable at `+0.004023` and the genuine triple term is `-0.000806`, an
order of magnitude smaller and stabilizing. The toy case is the opposite. Same
`kappa`, different mechanism. `E01` emits `mechanism_verdict` so this is decided
by the data on every future case. See `docs/RESULTS.md` R3.

## Read before extending

`docs/CLAIMS.md` records the relationship with the frozen TX3 envelope. TX3
measured finite re-equilibrated externalities on the matched IEEE-39 / three-GFL
benchmark and found them **immaterial** on the margin-setting mode (C3a and C3b
rejected, C3c unresolved, cycle/SCC claim unassessed). No TX3 number may be
reused to support a cycle, provenance or repair claim here.

`docs/METHODS.md` §10 records the open modelling decision between a frozen and a
re-equilibrated operating point. Every IEEE-39 number must say which one produced
it.
