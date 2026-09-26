# Gate 5 — held-out Monte Carlo: FAILED as frozen

Protocol `configs/ias2026/trackA_final_v1.yaml` (v1.1). 400 samples drawn,
283 accepted, 29.2 % rejected. Minimum MAC over every accepted sample and
configuration 0.811, above the 0.80 floor.

## The frozen hypotheses do not hold

| hypothesis | fraction | 95 % CI | required lower bound |
|---|---|---|---|
| H1 `d_alpha_Q > 0` — setpoint replacement degrades the branch | **0.106** | [0.074, 0.141] | 0.80 |
| H2 `d_alpha_V < 0` — voltage control mitigates | **0.099** | [0.067, 0.134] | 0.80 |
| H3 ordering unityPF > Qset > Vctrl | 0.099 | [0.067, 0.134] | — |

Median tracked-branch real part over accepted samples:

| configuration | median | 5th | 95th | unstable |
|---|---|---|---|---|
| reference | `-0.1253` | `-0.1359` | `-0.1117` | 0.0 % |
| Q setpoint | `-0.2267` | `-0.6987` | `+0.2599` | **9.2 %** |
| unity power factor | `-0.2286` | `-0.6845` | `+0.6121` | 21.9 % |
| voltage control | `-0.1116` | `-0.1329` | `-0.0829` | 0.0 % |

Median `d_alpha_Q = -0.103`. **At a typical operating point drawn from the frozen
distribution the four-replacement portfolio makes the tracked inter-area branch
MORE stable, not less.** The discovery operating point sits in roughly the
tenth percentile.

**The unconditional Track-A claim is not supported and is withdrawn.**

## Where the failure actually lives — post-hoc, not preregistered

Everything below was computed after seeing the held-out data. It is diagnosis of
a failed gate, it is labelled as such, and it does **not** license a claim. Any
of it that is to be asserted needs its own frozen protocol and its own held-out
campaign.

Rank correlation of `d_alpha_Q` with the sampled variables:

| variable | Spearman | p |
|---|---|---|
| PV availability | `-0.348` | `1.7e-09` |
| global active load | `-0.339` | `5.0e-09` |
| global reactive load | `+0.045` | `0.45` |

The marginal correlations are negative, but the effect is not monotone: it
concentrates in the corner where **both** are high, which is where the mechanism
was discovered.

| condition | n | H1 fraction | median `d_alpha_Q` | Q-setpoint unstable |
|---|---|---|---|---|
| all accepted | 283 | 0.106 | `-0.103` | 9.2 % |
| availability ≥ 0.90 | 102 | 0.157 | `-0.155` | 13.7 % |
| availability ≥ 0.95 | 50 | 0.240 | `-0.156` | 22.0 % |
| availability ≥ 0.90 and load ≥ 1.00 | 26 | 0.423 | `-0.280` | 38.5 % |
| **availability ≥ 0.95 and load ≥ 1.00** | **13** | **0.615** | **`+0.289`** | **61.5 %** |

The instability rate rises monotonically toward the discovery corner, from 9 % to
62 %. With thirteen samples in the last row this is a direction, not a number.

## The mitigation, conditional on the failure occurring

Among the 26 accepted samples where the Q-setpoint case is genuinely unstable:

| quantity | value |
|---|---|
| voltage control restores stability | **26 of 26** |
| samples still unstable under voltage control | **0** |
| median `d_alpha_V` | `-0.362` |
| `alpha_V < alpha_Q` | **26 of 26** |

So H2 failed as frozen for a formulation reason, not a physical one. As written,
H2 asks whether voltage control always lowers the branch; when the setpoint case
is already comfortably stable it instead brings the branch back toward the
reference, which counts as a failure without anything having gone wrong. The
mitigation claim was always conditional, and conditionally it holds in every
sample where it is asked. This too is post-hoc and needs a v2 protocol with H2
stated conditionally before it may be claimed.

## What survives and what does not

Withdrawn:

- any statement that replacing these four machines destabilizes the system at a
  generic operating point. It does not; it usually helps.

Untouched by this gate, because they are properties of a specified operating
point and were validated on their own terms:

- the order-4 audit at the discovery point (E15): orders 1 to 3 stable, exact
  unstable, `mu4 = +0.219`, minimum MAC 0.967;
- the exact 8x8 port representation and the invariant closure (E16, E17);
- the closure continuation `1.000 → 0.55-0.82 → 0.21-0.42 → 1.6e-07` (E23);
- the independent time-domain confirmation of the unstable case to
  `4.4e-04 Hz` and `1.7e-03 s^-1` (E25);
- the negative controls, including that the failure is not controller-caused
  (E14, E14B).

## Honest form of the claim now

> At a specified high-load, high-output operating condition, four individually
> and lower-order compatible SG to PV-GFL replacements operating on a reactive
> power setpoint drive a pre-existing inter-area branch unstable, and the
> irreducible fourth-order interaction is required for the crossing. Over a
> held-out sweep of operating conditions the same portfolio is unstable in about
> nine percent of cases and is otherwise stabilizing, so the phenomenon is a
> stress condition rather than a generic consequence of replacement.

That is narrower than the poster message assumed and it is what the evidence
supports.

## Recommended v2, for the user to decide

1. Restate the nominal case as a **stress condition**: PV availability at or near
   full output and active load at or above nominal, declared before sampling.
2. State H2 **conditionally**: among samples where the setpoint case is unstable,
   voltage control restores stability.
3. Sample within the stress condition, not across the whole operating envelope,
   and report the conditional probabilities with confidence intervals.
4. Fix the dispatch model: 29 % of samples were rejected because a machine
   slightly exceeded its rating once generation was scaled with load. A
   proportional dispatch respecting headroom would cut that.
