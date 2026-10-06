# Q4 — conditional closed-form PI target-root law

For a frozen design of every other device and a nonsingular `Delta_{-i}(s_b)`, removing the proportional/integral action of GFL `i` leaves a rank-one return term. The matrix determinant lemma gives the exact conditional root equation

`Kp_i g_p_i(s_b) + Ki_i g_I_i(s_b) = exp(s_b tau_i)`.

Writing `g_p=a+jb`, `g_I=c+jd`, and `s_b=-sigma+j omega` yields the real 2-by-2 system documented in the experiment prompt. If `ad-bc != 0`, its unique conditional gain pair follows by matrix inversion. This is an exact identity for a prescribed target `s_b`, the chosen frozen remainder, and the repository's delayed input channel. It does not say that the target is the rightmost pole or that the pair obeys hardware/design bounds.

Explicitly, for `g_p=a+jb`, `g_I=c+jd`, and `s_b=-sigma+j omega`,

`[a c; b d] [Kp_i; Ki_i] = exp(-sigma tau_i) [cos(omega tau_i); sin(omega tau_i)]`.

If `Delta_g=ad-bc != 0`,

`Kp_i = exp(-sigma tau_i) (d cos(omega tau_i)-c sin(omega tau_i))/Delta_g`,

`Ki_i = exp(-sigma tau_i) (a sin(omega tau_i)-b cos(omega tau_i))/Delta_g`.

The 45 checked combinations (five buses, three delays, and three target frequencies) solved the conditional equations with maximum gain-equation residual `2.29e-16` and full characteristic normalized residual `1.02e-16`. The root solve was initialized at the prescribed target; zero target-coordinate error therefore verifies the constructed conditional equation, not an independent nearest-pole search. The run found only 5/45 gain pairs inside both frozen gain intervals; all 45 are listed in `TABLE_Q04_PI_CLOSED_FORM_VALIDATION.csv`. Thus this experiment has not validated a full stability-boundary curve or a useful feasible tuning law.

Across this deliberately small correlated grid, Spearman correlation between `sigma_min(G_i)` and the Euclidean norm of the calculated gain pair was approximately `-0.951` (Pearson `-0.525`). Treat this only as exploratory description of these 45 points. The sample does not establish a general PI-authority predictor, independent statistical evidence, or causality. Conditioning, positivity, gain bounds, and whether the nearest pole is the target remain separate requirements.
