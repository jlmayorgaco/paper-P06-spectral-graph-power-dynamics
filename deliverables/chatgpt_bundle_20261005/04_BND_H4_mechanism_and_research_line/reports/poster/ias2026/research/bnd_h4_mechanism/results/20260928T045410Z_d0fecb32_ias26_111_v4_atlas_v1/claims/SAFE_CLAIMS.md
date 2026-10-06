# IAS26-111: claims supported within scope

Status: PASS_WITH_LIMITATIONS. This is a post-hoc descriptive reanalysis of the frozen synthetic IAS26-060 run. No power-flow, DAE, or eigenvalue solve was executed.

- Among 967 physically valid scenarios in the declared synthetic ensemble, 412 had no inclusion-minimal unstable V4 portfolio, 190 had H4={30,33,35,37} as the sole minimal blocker, and 365 had one or more alternative minimal triple blockers.
- Across the same 967 scenarios, all 4 singleton portfolios and all 6 pair portfolios were stable; there were no indeterminate portfolio classifications at tau_dec=1e-8.
- The minimum V4 blocker order was 3 for 365 scenarios, 4 for 190, and undefined for 412. This is not a full-network/V9 result.
- The observed minimal triples were {30,33,35} (365 scenario memberships), {30,33,37} (288), and {33,35,37} (137). All observed minimal triples include bus 33; this is descriptive and not evidence of causality.
- x_s=max over proper H4 subsets alpha_perp, y_s=alpha_perp(H4), and m_comp=min(y_s,-x_s) are reported per valid scenario. m_comp is a signed coordinate margin in the (x,y) plane, not a probability or a universal physical stability margin.
- Recomputed minimal masks and boundary coordinates reconcile to the frozen SCENARIO_METRICS table for all 967 valid IDs with zero mismatches. All 1,000 IDs remain visible; 33 physically invalid IDs were not replaced and were not admitted to the valid-cohort counts.

All rates are descriptive frequencies under this specific synthetic ensemble, not real-world probabilities.
