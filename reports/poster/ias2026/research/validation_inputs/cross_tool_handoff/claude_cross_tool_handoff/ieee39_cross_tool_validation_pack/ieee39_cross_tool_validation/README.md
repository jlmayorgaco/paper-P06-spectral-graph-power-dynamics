
# IEEE-39 cross-tool validation pack

Purpose: verify the *same physical case* with three independent stacks before using
the benchmark as evidence for the portfolio-composability theory.

Roles
-----
1. Pure Python project model
   - source of the frozen canonical case;
   - AC power flow and the project dynamic/port model.
2. pandapower
   - independent static AC/topology/Ybus/power-flow cross-check only.
3. ANDES
   - independent DAE/eigenvalue cross-check where the dynamic equations can be
     made equation-equivalent.
   - It must NOT be described as independent validation of the custom GFL if the
     same converter equations are not available.

Important
---------
Do NOT compare the stock `pandapower.networks.case39()` directly with the frozen
ANDES-derived benchmark and call the differences an error.  The shipped
pandapower case39 is the PYPOWER/MATPOWER static case and has a different solved
operating point.  The parity gate translates `configs/ias2026/ieee39_network.json`
into pandapower/PYPOWER format so all stacks solve the same physical case.

Recommended environment:
    ANDES 2.0.0
    pandapower 3.4.x
    numpy/scipy/pandas compatible with both

Commands
--------
From the research repository root:

    python /path/to/ieee39_cross_tool_validation/run_python_reference.py .
    python /path/to/ieee39_cross_tool_validation/run_pandapower_parity.py .
    python /path/to/ieee39_cross_tool_validation/run_andes_dynamic.py .
    python /path/to/ieee39_cross_tool_validation/compare_results.py .

Pass gates
----------
STATIC:
- Ybus max error <= 1e-9 pu
- Vm max error <= 1e-5 pu
- gauge-aligned Va max error <= 1e-5 rad
- generator P max error <= 1e-8 pu
- generator Q max error <= 1e-5 pu
- branch terminal P/Q max error <= 1e-3 MVA preferred; investigate before
  accepting anything >1e-2 MVA

DYNAMIC, equation-equivalent synchronous model:
- RHP count agrees for every tested subset
- critical-band alpha error <= 1e-4 1/s
- critical-band frequency error <= 1e-3 Hz
- whole electromechanical-band nearest-mode mismatch <= 1e-4 preferred

THEORY-SPECIFIC:
- determinant/port identity <= 1e-8
- the full-order interaction requirement must reproduce the frozen boundary
  anatomy; lower-order truncations must not be retrospectively tuned
- branch-sensitivity holdout should preserve sign and ranking when tested in an
  equation-equivalent independent model

Interpretation
--------------
If all gates pass, the strongest justified statement is:

    The network/equilibrium and equation-equivalent synchronous dynamics are
    independently reproduced, while the SG-to-custom-GFL portfolio conclusions
    remain validated by the documented project DAE rather than by ANDES.

This is enough to protect the mathematical framework from a bad static network
or an inconsistent synchronous benchmark without overclaiming converter
independence.
