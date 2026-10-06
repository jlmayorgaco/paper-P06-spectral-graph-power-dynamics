# CDW Africano / PV hosting benchmark plan (not run)

**Status: NOT RUN.** This plan does not assume that the Africano material is
available. Gate E17 (preregistration §2) decides, by searching for all of the
following:
- the source document;
- the feeder data;
- the weak-node metric definition (PVR or other);
- the PV placement rule;
- the voltage limits;
- the hosting-capacity method;
- the reported results.

If any of these is missing, the result is `docs/CDW_AFRICANO_MISSING_INPUTS.md`,
and Phases 18–21 (GOLD-E and GOLD-F) are BLOCKED. No generic feeder is ever
substituted and called a replication.

## Future sequence (only if the gate passes)

1. **Reproduce Africano exactly:** the feeder, loading and source data, with
   provenance frozen (paths and sha256).
2. **Reproduce the original static weak-node ranking** with the original metric
   and no reinterpretation.
3. **Reproduce the hosting-capacity result:** the original placement, limits and
   method.
4. **Compare with static sensitivities:**
   - `dV/dP`, `dV/dQ`, Thevenin impedance;
   - the distinction between frozen-state static sensitivity and the
     equilibrium-relaxed zero-frequency response.
   - It is **not** assumed that the Africano metric equals `K(0)`; that
     equivalence is unproved.
5. **Add a dynamic converter model.** A GFL model on the feeder is justified
   only if its parameters are documented for distribution-level inverters.
6. **Coalitions:** node combinations (pairs, triples) and minimal static failing
   coalitions for voltage and thermal constraints.
7. **Control-only optimization:** inverter Q(V) and PF, with no storage and no
   curtailment in the strongest variant.
8. **Topology-only optimization:** admissible switches, keeping radiality where
   required.
9. **Joint topology and control.**

## Guard rails

- Active-power delivery is kept fixed in the counterfactuals where possible.
- Line flows are not required to stay fixed.
- The objective is physical margins, not cycle scores.
- A static metric is never labelled "dynamic" without a dynamic model.
