# NEGATIVE RESULTS

(dated entries; nothing deleted)

## 2026-10-05 process corrections (kept, not deleted)
- Phase C, first pass: all-SG branch tracking up from rho=0.125 collapsed several modes onto the same Padé-catalogue root (duplicates, MAC ~0). Fix: a branch is kept only if its initial SG-subspace MAC >= 0.7 and its seed root is unused; the rerun lists no unmatched modes. PLL-family branches PLL2 and PLL5 had min MAC 0.76 / 0.64 in the rho continuation (possible local branch exchange between the 4.9 Hz cluster members); their continuation is therefore reported as a family, not branch by branch.
- Phase D, first pass (0.25 ms steps): the 4.649 Hz root and the 4.681 Hz root merged when tracked back to 40 ms (duplicate branches R0=R1), which hid one crossing. Fix: 0.05 ms back-steps plus a duplicate guard; corrected crossings are in TABLE_D01b.
- The 7-8 Hz family of the Python reduced model does NOT appear in the full model: the PLL/controller family is at 4.5-5.0 Hz (see Phase B/C/D). Frequency was not a required match.
- Phase E: local factors are regular but not large: min |L_i| = 0.011-0.034 at the critical roots (max 0.24-0.28); they never approach 0. The closure is collective (eig(Q) = -1), as required by definition at an exact root.
- Phase F: with the preregistered edge threshold (0.1 max |Q|) the interaction graph is almost strongly connected (SCC of 7-10 buses); localisation appears in the pair cycles, not in the SCC partition.

## 2026-10-05 engine v1 -> v2 (kept: derived/engine_v1/)
Engine v1 accepted every trust-region step (no acceptance test, last iterate returned). Its results are preserved in derived/engine_v1/: S0 Ki-only/Kp-only/joint all failed (worst root +0.69/+0.66/+0.67 s^-1), core4_ki failed, sens4_joint failed; all-ten nodal (Ki-only, joint) and uniform joint repaired in 2-4 iterations. Because v1 sometimes ended on a worse iterate (e.g. core scan k=6, joint S0 oscillation), AMENDMENT_04 introduces v2 (step acceptance, best-iterate return) and ALL sparse results are rerun with it; v2 results supersede v1 only as engine-quality, both are reported.
