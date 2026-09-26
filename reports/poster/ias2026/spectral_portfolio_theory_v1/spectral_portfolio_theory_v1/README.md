# Spectral portfolio certification — research package v1

Prepared 10 September 2026 for the Vancouver IAS research program.

## Files

- `TEORIA_Y_DEMOSTRACIONES.md`: Spanish technical note, literature context,
  eight main propositions/theorems, proofs, limits and experiment plan.
- `TEORIA_Y_DEMOSTRACIONES.pdf`: 24-page typeset copy of the note.
- `PROMPT_CLAUDE_CODE.md`: standalone implementation and validation contract for
  the existing repository. No permission to overwrite freezes or claim new
  IEEE results without running them.
- `theory_checks.py`: deterministic NumPy/SciPy/SymPy synthetic regression bank.
- `results/theory_check_results.json`: results from this session: 404 checks.
- `sources_manifest.json`: exact hashes of the uploaded evidence and literature
  identifiers used.

## Run synthetic checks

```bash
python theory_checks.py --out results
```

Dependencies: numpy, scipy, sympy. These checks do not need ANDES or an SDP
solver. They do NOT implement the Boolean SOS optimizer; the explicit 2x2
Boolean certificate is verified symbolically. The low-cost Lyapunov-subcube
search is implemented and compared with exhaustive enumeration on five small
synthetic affine families.

## Rebuild the PDF

Run `./build_pdf.sh` with pandoc, XeLaTeX, and the referenced system fonts installed. No font files are bundled.

## Important interpretation

The mathematical proofs establish conditional statements, not priority or
industrial validation. The proposed certificates may be conservative or fail
applicability on some SG–IBR port representations. `UNKNOWN` is never treated
as safe or unsafe. The package does not contain the user's full research
repository or rerun the IEEE-39/Kundur/68-bus experiments.

Two immediate audit concerns are explicit: arbitrary deletion of every
near-zero eigenvalue, and possible MW/MVA confusion in an older mitigation
table. Both require tracing the source model, not silently rewriting results.
