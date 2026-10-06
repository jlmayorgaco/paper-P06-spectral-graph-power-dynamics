"""Integrate the reviewed section into the existing open source, in place."""
from pathlib import Path
import hashlib
import json

OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[1]
target=ROOT/'experiments/theory_collective_damping_20261003/THEORY.tex'
text=target.read_text(encoding='utf-8')
backup=OUT/'THEORY_before_extension.tex'
if not backup.exists():
    backup.write_bytes(target.read_bytes())
start='% BEGIN NONLINEAR_FREQUENCY_LIMITS_20261004'
end='% END NONLINEAR_FREQUENCY_LIMITS_20261004'
section=(OUT/'FREQUENCY_LIMITS_SECTION.tex').read_text(encoding='utf-8')
block=start+'\n'+section+'\n'+end+'\n\n'
if start in text:
    a=text.index(start)
    b=text.index(end,a)+len(end)
    text=text[:a]+block.rstrip()+text[b:]
else:
    anchor=r'\section*{Reproducibility and research integrity}'
    assert text.count(anchor)==1
    text=text.replace(anchor,block+anchor)
paragraph=r"""
A separate nonlinear development in Section~\ref{sec:nonlinear-fixed-bus}
uses fixed measurement buses and arbitrary fixed events. It derives a
total-inertia obstruction from the two power-flow endpoints and proves a
sharp two-node threshold with an explicit feasible control above it.
This lossless swing result has a different model contract from the
detailed delayed IEEE--39 examples. An audit of the archived converter
energy sign and governor limits identifies why transfer to that full
model remains an explicit gate.
"""
introanchor='design. The ultimate maximum-replacement problem remains open.'
if paragraph.strip() not in text:
    assert text.count(introanchor)==1
    text=text.replace(introanchor,introanchor+'\n'+paragraph)
bibs=r"""
\bibitem{poollainertia} B. K. Poolla, S. Bolognani, and F. D\"orfler,
\emph{Optimal Placement of Virtual Inertia in Power Grids},
IEEE Transactions on Automatic Control,62(12):6209--6220,2017.
\url{https://doi.org/10.1109/TAC.2017.2703302}.
\bibitem{kohlhaas2026} S. Kohlhaas and P. Kotyczka,
\emph{Structured Modal-Energy Dissipation for Inter-Area Oscillations in
Low-Inertia Power Systems}, RE Grid Integration Week,2026, contribution111.
Official abstract consulted; full paper not obtained.
\url{https://coms.events/regridweek2026/data/abstracts/en/abstract_0094_0111_0653.html}.
\bibitem{jiangfs} Y. Jiang, W. Chen, Z. Lyu, X. Zhang, D. Wang, and S. Hara,
\emph{Frequency Shaping Control for Oscillation Damping in Weakly-Connected
Power Network: A Root Locus Method}, arXiv:2601.19665,2026.
\url{https://arxiv.org/abs/2601.19665}.
"""
if r'\bibitem{poollainertia}' not in text:
    text=text.replace(r'\end{thebibliography}',bibs+'\n'+r'\end{thebibliography}')
eol='\r\n' if b'\r\n' in backup.read_bytes() else '\n'
target.write_bytes(text.replace('\r\n','\n').replace('\n',eol).encode('utf-8'))
log={'target':str(target),'before_sha256':hashlib.sha256(backup.read_bytes()).hexdigest(),
 'after_sha256':hashlib.sha256(target.read_bytes()).hexdigest(),
 'section_sha256':hashlib.sha256(section.encode()).hexdigest(),
 'scope':'New reviewed nonlinear fixed-bus section, introductory pointer, three verified references. Existing results preserved.',
 'compiled':False}
(OUT/'MANUSCRIPT_INTEGRATION.json').write_text(json.dumps(log,indent=2),encoding='utf-8')
print(json.dumps(log,indent=2))
