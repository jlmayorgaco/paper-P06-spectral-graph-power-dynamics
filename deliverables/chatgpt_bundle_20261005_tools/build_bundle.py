"""Organise repository results into one ChatGPT-readable bundle (text-first, size-capped).
Usage: python build_bundle.py [--dry]"""
import os, sys, csv, hashlib, shutil, zipfile, pathlib, fnmatch
ROOT = pathlib.Path(__file__).resolve().parents[2]
OUT = ROOT / 'deliverables' / 'chatgpt_bundle_20261005'
DRY = '--dry' in sys.argv
CAP = {'.csv': 120e3, '.json': 120e3, '.toml': 200e3, '.md': 600e3, '.txt': 600e3, '.tex': 800e3,
       '.py': 120e3, '.jl': 120e3, '.png': 450e3, '.svg': 400e3, '.pdf': 2.5e6, '.bib': 300e3,
       '.yaml': 100e3, '.yml': 100e3, '.sha256': 50e3, '.ps1': 100e3, '.cff': 50e3}
SKIP_DIRS = {'.git', '__pycache__', 'build', 'node_modules', '.julia', '.venv', 'texmf-cache', '$qaDir', '.codex-remote-attachments',
             'chatgpt_bundle_20261005', 'chatgpt_bundle_20261005_tools'}  # never copy our own output (it lives under deliverables/)
SKIP_PAT = ['*_SENSORS.csv', '*.aux', '*.out', '*.nav', '*.snm', '*.toc', '*.synctex*', 'state_*.npz', '*_TRACE*.csv']
# category -> list of (source relative to ROOT, optional per-source depth note)
CATS = {
 '01_papers_series': ['reports/papers', 'reports/shared/bib'],
 '02_portfolios_incompatibility_TX3_TX4': ['reports/poster/ias2026/spectral_portfolio_theory_v1', 'results/pd39', 'artifacts/tx3', 'experiments/tx3',
     'docs/PD39_255PLUS1_FINAL_REPORT.md', 'docs/PD39_BLOCKER_ATLAS_FINAL_REPORT.md', 'docs/PD39_CONFIRMATORY_FINAL_REPORT.md', 'docs/PD39_RESULTS_REPORT.md',
     'docs/TX4_ROBUSTNESS_FINAL_REPORT.md', 'docs/TX4_FINAL_MODAL_SCOPE_REPORT.md', 'docs/TX4_BLIND_PREDICTION_FINAL_REPORT.md',
     'docs/TX4_CONTEXTUAL_RETURN_THEOREM.md', 'docs/TX4_IAS_FINAL_POSTER_CONTENT.md', 'docs/TX4_EXACT_P4_JULIA_CURRENT_STATE.md'],
 '03_weak_nodes_CDW_algebraic': ['reports/poster/ias2026/research/contextual_dynamic_weakness', 'reports/papers/cdw_contextual_dynamic_weakness',
     'reports/papers/tx2_when_is_weak_a_graph_object', 'temp/paper_ieee_conference_weak_elements'],
 '04_BND_H4_mechanism_and_research_line': ['reports/poster/ias2026/research/bnd_h4_mechanism', 'reports/poster/ias2026/research/README.md',
     'reports/poster/ias2026/research/docs', 'reports/poster/ias2026/research/theory', 'reports/poster/ias2026/research/results'],
 '05_codesign_campaigns_experiments': ['experiments'],
 '06_codesign_campaigns_reports_AtoQ': [f'reports/experiment_{x}' for x in 'A B C D E F_program F0 F1 F2 F3 F4 F5 F6 F7 G H J0 K M N P Q Q2 Q2B'.split()] +
     ['reports/analytic_iteration_20261001', 'reports/codesign_audit_20261001', 'reports/codesign_validation_20261001', 'reports/nonlinear_codesign_20261001',
      'reports/nonlinear_formulation_20261001', 'reports/retuning_boundary_20261002'],
 '07_poster_sources_and_outputs': ['reports/poster/ias2026/collective_interaction_20261003', 'reports/poster/ias2026/collective_interaction_20261004',
     'reports/poster/ias2026/sections', 'reports/poster/ias2026/main.tex', 'reports/poster/ias2026/generated', 'reports/poster/ias2026/SECTIONS.md',
     'reports/poster/ias2026/DESIGN_SYSTEM.md', 'reports/poster/ias2026/VISUAL_REVIEW.md', 'output/pdf', 'reports/conference'],
 '08_handoffs_and_state': ['reports/handoff', 'README.md', 'CITATION.cff', 'docs/architecture.md', 'docs/reproducibility.md', 'reports/poster/ias2026/AGENTS.md',
     'docs/PD39_RUNBOOK.md', 'docs/TX4_ROBUSTNESS_CHATGPT_HANDOFF.md', 'docs/TX4_BLIND_PREDICTION_CHATGPT_HANDOFF.md', 'docs/TX4_FINAL_MODAL_SCOPE_CHATGPT_HANDOFF.md',
     'docs/PD39_255PLUS1_CHATGPT_HANDOFF.md', 'docs/PD39_BLOCKER_ATLAS_CHATGPT_HANDOFF.md'],
 '09_model_source_code': ['src', 'Project.toml', 'pyproject.toml', 'Makefile', 'scripts'],
 '10_deliverable_reports_pdf': ['deliverables'],
}
def skipped(p):
    return any(fnmatch.fnmatch(p.name, pat) for pat in SKIP_PAT)
def walk(src):
    if src.is_file(): yield src; return
    for dp, dn, fn in os.walk(src):
        dn[:] = [d for d in dn if d not in SKIP_DIRS]
        for f in fn: yield pathlib.Path(dp) / f
rows, skippedrows, total = [], [], 0
if not DRY and OUT.exists(): shutil.rmtree(OUT)
seen = set()
for cat, srcs in CATS.items():
    for s in srcs:
        src = ROOT / s
        if not src.exists(): skippedrows.append((cat, s, 'MISSING_SOURCE', 0)); continue
        for p in walk(src):
            rel = p.relative_to(ROOT)
            if p in seen: continue
            ext = p.suffix.lower(); size = p.stat().st_size
            if skipped(p): skippedrows.append((cat, str(rel), 'pattern', size)); continue
            if ext not in CAP: skippedrows.append((cat, str(rel), 'ext', size)); continue
            if size > CAP[ext]: skippedrows.append((cat, str(rel), 'too_large', size)); continue
            seen.add(p); total += size
            rows.append((cat, str(rel), size))
            if not DRY:
                d = OUT / cat / rel; d.parent.mkdir(parents=True, exist_ok=True); shutil.copy2(p, d)
print(f'files={len(rows)} total_MB={total/1e6:.1f} skipped={len(skippedrows)}')
from collections import Counter
c = Counter(); [c.update({r[0]: r[2]}) for r in rows]
for k, v in c.items(): print(f'  {k}: {v/1e6:.1f} MB')
sk = Counter(r[2] for r in skippedrows); print('skip reasons', dict(sk))
if not DRY:
    with open(OUT / 'MANIFEST.csv', 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f); w.writerow(['category', 'repo_path', 'bytes', 'sha256'])
        for cat, rel, size in rows:
            w.writerow([cat, rel, size, hashlib.sha256((OUT / cat / rel).read_bytes()).hexdigest()])
    with open(OUT / 'SKIPPED_FILES.csv', 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f); w.writerow(['category', 'repo_path', 'reason', 'bytes']); w.writerows(skippedrows)

if not DRY:
    for f in pathlib.Path(__file__).parent.glob('0*.md'):
        shutil.copy2(f, OUT / f.name)
    with zipfile.ZipFile(ROOT/'deliverables'/'BND_IEEE39_ChatGPT_Bundle_20261005.zip','w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for fp in sorted(OUT.rglob('*')):
            if fp.is_file(): z.write(fp, 'BND_ChatGPT_Bundle/'+str(fp.relative_to(OUT)).replace(os.sep,'/'))
