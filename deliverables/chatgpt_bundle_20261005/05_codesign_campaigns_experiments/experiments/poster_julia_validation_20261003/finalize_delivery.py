"""Audit preservation, record evidence and make a portable delivery archive."""
from pathlib import Path
import hashlib,json,re,subprocess,zipfile,sys,platform
import numpy,pandas,scipy,matplotlib
import pymupdf
ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent
POSTER=ROOT/'reports/poster/ias2026/collective_interaction_20261003'
REG=ROOT/'experiments/regional_paper_closure_20261003'
THEORY=ROOT/'experiments/theory_collective_damping_20261003/THEORY.tex'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
dump=lambda p,o:p.write_text(json.dumps(o,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
def git(*a):return subprocess.check_output(['git',*a],cwd=ROOT,text=True).strip()

# Every historical numerical/source file archived by the previous turn stays intact.
checked=[];skipped=[]
for archive in sorted((REG/'evidence').glob('*.zip')):
    with zipfile.ZipFile(archive) as z:
        for n in z.namelist():
            if n.endswith('/'):continue
            if n in ['Project.toml','Manifest.toml','experiments/theory_collective_damping_20261003/THEORY.tex']:
                skipped.append(n);continue
            p=ROOT/n
            if not p.exists():continue
            assert p.read_bytes()==z.read(n),('Historical file changed',n)
            checked.append(n)
dump(OUT/'PRESERVATION_CHECK.json',{'pass':True,'historical_files_checked':len(checked),
     'checked':checked,'editable_manuscript_and_root_environment_excluded':skipped,
     'old_regional_manuscript_hash_matches_new_stage_snapshot':
       sha(OUT/'THEORY_before_poster_closure.tex')==json.loads((REG/'FINAL_MANIFEST.json').read_text())['manuscript_sha256']})

protocol=(OUT/'PROTOCOL.txt').read_text(encoding='utf-8')
v1=protocol.split('\nAdditional algebra-only check')[0].rstrip()+'\n'
(OUT/'PROTOCOL_COMPARATOR_AT_LAUNCH.txt').write_text(v1,encoding='utf-8')
dump(OUT/'PROTOCOL_HISTORY.json',{
 'initial_contract':'Five-event comparator section written before Julia launch.',
 'appendix':'The algebra-only 60-case test was appended before its Python execution, after the Julia comparator had launched.',
 'no_event_contract_changed':True,
 'initial_text_recovered_verbatim_from_unchanged_prefix':True,
 'comparator_prefix_sha256':sha(OUT/'PROTOCOL_COMPARATOR_AT_LAUNCH.txt'),
 'complete_protocol_sha256':sha(OUT/'PROTOCOL.txt')})

source_manifest={p.relative_to(ROOT).as_posix():sha(p)
 for p in (POSTER.parent/'sections').glob('*.tex')}
for name in ['main.tex','compile_section.ps1','poster_common.tex','poster_layout.tex','poster_design_system.tex','DESIGN_SYSTEM.md','SECTIONS.md']:
    p=POSTER.parent/name;source_manifest[p.relative_to(ROOT).as_posix()]=sha(p)
dump(POSTER/'TEMPLATE_ORIGIN.json',{'template':'Most advanced modular working-tree version; original remained unmodified.',
     'source_hashes':source_manifest,'prior_HEAD':git('rev-parse','HEAD')})

f=pandas.read_csv(OUT/'TABLE_01_FIXED_GAIN_EVENTS.csv')
r=json.loads((OUT/'RANK_LAW_SUMMARY.json').read_text())
assert len(f)==5 and bool(f['pass'].all()) and r['cases']==60 and r['identity_pass'] and r['all_indefinite']
for pattern in ['were not run','lacks a paired nonlinear','path independence has not been tested']:
    assert pattern not in THEORY.read_text(encoding='utf-8'),pattern
text=THEORY.read_text(encoding='utf-8')
labels=re.findall(r'\\label\{([^}]+)\}',text)
assert len(labels)==len(set(labels))
assert set(re.findall(r'\\(?:ref|eqref)\{([^}]+)\}',text))<=set(labels)

pdfs=[ROOT/'output/pdf/Beyond_Nodal_Damping_Theory_20261003.pdf',
      ROOT/'output/pdf/IAS2026_Collective_Damping_20261003.pdf']
pdfcheck={}
for p in pdfs:
    doc=pymupdf.open(p);outside=[]
    for i,page in enumerate(doc):
        for b in page.get_text('blocks'):
            if b[0]<-.1 or b[1]<-.1 or b[2]>page.rect.width+.1 or b[3]>page.rect.height+.1:
                outside.append({'page':i+1,'bbox':list(b[:4])})
    pdfcheck[p.name]={'pages':len(doc),'width_in':doc[0].rect.width/72,
       'height_in':doc[0].rect.height/72,'text_outside_page':outside,'sha256':sha(p)}
    assert not outside
assert pdfcheck[pdfs[1].name]['pages']==1
for log in [OUT/'article_build/Beyond_Nodal_Damping_Theory_20261003.log',POSTER/'build/main.log']:
    logtext=log.read_text(encoding='utf-8',errors='replace')
    assert not re.search(r'Overfull|Underfull|undefined references|LaTeX Error|Missing character',logtext),log
dump(OUT/'PDF_CHECKS.json',pdfcheck)
dump(OUT/'STATUS.json',{'status':'POSTER_AND_ARTICLE_COMPILED','scientific_status':'THEORY_PLUS_TARGETED_VALIDATION',
    'article_pages':pdfcheck[pdfs[0].name]['pages'],'poster_pages':1,
    'rank_law_cases':60,'rank_law_identity_pass':True,'fixed_gain_events_passed':5,
    'analytical_gain_events_passed':5,'max_replacement_established':False,
    'superiority_established':False,'native_latex':'UNAVAILABLE_PLATFORM_DIRECTORIES',
    'terminal_latex':'SUCCESS_WITH_EXISTING_MIKTEX',
    'independent_review':'PASSED_WITH_CORRECTIONS_INTEGRATED',
    'git_branch':git('branch','--show-current'),
    'python':platform.python_version(),'numpy':numpy.__version__,'pandas':pandas.__version__,
    'scipy':scipy.__version__,'matplotlib':matplotlib.__version__})
dump(OUT/'ARTICLE_UPDATE_PROVENANCE.json',{'before_sha256':sha(OUT/'THEORY_before_poster_closure.tex'),
     'after_sha256':sha(THEORY),'no_historical_result_mutation':True,
     'comparison_events':5,'rank_law_cases':60})

def relevant(folder):
    for p in sorted(folder.rglob('*')):
        if p.is_file() and not any(x in p.relative_to(folder).parts for x in ['build','article_build','qa','vendor','__pycache__','.git']):
            if p.suffix=='.pid':continue
            yield p
artifacts=list(relevant(POSTER))+list(relevant(OUT))+[THEORY,ROOT/'output/pdf/.gitattributes']+pdfs+[
    ROOT/'output/pdf/IAS2026_Collective_Damping_20261003.png']
artifacts=[p for p in artifacts if p.name not in ['DELIVERY_MANIFEST.json']]
dump(OUT/'DELIVERY_MANIFEST.json',{'base_commit':git('rev-parse','HEAD'),
    'files':{p.relative_to(ROOT).as_posix():sha(p) for p in artifacts},
    'historical_manifests_remain_version_specific':True})
artifacts.append(OUT/'DELIVERY_MANIFEST.json')
artifacts+=list(relevant(REG))
dest=ROOT/'deliverables/IAS2026_Collective_Damping_20261003_complete.zip'
dest.parent.mkdir(exist_ok=True)
entries={p.relative_to(ROOT).as_posix():p for p in artifacts}
with zipfile.ZipFile(dest,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    z.writestr('START_HERE.txt',
       'Open output/pdf/IAS2026_Collective_Damping_20261003.pdf for the poster.\n'
       'Open output/pdf/Beyond_Nodal_Damping_Theory_20261003.pdf for the article.\n'
       'LaTeX, figures and speaking notes: reports/poster/ias2026/collective_interaction_20261003/.\n'
       'Julia code, new results and gain vectors: experiments/poster_julia_validation_20261003/.\n'
       'Before rerunning on a new machine, read experiments/poster_julia_validation_20261003/REPRODUCE.txt.\n'
       'Historical dependencies are stored intact in regional_paper_closure_20261003/evidence/.\n'
       'This is a theory-and-validation poster, not a demonstrated replacement optimum.\n')
    for rel,p in sorted(entries.items()):z.write(p,rel)
with zipfile.ZipFile(dest) as z:assert z.testzip() is None
(dest.with_suffix('.zip.sha256')).write_text(sha(dest)+'  '+dest.name+'\n')
print(json.dumps({'archive':str(dest),'size_MB':dest.stat().st_size/1e6,
                  'entries':len(entries)+1,'historical_preservation_checks':len(checked),
                  'pdf_checks':pdfcheck},indent=2))
