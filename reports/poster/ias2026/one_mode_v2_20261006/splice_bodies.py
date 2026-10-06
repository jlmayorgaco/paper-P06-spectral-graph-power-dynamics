"""Replace the content block of each section file with bodies/NN.tex (keeps the file's preview/compile scaffold)."""
import pathlib
H = pathlib.Path(__file__).resolve().parent
MAP = {'01': '01_motivation_question', '02': '02_flagship_ieee39', '03': '03_target_vs_path', '04': '04_network_closure',
       '05': '05_policy_atlas', '06': '06_synthetic_envelope', '07': '07_retuning', '08': '08_safe_paths', '09': '09_evidence_scope', 'strip': 'result_strip', 'footer': 'footer'}
BS = chr(92)
head = BS + 'long' + BS + 'def' + BS + 'PosterSectionContent{%\n'
tail = '}%\n' + BS + 'ifdefined' + BS + 'PosterFull'
for k, name in MAP.items():
    if not (H / 'bodies' / f'{k}.tex').exists(): continue
    p = H / 'sections' / f'{name}.tex'
    s = p.read_text(encoding='utf-8')
    a = s.index(head) + len(head)
    b = s.index(tail)
    body = (H / 'bodies' / f'{k}.tex').read_text(encoding='utf-8')
    p.write_text(s[:a] + body + s[b:], encoding='utf-8', newline='\n')
    print('spliced', name)
