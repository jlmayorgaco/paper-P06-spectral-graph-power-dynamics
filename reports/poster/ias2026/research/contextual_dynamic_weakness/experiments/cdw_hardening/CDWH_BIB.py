# ruff: noqa: E501
"""Build references.bib for the CDW paper from Crossref metadata of verified DOIs (no personal data
is sent). Book/chapter entries without a single DOI and the arXiv preprint are written by hand
from the verified gap-matrix rows."""

from __future__ import annotations

import _hinfra as HI  # noqa: I001

import json
import subprocess
import time

PAPER = HI.I.REPO / "reports" / "papers" / "cdw_contextual_dynamic_weakness"
DOIS = {
    "wu2018sdscr": "10.1109/TSTE.2017.2764871", "dong2019gscr": "10.1109/TPWRS.2018.2875305", "zhou2023heterogeneous": "10.1109/TPWRS.2022.3183005",
    "liu2024goscr": "10.1109/TPWRS.2023.3340158", "ma2024strength": "10.1109/TPWRD.2024.3432582", "henderson2024gsim": "10.1109/TPWRD.2022.3233455",
    "lamrani2025fagim": "10.1109/PowerTech59965.2025.11180496", "ebrahimzadeh2018busPF": "10.1109/TPEL.2018.2803846",
    "zhu2022greybox": "10.1109/TPWRS.2021.3088345", "zhu2023rootcause": "10.1109/TPWRS.2022.3179143", "coletta2016braess": "10.1103/PhysRevE.93.032222",
    "schafer2022braess": "10.1038/s41467-022-32917-6", "song2018network": "10.1109/TCNS.2017.2654162", "li2023intrinsic": "10.1109/TCSI.2023.3303644",
    "markovic2021lowinertia": "10.1109/TPWRS.2021.3061434", "stanojev2023vim": "10.1109/TPWRS.2022.3187789", "smed1993feasible": "10.1109/59.260827",
    "nam2000eigsens": "10.1109/59.852145", "mendozaarmenta2016redispatch": "10.1109/TPWRS.2015.2485519", "li2019opsens": "10.1016/j.epsr.2019.04.037",
    "saric2015rapid": "10.1109/TPWRS.2014.2342494", "li2018switching": "10.1109/TSG.2017.2656885", "summers2016submodularity": "10.1109/TCNS.2015.2453711",
    "olshevsky2018nonsupermodularity": "10.1109/TCNS.2017.2691463", "topkis1978": "10.1287/opre.26.2.305", "lovasz1983": "10.1007/978-3-642-68874-4_10",
    "bach2013": "10.1561/2200000039", "pourbeik2017generic": "10.1109/TEC.2016.2639050", "cui2021andes": "10.1109/TPWRS.2020.3017019",
    "athay1979": "10.1109/TPAS.1979.319407", "yang2021gfm": "10.1109/TPWRS.2020.3042741", "xin2025howmany": "10.1109/TPWRS.2024.3393877",
    "khaji2017switching": "10.1016/j.ijepes.2016.10.011", "huang2022gridstructure": "10.1016/j.ifacol.2022.07.270",
}
MANUAL = r"""
@book{sauer1998,
  author = {Sauer, Peter W. and Pai, M. A.},
  title = {Power System Dynamics and Stability},
  publisher = {Prentice Hall}, address = {Upper Saddle River, NJ}, year = {1998}
}
@book{fujishige2005,
  author = {Fujishige, Satoru},
  title = {Submodular Functions and Optimization},
  edition = {2nd}, series = {Annals of Discrete Mathematics}, volume = {58},
  publisher = {Elsevier}, address = {Amsterdam}, year = {2005}
}
@misc{joswigjones2026sens,
  author = {Joswig-Jones, T. and Dong, W. and Tan, B. and others},
  title = {Sensitivity-Based System Strength Assessment: Mapping Power Flow and Network Topology Perturbations to System Eigenvalues},
  howpublished = {arXiv:2607.28764 (preprint, not peer reviewed)}, year = {2026}
}
@article{xin2016gscr,
  author = {Xin, Huanhai and Dong, Wei and Yuan, Xiaoming and Gan, Deqiang and Wang, Kang and Xie, Huan},
  title = {Generalized Short Circuit Ratio for Multi Power Electronic Based Devices Infeed to Power Systems},
  journal = {Proceedings of the CSEE},
  year = {2016}, volume = {36}, number = {22}, pages = {6013--6027},
  doi = {10.13334/j.0258-8013.pcsee.161682}
}
"""


def tex(s):
    import html

    s = html.unescape(s or "").replace("®", "").replace("’", "'")
    return s.replace("&", r"\&").replace("%", r"\%").replace("_", r"\_").replace("#", r"\#")


def fetch(doi):
    out = subprocess.run(["curl", "-s", "-m", "30", f"https://api.crossref.org/works/{doi}"], capture_output=True, text=True, encoding="utf-8")
    return json.loads(out.stdout)["message"]


def entry(key, m):
    typ = m.get("type")
    auth = " and ".join(f"{tex(a.get('family', ''))}, {tex(a.get('given', ''))}".strip(", ") for a in m.get("author", []))
    year = (m.get("published-print") or m.get("published-online") or m.get("issued"))["date-parts"][0][0]
    title = tex(m["title"][0]).replace("\n", " ")
    cont = tex((m.get("container-title") or [""])[0])
    f = {"author": auth, "title": "{{" + title + "}}", "year": str(year), "doi": m["DOI"]}
    if m.get("volume"):
        f["volume"] = m["volume"]
    if m.get("issue"):
        f["number"] = m["issue"]
    if m.get("page"):
        f["pages"] = m["page"].replace("-", "--")
    if typ == "journal-article":
        kind, f["journal"] = "article", cont
    elif typ == "book-chapter":
        kind, f["booktitle"] = "incollection", cont
        f["publisher"] = tex(m.get("publisher", ""))
    else:
        kind, f["booktitle"] = "inproceedings", cont
    body = ",\n".join(f"  {k} = {{{v}}}" if k != "title" else f"  {k} = {v}" for k, v in f.items())
    return f"@{kind}{{{key},\n{body}\n}}\n"


def main():
    parts = [MANUAL.strip() + "\n"]
    for key, doi in DOIS.items():
        m = fetch(doi)
        parts.append(entry(key, m))
        time.sleep(0.3)
    (PAPER / "references.bib").write_text("\n".join(parts), encoding="utf-8")
    print("entries", len(DOIS) + 4)


if __name__ == "__main__":
    main()
