from pathlib import Path
import re
from xml.sax.saxutils import escape
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT, TA_JUSTIFY, TA_CENTER
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import BaseDocTemplate, PageTemplate, Frame, Paragraph, Spacer, PageBreak, NextPageTemplate, Table, TableStyle, Image

ROOT=Path(__file__).resolve().parents[3]
SRC=ROOT/"reports/handoff/ai_migration_20261002/HANDOFF_MASTER.md"
OUT=ROOT/"output/pdf/IEEE39_BND_PROJECT_HANDOFF_20261002.pdf"
FONT=Path("C:/Windows/Fonts")
for n,f in [("Arial","arial.ttf"),("Arial-Bold","arialbd.ttf"),("Arial-Italic","ariali.ttf"),("Arial-BoldItalic","arialbi.ttf")]:
    pdfmetrics.registerFont(TTFont(n,str(FONT/f)))
pdfmetrics.registerFontFamily("Arial",normal="Arial",bold="Arial-Bold",italic="Arial-Italic",boldItalic="Arial-BoldItalic")
NAVY=colors.HexColor("#12304A"); BLUE=colors.HexColor("#1D6485"); TEAL=colors.HexColor("#147D78")
INK=colors.HexColor("#202B35"); MUTED=colors.HexColor("#5D6B75"); PALE=colors.HexColor("#EEF4F7"); RULE=colors.HexColor("#D5E0E6")
ss=getSampleStyleSheet()
def sty(name,**kw): ss.add(ParagraphStyle(name=name,**kw))
sty("CoverTitle",fontName="Arial-Bold",fontSize=22,leading=26,textColor=NAVY,spaceAfter=8)
sty("CoverSub",fontName="Arial",fontSize=9,leading=13,textColor=MUTED,spaceAfter=12)
sty("CoverLabel",fontName="Arial-Bold",fontSize=8,leading=10,textColor=BLUE,spaceBefore=5,spaceAfter=4)
sty("CoverBody",fontName="Arial",fontSize=8.6,leading=12,textColor=INK,alignment=TA_JUSTIFY,spaceAfter=6)
sty("Body",fontName="Arial",fontSize=7.1,leading=9.1,textColor=INK,alignment=TA_JUSTIFY,spaceAfter=3.5)
sty("BulletCustom",fontName="Arial",fontSize=7,leading=9,textColor=INK,leftIndent=9,firstLineIndent=-7,spaceAfter=2)
sty("H1",fontName="Arial-Bold",fontSize=11.5,leading=13.5,textColor=NAVY,spaceBefore=7,spaceAfter=4,keepWithNext=True)
sty("H2",fontName="Arial-Bold",fontSize=8.7,leading=10.5,textColor=BLUE,spaceBefore=5,spaceAfter=3,keepWithNext=True)
sty("H3",fontName="Arial-Bold",fontSize=7.7,leading=9.2,textColor=TEAL,spaceBefore=4,spaceAfter=2,keepWithNext=True)
sty("Math",fontName="Arial",fontSize=6.7,leading=8.5,textColor=NAVY,leftIndent=5,rightIndent=3,borderColor=RULE,borderWidth=.4,borderPadding=4,backColor=PALE,spaceBefore=3,spaceAfter=5)
sty("Cell",fontName="Arial",fontSize=5.7,leading=7.1,textColor=INK)
sty("Head",fontName="Arial-Bold",fontSize=5.7,leading=7.1,textColor=colors.white)
sty("Caption",fontName="Arial-Italic",fontSize=6,leading=7.5,textColor=MUTED,alignment=TA_CENTER,spaceBefore=2,spaceAfter=4)
sty("Small",fontName="Arial",fontSize=6.3,leading=8,textColor=MUTED,spaceAfter=3)

def plain_math(s):
    s=s.replace(r"\(|\Delta f|\)","absolute delta f")
    # Replace simple LaTeX fractions with readable quotient notation.
    while r"\frac" in s:
        p=s.find(r"\frac"); i=p+5
        while i<len(s) and s[i].isspace(): i+=1
        def group(j):
            depth=0
            for k in range(j,len(s)):
                if s[k]=="{": depth+=1
                elif s[k]=="}":
                    depth-=1
                    if depth==0:return s[j+1:k],k+1
            return s[j+1:],len(s)
        if i>=len(s) or s[i]!="{": break
        a,j=group(i)
        while j<len(s) and s[j].isspace():j+=1
        if j>=len(s) or s[j]!="{":break
        b,k=group(j); s=s[:p]+"("+a+")/("+b+")"+s[k:]
    for a,b in [(r"\dot{x}","dx/dt"),(r"\dot{S}","dS/dt"),(r"\mathcal{W}","W"),(r"\mathcal{K}","K"),(r"\mathcal{B}","B"),(r"\mathbb{R}","R"),(r"\nabla","grad "),(r"\lambda","lambda"),(r"\rho","rho"),(r"\kappa","kappa"),(r"\theta","theta"),(r"\eta","eta"),(r"\tau","tau"),(r"\phi","phi"),(r"\Delta","Delta"),(r"\Pi","Pi"),(r"\nu","nu"),(r"\cdot","*"),(r"\times","x"),(r"\le","<="),(r"\ge",">="),(r"\in"," in "),(r"\neq","!="),(r"\|","||"),(r"\partial","d"),(r"\max","max"),(r"\min","min"),(r"\sup","sup"),(r"\sum","sum"),(r"\quad"," "), (r"\qquad"," "), (r"\left",""),(r"\right",""),(r"\mathrm",""),(r"\text",""),(r"\operatorname",""),(r"\begin{aligned}",""),(r"\end{aligned}",""),(r"\\"," ; "),("&"," ")]:
        s=s.replace(a,b)
    s=s.replace(r"\{","{").replace(r"\}","}").replace("{","").replace("}","").replace("~"," ")
    s=re.sub(r"\\[A-Za-z]+","",s)
    return re.sub(r"\s+"," ",s).strip()

def markup(s):
    s=s.replace(r"\(","MATHOPEN").replace(r"\)","MATHCLOSE")
    s=re.sub(r"MATHOPEN(.*?)MATHCLOSE",lambda m:plain_math(m.group(1)),s)
    s=escape(s)
    s=re.sub(r"\*\*(.+?)\*\*",r"<b>\1</b>",s)
    s=re.sub(r"(?<!\*)\*([^*]+)\*(?!\*)",r"<i>\1</i>",s)
    return s

def table_from(lines):
    rows=[]
    for line in lines:
        line=line.strip().strip("|").replace(r"\(|\Delta f|\)","absolute delta f")
        cells=[c.strip() for c in line.split("|")]
        if cells and all(re.fullmatch(r":?-{2,}:?",c or "-") for c in cells):continue
        rows.append(cells)
    if not rows:return None
    n=max(map(len,rows)); out=[]
    for ri,row in enumerate(rows):
        row+=[""]*(n-len(row)); style=ss["Head"] if ri==0 else ss["Cell"]
        out.append([Paragraph(markup(c),style) for c in row])
    w=256; widths=([w*.32,w*.68] if n==2 else [w*.27,w*.38,w*.35] if n==3 else [w/n]*n)
    t=Table(out,colWidths=widths,repeatRows=1,hAlign="LEFT",splitByRow=1)
    t.setStyle(TableStyle([("BACKGROUND",(0,0),(-1,0),NAVY),("ROWBACKGROUNDS",(0,1),(-1,-1),[colors.white,PALE]),("GRID",(0,0),(-1,-1),.35,RULE),("VALIGN",(0,0),(-1,-1),"TOP"),("LEFTPADDING",(0,0),(-1,-1),3),("RIGHTPADDING",(0,0),(-1,-1),3),("TOPPADDING",(0,0),(-1,-1),3),("BOTTOMPADDING",(0,0),(-1,-1),3)]))
    return t

def get_cover(lines):
    title=next((x[2:].strip() for x in lines if x.startswith("# ")),"Research handoff")
    meta_line=next((x for x in lines if x.startswith("**Internal project-transfer")),"")
    meta_match=re.search(r"\*\*(.*?)\*\*",meta_line)
    meta=meta_match.group(1) if meta_match else meta_line
    i=next(k for k,x in enumerate(lines) if x.strip()=="## Abstract")
    ab=[]
    for x in lines[i+1:]:
        if x.startswith("**Index terms:**") or x.startswith("## "):break
        if x.strip():ab.append(x.strip())
    j=next(k for k,x in enumerate(lines) if x.startswith("## 1. "))
    return title,meta," ".join(ab),lines[j:]

def story():
    lines=SRC.read_text(encoding="utf-8-sig").splitlines()
    title,meta,abstract,body=get_cover(lines)
    s=[Spacer(1,4),Paragraph(escape(title),ss["CoverTitle"]),Paragraph(escape(meta),ss["CoverSub"]),Paragraph("ABSTRACT",ss["CoverLabel"]),Paragraph(markup(abstract),ss["CoverBody"])]
    box=Table([[Paragraph("<b>Resultado físico más reciente</b>",ss["Body"])],[Paragraph("A 90.7854% de GFL fijo, el PLL co-diseñado redujo los seis picos de desviación de frecuencia y cambió los incumplimientos exploratorios de 2/6 a 0/6. El máximo RoCoF aumentó 37.7%. Es una comparación finita, no una frontera de capacidad.",ss["Body"])]],colWidths=[7.15*inch])
    box.setStyle(TableStyle([("BACKGROUND",(0,0),(-1,-1),PALE),("BOX",(0,0),(-1,-1),.8,BLUE),("LINEBELOW",(0,0),(-1,0),.4,RULE),("LEFTPADDING",(0,0),(-1,-1),9),("RIGHTPADDING",(0,0),(-1,-1),9),("TOPPADDING",(0,0),(-1,-1),7),("BOTTOMPADDING",(0,0),(-1,-1),7)]))
    s += [Spacer(1,4),box,Spacer(1,8),Paragraph("EVIDENCE BOUNDARY",ss["CoverLabel"]),Paragraph("This dossier preserves model mismatches, failed optimization and TDS gates, incomplete robustness/globality, a failed graph modal-separation test, and missing hard current/DC limits. The 90.7854%, 90.0470%, and 91.885% candidates come from different runs and contracts. No global optimum, universal robustness, or established novelty is claimed.",ss["CoverBody"]),Paragraph("HOW TO USE THIS PACKET",ss["CoverLabel"]),Paragraph("Start with START_HERE_FOR_NEXT_AI.md, then read this dossier and ARTIFACT_MAP.json. Preserve the working tree. Verify every numerical claim against raw CSV/TOML outputs and source scripts before extending it.",ss["CoverBody"]),Paragraph("Repository snapshot 2026-10-02 | Internal continuity record | IEEE-style two-column report",ss["Small"]),NextPageTemplate("body"),PageBreak()]
    i=0; para=[]
    def flush():
        nonlocal para
        if para:
            text=" ".join(x.strip() for x in para)
            if text:s.append(Paragraph(markup(text),ss["Body"]))
            para=[]
    while i<len(body):
        x=body[i]
        if not x.strip():flush();i+=1;continue
        if x.startswith("|"):
            flush(); rows=[]
            while i<len(body) and body[i].startswith("|"):rows.append(body[i]);i+=1
            t=table_from(rows)
            if t:s.extend([Spacer(1,2),t,Spacer(1,4)])
            continue
        if x.startswith("!["):
            flush(); m=re.search(r"\]\(([^)]+)\)",x)
            if m:
                p=ROOT/m.group(1)
                if p.exists():
                    im=Image(str(p)); scale=min(248/im.imageWidth,150/im.imageHeight)
                    im.drawWidth*=scale;im.drawHeight*=scale;s.extend([im,Paragraph("Fixed-mix PLL retuning tradeoff across six independent IEEE-39 load-step tests.",ss["Caption"])])
            i+=1;continue
        if x.startswith("## "):flush();s.append(Paragraph(markup(x[3:]),ss["H1"]));i+=1;continue
        if x.startswith("### "):flush();s.append(Paragraph(markup(x[4:]),ss["H2"]));i+=1;continue
        if x.startswith("#### "):flush();s.append(Paragraph(markup(x[5:]),ss["H3"]));i+=1;continue
        if x.startswith(r"\["):
            flush(); math=x[2:]; i+=1
            while r"\]" not in math and i<len(body):
                math+=" "+body[i];i+=1
            math=math.split(r"\]")[0];s.append(Paragraph(escape(plain_math(math)),ss["Math"]));continue
        if re.match(r"^\s*[-*]\s+",x):
            flush(); item=re.sub(r"^\s*[-*]\s+","",x);s.append(Paragraph("&#8226; "+markup(item),ss["BulletCustom"]));i+=1;continue
        m=re.match(r"^\s*(\d+)\.\s+(.*)$",x)
        if m:
            flush();s.append(Paragraph("<b>"+m.group(1)+".</b> "+markup(m.group(2)),ss["BulletCustom"]));i+=1;continue
        para.append(x);i+=1
    flush();return s

class Doc(BaseDocTemplate):
    def __init__(self,path):
        super().__init__(path,pagesize=letter,leftMargin=.48*inch,rightMargin=.48*inch,topMargin=.52*inch,bottomMargin=.5*inch,title="IEEE-Style Research Handoff: IEEE-39 SG-to-GFL and Beyond Nodal Damping",author="Research project continuity record")
        w,h=letter;x0=self.leftMargin;x1=w-self.rightMargin;y0=self.bottomMargin;y1=h-self.topMargin
        cover=Frame(x0,y0,x1-x0,y1-y0,id="cover",leftPadding=0,rightPadding=0,topPadding=0,bottomPadding=0)
        cw=(x1-x0-16)/2
        left=Frame(x0,y0,cw,y1-y0,id="left",leftPadding=0,rightPadding=0,topPadding=0,bottomPadding=0)
        right=Frame(x0+cw+16,y0,cw,y1-y0,id="right",leftPadding=0,rightPadding=0,topPadding=0,bottomPadding=0)
        self.addPageTemplates([PageTemplate(id="cover",frames=[cover],onPage=self.decorate),PageTemplate(id="body",frames=[left,right],onPage=self.decorate)])
    def decorate(self,c,doc):
        w,h=letter;c.saveState();c.setStrokeColor(RULE);c.setLineWidth(.5)
        if doc.page>1:
            c.line(self.leftMargin,h-.36*inch,w-self.rightMargin,h-.36*inch);c.setFont("Arial-Bold",7);c.setFillColor(NAVY);c.drawString(self.leftMargin,h-.28*inch,"PROJECT HANDOFF | IEEE-39 SG-to-GFL / BEYOND NODAL DAMPING")
        c.line(self.leftMargin,.36*inch,w-self.rightMargin,.36*inch);c.setFont("Arial",6.5);c.setFillColor(MUTED);c.drawString(self.leftMargin,.22*inch,"Internal research-transfer dossier | 2026-10-02");c.drawRightString(w-self.rightMargin,.22*inch,str(doc.page));c.restoreState()

if __name__=="__main__":
    OUT.parent.mkdir(parents=True,exist_ok=True);Doc(str(OUT)).build(story());print(f"Wrote {OUT} ({OUT.stat().st_size} bytes)")


