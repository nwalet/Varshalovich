#!/usr/bin/env python3
r"""Chapter 12 diagram/equation triad-consistency check (prototype).

Idea (per NW): a 3nj symbol in an equation imposes triangle conditions on
triples of its arguments; in the accompanying diagram each such triple must
appear as three labelled lines meeting at one node.  So:

  equation  --(extract \sixj/\ninej/\threej)-->  set of triads (triples)
  figure    --(macros + raw \draw edges)------->  graph; node label-sets
  check     : every equation triad == some node's 3-line label-set

3nj triads:  6j {abc;def} -> (abc)(aef)(dbf)(dec);
             9j rows+cols; 3j -> (abc).
Figure graph: geometric macros (\dsixjsq, \dsixjtri, \dsixjr/l, \dninehex,
\dhexflat, \dsq) expand to labelled edges on their *exposed* coordinate names,
so caller-drawn raw diameters (\draw (hT)--(hB)...) join on the same nodes;
theta bubbles (\ctheta[v], \dtheta[v]) are self-contained 2-node/3-line triads.

Scoping: geometric macros place their exposed coordinates at FIXED names (sqTL,
hT, ...), so a figure with several sub-diagrams reuses them, redefined inside
each \begin{scope}[xshift=..].  fig_triads() therefore builds the incidence
graph per scope (plus the picture-level connective draws) and unions the triads,
so sibling sub-diagrams don't cross-contaminate.  Pairing: Chap12 is laid out
"<equation(s)> then <figure>" repeatedly, so a figure's symbols are exactly the
3nj's in the text since the PREVIOUS figure (pairs_from_chapter).

OUTPUT is a triage list, not a pass/fail oracle: PASS = every equation triad
found as a diagram node (trustworthy).  For the rest we print how many triads
were located vs total and rank by fraction missing:
  * REVIEW  (few missing, most located) -- strong candidate for a real
    equation/diagram label mismatch; check these first against the scan.
  * (many missing) -- usually the whole sub-diagram uses a drawing style this
    parser does not model (3jm open-tree nodes, nested 12j blocks, ...); low
    priority.  Modelled well: 6j-square (\dsixjsq/\dsq), 9j-hex (\dninehex/
    \dhexflat), theta bubbles.

KNOWN residual false positives (all verified harmless on the 2026-08-30 pass):
a missing triad that CONTAINS a summed index (X/Y/Z, shown lower-cased x/y/z)
usually sits at a vertex where two sub-diagrams are glued along that shared line
-- an incidence that spans two \begin{scope}s and so is invisible to the
per-scope graph (e.g. the two hexagons of a 12j: fig12_2_11/12/33c/45).  A
missing triad with only NON-summed labels stays the strong real candidate.
"""
import re, os, sys
GREEK=r'\\(?:alpha|beta|gamma|delta|epsilon|varepsilon|eta|xi|kappa|sigma|psi|rho|theta|tau|lambda|mu|nu|zeta|omega|phi|chi|varphi)'
def jlabel(lab):
    lab=re.sub(r'\{|\}','',lab).strip().strip('$')
    lab=re.sub(GREEK,'',lab); lab=re.sub(r"[\\,\s']",'',lab)
    return lab.lower()                 # summed indices are Capitalised in diagrams
def sixj(a,b,c,d,e,f): return [{a,b,c},{a,e,f},{d,b,f},{d,e,c}]
def ninej(a,b,c,d,e,f,g,h,i): return [{a,b,c},{d,e,f},{g,h,i},{a,d,g},{b,e,h},{c,f,i}]
def threej(a,b,c,*_): return [{a,b,c}]
def eq_triads(t):
    out=[]
    for tag,fn,n in (('6j',sixj,6),('9j',ninej,9),('3j',threej,6)):
        for m in re.finditer(r'\\'+{'6j':'sixj','9j':'ninej','3j':'threej'}[tag]+r'((?:\{[^{}]*\}){'+str(n)+'})',t):
            out.append((tag,fn(*[jlabel(x) for x in re.findall(r'\{([^{}]*)\}',m.group(1))])))
    return out
def macro_edges(name,a):
    if name=='dsixjsq': tl,tr,bl,br,hd,vd=a[2:8]; return [('ql','qt',tl),('qt','qr',tr),('ql','qb',bl),('qb','qr',br),('ql','qr',hd),('qt','qb',vd)]
    if name=='dsixjtri': eL,eR,eB,sT,sBL,sBR=a[2:8]; return [('vL','vT',eL),('vT','vR',eR),('vL','vR',eB),('ct','vT',sT),('ct','vL',sBL),('ct','vR',sBR)]
    if name in ('dsixjr','dsixjl'): l,t,b,mt,mb,mr=a[0:6]; return [('P','Q',l),('P','V',t),('Q','V',b),('c','P',mt),('c','Q',mb),('c','V',mr)]
    if name=='dninehex': e=a[2:8]; V=['hT','hUL','hLL','hB','hLR','hUR']; return [(V[i],V[(i+1)%6],e[i]) for i in range(6)]
    if name=='dhexflat': e=a[2:8]; return [('v2','v1',e[0]),('v1','v0',e[1]),('v0','v5',e[2]),('v5','v4',e[3]),('v4','v3',e[4]),('v3','v2',e[5])]
    if name=='dsq': top,bot,left,right=a[2:6]; return [('sqTL','sqTR',top),('sqBL','sqBR',bot),('sqTL','sqBL',left),('sqTR','sqBR',right)]
    return None
def fig_edges(t):
    edges=[]
    for m in re.finditer(r'\\(d[a-z]+)(?:\[[^\]]*\])?((?:\{[^{}]*\})+)',t):
        e=macro_edges(m.group(1),re.findall(r'\{([^{}]*)\}',m.group(2)))
        if e: edges+=[(x,y,jlabel(l)) for x,y,l in e]
    for i,m in enumerate(re.finditer(r'\\(cthetav|ctheta|dthetav|dtheta)(?:\[[^\]]*\])?((?:\{[^{}]*\})+)',t)):
        a=re.findall(r'\{([^{}]*)\}',m.group(2))
        if len(a)>=5:
            for l in a[2:5]: edges.append((f'th{i}a',f'th{i}b',jlabel(l)))
    for m in re.finditer(r'\\draw\[[^\]]*\]\s*\(([^)]*)\)--\(([^)]*)\)[^;]*?node[^;]*?\{\$?([^${}]*)\$?\}',t):
        A=m.group(1).split('+')[0].strip(); B=m.group(2).split('+')[0].strip()
        edges.append((A,B,jlabel(m.group(3))))
    return edges
def node_triads(edges):
    from collections import defaultdict
    inc=defaultdict(set)
    for A,B,l in edges:
        if l: inc[A].add(l); inc[B].add(l)
    return [s for s in inc.values() if len(s)>=3]

def fig_triads(figtext):
    # A geometric macro (\dsq, \dninehex, ...) always places its exposed
    # coordinates (sqTL, hT, ...) at the SAME names.  When a figure draws several
    # sub-diagrams, each sits in its own \begin{scope}[xshift=..] and REDEFINES
    # those names, drawing its labelled edges immediately -- so within one scope
    # the incidences are right, but unioning names across sibling scopes scrambles
    # them.  Analyse each scope together with the picture-level (outside-scope)
    # connective draws, but never with a sibling scope.  Union the triads.
    pics=re.findall(r'\\begin\{tikzpicture\}.*?\\end\{tikzpicture\}',figtext,re.S) or [figtext]
    triads=[]
    for pic in pics:
        scopes=re.findall(r'\\begin\{scope\}.*?\\end\{scope\}',pic,re.S)
        if scopes:
            base=re.sub(r'\\begin\{scope\}.*?\\end\{scope\}','',pic,flags=re.S)
            units=[base+s for s in scopes]
        else:
            units=[pic]
        for u in units:
            triads+=node_triads(fig_edges(u))
    return triads

def pairs_from_chapter(chap='Chap12.tex'):
    # Chapter 12 lays out "<equation(s)> then <figure>" repeatedly.  A figure's
    # symbols are exactly those in the text since the PREVIOUS figure, so bound
    # the equation-region that way -- this keeps each figure from being charged
    # with a neighbour's 3nj symbol, and naturally gathers all symbols (LHS+RHS)
    # that the figure's sub-diagrams are meant to illustrate.
    src=open(chap).read(); out=[]; prev=0
    for m in re.finditer(r'\\centerline\{\\includegraphics(?:\[[^\]]*\])?\{(fig12[^}]*)\}\}',src):
        out.append((m.group(1), src[prev:m.start()])); prev=m.end()
    return out

if __name__=='__main__':
    only=sys.argv[1] if len(sys.argv)>1 else None
    npass=0; rows=[]
    for fig,eqtext in pairs_from_chapter():
        if only and only not in fig: continue
        fp=f'tikz_files/{fig}.tex'
        if not os.path.exists(fp): continue
        ets=eq_triads(eqtext)
        if not ets: continue
        diat=fig_triads(open(fp).read())
        want=[tr for _,triads in ets for tr in triads if len(tr)>=3]
        miss=sorted({tuple(sorted(tr)) for tr in want if not any(tr<=d for d in diat)})
        if not miss:
            npass+=1
        else:
            rows.append((len(miss)/len(want), len(want)-len(miss), len(want), fig, miss))
    # few-missing (strong candidates) first
    rows.sort()
    print(f"{npass} PASS  |  {len(rows)} to review  (Chap12; sorted best-candidate first)\n")
    for frac,got,tot,fig,miss in rows:
        tag='REVIEW ' if frac<=0.2 else '       '   # >=80% located -> a couple of anomalies stand out
        print(f"{tag}{fig:<13} {got:>2}/{tot:<2} triads located; missing: "
              +', '.join('('+' '.join(t)+')' for t in miss))
