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

STATUS: proof-of-concept.  Cleanly validates single-diagram figures (44/68 in
Chap12 at time of writing, incl. every standard 6j-square / 9j-hexagon).
Known gaps (cause false "missing", NOT label errors): (a) \dsq squares whose
triad-forming diagonals are raw \draw'd with coordinate arithmetic the regex
skips; (b) 3jm "open" tree nodes; (c) multi-part derivation figures where one
equation carries several symbols mapping to several sub-diagrams -- needs
per-subpicture pairing.  Use PASS as trustworthy; investigate PARTIAL by eye.
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
    if name=='dsq': top,bot,left,right=a[0:4]; return [('sqTL','sqTR',top),('sqBL','sqBR',bot),('sqTL','sqBL',left),('sqTR','sqBR',right)]
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

def pairs_from_chapter(chap='Chap12.tex'):
    src=open(chap).read(); out=[]
    for m in re.finditer(r'\\centerline\{\\includegraphics(?:\[[^\]]*\])?\{(fig12[^}]*)\}\}',src):
        eqs=list(re.finditer(r'\\begin\{(equation|align|gather)\*?\}(.*?)\\end\{\1\*?\}',src[:m.start()],re.S))
        out.append((m.group(1), eqs[-1].group(2) if eqs else ''))
    return out

if __name__=='__main__':
    only=sys.argv[1] if len(sys.argv)>1 else None
    npass=npart=0
    for fig,eqtext in pairs_from_chapter():
        if only and only not in fig: continue
        fp=f'tikz_files/{fig}.tex'
        if not os.path.exists(fp): continue
        ets=eq_triads(eqtext)
        if not ets: continue
        diat=node_triads(fig_edges(open(fp).read()))
        miss=[(k,sorted(tr)) for k,triads in ets for tr in triads if len(tr)>=3 and not any(tr<=d for d in diat)]
        if miss:
            npart+=1; print(f"PARTIAL {fig}: {len(miss)} triad(s) not located: {sorted(set(map(lambda x:tuple(x[1]),miss)))}")
        else:
            npass+=1; print(f"PASS    {fig}")
    print(f"\n{npass} PASS, {npart} PARTIAL (investigate parser-gap vs real mismatch by eye)")
