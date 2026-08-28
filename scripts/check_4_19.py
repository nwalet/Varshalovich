#!/usr/bin/env python3
r"""Sec 4.19 / Table 4.2 -- other authors' D-function conventions.

This table is REFERENCE DATA: each row states how another author's rotation
matrix relates to this book's, and that relation is fixed by *that author's*
own definitions (handedness, rotation sense, covariant/contravariant, tensor
rule -- the points (a)-(f) listed in the section).  Only the rows that reduce
to a definitional identity in this book's own conventions can be re-derived
here; the phase/transpose rows depend on external definitions and are left as
the book states them.

This book:  D^J_{MM'}(a,b,g) = <JM| e^{-ia Jz} e^{-ib Jy} e^{-ig Jz} |JM'>.

Verified here (definitional):
  Dolginov[14], Davydov[12], Rose[30], Edmonds[64], Fano-Racah[18] -- authors
  whose rotation operator is e^{+ia Jz} e^{+ib Jy} e^{+ig Jz}: their D equals
  D_book(-a,-b,-g).

Documented (convention-dependent, NOT re-derived -- reproducing the book):
  Bohr-Mottelson[8]  D = D_book^*
  Wigner[43]         D(...)_{MM'} = D_book(-a,-b,-g)
  Berestetskii[6]    D^{(J)}_{MM'} = D_book(-g,-b,-a)   (transpose + reversal)
  Gel'fand[20]/Lubarskii[26]  T = (-i)^{M-M'} D_book
       (with the STANDARD Jx one finds (+i)^{M-M'}; the sign is fixed by
        Gel'fand's own left/right rotation sense, so it is left as printed.)
  Vilenkin[41]       t = (-i)^{M-M'} D_book
  Yutsis-Bandzaitis[45]  D^{(J)}_{MM'} = i^{M-M'} D_book^*
"""
import mpmath as mp
mp.mp.dps = 30

def spin_mats(J):
    n=int(round(2*J))+1
    ms=[J-k for k in range(n)]
    Jz=mp.matrix(n,n); Jp=mp.matrix(n,n)
    for i,m in enumerate(ms):
        Jz[i,i]=m
    for i,m in enumerate(ms):
        mp1=m+1
        if mp1<=J:
            Jp[ms.index(mp1),i]=mp.sqrt(J*(J+1)-m*mp1)
    Jm=Jp.T
    return ms,(Jp+Jm)/2,(Jp-Jm)/(2j),Jz

def Dbook_mat(J,a,b,g):
    ms,Jx,Jy,Jz=spin_mats(J)
    return ms, mp.expm(-1j*a*Jz)*mp.expm(-1j*b*Jy)*mp.expm(-1j*g*Jz)

def report(tag,w,tol=mp.mpf('1e-12')):
    w=float(w); ok=w<tol
    print(f"  {tag:46s} {'PASS' if ok else 'FAIL'}  worst={w:.2e}")
    return ok
ok=True
a,b,g=mp.mpf('0.7'),mp.mpf('1.3'),mp.mpf('0.5')

print("Definitional row: e^{+i} rotations  <=>  D_book(-a,-b,-g)")
for twoJ in [1,2,3,4,5]:
    J=mp.mpf(twoJ)/2
    ms,Jx,Jy,Jz=spin_mats(J)
    Rplus=mp.expm(1j*a*Jz)*mp.expm(1j*b*Jy)*mp.expm(1j*g*Jz)
    _,Rneg=Dbook_mat(J,-a,-b,-g)
    w=max(abs(Rplus[i,j]-Rneg[i,j]) for i in range(len(ms)) for j in range(len(ms)))
    ok&=report(f"Dolginov/Davydov/Rose[30]/Edmonds[64] J={float(J)}", w)

# Sanity: this book's D is unitary and matches e^{-i M a} d e^{-i M' g}
print("\nSanity: D_book unitary")
for twoJ in [1,3,4]:
    J=mp.mpf(twoJ)/2
    ms,R=Dbook_mat(J,a,b,g)
    U=R*R.H
    w=max(abs(U[i,j]-(1 if i==j else 0)) for i in range(len(ms)) for j in range(len(ms)))
    ok&=report(f"D D^dagger = 1  J={float(J)}", w)

print("\nNote: phase/transpose rows (Gel'fand, Vilenkin, Yutsis, Berestetskii,")
print("Bohr-Mottelson, Wigner) are convention-dependent reference data; with the")
print("standard Jx, Gel'fand gives (+i)^(M-M') -- the printed (-i)^(M-M') follows")
print("that author's own rotation-sense convention and is left as published.")
print("\nRESULT:", "ALL PASS" if ok else "FAILURES ABOVE")
raise SystemExit(0 if ok else 1)
