#!/usr/bin/env python3
r"""Checks for Secs 5.4 (symmetry), 5.7 (recursions), 5.8 (differential
relations) of Chapter 5, VMK.  Reference Y = mpmath.spherharm.
"""
import math, cmath
import mpmath as mp
mp.mp.dps = 30
TOL = 1e-8
TH = [0.4, 0.9, 1.3, 1.9, 2.6]
PH = [0.3, 1.8, 4.0]

def Y(l, m, th, ph):
    if l < 0 or abs(m) > l: return 0j
    return complex(mp.spherharm(l, m, th, ph))
def sq(x): return math.sqrt(x) if x >= 0 else float('nan')
def report(tag, w, tol=TOL):
    ok = w < tol; print(f"  {tag:34s} {'PASS' if ok else 'FAIL'}  worst={w:.2e}"); return ok

def main():
    ok = True
    print("Sec 5.4 symmetry")
    # 5.4.1 Y* = Y(th,-ph) = (-1)^m Y_{l,-m}
    w = 0.0
    for l in range(0,5):
        for m in range(-l,l+1):
            for th in TH:
                for ph in PH:
                    c = Y(l,m,th,ph).conjugate()
                    w = max(w, abs(c-Y(l,m,th,-ph)), abs(c-(-1)**m*Y(l,-m,th,ph)))
    ok &= report("5.4.1 conjugation", w)
    # 5.4.2 Y_{l,-m} = (-1)^m Y(th,-ph) = (-1)^m e^{-2imph} Y
    w = 0.0
    for l in range(0,5):
        for m in range(-l,l+1):
            for th in TH:
                for ph in PH:
                    w = max(w, abs(Y(l,-m,th,ph)-(-1)**m*Y(l,m,th,-ph)),
                            abs(Y(l,-m,th,ph)-(-1)**m*cmath.exp(-2j*m*ph)*Y(l,m,th,ph)))
    ok &= report("5.4.2 sign reversal of m", w)
    # 5.4.4-6, 5.4.7-9
    def chk(tag, f, lmin=0):
        w = 0.0
        for l in range(lmin,6):
            for m in range(-l,l+1):
                for th in TH:
                    for ph in PH:
                        try: v = f(l,m,th,ph)
                        except (ValueError, ZeroDivisionError): continue
                        if v is None or not (math.isfinite(v.real) and math.isfinite(v.imag)): continue
                        w = max(w, abs(v))
        return report(tag, w)
    ok &= chk("5.4.4 Y(pi-th,ph)", lambda l,m,th,ph: Y(l,m,math.pi-th,ph)-(-1)**(l+m)*Y(l,m,th,ph))
    ok &= chk("5.4.5 Y(th,pi+ph)", lambda l,m,th,ph: Y(l,m,th,math.pi+ph)-(-1)**m*Y(l,m,th,ph))
    ok &= chk("5.4.6 Y(pi-th,pi+ph)", lambda l,m,th,ph: Y(l,m,math.pi-th,math.pi+ph)-(-1)**l*Y(l,m,th,ph))
    ok &= chk("5.4.8 Y(th,-ph)", lambda l,m,th,ph: Y(l,m,th,-ph)-(-1)**m*Y(l,-m,th,ph))
    # 5.4.7 Y(-th,ph)=(-1)^m Y and 5.4.9 Y(-th,-ph)=Y_{l,-m} use the book's
    # continuation of Y to negative theta via the (sin th)^m branch; mpmath's
    # spherharm continues cos-evenly (P_l^m(cos th)), so these are correct under
    # VMK's convention but not reproducible with spherharm. (Not checked here.)
    print("  5.4.7/5.4.9  negative-theta continuation -- convention (not checked)")

    print("\nSec 5.7 recursions")
    ok &= chk("5.7.1 -2m cot th Y", lambda l,m,th,ph:
        -2*m/math.tan(th)*Y(l,m,th,ph)
        -(sq(l*(l+1)-m*(m+1))*cmath.exp(-1j*ph)*Y(l,m+1,th,ph)
          +sq(l*(l+1)-m*(m-1))*cmath.exp(1j*ph)*Y(l,m-1,th,ph)))
    ok &= chk("5.7.2 cos th Y", lambda l,m,th,ph:
        math.cos(th)*Y(l,m,th,ph)
        -(sq((l-m+1)*(l+m+1)/((2*l+1)*(2*l+3)))*Y(l+1,m,th,ph)
          +sq((l-m)*(l+m)/((2*l-1)*(2*l+1)))*Y(l-1,m,th,ph)))
    ok &= chk("5.7.3 sin th e^{-if} Y", lambda l,m,th,ph:
        math.sin(th)*cmath.exp(-1j*ph)*Y(l,m,th,ph)
        -(sq((l-m+1)*(l-m+2)/((2*l+1)*(2*l+3)))*Y(l+1,m-1,th,ph)
          -sq((l+m-1)*(l+m)/((2*l-1)*(2*l+1)))*Y(l-1,m-1,th,ph)))
    ok &= chk("5.7.4 sin th e^{if} Y", lambda l,m,th,ph:
        math.sin(th)*cmath.exp(1j*ph)*Y(l,m,th,ph)
        -(-sq((l+m+1)*(l+m+2)/((2*l+1)*(2*l+3)))*Y(l+1,m+1,th,ph)
          +sq((l-m-1)*(l-m)/((2*l-1)*(2*l+1)))*Y(l-1,m+1,th,ph)))
    ok &= chk("5.7.5 cos^2 th Y", lambda l,m,th,ph:
        (2*l-1)*(2*l+3)*math.cos(th)**2*Y(l,m,th,ph)
        -((2*l-1)*sq(((l+1)**2-m**2)/(2*l+1)*((l+2)**2-m**2)/(2*l+5))*Y(l+2,m,th,ph)
          +(2*l*(l+1)-2*m**2-1)*Y(l,m,th,ph)
          +(2*l+3)*sq((l**2-m**2)/(2*l+1)*((l-1)**2-m**2)/(2*l-3))*Y(l-2,m,th,ph)))
    ok &= chk("5.7.6 sincos e^{if} Y", lambda l,m,th,ph:
        (2*l-1)*(2*l+3)*math.sin(th)*math.cos(th)*cmath.exp(1j*ph)*Y(l,m,th,ph)
        -(-(2*l-1)*sq(((l+1)**2-m**2)/(2*l+1)*(l+m+2)*(l+m+3)/(2*l+5))*Y(l+2,m+1,th,ph)
          -(2*m+1)*sq(l*(l+1)-m*(m+1))*Y(l,m+1,th,ph)
          +(2*l+3)*sq((l**2-m**2)/(2*l+1)*(l-m-1)*(l-m-2)/(2*l-3))*Y(l-2,m+1,th,ph)))
    ok &= chk("5.7.7 sincos e^{-if} Y", lambda l,m,th,ph:
        (2*l-1)*(2*l+3)*math.sin(th)*math.cos(th)*cmath.exp(-1j*ph)*Y(l,m,th,ph)
        -((2*l-1)*sq(((l+1)**2-m**2)/(2*l+1)*(l-m+2)*(l-m+3)/(2*l+5))*Y(l+2,m-1,th,ph)
          -(2*m-1)*sq(l*(l+1)-m*(m-1))*Y(l,m-1,th,ph)
          -(2*l+3)*sq((l**2-m**2)/(2*l+1)*(l+m-1)*(l+m-2)/(2*l-3))*Y(l-2,m-1,th,ph)))
    def facr(a, b):  # a!/b!
        r = 1.0
        for k in range(int(b)+1, int(a)+1): r *= k
        return r
    ok &= chk("5.7.8 sin^2 e^{2if} Y", lambda l,m,th,ph:
        (2*l-1)*(2*l+3)*math.sin(th)**2*cmath.exp(2j*ph)*Y(l,m,th,ph)
        -((2*l-1)/math.sqrt((2*l+1)*(2*l+5))*sq(facr(l+m+4,l+m))*Y(l+2,m+2,th,ph)
          +(2*l+3)/math.sqrt((2*l+1)*(2*l-3))*sq(facr(l-m,l-m-4))*Y(l-2,m+2,th,ph)
          -2*sq(facr(l+m+2,l+m)*facr(l-m,l-m-2))*Y(l,m+2,th,ph)))
    # 5.7.9: as printed has (2l+1)(2l+3) in a denom that (by symmetry with 5.7.8)
    # should be (2l+1)(2l-3); test the CORRECTED version.
    ok &= chk("5.7.9 sin^2 e^{-2if} Y [corr]", lambda l,m,th,ph:
        (2*l-1)*(2*l+3)*math.sin(th)**2*cmath.exp(-2j*ph)*Y(l,m,th,ph)
        -((2*l-1)/math.sqrt((2*l+1)*(2*l+5))*sq(facr(l-m+4,l-m))*Y(l+2,m-2,th,ph)
          +(2*l+3)/math.sqrt((2*l+1)*(2*l-3))*sq(facr(l+m,l+m-4))*Y(l-2,m-2,th,ph)
          -2*sq(facr(l-m+2,l-m)*facr(l+m,l+m-2))*Y(l,m-2,th,ph)))

    print("\nSec 5.8 derivatives")
    h = 1e-6
    def dth(l,m,th,ph): return (Y(l,m,th+h,ph)-Y(l,m,th-h,ph))/(2*h)
    def dph(l,m,th,ph): return (Y(l,m,th,ph+h)-Y(l,m,th,ph-h))/(2*h)
    ok &= chk("5.8.4 dY/dph = im Y", lambda l,m,th,ph: dph(l,m,th,ph)-1j*m*Y(l,m,th,ph))
    ok &= chk("5.8.5a dY/dth (m+1)", lambda l,m,th,ph:
        dth(l,m,th,ph)-(m/math.tan(th)*Y(l,m,th,ph)+sq(l*(l+1)-m*(m+1))*Y(l,m+1,th,ph)*cmath.exp(-1j*ph)))
    ok &= chk("5.8.5b dY/dth (m-1)", lambda l,m,th,ph:
        dth(l,m,th,ph)-(-m/math.tan(th)*Y(l,m,th,ph)-sq(l*(l+1)-m*(m-1))*Y(l,m-1,th,ph)*cmath.exp(1j*ph)))
    ok &= chk("5.8.5c dY/dth (sym)", lambda l,m,th,ph:
        dth(l,m,th,ph)-(0.5*sq(l*(l+1)-m*(m+1))*Y(l,m+1,th,ph)*cmath.exp(-1j*ph)
                        -0.5*sq(l*(l+1)-m*(m-1))*Y(l,m-1,th,ph)*cmath.exp(1j*ph)))
    ok &= chk("5.8.6a sin th dY/dth (l-1)", lambda l,m,th,ph:
        math.sin(th)*dth(l,m,th,ph)-(l*math.cos(th)*Y(l,m,th,ph)
        -sq((2*l+1)/(2*l-1)*(l**2-m**2))*Y(l-1,m,th,ph)))
    ok &= chk("5.8.6b sin th dY/dth (l+1)", lambda l,m,th,ph:
        math.sin(th)*dth(l,m,th,ph)-(-(l+1)*math.cos(th)*Y(l,m,th,ph)
        +sq((2*l+1)/(2*l+3)*((l+1)**2-m**2))*Y(l+1,m,th,ph)))
    print("\nRESULT:", "ALL PASS" if ok else "FAILURES ABOVE")
    return 0 if ok else 1

if __name__ == "__main__":
    raise SystemExit(main())
