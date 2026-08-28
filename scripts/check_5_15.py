#!/usr/bin/env python3
r"""Verify Sec 5.15 zero/extremum formulas vs mpmath.spherharm.
5.15.2 gives cos^2(theta_alpha) for zeros of Y_lm; 5.15.4 gives
cos^2(theta_beta) for zeros of dY_lm/dtheta.  We plug each predicted cos^2 value
back and confirm Y (resp. dY/dtheta) vanishes there.  Also checks the zero
COUNTS stated in 5.15.1.
"""
import math
import mpmath as mp
mp.mp.dps = 30
pi = mp.pi
def Y(l,m,th): return mp.spherharm(l,m,th,0)
def dY(l,m,th):
    h=mp.mpf('1e-8')
    return (Y(l,m,th+h)-Y(l,m,th-h))/(2*h)
def report(tag,w,tol=mp.mpf('1e-7')):
    w=float(w); ok=w<tol; print(f"  {tag:44s} {'PASS' if ok else 'FAIL'}  worst={w:.2e}"); return ok
ok=True

def check_zeros(tag, cos2_list_fn, func, mrel):
    """cos2_list_fn(l)->list of cos^2 theta in (0,1) (interior, excl 0 and 1);
    func = Y or dY; mrel: m = l - mrel."""
    w=mp.mpf(0)
    for l in range(mrel+ (1 if func is dY else 1), 12):
        m=l-mrel
        if m<1: continue
        for v in cos2_list_fn(l):
            if v<=0 or v>=1: continue
            th=mp.acos(mp.sqrt(v))
            # normalize by the harmonic's scale so "zero" is meaningful
            scale=max(abs(func(l,m,mp.mpf('0.5'))), abs(func(l,m,mp.mpf('1.2'))), mp.mpf('1e-6'))
            w=max(w, abs(func(l,m,th))/scale)
    return report(tag,w)

print("Sec 5.15.2  zeros of Y_lm  (cos^2 theta_alpha)")
# m=+-(l-2): 1/(2l-1)
ok&=check_zeros("5.15.2 m=l-2: 1/(2l-1)", lambda l:[mp.mpf(1)/(2*l-1)], Y, 2)
# m=+-(l-3): 3/(2l-1)  (plus 0 boundary)
ok&=check_zeros("5.15.2 m=l-3: 3/(2l-1)", lambda l:[mp.mpf(3)/(2*l-1)], Y, 3)
# m=+-(l-4): [3 +- 2 sqrt(3(l-2)/(2l-3))]/(2l-1)
def a4(l):
    d=2*mp.sqrt(3*(l-2)/mp.mpf(2*l-3))
    return [(3+d)/(2*l-1),(3-d)/(2*l-1)]
ok&=check_zeros("5.15.2 m=l-4: [3+-2 sqrt(3(l-2)/(2l-3))]/(2l-1)", a4, Y, 4)
# m=+-(l-5): [5 +- 2 sqrt(5(l-3)/(2l-3))]/(2l-1)  (plus 0)
def a5(l):
    d=2*mp.sqrt(5*(l-3)/mp.mpf(2*l-3))
    return [(5+d)/(2*l-1),(5-d)/(2*l-1)]
ok&=check_zeros("5.15.2 m=l-5: [5+-2 sqrt(5(l-3)/(2l-3))]/(2l-1)", a5, Y, 5)

print("\nSec 5.15.4  zeros of dY/dtheta  (cos^2 theta_beta)")
# m=+-(l-1): 1/l
ok&=check_zeros("5.15.4 m=l-1: 1/l", lambda l:[mp.mpf(1)/l], dY, 1)
# m=+-(l-2): (5l-4)/(l(2l-1))
ok&=check_zeros("5.15.4 m=l-2: (5l-4)/(l(2l-1))", lambda l:[mp.mpf(5*l-4)/(l*(2*l-1))], dY, 2)
# m=+-(l-3): [9(l-1) +- sqrt(3(19l^2-50l+27))]/(2l(2l-1))
def b3(l):
    d=mp.sqrt(3*(19*l*l-50*l+27))
    return [(9*(l-1)+d)/(2*l*(2*l-1)),(9*(l-1)-d)/(2*l*(2*l-1))]
ok&=check_zeros("5.15.4 m=l-3: [9(l-1)+-sqrt(3(19l^2-50l+27))]/(2l(2l-1))", b3, dY, 3)
# m=+-(l-4): [7l-8 +- 2 sqrt((11l^3-62l^2+104l-48)/(2l-3))]/(l(2l-1))
def b4(l):
    d=2*mp.sqrt((11*l**3-62*l*l+104*l-48)/mp.mpf(2*l-3))
    return [(7*l-8+d)/(l*(2*l-1)),(7*l-8-d)/(l*(2*l-1))]
ok&=check_zeros("5.15.4 m=l-4: [7l-8+-2 sqrt(...)]/(l(2l-1))", b4, dY, 4)

print("\nSec 5.15.1  zero counts")
# m!=0: Y_lm has (l-|m|) interior zeros in (0,pi)
def count_sign_changes(f,l,m,N=4000):
    prev=None; c=0
    for i in range(1,N):
        th=pi*i/N
        v=f(l,m,th)
        s=1 if v.real>0 else -1
        if prev is not None and s!=prev: c+=1
        prev=s
    return c
w=0
for l in range(1,7):
    for m in range(1,l+1):
        c=count_sign_changes(lambda l,m,th:Y(l,m,th),l,m)
        if c!=l-m: w=max(w,1); print(f"    Y count mismatch l={l} m={m}: got {c}, expect {l-m}")
report("5.15.1 Y_lm interior zero count = l-|m|", w, tol=0.5)
ok&= (w==0)

print("\nRESULT:", "ALL PASS" if ok else "FAILURES ABOVE")
raise SystemExit(0 if ok else 1)
