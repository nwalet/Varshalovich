#!/usr/bin/env python3
r"""Verify Sec 5.14 (Y and dY/dtheta at special theta = 0, pi/2, pi, n*pi)
vs mpmath.spherharm.  Derivatives via high-precision finite differences:
one-sided 2nd-order at the poles, centered at pi/2.
"""
import mpmath as mp
mp.mp.dps = 40
TOL = mp.mpf('1e-9')
PH = [mp.mpf('0.3'), mp.mpf('1.8')]
pi = mp.pi

def Y(l, m, th, ph):
    if abs(m) > l: return mp.mpc(0)
    return mp.spherharm(l, m, th, ph)
def dfac(n):
    n = int(round(n))
    if n <= 0: return mp.mpf(1)
    r = mp.mpf(1)
    while n > 1: r *= n; n -= 2
    return r
def report(tag, w):
    ok = w < TOL; print(f"  {tag:34s} {'PASS' if ok else 'FAIL'}  worst={mp.nstr(w,3)}"); return ok

# high-precision derivatives
H = mp.mpf('1e-6')
def dth_c(l, m, th, ph):  # centered (interior)
    return (Y(l,m,th+H,ph) - Y(l,m,th-H,ph)) / (2*H)
def dth_fwd(l, m, th, ph):  # one-sided 2nd order forward
    return (-3*Y(l,m,th,ph) + 4*Y(l,m,th+H,ph) - Y(l,m,th+2*H,ph)) / (2*H)
def dth_bwd(l, m, th, ph):  # one-sided 2nd order backward
    return (3*Y(l,m,th,ph) - 4*Y(l,m,th-H,ph) + Y(l,m,th-2*H,ph)) / (2*H)

ok = True
print("Sec 5.14  Y at special theta")
# 5.14.1  Y(0) = d_{m0} sqrt((2l+1)/4pi)
w = mp.mpf(0)
for l in range(0,7):
    for m in range(-l,l+1):
        for ph in PH:
            rhs = (mp.sqrt((2*l+1)/(4*pi)) if m==0 else mp.mpf(0))
            w = max(w, abs(Y(l,m,0,ph)-rhs))
ok &= report("5.14.1 Y(0)", w)

# 5.14.2  Y(pi/2): l+m even -> (-1)^{(l+m)/2} e^{imph} sqrt((2l+1)/4pi (l+m-1)!!/(l+m)!! (l-m-1)!!/(l-m)!!); else 0
w = mp.mpf(0)
for l in range(0,7):
    for m in range(-l,l+1):
        for ph in PH:
            if (l+m) % 2 == 0:
                rhs = ((-1)**((l+m)//2)*mp.e**(1j*m*ph)
                       * mp.sqrt((2*l+1)/(4*pi)*dfac(l+m-1)/dfac(l+m)*dfac(l-m-1)/dfac(l-m)))
            else:
                rhs = mp.mpf(0)
            w = max(w, abs(Y(l,m,pi/2,ph)-rhs))
ok &= report("5.14.2 Y(pi/2)", w)

# 5.14.3  Y(pi) = d_{m0} (-1)^l sqrt((2l+1)/4pi)
w = mp.mpf(0)
for l in range(0,7):
    for m in range(-l,l+1):
        for ph in PH:
            rhs = ((-1)**l*mp.sqrt((2*l+1)/(4*pi)) if m==0 else mp.mpf(0))
            w = max(w, abs(Y(l,m,pi,ph)-rhs))
ok &= report("5.14.3 Y(pi)", w)

# 5.14.4  Y(+-n pi) = d_{m0} (-1)^{nl} sqrt((2l+1)/4pi)   (m=0 -> phi-indep, real)
w = mp.mpf(0)
for n in (1,2,3):
    for l in range(0,6):
        for ph in PH:
            rhs = (-1)**(n*l)*mp.sqrt((2*l+1)/(4*pi))
            w = max(w, abs(Y(l,0,n*pi,ph)-rhs))
ok &= report("5.14.4 Y(n*pi), m=0", w)

print("\nSec 5.14  dY/dtheta at special theta")
# 5.14.5  dY|_0 = (d_{m,-1}e^{-iph} - d_{m,1}e^{iph}) sqrt(l(l+1)(2l+1)/16pi)
w = mp.mpf(0)
for l in range(1,7):
    for m in range(-l,l+1):
        for ph in PH:
            k = (1 if m==-1 else 0)*mp.e**(-1j*ph) - (1 if m==1 else 0)*mp.e**(1j*ph)
            rhs = k*mp.sqrt(l*(l+1)*(2*l+1)/(16*pi))
            w = max(w, abs(dth_fwd(l,m,0,ph)-rhs))
ok &= report("5.14.5 dY(0)", w)

# 5.14.6  dY|_{pi/2}: l+m odd -> (-1)^{l+(m+1)/2} e^{imph} sqrt((2l+1)/4pi (l+m)!!(l-m)!!/((l+m-1)!!(l-m-1)!!)); else 0
w = mp.mpf(0)
for l in range(0,7):
    for m in range(-l,l+1):
        for ph in PH:
            if (l+m) % 2 == 1:
                # CORRECTED phase (-1)^{(l+m+1)/2} (printed as (-1)^{l+(m+1)/2}, misprint)
                rhs = ((-1)**((l+m+1)//2)*mp.e**(1j*m*ph)
                       * mp.sqrt((2*l+1)/(4*pi)*dfac(l+m)*dfac(l-m)/(dfac(l+m-1)*dfac(l-m-1))))
            else:
                rhs = mp.mpf(0)
            w = max(w, abs(dth_c(l,m,pi/2,ph)-rhs))
ok &= report("5.14.6 dY(pi/2)", w)

# 5.14.7  dY|_pi = (-1)^l (d_{m,-1}e^{-iph} - d_{m,1}e^{iph}) sqrt(l(l+1)(2l+1)/16pi)
w = mp.mpf(0)
for l in range(1,7):
    for m in range(-l,l+1):
        for ph in PH:
            k = (1 if m==-1 else 0)*mp.e**(-1j*ph) - (1 if m==1 else 0)*mp.e**(1j*ph)
            rhs = (-1)**l*k*mp.sqrt(l*(l+1)*(2*l+1)/(16*pi))
            w = max(w, abs(dth_bwd(l,m,pi,ph)-rhs))
ok &= report("5.14.7 dY(pi)", w)

# 5.14.8  dY|_{+-n pi} = (-1)^{nl}(d_{m,-1}e^{-iph} - d_{m,1}e^{iph}) sqrt(l(l+1)(2l+1)/16pi)
# spherharm does not follow VMK's analytic continuation past the poles, so instead
# evaluate dY via the (already-verified) recurrence 5.8.5c at th=n*pi, where
# Y_{l,m'}(n*pi)=0 for m'!=0 and Y_{l,0}(n*pi)=(-1)^{nl} sqrt((2l+1)/4pi).
def dth_rec(l, m, th, ph):
    return (mp.mpf('0.5')*mp.sqrt(l*(l+1)-m*(m+1))*Y(l,m+1,th,ph)*mp.e**(-1j*ph)
            - mp.mpf('0.5')*mp.sqrt(l*(l+1)-m*(m-1))*Y(l,m-1,th,ph)*mp.e**(1j*ph))
w = mp.mpf(0)
for n in (1,2,3):
    for l in range(1,6):
        for m in range(-l,l+1):
            for ph in PH:
                k = (1 if m==-1 else 0)*mp.e**(-1j*ph) - (1 if m==1 else 0)*mp.e**(1j*ph)
                rhs = (-1)**(n*l)*k*mp.sqrt(l*(l+1)*(2*l+1)/(16*pi))
                w = max(w, abs(dth_rec(l,m,n*pi,ph)-rhs))
ok &= report("5.14.8 dY(n*pi) [via 5.8.5c]", w)

print("\nRESULT:", "ALL PASS" if ok else "FAILURES ABOVE")
raise SystemExit(0 if ok else 1)
