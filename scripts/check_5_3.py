#!/usr/bin/env python3
r"""Verify Sec 5.3 (integral representations of Y_lm) vs mpmath.spherharm.
Integrals evaluated with mpmath.quad / quadosc.  Each Y_lm(th,ph)=e^{i m ph}
A_lm(th); most reps carry the e^{i m ph} explicitly, so we compare to
spherharm(l,m,th,ph) directly.
"""
import mpmath as mp
mp.mp.dps = 25
pi = mp.pi
def Y(l,m,th,ph): return mp.spherharm(l,m,th,ph)
def fac(n): return mp.factorial(int(round(n)))
def dfac(n):
    n=int(round(n))
    if n<=0: return mp.mpf(1)
    r=mp.mpf(1)
    while n>1: r*=n; n-=2
    return r
def report(tag,w,tol=mp.mpf('1e-10')):
    ok=w<tol; print(f"  {tag:36s} {'PASS' if ok else 'FAIL'}  worst={mp.nstr(w,3)}"); return ok

TH = [mp.mpf('0.5'), mp.mpf('1.0'), mp.mpf('2.0')]   # last > pi/2
THlo = [mp.mpf('0.5'), mp.mpf('1.0')]                # < pi/2 (for e^{-k cos} etc.)
PH = mp.mpf('0.7')
ok = True

print("Sec 5.3.1 indefinite integrals")
# 5.3.2 (m=0): P_l Mehler-Dirichlet, cos form.  (m>=1 is a finite-part integral.)
w = mp.mpf(0)
for l in range(0,6):
    for th in TH:
        integ = mp.quad(lambda ps: mp.cos((2*l+1)*ps/2)/mp.sqrt(mp.cos(ps)-mp.cos(th)), [0, th])
        rhs = mp.sqrt(2)/pi*mp.sqrt((2*l+1)/(4*pi))*integ
        w = max(w, abs(rhs - Y(l,0,th,0)))
ok &= report("5.3.2 Mehler-Dirichlet cos (m=0)", w)
# 5.3.3 (m=0): sin form, integral th..pi
w = mp.mpf(0)
for l in range(0,6):
    for th in TH:
        integ = mp.quad(lambda ps: mp.sin((2*l+1)*ps/2)/mp.sqrt(mp.cos(th)-mp.cos(ps)), [th, pi])
        rhs = mp.sqrt(2)/pi*mp.sqrt((2*l+1)/(4*pi))*integ
        w = max(w, abs(rhs - Y(l,0,th,0)))
ok &= report("5.3.3 Mehler-Dirichlet sin (m=0)", w)
# 5.3.4  (m>0): int_{cos th}^1 P_l(x)(x-cos th)^{m-1} dx
w = mp.mpf(0)
for l in range(1,6):
    for m in range(1,l+1):
        for th in TH:
            c = mp.cos(th)
            integ = mp.quad(lambda x,l=l,m=m,c=c: mp.legendre(l,x)*(x-c)**(m-1), [c, 1])
            rhs = ((-1)**m*mp.e**(1j*m*PH)*mp.sqrt((2*l+1)/(4*pi)*fac(l+m)/fac(l-m))
                   /(fac(m-1)*mp.sin(th)**m)*integ)
            w = max(w, abs(rhs - Y(l,m,th,PH)))
ok &= report("5.3.4 P_l (cos-cos)^{m-1}", w)

print("\nSec 5.3.2 definite integrals")
# 5.3.5  (+i)^m/pi e sqrt((2l+1)/4pi (l+m)!(l-m)!)/l! int_0^pi [cos +i sin cos ps]^l cos(m ps) dps
for sgn in (+1,-1):
    w = mp.mpf(0)
    for l in range(0,5):
        for m in range(0,l+1):
            for th in TH:
                integ = mp.quad(lambda ps,l=l,m=m,th=th:
                    (mp.cos(th)+sgn*1j*mp.sin(th)*mp.cos(ps))**l*mp.cos(m*ps), [0, pi])
                rhs = ((sgn*1j)**m/pi*mp.e**(1j*m*PH)
                       *mp.sqrt((2*l+1)/(4*pi)*fac(l+m)*fac(l-m))/fac(l)*integ)
                w = max(w, abs(rhs - Y(l,m,th,PH)))
    ok &= report(f"5.3.5 [cos+-i sin cos]^l  (sgn={sgn:+d})", w)
# 5.3.6  (+i)^m/pi e sqrt((2l+1)/(4pi(l+m)!(l-m)!)) l! int_0^pi cos(m ps)/(cos -+i sin cos)^{l+1} dps
for sgn in (+1,-1):
    w = mp.mpf(0)
    for l in range(0,5):
        for m in range(0,l+1):
            for th in THlo:
                integ = mp.quad(lambda ps,l=l,m=m,th=th:
                    mp.cos(m*ps)/(mp.cos(th)-sgn*1j*mp.sin(th)*mp.cos(ps))**(l+1), [0, pi])
                rhs = ((sgn*1j)**m/pi*mp.e**(1j*m*PH)
                       *mp.sqrt((2*l+1)/(4*pi*fac(l+m)*fac(l-m)))*fac(l)*integ)
                w = max(w, abs(rhs - Y(l,m,th,PH)))
    ok &= report(f"5.3.6 1/(cos-+i sin cos)^{{l+1}} (sgn={sgn:+d})", w)
# 5.3.7  (+i)^m/2pi sqrt((2l+1)/4pi (l+m)!(l-m)!)/l! int_0^2pi [cos +i sin cos(ps-ph)]^l e^{im ps} dps
for sgn in (+1,-1):
    w = mp.mpf(0)
    for l in range(0,5):
        for m in range(-l,l+1):
            for th in TH:
                integ = mp.quad(lambda ps,l=l,m=m,th=th:
                    (mp.cos(th)+sgn*1j*mp.sin(th)*mp.cos(ps-PH))**l*mp.e**(1j*m*ps), [0, 2*pi])
                rhs = (sgn*1j)**m/(2*pi)*mp.sqrt((2*l+1)/(4*pi)*fac(l+m)*fac(l-m))/fac(l)*integ
                w = max(w, abs(rhs - Y(l,m,th,PH)))
    ok &= report(f"5.3.7 2pi form (sgn={sgn:+d})", w)
# 5.3.8  (+i)^m/2pi sqrt((2l+1)/(4pi(l+m)!(l-m)!)) l! int_0^2pi [cos -+i sin cos(ps-ph)]^{-l-1} e^{im ps} dps
for sgn in (+1,-1):
    w = mp.mpf(0)
    for l in range(0,5):
        for m in range(-l,l+1):
            for th in THlo:
                integ = mp.quad(lambda ps,l=l,m=m,th=th:
                    (mp.cos(th)-sgn*1j*mp.sin(th)*mp.cos(ps-PH))**(-l-1)*mp.e**(1j*m*ps), [0, 2*pi])
                rhs = (sgn*1j)**m/(2*pi)*mp.sqrt((2*l+1)/(4*pi*fac(l+m)*fac(l-m)))*fac(l)*integ
                w = max(w, abs(rhs - Y(l,m,th,PH)))
    ok &= report(f"5.3.8 2pi inverse (sgn={sgn:+d})", w)
# 5.3.9  (-1)^m/pi e sqrt((2l+1)/4pi (l+m)!/(l-m)!) (sin)^m/(2m-1)!! int_0^pi (cos+i sin cos ps)^{l-m}(sin ps)^{2m} dps
for sgn in (+1,-1):
    w = mp.mpf(0)
    for l in range(0,5):
        for m in range(0,l+1):
            for th in TH:
                integ = mp.quad(lambda ps,l=l,m=m,th=th:
                    (mp.cos(th)+sgn*1j*mp.sin(th)*mp.cos(ps))**(l-m)*mp.sin(ps)**(2*m), [0, pi])
                rhs = ((-1)**m/pi*mp.e**(1j*m*PH)*mp.sqrt((2*l+1)/(4*pi)*fac(l+m)/fac(l-m))
                       *mp.sin(th)**m/dfac(2*m-1)*integ)
                w = max(w, abs(rhs - Y(l,m,th,PH)))
    ok &= report(f"5.3.9 (sin ps)^2m [corr (l+m)!] (sgn={sgn:+d})", w)
# 5.3.10  (-1)^m/pi e sqrt((2l+1)/4pi (l+m)!/(l-m)!) (sin)^m/(2m-1)!! int_0^pi (sin chi)^{2m}/(cos -+i sin cos chi)^{l+m+1} dchi
for sgn in (+1,-1):
    w = mp.mpf(0)
    for l in range(0,5):
        for m in range(0,l+1):
            for th in THlo:
                integ = mp.quad(lambda ch,l=l,m=m,th=th:
                    mp.sin(ch)**(2*m)/(mp.cos(th)-sgn*1j*mp.sin(th)*mp.cos(ch))**(l+m+1), [0, pi])
                rhs = ((-1)**m/pi*mp.e**(1j*m*PH)*mp.sqrt((2*l+1)/(4*pi)*fac(l+m)/fac(l-m))
                       *mp.sin(th)**m/dfac(2*m-1)*integ)
                w = max(w, abs(rhs - Y(l,m,th,PH)))
    ok &= report(f"5.3.10 1/(..)^{{l+m+1}} [corr] (sgn={sgn:+d})", w)

print("\nSec 5.3.3 improper integrals")
# 5.3.12  i^{m+1}/pi e sqrt((2l+1)/(4pi(l+m)!(l-m)!)) l! {(-1)^m int cosh(mt)/(cos+i sin cosh)^{l+1} - int cosh(mt)/(cos-i sin cosh)^{l+1}}
w = mp.mpf(0)
for l in range(1,5):
    for m in range(0,l+1):
        for th in TH:
            I1 = mp.quad(lambda t,l=l,m=m,th=th: mp.cosh(m*t)/(mp.cos(th)+1j*mp.sin(th)*mp.cosh(t))**(l+1), [0, mp.inf])
            I2 = mp.quad(lambda t,l=l,m=m,th=th: mp.cosh(m*t)/(mp.cos(th)-1j*mp.sin(th)*mp.cosh(t))**(l+1), [0, mp.inf])
            rhs = (1j**(m+1)/pi*mp.e**(1j*m*PH)*mp.sqrt((2*l+1)/(4*pi*fac(l+m)*fac(l-m)))*fac(l)
                   *((-1)**m*I1 - I2))
            w = max(w, abs(rhs - Y(l,m,th,PH)))
ok &= report("5.3.12 cosh(mt) improper", w)
# 5.3.13  i/pi e sqrt((2l+1)/4pi (l+m)!/(l-m)!) (sin)^m/(2m-1)!! {int (sinh)^2m/(cos+i sin cosh)^{l+m+1} - int (sinh)^2m/(cos-i sin cosh)^{l+m+1}}
w = mp.mpf(0)
for l in range(1,5):
    for m in range(0,l+1):
        for th in TH:
            I1 = mp.quad(lambda t,l=l,m=m,th=th: mp.sinh(t)**(2*m)/(mp.cos(th)+1j*mp.sin(th)*mp.cosh(t))**(l+m+1), [0, mp.inf])
            I2 = mp.quad(lambda t,l=l,m=m,th=th: mp.sinh(t)**(2*m)/(mp.cos(th)-1j*mp.sin(th)*mp.cosh(t))**(l+m+1), [0, mp.inf])
            rhs = (1j/pi*mp.e**(1j*m*PH)*mp.sqrt((2*l+1)/(4*pi)*fac(l+m)/fac(l-m))
                   *mp.sin(th)**m/dfac(2*m-1)*(I1 - I2))
            w = max(w, abs(rhs - Y(l,m,th,PH)))
ok &= report("5.3.13 (sinh)^2m improper", w)
# 5.3.14  (-1)^m e sqrt((2l+1)/(4pi(l+m)!(l-m)!)) int_0^inf e^{-k cos} J_m(k sin) k^l dk   (cos>0)
w = mp.mpf(0)
for l in range(0,5):
    for m in range(0,l+1):
        for th in THlo:
            integ = mp.quad(lambda k,l=l,m=m,th=th: mp.e**(-k*mp.cos(th))*mp.besselj(m,k*mp.sin(th))*k**l, [0, mp.inf])
            rhs = (-1)**m*mp.e**(1j*m*PH)*mp.sqrt((2*l+1)/(4*pi*fac(l+m)*fac(l-m)))*integ
            w = max(w, abs(rhs - Y(l,m,th,PH)))
ok &= report("5.3.14 e^{-k cos} J_m k^l", w)
# 5.3.15  |Y_lm|^2 = (2l+1)/4pi int_0^inf [J_m(t sin/2)]^2 J_{2l+1}(t) dt
# (oscillatory quadrature -- slow; pass "skip15" as argv to omit)
import sys
if 'skip15' not in sys.argv:
    w = mp.mpf(0)
    for l in range(0,4):
        for m in range(0,l+1):
            for th in TH:                 # valid for all theta (incl. > pi/2)
                integ = mp.quadosc(lambda t,l=l,m=m,th=th: mp.besselj(m,t*mp.sin(th)/2)**2*mp.besselj(2*l+1,t),
                                   [0, mp.inf], period=2*pi)
                rhs = (2*l+1)/(4*pi)*integ
                w = max(w, abs(rhs - abs(Y(l,m,th,PH))**2))
    ok &= report("5.3.15 |Y|^2 Bessel integral", w, tol=mp.mpf('1e-6'))

print("\nRESULT:", "ALL PASS" if ok else "FAILURES ABOVE")
raise SystemExit(0 if ok else 1)
