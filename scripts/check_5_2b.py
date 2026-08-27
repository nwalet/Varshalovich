#!/usr/bin/env python3
r"""Verify the remaining Sec 5.2 forms (scan-read, printed pp.135-138) vs
mpmath.spherharm, before de-garbling the .tex.
  5.2.15/16 (trig th/2 multline), 5.2.18/19/21/22 (trig th),
  5.2.29/30/32 (hypergeom th), 5.2.35/36 (exponential), 5.2.38 (Jacobi).
"""
import math, cmath
import mpmath as mp
from scipy.special import eval_jacobi
mp.mp.dps = 30
TOL = 1e-8
TH = [0.4, 0.9, 1.3, 1.9, 2.6]
PH = [0.0, 1.8]

def Y(l, m, th, ph):
    return complex(mp.spherharm(l, m, th, ph)) if abs(m) <= l else 0j
def fac(n):
    n = int(round(n)); return math.factorial(n) if n >= 0 else math.inf
def dfac(n):
    n = int(round(n))
    if n < -1: return math.inf
    r = 1.0
    while n > 1: r *= n; n -= 2
    return r
def xi(m): return (-1)**m if m > 0 else 1
def F(a, b, c, z): return complex(mp.hyp2f1(a, b, c, z))
def report(tag, w, tol=TOL): print(f"  {tag:34s} {'PASS' if w<tol else 'FAIL'}  worst={w:.2e}"); return w < tol
def sweep(tag, form, ls=range(0,6), mpos=False, ths=TH, tol=TOL):
    w = 0.0
    for l in ls:
        for m in (range(0,l+1) if mpos else range(-l,l+1)):
            for th in ths:
                for ph in PH:
                    try: v = form(l, m, th, ph)
                    except (ValueError, ZeroDivisionError, OverflowError): continue
                    if not (math.isfinite(v.real) and math.isfinite(v.imag)): continue
                    w = max(w, abs(v - Y(l, m, th, ph)))
    return report(tag, w, tol)

ok = True
e = cmath.exp
# 5.2.15  xi e sqrt((2l+1)/4pi) sqrt((l+m)!(l-m)!l!|m|!) (cos th/2)^{2l}
#          sum_s (-1)^s (tan th/2)^{2s+|m|}/(s!(s+|m|)!(l-s)!(l-|m|-s)!)
def f15(l, m, th, ph):
    am = abs(m); c2 = math.cos(th/2)
    tot = 0.0
    for s in range(0, l-am+1):
        if fac(l-am-s) == math.inf or fac(l-s) == math.inf: continue
        tot += (-1)**s*math.tan(th/2)**(2*s+am)/(fac(s)*fac(s+am)*fac(l-s)*fac(l-am-s))
    return (xi(m)*e(1j*m*ph)*math.sqrt((2*l+1)/(4*math.pi))
            * math.sqrt(fac(l+m)*fac(l-m)*fac(l)*fac(am))*c2**(2*l)*tot)
ok &= sweep("5.2.15 tan^{2s+|m|} multline", f15)

# 5.2.16  xi (-1)^{l-m} e sqrt((2l+1)/4pi) sqrt((l+m)!(l-m)!l!|m|!) (sin th/2)^{2l}
#          sum_s (-1)^s (cot th/2)^{2s+|m|}/(...)
def f16(l, m, th, ph):
    am = abs(m); s2 = math.sin(th/2)
    tot = 0.0
    for s in range(0, l-am+1):
        if fac(l-am-s) == math.inf or fac(l-s) == math.inf: continue
        tot += (-1)**s*(1/math.tan(th/2))**(2*s+am)/(fac(s)*fac(s+am)*fac(l-s)*fac(l-am-s))
    return (xi(m)*(-1)**(l-m)*e(1j*m*ph)*math.sqrt((2*l+1)/(4*math.pi))
            * math.sqrt(fac(l+m)*fac(l-m)*fac(l)*fac(am))*s2**(2*l)*tot)
ok &= sweep("5.2.16 cot^{2s+|m|} multline", f16)

# 5.2.18  e sqrt((2l+1)/4pi/((l+m)!(l-m)!)) (sin th)^l {even/odd}
def f18(l, m, th, ph):
    am = abs(m); s = math.sin(th)
    pre = e(1j*m*ph)*math.sqrt((2*l+1)/(4*math.pi*fac(l+m)*fac(l-m)))*s**l
    if (l-m) % 2 == 0:
        tot = 0.0
        for k in range(0, l-am+1, 2):
            tot += ((-1)**((l+m-k)//2)*dfac(l+m)*dfac(l-m)/(dfac(l+m-k)*dfac(l-m-k))
                    * dfac(2*l-k-1)/(fac(k)*s**k))
        return pre*tot
    else:
        tot = 0.0
        for k in range(1, l-am+1, 2):
            tot += ((-1)**((l+m-k)//2)*dfac(l+m-1)*dfac(l-m-1)/(dfac(l+m-k)*dfac(l-m-k))
                    * dfac(2*l-k)/(dfac(k-1)*s**k))
        return pre*math.cos(th)*tot
ok &= sweep("5.2.18 (sin th)^l series", f18)

# 5.2.21  e sqrt((2l+1)/4pi (l+m)!(l-m)!) (sin th)^l sum_s (-1)^{(l+m-s)/2}
#          /((l+m-s)!!(l-m-s)!!) (cot th)^s/s!   (l+m-s even)
def f21(l, m, th, ph):
    s0 = math.sin(th)
    tot = 0.0
    for k in range(0, l+1):
        if (l+m-k) % 2 != 0: continue
        if dfac(l+m-k) == math.inf or dfac(l-m-k) == math.inf: continue
        tot += (-1)**((l+m-k)//2)/(dfac(l+m-k)*dfac(l-m-k))*(1/math.tan(th))**k/fac(k)
    return e(1j*m*ph)*math.sqrt((2*l+1)/(4*math.pi)*fac(l+m)*fac(l-m))*s0**l*tot
ok &= sweep("5.2.21 (cot th)^s series", f21, ths=[0.4,0.9,1.3])

# 5.2.22  |Y|^2 = (2l+1)/4pi sum_{s=|m|,|m|+1..}^l (-1)^{s+m} (l+s)!/(l-s)! (2s-1)!!/(2s)!! (sin th)^{2s}/((s-m)!(s+m)!)
def f22(l, m, th):
    tot = 0.0
    for s in range(abs(m), l+1):
        tot += ((-1)**(s+m)*fac(l+s)/fac(l-s)*dfac(2*s-1)/dfac(2*s)
                * math.sin(th)**(2*s)/(fac(s-m)*fac(s+m)))
    return (2*l+1)/(4*math.pi)*tot
w = 0.0
for l in range(0,5):
    for m in range(-l,l+1):
        for th in TH:
            w = max(w, abs(f22(l,m,th) - abs(Y(l,m,th,0.7))**2))
ok &= report("5.2.22 |Y|^2", w)

# 5.2.29  xi e sqrt((2l+1)/4pi (l+|m|)!/(l-|m|)!) (sin th)^|m|/(2^|m| |m|!) {even/odd}
def f29(l, m, th, ph):
    am = abs(m)
    pre = (xi(m)*e(1j*m*ph)*math.sqrt((2*l+1)/(4*math.pi)*fac(l+am)/fac(l-am))
           * math.sin(th)**am/(2**am*fac(am)))
    if (l+m) % 2 == 0:
        return pre*F(-(l-am)/2, (l+am+1)/2, am+1, math.sin(th)**2)
    else:
        return pre*math.cos(th)*F(-(l-am-1)/2, (l+am+2)/2, am+1, math.sin(th)**2)
ok &= sweep("5.2.29 2F1(sin^2 th)", f29)

# 5.2.30  e sqrt((2l+1)/4pi/((l+m)!(l-m)!)) (2l-1)!! (sin th)^l {even/odd} 1/sin^2 th
def f30(l, m, th, ph):
    pre = e(1j*m*ph)*math.sqrt((2*l+1)/(4*math.pi*fac(l+m)*fac(l-m)))*dfac(2*l-1)*math.sin(th)**l
    z = 1/math.sin(th)**2
    if (l+m) % 2 == 0:
        return pre*(-1)**((l+m)//2)*F(-(l+m)/2, -(l-m)/2, -(2*l-1)/2, z)
    else:
        return pre*(-1)**((l+m-1)//2)*(1/math.tan(th))*F(-(l+m-1)/2, -(l-m-1)/2, -(2*l-1)/2, z)
ok &= sweep("5.2.30 2F1(1/sin^2 th)", f30, ths=[0.4,0.9,1.3,1.9])

# 5.2.32  (-1)^m e sqrt((2l+1)/4pi/((l+m)!(l-m)!)) (2l-1)!! (cos th)^l (tan th)^m F(..;1/cos^2 th)
def f32(l, m, th, ph):
    return ((-1)**m*e(1j*m*ph)*math.sqrt((2*l+1)/(4*math.pi*fac(l+m)*fac(l-m)))
            * dfac(2*l-1)*math.cos(th)**l*math.tan(th)**m
            * F(-(l-m)/2, -(l-m-1)/2, -(2*l-1)/2, 1/math.cos(th)**2))
ok &= sweep("5.2.32 2F1(1/cos^2 th) m>=0", f32, mpos=True, ths=[0.4,0.9,1.3])

# 5.2.35  -i e/pi sqrt((2l+1)/4pi (l+m)!(l-m)!) 2^{l+m+1}(sin th)^m/(2l+1)!! {exp F - exp F}
def f35(l, m, th, ph):
    pre = (-1j*e(1j*m*ph)/math.pi*math.sqrt((2*l+1)/(4*math.pi)*fac(l+m)*fac(l-m))
           * 2**(l+m+1)*math.sin(th)**m/dfac(2*l+1))
    A = e(-1j*(l+m+1)*th)*F(m+0.5, l+m+1, l+1.5, e(-2j*th))
    B = e(1j*(l+m+1)*th)*F(m+0.5, l+m+1, l+1.5, e(2j*th))
    return pre*(A - B)
ok &= sweep("5.2.35 exp 2F1(e^{-+2i th})", f35, ls=range(0,5), tol=1e-6)

# 5.2.38  xi e/(2^|m| l!) sqrt((2l+1)/4pi (l+m)!(l-m)!) (sin th)^|m| P^{(|m|,|m|)}_{l-|m|}(cos th)
def f38(l, m, th, ph):
    am = abs(m)
    return (xi(m)*e(1j*m*ph)/(2**am*fac(l))
            * math.sqrt((2*l+1)/(4*math.pi)*fac(l+m)*fac(l-m))
            * math.sin(th)**am*eval_jacobi(l-am, am, am, math.cos(th)))
ok &= sweep("5.2.38 Jacobi P^{(|m|,|m|)}", f38)

print("\nRESULT:", "ALL PASS" if ok else "FAILURES ABOVE")
raise SystemExit(0 if ok else 1)
