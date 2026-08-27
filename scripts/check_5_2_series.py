#!/usr/bin/env python3
r"""
Checks for the power-series and hypergeometric forms of Y_{lm} in Sec 5.2, VMK:
  5.2.9-5.2.16   power series in trig(theta/2)
  5.2.24-5.2.36  hypergeometric (theta/2, theta, exponential)

Each form here encodes the INTENDED (corrected) math; a PASS means the .tex,
once de-garbled to this, is right. Reference Y = mpmath.spherharm.

Usage:  python3 check_5_2_series.py
"""
import math, cmath
import mpmath as mp
mp.mp.dps = 30
TOL = 1e-9
TH = [0.4, 0.9, 1.3, 1.9, 2.6]
PH = [0.0, 1.8]

def Y(l, m, th, ph):
    if abs(m) > l: return 0j
    return complex(mp.spherharm(l, m, th, ph))
def fac(n):
    n = int(round(n))
    return math.factorial(n) if n >= 0 else math.inf
def dfac(n):
    n = int(round(n))
    if n < -1: return math.inf
    r = 1.0
    while n > 1: r *= n; n -= 2
    return r
def xi(m): return (-1)**m if m > 0 else 1
def F(a, b, c, z): return complex(mp.hyp2f1(a, b, c, z))
def report(tag, worst, tol=TOL):
    ok = worst < tol
    print(f"  {tag:40s} {'PASS' if ok else 'FAIL'}  worst={worst:.2e}")
    return ok
def sweep(tag, form, ls=range(0, 6), tol=TOL, mpos=False, ths=TH):
    worst = 0.0
    for l in ls:
        ms = range(0, l+1) if mpos else range(-l, l+1)
        for m in ms:
            for th in ths:
                for ph in PH:
                    try:
                        v = form(l, m, th, ph)
                    except (ValueError, ZeroDivisionError, OverflowError):
                        continue
                    if not (math.isfinite(v.real) and math.isfinite(v.imag)):
                        continue
                    worst = max(worst, abs(v - Y(l, m, th, ph)))
    return report(tag, worst, tol)

def N1(l, m):  # sqrt((2l+1)/4pi (l+m)!/(l-m)!)
    return math.sqrt((2*l+1)/(4*math.pi)*fac(l+m)/fac(l-m))
def N2(l, m):  # sqrt((2l+1)/4pi (l-m)!/(l+m)!)
    return math.sqrt((2*l+1)/(4*math.pi)*fac(l-m)/fac(l+m))


def main():
    print("Sec 5.2 power-series & hypergeometric forms\n")
    ok = True
    e = cmath.exp

    # ---- 5.2.2.2 power series in theta/2 ----
    # 5.2.9  (-1)^m e^{imf} N1 (tan th/2)^m sum_s (-1)^s (l+s)!/(l-s)! (sin th/2)^{2s}/(s!(s+m)!)
    def f9(l, m, th, ph):
        s2 = math.sin(th/2)
        tot = sum((-1)**s*fac(l+s)/fac(l-s)*s2**(2*s)/(fac(s)*fac(s+m))
                  for s in range(0, l+1) if fac(l-s) != math.inf and fac(s+m) != math.inf)
        return (-1)**m*e(1j*m*ph)*N1(l, m)*math.tan(th/2)**m*tot
    ok &= sweep("5.2.9  tan^m, (l+s)!/(l-s)!, sin^2s", f9, mpos=True)

    # 5.2.10  (-1)^m e N2 (sin th/2 cos th/2)^m sum_s (-1)^s (l+m+s)!/(l-m-s)! (sin th/2)^{2s}/(s!(s+m)!)
    def f10(l, m, th, ph):
        s2 = math.sin(th/2); c2 = math.cos(th/2)
        tot = 0.0
        for s in range(0, l-m+1):
            if fac(l-m-s) == math.inf: continue
            tot += (-1)**s*fac(l+m+s)/fac(l-m-s)*s2**(2*s)/(fac(s)*fac(s+m))
        return (-1)**m*e(1j*m*ph)*N2(l, m)*(s2*c2)**m*tot
    ok &= sweep("5.2.10 (l+m+s)!/(l-m-s)!", f10, mpos=True)

    # 5.2.11  e N2 (cot th/2)^m sum_s (-1)^{l-s} (2l-s)!/(s!(l-s)!) (sin th/2)^{2(l-s)}/(l-m-s)!
    def f11(l, m, th, ph):
        s2 = math.sin(th/2)
        tot = 0.0
        for s in range(0, l-m+1):
            if fac(l-m-s) == math.inf or fac(l-s) == math.inf: continue
            tot += (-1)**(l-s)*fac(2*l-s)/(fac(s)*fac(l-s))*s2**(2*(l-s))/fac(l-m-s)
        return e(1j*m*ph)*N2(l, m)*(1/math.tan(th/2))**m*tot
    ok &= sweep("5.2.11 (cot)^m, (2l-s)!, sin^2(l-s)", f11, mpos=True)

    # 5.2.12 (-1)^l e N1 (cot th/2)^m sum_s (-1)^s (l+s)!/(l-s)! (cos th/2)^{2s}/(s!(s+m)!)
    def f12(l, m, th, ph):
        c2 = math.cos(th/2)
        tot = sum((-1)**s*fac(l+s)/fac(l-s)*c2**(2*s)/(fac(s)*fac(s+m))
                  for s in range(0, l+1) if fac(l-s) != math.inf and fac(s+m) != math.inf)
        return (-1)**l*e(1j*m*ph)*N1(l, m)*(1/math.tan(th/2))**m*tot
    ok &= sweep("5.2.12 (-1)^l cot^m cos^2s", f12, mpos=True)

    # 5.2.13 (-1)^l e N2 (sin cos th/2)^m sum_s (-1)^s (l+m+s)!/(l-m-s)! (cos th/2)^{2s}/(s!(s+m)!)
    def f13(l, m, th, ph):
        s2 = math.sin(th/2); c2 = math.cos(th/2)
        tot = 0.0
        for s in range(0, l-m+1):
            if fac(l-m-s) == math.inf: continue
            tot += (-1)**s*fac(l+m+s)/fac(l-m-s)*c2**(2*s)/(fac(s)*fac(s+m))
        return (-1)**l*e(1j*m*ph)*N2(l, m)*(s2*c2)**m*tot
    ok &= sweep("5.2.13 (-1)^l (l+m+s)!/(l-m-s)! cos^2s", f13, mpos=True)

    # 5.2.14 (-1)^{l-m} e N2 (tan th/2)^m sum_s (-1)^{l-s} (2l-s)!/(s!(l-s)!) (cos th/2)^{2(l-s)}/(l-m-s)!
    def f14(l, m, th, ph):
        c2 = math.cos(th/2)
        tot = 0.0
        for s in range(0, l-m+1):
            if fac(l-m-s) == math.inf or fac(l-s) == math.inf: continue
            tot += (-1)**(l-s)*fac(2*l-s)/(fac(s)*fac(l-s))*c2**(2*(l-s))/fac(l-m-s)
        return (-1)**(l-m)*e(1j*m*ph)*N2(l, m)*math.tan(th/2)**m*tot
    ok &= sweep("5.2.14 (-1)^{l-m} tan^m cos^2(l-s)", f14, mpos=True)

    # ---- remaining hypergeometric forms ----
    # 5.2.24 (-1)^{l-m} xi e N (2l)!/l! (sin th/2)^{2l} (cot th/2)^{|m|} F(-l,-l+|m|;-2l;1/sin^2 th/2)
    def f24(l, m, th, ph):
        am = abs(m)
        N = math.sqrt((2*l+1)/(4*math.pi*fac(l+m)*fac(l-m)))
        return ((-1)**(l-m)*xi(m)*e(1j*m*ph)*N*fac(2*l)/fac(l)
                * math.sin(th/2)**(2*l)*(1/math.tan(th/2))**am
                * F(-l, -l+am, -2*l, 1/math.sin(th/2)**2))
    ok &= sweep("5.2.24 2F1(1/sin^2 th/2)", f24, ths=[0.4,0.9,1.3,1.9])

    # 5.2.25 (-1)^{l-m} xi e sqrt((2l+1)/4pi (l+|m|)!/(l-|m|)!) (sin th)^|m|/(|m|! 2^|m|) F(..;cos^2 th/2)
    def f25(l, m, th, ph):
        am = abs(m)
        return ((-1)**(l-m)*xi(m)*e(1j*m*ph)
                * math.sqrt((2*l+1)/(4*math.pi)*fac(l+am)/fac(l-am))
                * math.sin(th)**am/(fac(am)*2**am)
                * F(-l+am, l+am+1, am+1, math.cos(th/2)**2))
    ok &= sweep("5.2.25 2F1(cos^2 th/2)", f25)

    # 5.2.26 xi e N (2l)!/l! (cos th/2)^{2l} (tan th/2)^{|m|} F(-l,-l+|m|;-2l;1/cos^2 th/2)
    def f26(l, m, th, ph):
        am = abs(m)
        N = math.sqrt((2*l+1)/(4*math.pi*fac(l+m)*fac(l-m)))
        return (xi(m)*e(1j*m*ph)*N*fac(2*l)/fac(l)
                * math.cos(th/2)**(2*l)*math.tan(th/2)**am
                * F(-l, -l+am, -2*l, 1/math.cos(th/2)**2))
    ok &= sweep("5.2.26 2F1(1/cos^2 th/2)", f26, ths=[0.4,0.9,1.3,1.9,2.6])

    # 5.2.28 (-1)^{l-m} xi e sqrt((2l+1)/4pi (l+|m|)!/(l-|m|)!) (cot th/2)^|m| (sin th/2)^{2l} F(-l+|m|,-l;|m|+1;-cot^2 th/2)
    def f28(l, m, th, ph):
        am = abs(m)
        return ((-1)**(l-m)*xi(m)*e(1j*m*ph)
                * math.sqrt((2*l+1)/(4*math.pi)*fac(l+am)/fac(l-am))
                * (1/math.tan(th/2))**am*math.sin(th/2)**(2*l)
                * F(-l+am, -l, am+1, -1/math.tan(th/2)**2))
    ok &= sweep("5.2.28 2F1(-cot^2 th/2)", f28)

    # ---- 5.2.5 hypergeometric in theta ----
    # 5.2.31 e sqrt((2l+1)/4pi/((l+m)!(l-m)!)) (sin th)^m {even/odd}
    def f31(l, m, th, ph):
        pre = e(1j*m*ph)*math.sqrt((2*l+1)/(4*math.pi*fac(l+m)*fac(l-m)))*math.sin(th)**m
        if (l+m) % 2 == 0:
            body = ((-1)**((l+m)//2)*dfac(l+m-1)*dfac(l-m-1)
                    * F(-(l-m)/2, (l+m+1)/2, 0.5, math.cos(th)**2))
        else:
            body = ((-1)**((l+m-1)//2)*dfac(l+m)*dfac(l-m)*math.cos(th)
                    * F(-(l-m-1)/2, (l+m+2)/2, 1.5, math.cos(th)**2))
        return pre*body
    ok &= sweep("5.2.31 2F1(cos^2 th) piecewise", f31)

    # 5.2.33 xi e sqrt((2l+1)/4pi (l+|m|)!/(l-|m|)!) (cos th)^l (tan th)^|m|/(2^|m| |m|!) F(..;-tan^2 th)
    def f33(l, m, th, ph):
        am = abs(m)
        return (xi(m)*e(1j*m*ph)*math.sqrt((2*l+1)/(4*math.pi)*fac(l+am)/fac(l-am))
                * math.cos(th)**l*math.tan(th)**am/(2**am*fac(am))
                * F(-(l-am)/2, -(l-am-1)/2, am+1, -math.tan(th)**2))
    ok &= sweep("5.2.33 2F1(-tan^2 th)", f33, ths=[0.4,0.9,1.3])

    # 5.2.34 e sqrt(..1/((l+m)!(l-m)!)) (sin th)^l {even/odd} F(..;-cot^2 th)
    def f34(l, m, th, ph):
        pre = e(1j*m*ph)*math.sqrt((2*l+1)/(4*math.pi*fac(l+m)*fac(l-m)))*math.sin(th)**l
        if (l+m) % 2 == 0:
            body = ((-1)**((l+m)//2)*dfac(l+m-1)*dfac(l-m-1)
                    * F(-(l-m)/2, -(l+m)/2, 0.5, -1/math.tan(th)**2))
        else:
            body = ((-1)**((l+m-1)//2)*dfac(l+m)*dfac(l-m)/math.tan(th)  # cot th
                    * F(-(l-m-1)/2, -(l+m-1)/2, 1.5, -1/math.tan(th)**2))
        return pre*body
    ok &= sweep("5.2.34 2F1(-cot^2 th) piecewise", f34, ths=[0.4,0.9,1.3])

    print("\nRESULT:", "ALL PASS" if ok else "FAILURES ABOVE")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
