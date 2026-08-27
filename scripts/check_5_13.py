#!/usr/bin/env python3
r"""Verify Sec 5.13 live forms vs mpmath.spherharm.
  5.13.2 (eq:5:13:1..5): Y_{lm} for |m|=0..4 via Legendre P_l.
  5.13.3 (eq:5:13:6..11): Y_{lm} for |m|=l..l-5 via trig.
Reference Y = mpmath.spherharm (VMK = Condon-Shortley).
"""
import math, cmath
import mpmath as mp
mp.mp.dps = 30
TOL = 1e-10
TH = [0.4, 0.9, 1.3, 1.9, 2.6]
PH = [0.3, 1.8, 4.0]

def Y(l, m, th, ph):
    if l < 0 or abs(m) > l: return 0j
    return complex(mp.spherharm(l, m, th, ph))
def P(l, x): return float(mp.legendre(l, x))
def fac(n): return math.factorial(int(round(n)))
def dfac(n):
    n = int(round(n))
    if n <= 0: return 1.0
    r = 1.0
    while n > 1: r *= n; n -= 2
    return r
def report(tag, w, tol=TOL):
    ok = w < tol; print(f"  {tag:30s} {'PASS' if ok else 'FAIL'}  worst={w:.2e}"); return ok

def sweep(tag, form, lm_pairs, ths=TH, tol=TOL):
    w = 0.0
    for (l, m) in lm_pairs:
        for th in ths:
            for ph in PH:
                try: v = form(l, m, th, ph)
                except (ValueError, ZeroDivisionError, OverflowError): continue
                if not (math.isfinite(v.real) and math.isfinite(v.imag)): continue
                w = max(w, abs(v - Y(l, m, th, ph)))
    return report(tag, w, tol)

e = cmath.exp
ok = True
print("Sec 5.13.2  |m|=0..4 via Legendre P_l")

# 5.13.1  Y_{l0} = sqrt((2l+1)/4pi) P_l(cos th)
ok &= sweep("5.13.1 Y_{l0}",
    lambda l, m, th, ph: math.sqrt((2*l+1)/(4*math.pi))*P(l, math.cos(th)),
    [(l, 0) for l in range(0, 8)])

# 5.13.2  Y_{l,+-1} = -+ e^{+-i ph}/sin th sqrt(l(l+1)/(4pi(2l+1)))[P_{l-1}-P_{l+1}]
def f2(l, m, th, ph):
    s = 1 if m > 0 else -1  # sign of m; leading -+ means -sign(m)
    return (-s*e(1j*m*ph)/math.sin(th)
            * math.sqrt(l*(l+1)/(4*math.pi*(2*l+1)))
            * (P(l-1, math.cos(th)) - P(l+1, math.cos(th))))
ok &= sweep("5.13.2 Y_{l,+-1}", f2,
    [(l, m) for l in range(1, 8) for m in (-1, 1)])

# 5.13.3  Y_{l,+-2} = e^{+-i2ph}/sin^2 th sqrt((l-1)l(l+1)(l+2)/(4pi(2l+1)))
#   [P_{l-2}/(2l-1) - 2(2l+1)/((2l-1)(2l+3)) P_l + P_{l+2}/(2l+3)]
def f3(l, m, th, ph):
    c = math.cos(th)
    return (e(1j*m*ph)/math.sin(th)**2
            * math.sqrt((l-1)*l*(l+1)*(l+2)/(4*math.pi*(2*l+1)))
            * (P(l-2, c)/(2*l-1)
               - 2*(2*l+1)/((2*l-1)*(2*l+3))*P(l, c)
               + P(l+2, c)/(2*l+3)))
ok &= sweep("5.13.3 Y_{l,+-2}", f3,
    [(l, m) for l in range(2, 8) for m in (-2, 2)])

# 5.13.4  Y_{l,+-3} = -+ e^{+-i3ph}/sin^3 th sqrt((l+3)!/(l-3)! /(4pi(2l+1)))
#   [P_{l-3}/((2l-3)(2l-1)) - 3 P_{l-1}/((2l-3)(2l+3))
#    + 3 P_{l+1}/((2l-1)(2l+5)) - P_{l+3}/((2l+3)(2l+5))]
def f4(l, m, th, ph):
    s = 1 if m > 0 else -1
    c = math.cos(th)
    return (-s*e(1j*m*ph)/math.sin(th)**3
            * math.sqrt(fac(l+3)/fac(l-3)/(4*math.pi*(2*l+1)))
            * (P(l-3, c)/((2*l-3)*(2*l-1))
               - 3*P(l-1, c)/((2*l-3)*(2*l+3))
               + 3*P(l+1, c)/((2*l-1)*(2*l+5))
               - P(l+3, c)/((2*l+3)*(2*l+5))))
ok &= sweep("5.13.4 Y_{l,+-3}", f4,
    [(l, m) for l in range(3, 8) for m in (-3, 3)])

# 5.13.5  Y_{l,+-4} = e^{+-i4ph}/sin^4 th sqrt((l+4)!/(l-4)! /(4pi(2l+1)))
#   [P_{l-4}/((2l-5)(2l-3)(2l-1)) - 4 P_{l-2}/((2l-5)(2l-1)(2l+3))
#    + 6(2l+1) P_l/((2l-3)(2l-1)(2l+3)(2l+5))
#    - 4 P_{l+2}/((2l-1)(2l+3)(2l+7)) + P_{l+4}/((2l+3)(2l+5)(2l+7))]
def f5(l, m, th, ph):
    c = math.cos(th)
    return (e(1j*m*ph)/math.sin(th)**4
            * math.sqrt(fac(l+4)/fac(l-4)/(4*math.pi*(2*l+1)))
            * (P(l-4, c)/((2*l-5)*(2*l-3)*(2*l-1))
               - 4*P(l-2, c)/((2*l-5)*(2*l-1)*(2*l+3))
               + 6*(2*l+1)*P(l, c)/((2*l-3)*(2*l-1)*(2*l+3)*(2*l+5))
               - 4*P(l+2, c)/((2*l-1)*(2*l+3)*(2*l+7))
               + P(l+4, c)/((2*l+3)*(2*l+5)*(2*l+7))))
ok &= sweep("5.13.5 Y_{l,+-4}", f5,
    [(l, m) for l in range(4, 9) for m in (-4, 4)])

print("\nSec 5.13.3  |m|=l..l-5 via trig")
# 5.13.6  Y_{l,+-l} = (-+1)^l e^{+-il ph} sqrt((2l+1)!!/(4pi(2l)!!)) (sin th)^l
def g6(l, m, th, ph):
    s = 1 if m > 0 else -1
    return ((-s)**l*e(1j*m*ph)*math.sqrt(dfac(2*l+1)/(4*math.pi*dfac(2*l)))
            * math.sin(th)**l)
ok &= sweep("5.13.6 Y_{l,+-l}", g6,
    [(l, s*l) for l in range(1, 9) for s in (-1, 1)])

# 5.13.7  Y_{l,+-(l-1)} = (-+1)^{l-1} e^{+-i(l-1)ph} sqrt((2l+1)!!/(4pi(2l-2)!!)) cos th (sin th)^{l-1}
def g7(l, m, th, ph):
    s = 1 if m > 0 else -1
    return ((-s)**(l-1)*e(1j*m*ph)*math.sqrt(dfac(2*l+1)/(4*math.pi*dfac(2*l-2)))
            * math.cos(th)*math.sin(th)**(l-1))
ok &= sweep("5.13.7 Y_{l,+-(l-1)}", g7,
    [(l, s*(l-1)) for l in range(1, 9) for s in (-1, 1)])

# 5.13.8  Y_{l,+-(l-2)} = (-+1)^{l-2} e^{...} sqrt((2l+1)/8pi (2l-3)!!/(2l-2)!!) (sin th)^{l-2}[(2l-1)cos^2-1]
def g8(l, m, th, ph):
    s = 1 if m > 0 else -1
    return ((-s)**(l-2)*e(1j*m*ph)
            * math.sqrt((2*l+1)/(8*math.pi)*dfac(2*l-3)/dfac(2*l-2))
            * math.sin(th)**(l-2)*((2*l-1)*math.cos(th)**2 - 1))
ok &= sweep("5.13.8 Y_{l,+-(l-2)}", g8,
    [(l, s*(l-2)) for l in range(2, 9) for s in (-1, 1)])

# 5.13.9  Y_{l,+-(l-3)} = (-+1)^{l-3} e^{...} sqrt((2l+1)/24pi (2l-3)!!/(2l-4)!!)(sin)^{l-3}cos[(2l-1)cos^2-3]
def g9(l, m, th, ph):
    s = 1 if m > 0 else -1
    return ((-s)**(l-3)*e(1j*m*ph)
            * math.sqrt((2*l+1)/(24*math.pi)*dfac(2*l-3)/dfac(2*l-4))
            * math.sin(th)**(l-3)*math.cos(th)*((2*l-1)*math.cos(th)**2 - 3))
ok &= sweep("5.13.9 Y_{l,+-(l-3)}", g9,
    [(l, s*(l-3)) for l in range(3, 9) for s in (-1, 1)])

# 5.13.10  Y_{l,+-(l-4)} = (-+1)^{l-4} e sqrt((2l+1)/96pi (2l-5)!!/(2l-4)!!)(sin)^{l-4}
#   [(2l-1)(2l-3)cos^4 - 6(2l-3)cos^2 + 3]
def g10(l, m, th, ph):
    s = 1 if m > 0 else -1
    c = math.cos(th)
    return ((-s)**(l-4)*e(1j*m*ph)
            * math.sqrt((2*l+1)/(96*math.pi)*dfac(2*l-5)/dfac(2*l-4))
            * math.sin(th)**(l-4)
            * ((2*l-1)*(2*l-3)*c**4 - 6*(2*l-3)*c**2 + 3))
ok &= sweep("5.13.10 Y_{l,+-(l-4)}", g10,
    [(l, s*(l-4)) for l in range(4, 9) for s in (-1, 1)])

# 5.13.11  Y_{l,+-(l-5)} = (-+1)^{l-5} e sqrt((2l+1)/480pi (2l-5)!!/(2l-6)!!)(sin)^{l-5}cos
#   [(2l-1)(2l-3)cos^4 - 10(2l-3)cos^2 + 15]
def g11(l, m, th, ph):
    s = 1 if m > 0 else -1
    c = math.cos(th)
    return ((-s)**(l-5)*e(1j*m*ph)
            * math.sqrt((2*l+1)/(480*math.pi)*dfac(2*l-5)/dfac(2*l-6))
            * math.sin(th)**(l-5)*c
            * ((2*l-1)*(2*l-3)*c**4 - 10*(2*l-3)*c**2 + 15))
ok &= sweep("5.13.11 Y_{l,+-(l-5)}", g11,
    [(l, s*(l-5)) for l in range(5, 10) for s in (-1, 1)])

print("\nRESULT:", "ALL PASS" if ok else "FAILURES ABOVE")
raise SystemExit(0 if ok else 1)
