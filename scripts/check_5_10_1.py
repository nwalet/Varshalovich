#!/usr/bin/env python3
r"""Verify Sec 5.10.1 (sums over m at fixed l) vs mpmath.spherharm."""
import mpmath as mp
mp.mp.dps = 30
TOL = mp.mpf('1e-12')
TH = [mp.mpf('0.4'), mp.mpf('0.9'), mp.mpf('1.3'), mp.mpf('1.9'), mp.mpf('2.6')]
PH = [mp.mpf('0.3'), mp.mpf('1.8'), mp.mpf('4.0')]
pi = mp.pi
def Y(l, m, th, ph):
    return mp.spherharm(l, m, th, ph) if abs(m) <= l else mp.mpc(0)
def fac(n): return mp.factorial(int(round(n)))
def report(tag, w): ok = w < TOL; print(f"  {tag:30s} {'PASS' if ok else 'FAIL'}  worst={mp.nstr(w,3)}"); return ok

ok = True
# 5.10.1  sum |Y|^2 = (2l+1)/4pi
w = mp.mpf(0)
for l in range(0,7):
    for th in TH:
        for ph in PH:
            s = sum(abs(Y(l,m,th,ph))**2 for m in range(-l,l+1))
            w = max(w, abs(s - (2*l+1)/(4*pi)))
ok &= report("5.10.1 sum |Y|^2", w)
# 5.10.2  sum m|Y|^2 = 0
w = mp.mpf(0)
for l in range(0,7):
    for th in TH:
        for ph in PH:
            s = sum(m*abs(Y(l,m,th,ph))**2 for m in range(-l,l+1))
            w = max(w, abs(s))
ok &= report("5.10.2 sum m|Y|^2", w)
# 5.10.3  sum m^2|Y|^2 = l(l+1)(2l+1)/8pi sin^2 th   [printed lower limit -1 -> -l]
w = mp.mpf(0)
for l in range(0,7):
    for th in TH:
        for ph in PH:
            s = sum(m*m*abs(Y(l,m,th,ph))**2 for m in range(-l,l+1))
            rhs = l*(l+1)*(2*l+1)/(8*pi)*mp.sin(th)**2
            w = max(w, abs(s - rhs))
ok &= report("5.10.3 sum m^2|Y|^2", w)
# 5.10.4  sum sqrt((l^2-m^2)((l+1)^2-m^2)/((2l-1)(2l+3))) Y_{l-1,m} Y*_{l+1,m}
#          = l(l+1)/8pi (3cos^2-1)
w = mp.mpf(0)
for l in range(1,7):
    for th in TH:
        for ph in PH:
            s = mp.mpc(0)
            for m in range(-l,l+1):
                num = (l*l-m*m)*((l+1)**2-m*m)
                s += mp.sqrt(mp.mpf(num)/((2*l-1)*(2*l+3)))*Y(l-1,m,th,ph)*mp.conj(Y(l+1,m,th,ph))
            rhs = l*(l+1)/(8*pi)*(3*mp.cos(th)**2-1)
            w = max(w, abs(s - rhs))
ok &= report("5.10.4 sum Y_{l-1}Y*_{l+1}", w)
# 5.10.5  sum (-+i)^m/sqrt((l-m)!(l+m)!) Y_lm = sqrt((2l+1)/4pi)(cos +- i sin cos ph)^l / l!
w = mp.mpf(0)
for sgn in (+1,-1):  # -+ : upper sign uses -i, lower uses +i
    for l in range(0,7):
        for th in TH:
            for ph in PH:
                s = sum((-sgn*1j)**m/mp.sqrt(fac(l-m)*fac(l+m))*Y(l,m,th,ph) for m in range(-l,l+1))
                rhs = mp.sqrt((2*l+1)/(4*pi))*(mp.cos(th)+sgn*1j*mp.sin(th)*mp.cos(ph))**l/fac(l)
                w = max(w, abs(s - rhs))
ok &= report("5.10.5 sum (-+i)^m Y", w)

print("\nRESULT:", "ALL PASS" if ok else "FAILURES ABOVE")
raise SystemExit(0 if ok else 1)
