#!/usr/bin/env python3
r"""Verify Sec 5.10.2 (sums over l at fixed m>=0) vs mpmath.spherharm.
Infinite sums truncated at l<=LMAX (|t|<1 => geometric convergence);
convergence cross-checked implicitly by the tight tolerance.
Notation resolved from scan (PDF p.163-164):
  j(t) -> j_l(t) spherical Bessel (1st kind); z_l -> a 2nd solution (use y_l).
  o_1(;m+1;.) -> 0F1(;m+1;.).
"""
import mpmath as mp
mp.mp.dps = 30
TOL = mp.mpf('1e-10')
pi = mp.pi
def Y(l, m, th, ph):
    return mp.spherharm(l, m, th, ph) if abs(m) <= l else mp.mpc(0)
def fac(n): return mp.factorial(int(round(n)))
def dfac(n):
    n = int(round(n))
    if n <= 0: return mp.mpf(1)
    r = mp.mpf(1)
    while n > 1: r *= n; n -= 2
    return r
def jl(l, x):  # spherical Bessel 1st kind
    x = mp.mpf(x)
    if x == 0: return mp.mpf(1) if l == 0 else mp.mpf(0)
    return mp.sqrt(pi/(2*x))*mp.besselj(l+mp.mpf(1)/2, x)
def yl(l, x):  # spherical Bessel 2nd kind (Neumann)
    x = mp.mpf(x)
    return mp.sqrt(pi/(2*x))*mp.bessely(l+mp.mpf(1)/2, x)
def F01(b, z): return mp.hyp0f1(b, z)
def F21(a, b, c, z): return mp.hyp2f1(a, b, c, z)
def report(tag, w): ok = w < TOL; print(f"  {tag:34s} {'PASS' if ok else 'FAIL'}  worst={mp.nstr(w,3)}"); return ok

LMAX = 90
TH = [mp.mpf('0.4'), mp.mpf('1.3'), mp.mpf('2.2')]
PH = [mp.mpf('0.3'), mp.mpf('1.8')]
ok = True

# 5.10.7  sum_{l=m}^inf sqrt(4pi/((2l+1)(l-m)!(l+m)!)) (1/l!) t^{l-m} Y_lm
#  = (-sin e^{if})^m/(2^m (m!)^2) 0F1(;m+1;t cos^2(th/2)) 0F1(;m+1;-t sin^2(th/2))
w = mp.mpf(0)
for m in range(0,4):
    for t in (mp.mpf('0.3'), mp.mpf('0.6')):
        for th in TH:
            for ph in PH:
                s = mp.mpc(0)
                for l in range(m, LMAX+1):
                    s += (mp.sqrt(4*pi/((2*l+1)*fac(l-m)*fac(l+m)))/fac(l)
                          * t**(l-m)*Y(l,m,th,ph))
                rhs = ((-mp.sin(th)*mp.e**(1j*ph))**m/(2**m*fac(m)**2)
                       * F01(m+1, t*mp.cos(th/2)**2)*F01(m+1, -t*mp.sin(th/2)**2))
                w = max(w, abs(s-rhs))
ok &= report("5.10.7 product of 0F1", w)

# 5.10.8  sum (n+l-m)! t^{l-m}/sqrt((l+m)!(l-m)!) sqrt(4pi/(2l+1)) Y_lm
#  = n!/(2^m m!) (-sin e)^m/(1-t cos)^{n+1} 2F1((n+1)/2, n/2+1; m+1; -(t sin/(1-t cos))^2)
w = mp.mpf(0)
for n in range(0,3):
    for m in range(0,4):
        for t in (mp.mpf('0.3'), mp.mpf('0.55')):
            for th in TH:
                for ph in PH:
                    s = mp.mpc(0)
                    for l in range(m, LMAX+1):
                        s += (fac(n+l-m)*t**(l-m)/mp.sqrt(fac(l+m)*fac(l-m))
                              * mp.sqrt(4*pi/(2*l+1))*Y(l,m,th,ph))
                    den = 1-t*mp.cos(th)
                    rhs = (fac(n)/(2**m*fac(m))*(-mp.sin(th)*mp.e**(1j*ph))**m/den**(n+1)
                           * F21((n+1)/mp.mpf(2), mp.mpf(n)/2+1, m+1,
                                 -(t*mp.sin(th)/den)**2))
                    w = max(w, abs(s-rhs))
ok &= report("5.10.8 2F1 generating fn", w)

# 5.10.9  sum sqrt(4pi/((2l+1)(l-m)!(l+m)!)) (l+m-n)!(l-m+n-1)!/l! t^{l-m} Y_lm
#  = (2m-n)!(n-1)!/(2^m (m!)^2) (-sin e)^m F(n,2m+1-n;m+1;(1-t-R)/2) F(n,2m+1-n;m+1;(1+t-R)/2)
#  requires n>=1 and 2m>=n and l+m>=n (else (l+m-n)! ill-defined for small l)
w = mp.mpf(0)
for n in (1,2):
    for m in range(0,4):
        if 2*m < n: continue
        for t in (mp.mpf('0.3'), mp.mpf('0.5')):
            for th in TH:
                for ph in PH:
                    s = mp.mpc(0)
                    for l in range(m, LMAX+1):
                        if l+m-n < 0: continue
                        s += (mp.sqrt(4*pi/((2*l+1)*fac(l-m)*fac(l+m)))
                              * fac(l+m-n)*fac(l-m+n-1)/fac(l)*t**(l-m)*Y(l,m,th,ph))
                    R = mp.sqrt(1-2*t*mp.cos(th)+t*t)
                    rhs = (fac(2*m-n)*fac(n-1)/(2**m*fac(m)**2)*(-mp.sin(th)*mp.e**(1j*ph))**m
                           * F21(n, 2*m+1-n, m+1, (1-t-R)/2)
                           * F21(n, 2*m+1-n, m+1, (1+t-R)/2))
                    w = max(w, abs(s-rhs))
ok &= report("5.10.9 product of 2F1", w)

# 5.10.10  sum sqrt(4pi(l-m)!/((2l+1)(l+m)!)) L_{l-m}^{2m}(y) t^{l-m} Y_lm
#  = 1/(2^m m!) (-sin e)^m/R^{2m+1} exp[-y t(cos-t)/R^2] 0F1(;m+1; -y^2 t^2 sin^2/(4R^4))
w = mp.mpf(0)
for m in range(0,4):
    for y in (mp.mpf('0.7'), mp.mpf('1.5')):
        for t in (mp.mpf('0.3'), mp.mpf('0.55')):
            for th in TH:
                for ph in PH:
                    s = mp.mpc(0)
                    for l in range(m, LMAX+1):
                        s += (mp.sqrt(4*pi*fac(l-m)/((2*l+1)*fac(l+m)))
                              * mp.laguerre(l-m, 2*m, y)*t**(l-m)*Y(l,m,th,ph))
                    R = mp.sqrt(1-2*t*mp.cos(th)+t*t)
                    rhs = (1/(2**m*fac(m))*(-mp.sin(th)*mp.e**(1j*ph))**m/R**(2*m+1)
                           * mp.e**(-y*t*(mp.cos(th)-t)/R**2)
                           * F01(m+1, -y*y*t*t*mp.sin(th)**2/(4*R**4)))
                    w = max(w, abs(s-rhs))
ok &= report("5.10.10 Laguerre gen fn", w)

# 5.10.11  sum i^{l-m} sqrt(4pi(2l+1)(l+m)!/(l-m)!) j_l(t) Y_lm = (-t sin e)^m e^{it cos}
w = mp.mpf(0)
for m in range(0,4):
    for t in (mp.mpf('1.0'), mp.mpf('2.5')):
        for th in TH:
            for ph in PH:
                s = mp.mpc(0)
                for l in range(m, LMAX+1):
                    s += ((1j)**(l-m)*mp.sqrt(4*pi*(2*l+1)*fac(l+m)/fac(l-m))
                          * jl(l,t)*Y(l,m,th,ph))
                rhs = (-t*mp.sin(th)*mp.e**(1j*ph))**m*mp.e**(1j*t*mp.cos(th))
                w = max(w, abs(s-rhs))
ok &= report("5.10.11 Rayleigh (full)", w)

# 5.10.12  sum_{l=m,m+2,..} i^{l-m} sqrt(...) j_l(t) Y_lm = (-t sin e)^m cos(t cos)
w = mp.mpf(0)
for m in range(0,4):
    for t in (mp.mpf('1.0'), mp.mpf('2.5')):
        for th in TH:
            for ph in PH:
                s = mp.mpc(0)
                for l in range(m, LMAX+1, 2):
                    s += ((1j)**(l-m)*mp.sqrt(4*pi*(2*l+1)*fac(l+m)/fac(l-m))
                          * jl(l,t)*Y(l,m,th,ph))
                rhs = (-t*mp.sin(th)*mp.e**(1j*ph))**m*mp.cos(t*mp.cos(th))
                w = max(w, abs(s-rhs))
ok &= report("5.10.12 Rayleigh (even)", w)

# 5.10.13  sum_{l=m+1,m+3,..} i^{l-m-1} sqrt(...) j_l(t) Y_lm = (-t sin e)^m sin(t cos)
w = mp.mpf(0)
for m in range(0,4):
    for t in (mp.mpf('1.0'), mp.mpf('2.5')):
        for th in TH:
            for ph in PH:
                s = mp.mpc(0)
                for l in range(m+1, LMAX+1, 2):
                    s += ((1j)**(l-m-1)*mp.sqrt(4*pi*(2*l+1)*fac(l+m)/fac(l-m))
                          * jl(l,t)*Y(l,m,th,ph))
                rhs = (-t*mp.sin(th)*mp.e**(1j*ph))**m*mp.sin(t*mp.cos(th))
                w = max(w, abs(s-rhs))
ok &= report("5.10.13 Rayleigh (odd)", w)

# 5.10.15  sum sqrt(4pi(2l+1)(l+m)!/(l-m)!) j_l(x) z_l(y) Y_lm
#  = (-xy sin e/D)^m z_m(D),  D=sqrt(x^2+y^2-2xy cos),  z=y_l, need x<y for convergence
w = mp.mpf(0)
for m in range(0,4):
    for (x,y) in ((mp.mpf('0.5'), mp.mpf('2.0')), (mp.mpf('0.8'), mp.mpf('3.0'))):
        for th in TH:
            for ph in PH:
                s = mp.mpc(0)
                for l in range(m, LMAX+1):
                    s += (mp.sqrt(4*pi*(2*l+1)*fac(l+m)/fac(l-m))
                          * jl(l,x)*yl(l,y)*Y(l,m,th,ph))
                D = mp.sqrt(x*x+y*y-2*x*y*mp.cos(th))
                rhs = (-x*y*mp.sin(th)*mp.e**(1j*ph)/D)**m*yl(m, D)
                w = max(w, abs(s-rhs))
ok &= report("5.10.15 Bessel addition (z=y_l)", w)

# 5.10.14  x=y=t boundary case of 5.10.15.  The direct sum diverges at x=y
# (terms ~ l^{m-1/2}); the identity holds as the limit x->y^-.  Since 5.10.15's
# LHS=RHS is verified above, confirm 5.10.14 by the exact identity
# RHS_{5.10.14}(t) == RHS_{5.10.15}(x=y=t) to machine precision.
w = mp.mpf(0)
for m in range(0,4):
    for t in (mp.mpf('0.6'), mp.mpf('1.2'), mp.mpf('2.0')):
        for th in TH:
            for ph in PH:
                D = mp.sqrt(2*t*t - 2*t*t*mp.cos(th))          # x=y=t in 5.10.15
                rhs15 = (-t*t*mp.sin(th)*mp.e**(1j*ph)/D)**m*yl(m, D)
                rhs14 = (-t*mp.cos(th/2)*mp.e**(1j*ph))**m*yl(m, 2*t*mp.sin(th/2))
                w = max(w, abs(rhs14 - rhs15))
ok &= report("5.10.14 = 5.10.15|_{x=y=t}", w)


# 5.10.16  sum i^{l-m} j_l(t) Y_lm(th1,ph1) Y*_lm(th2,ph2)
#  = 1/4pi J_m(t sin1 sin2) e^{it cos1 cos2} e^{im(ph1-ph2)}
w = mp.mpf(0)
for m in range(0,4):
    for t in (mp.mpf('1.0'), mp.mpf('2.5')):
        for (t1,p1,t2,p2) in ((mp.mpf('0.7'),mp.mpf('0.3'),mp.mpf('1.9'),mp.mpf('1.2')),
                              (mp.mpf('1.1'),mp.mpf('2.0'),mp.mpf('0.5'),mp.mpf('4.0'))):
            s = mp.mpc(0)
            for l in range(m, LMAX+1):
                s += (1j)**(l-m)*jl(l,t)*Y(l,m,t1,p1)*mp.conj(Y(l,m,t2,p2))
            rhs = (1/(4*pi)*mp.besselj(m, t*mp.sin(t1)*mp.sin(t2))
                   * mp.e**(1j*t*mp.cos(t1)*mp.cos(t2))*mp.e**(1j*m*(p1-p2)))
            w = max(w, abs(s-rhs))
ok &= report("5.10.16 bilinear (J_m)", w)

print("\nRESULT:", "ALL PASS" if ok else "FAILURES ABOVE")
raise SystemExit(0 if ok else 1)
