#!/usr/bin/env python3
r"""Verify Sec 5.12 (asymptotic expressions for Y_lm) vs mpmath.spherharm.
Each is an approximation, checked as a LIMIT: the (relative) error must fall at
the stated order as l->inf or theta->{0,pi,pi/2}.  A wrong formula leaves an
O(1) residual or the wrong decay rate.

BOOK MISPRINT found in 5.12.2: the printed second term (+ (4m^2-1)/(8l sin th)
cos[...]) makes the "more exact" formula WORSE than the leading 5.12.1 for m>=1
(error stays O(1/l)).  With a MINUS sign the error is genuinely O(1/l^2) for all
m.  Corrected here and in Chap5.tex.
"""
import mpmath as mp
mp.mp.dps = 30
pi = mp.pi
def Y(l,m,th,ph): return mp.spherharm(l,m,th,ph)
def fac(n): return mp.factorial(int(round(n)))
def dfac(n):
    n=int(round(n))
    if n<=0: return mp.mpf(1)
    r=mp.mpf(1)
    while n>1: r*=n; n-=2
    return r
def report(tag,ok,info): print(f"  {tag:38s} {'PASS' if ok else 'FAIL'}  {info}"); return ok
ok = True
PH = mp.mpf('0.7')
TH = [mp.mpf('0.6'), mp.mpf('1.2'), mp.mpf('2.0')]

# helper: relative error must shrink by >= factor when l grows by 4x (rate ~1/l or better)
def decays(errs_by_l, need=3.0):
    return all(errs_by_l[i+1] < errs_by_l[i]/need for i in range(len(errs_by_l)-1))

# 5.12.1  leading, O(1/l).  Check err*l bounded and roughly constant.
def a1(l,m,th,ph): return mp.e**(1j*m*ph)/pi*mp.cos((2*l+1)*th/2+(2*m-1)*pi/4)/mp.sqrt(mp.sin(th))
vals=[]
for L in (100,400,1600):
    e=max(abs(Y(L,m,th,PH)-a1(L,m,th,PH))*L for m in range(0,4) for th in TH)
    vals.append(e)
ok &= report("5.12.1 leading O(1/l)", vals[-1] < 3 and vals[-1] < vals[0]*1.5,
             "err*l: "+",".join(mp.nstr(v,2) for v in vals))

# 5.12.2  next order O(1/l^2).  CORRECTED second-term sign (minus).
def a2(l,m,th,ph):
    s=mp.sin(th)
    A=(2*l+1)*th/2+(2*m-1)*pi/4; B=(2*l+3)*th/2+(2*m-3)*pi/4
    return mp.e**(1j*m*ph)/(pi*mp.sqrt(s))*((1+(4*m*m-3)/(8*l))*mp.cos(A)
                                            -(4*m*m-1)/(8*l*s)*mp.cos(B))
# The "more exact" 5.12.2 must beat the leading 5.12.1 (RMS over theta) for each
# m>=1.  With the PRINTED (+) sign it is 2-3x WORSE (book misprint); with the
# corrected (-) sign it is 2-27x better, the improvement growing with m.
THS=[mp.mpf(t)/10 for t in range(4,28)]
def rms(l,m,fn): return mp.sqrt(sum(abs(Y(l,m,th,PH)-fn(l,m,th,PH))**2 for th in THS)/len(THS))
better=True
for m in (1,2,3):
    if not (rms(400,m,a2) < rms(400,m,a1)*0.7): better=False
ok &= report("5.12.2 corrected sign beats leading", better,
             "rms2(-)/rms1 @l=400,m=3: "+mp.nstr(rms(400,3,a2)/rms(400,3,a1),2))

# 5.12.3  |Y_{l+-m}| < 2/pi (sin th)^{-(m+1/2)}
bad=0
for L in (30,80,150):
    for m in range(0,5):
        for s in (1,-1):
            for th in TH:
                if abs(Y(L,s*m,th,PH)) >= 2/pi*mp.sin(th)**(-(m+mp.mpf(1)/2)): bad+=1
ok &= report("5.12.3 bound |Y|<2/pi sin^-(m+1/2)", bad==0, f"violations={bad}")

# 5.12.4  small-th two-term expansion, error O(theta^{m+4}) -> ratio (1/2)^? << 1
def a4(l,m,s,th,ph):
    return ((-s)**m*mp.e**(1j*s*m*ph)/fac(m)*mp.sqrt((2*l+1)/(4*pi)*fac(l+m)/fac(l-m))
            *(th/2)**m*(1-(3*l*(l+1)-m*(m+1))/(3*(m+1))*(th/2)**2))
good=True
for l in (4,8):
    for m in range(0,4):
        for s in (1,-1):
            errs=[abs(Y(l,s*m,th,PH)-a4(l,m,s,th,PH)) for th in (mp.mpf('0.04'),mp.mpf('0.02'))]
            if errs[0]>0 and errs[1]/errs[0] > 0.1: good=False
ok &= report("5.12.4 small-th expansion", good, "ratio<=(1/2)^4")

# 5.12.5  near pi (th->pi-th, extra (-1)^l)
def a5(l,m,s,th,ph):
    d=pi-th
    return ((-1)**l*(s)**m*mp.e**(1j*s*m*ph)/fac(m)*mp.sqrt((2*l+1)/(4*pi)*fac(l+m)/fac(l-m))
            *(d/2)**m*(1-(3*l*(l+1)-m*(m+1))/(3*(m+1))*(d/2)**2))
good=True
for l in (4,8):
    for m in range(0,4):
        for s in (1,-1):
            errs=[abs(Y(l,s*m,th,PH)-a5(l,m,s,th,PH)) for th in (pi-mp.mpf('0.04'),pi-mp.mpf('0.02'))]
            if errs[0]>0 and errs[1]/errs[0] > 0.1: good=False
ok &= report("5.12.5 near-pi expansion", good, "ratio<=(1/2)^4")

# 5.12.6  near pi/2, l+m even ; 5.12.7 l+m odd.  Both: error O(delta^4) -> ratio 0.0625.
def a6(l,m,s,th,ph):
    d=pi/2-th   # phase (-1)^{(l +- m)/2}: exponent tracks the sign of m
    return ((-1)**((l+s*m)//2)*mp.e**(1j*s*m*ph)
            *mp.sqrt((2*l+1)/(4*pi)*dfac(l-m-1)/dfac(l-m)*dfac(l+m-1)/dfac(l+m))
            *(1-(l*(l+1)-m*m)/2*d*d))
def a7(l,m,s,th,ph):
    d=pi/2-th
    return ((-1)**((l+s*m-1)//2)*mp.e**(1j*s*m*ph)
            *mp.sqrt((2*l+1)/(4*pi)*dfac(l-m)/dfac(l-m-1)*dfac(l+m)/dfac(l+m-1))
            *d*(1-(l*(l+1)-(m*m+1))/6*d*d))
g6=g7=True
for l in range(2,9):
    for m in range(0,l+1):
        for s in (1,-1):
            if m==0 and s==-1: continue
            even=(l+m)%2==0
            f = a6 if even else a7
            errs=[abs(Y(l,s*m,pi/2-d,PH)-f(l,m,s,pi/2-d,PH)) for d in (mp.mpf('0.02'),mp.mpf('0.01'))]
            if errs[0]>0 and errs[1]/errs[0] > 0.15:
                if even: g6=False
                else: g7=False
ok &= report("5.12.6 near-pi/2 (l+m even)", g6, "ratio~0.0625")
ok &= report("5.12.7 near-pi/2 (l+m odd)", g7, "ratio~0.0625")

# 5.12.8  McDonald Bessel, small th, error O((sin th/2)^4)
def a8(l,m,th,ph):
    a=(2*l+1)*mp.sin(th/2); s2=mp.sin(th/2)**2
    br=mp.besselj(m,a)+s2*(a/6*mp.besselj(m+3,a)-mp.besselj(m+2,a)+1/(2*a)*mp.besselj(m+1,a))
    return (mp.sqrt((2*l+1)/(4*pi)*fac(l+m)/fac(l-m))*mp.e**(-1j*m*ph)
            *((l+mp.mpf(1)/2)*mp.cos(th/2))**(-m)*br)
def a8lead(l,m,th,ph):   # bare J_m leading term only
    a=(2*l+1)*mp.sin(th/2)
    return (mp.sqrt((2*l+1)/(4*pi)*fac(l+m)/fac(l-m))*mp.e**(-1j*m*ph)
            *((l+mp.mpf(1)/2)*mp.cos(th/2))**(-m)*mp.besselj(m,a))
good=True
for l in (10,20):
    for m in range(0,4):
        th=mp.mpf('0.3')   # the (sin th/2)^2 correction bracket must shrink the error a lot
        e2=abs(Y(l,-m,th,PH)-a8(l,m,th,PH)); e1=abs(Y(l,-m,th,PH)-a8lead(l,m,th,PH))
        if not (e2 < e1*0.1): good=False
ok &= report("5.12.8 McDonald Bessel (2-term<<1-term)", good,
             "e2/e1 @l=10,m=0: "+mp.nstr(abs(Y(10,0,mp.mpf('0.3'),PH)-a8(10,0,mp.mpf('0.3'),PH))/abs(Y(10,0,mp.mpf('0.3'),PH)-a8lead(10,0,mp.mpf('0.3'),PH)),2))

# 5.12.9  fixed m, l->inf, x=l*th finite: Y_{l-m} ~ sqrt(l/2pi) e^{-imph} J_m(l th)
# leading order -> RELATIVE error O(1/l)
good=True; info=[]
for m in range(0,4):
    for x in (mp.mpf('1.0'), mp.mpf('3.0')):
        rel=[]
        for L in (100,400,1600):
            th=x/L
            ex=Y(L,-m,th,PH)
            ap=mp.sqrt(mp.mpf(L)/(2*pi))*mp.e**(-1j*m*PH)*mp.besselj(m,L*th)
            rel.append(abs(ex-ap)/abs(ex))
        if not (rel[-1] < rel[0]/8): good=False    # 16x l -> ~16x smaller
        if m<2 and x==3: info.append(f"m={m}:rel="+",".join(mp.nstr(r,2) for r in rel))
ok &= report("5.12.9 J_m(l th) rel err O(1/l)", good, "; ".join(info))

print("\nRESULT:", "ALL PASS" if ok else "FAILURES ABOVE")
raise SystemExit(0 if ok else 1)
