#!/usr/bin/env python3
r"""Verify the (commented-out) explicit l<=5 table in Sec 5.13.1 vs
mpmath.spherharm.  Both the primary trig form (LHS of the '=...=') and the
multiple-angle RHS are checked.  OCR garbles in labels are corrected here
(Y_{5+9}->Y_{5+3}, Y_{delta+1}->Y_{5+1}, sin phi->sin vartheta, etc).
"""
import math, cmath
import mpmath as mp
mp.mp.dps = 30
TOL = 1e-12
TH = [0.4, 0.9, 1.3, 1.9, 2.6]
PH = [0.3, 1.8]

def Y(l, m, th, ph):
    return complex(mp.spherharm(l, m, th, ph))
e = cmath.exp
def chk(tag, l, m, prim, mult):
    w = 0.0
    for th in TH:
        for ph in PH:
            ref = Y(l, m, th, ph)
            w = max(w, abs(prim(th, ph) - ref), abs(mult(th, ph) - ref))
    ok = w < TOL
    print(f"  Y_{{{l:>1},{m:>+2}}}  {'PASS' if ok else 'FAIL'}  worst={w:.2e}   ({tag})")
    return ok

ok = True
s3, s5, s7, s11 = math.sqrt, None, None, None
sqrt = math.sqrt
sin = math.sin; cos = math.cos; pi = math.pi

# l=0,1
ok &= chk("00", 0, 0, lambda t,p: 1/(2*sqrt(pi)), lambda t,p: 1/(2*sqrt(pi)))
ok &= chk("1+1",1, 1, lambda t,p:-0.5*sqrt(3/(2*pi))*sin(t)*e(1j*p),
                       lambda t,p:-0.5*sqrt(3/(2*pi))*sin(t)*e(1j*p))
ok &= chk("10", 1, 0, lambda t,p:0.5*sqrt(3/pi)*cos(t), lambda t,p:0.5*sqrt(3/pi)*cos(t))
ok &= chk("1-1",1,-1, lambda t,p:0.5*sqrt(3/(2*pi))*sin(t)*e(-1j*p),
                       lambda t,p:0.5*sqrt(3/(2*pi))*sin(t)*e(-1j*p))
# l=2
ok &= chk("2+2",2, 2, lambda t,p:0.25*sqrt(15/(2*pi))*sin(t)**2*e(2j*p),
                       lambda t,p:1/8*sqrt(15/(2*pi))*(1-cos(2*t))*e(2j*p))
ok &= chk("2+1",2, 1, lambda t,p:-0.5*sqrt(15/(2*pi))*cos(t)*sin(t)*e(1j*p),
                       lambda t,p:-0.25*sqrt(15/(2*pi))*sin(2*t)*e(1j*p))
ok &= chk("20", 2, 0, lambda t,p:0.25*sqrt(5/pi)*(3*cos(t)**2-1),
                       lambda t,p:1/8*sqrt(5/pi)*(1+3*cos(2*t)))
ok &= chk("2-1",2,-1, lambda t,p:0.5*sqrt(15/(2*pi))*cos(t)*sin(t)*e(-1j*p),
                       lambda t,p:0.25*sqrt(15/(2*pi))*sin(2*t)*e(-1j*p))
ok &= chk("2-2",2,-2, lambda t,p:0.25*sqrt(15/(2*pi))*sin(t)**2*e(-2j*p),
                       lambda t,p:1/8*sqrt(15/(2*pi))*(1-cos(2*t))*e(-2j*p))
# l=3
ok &= chk("3+3",3, 3, lambda t,p:-1/8*sqrt(35/pi)*sin(t)**3*e(3j*p),
                       lambda t,p:-1/32*sqrt(35/pi)*(3*sin(t)-sin(3*t))*e(3j*p))
ok &= chk("3+2",3, 2, lambda t,p:0.25*sqrt(105/(2*pi))*cos(t)*sin(t)**2*e(2j*p),
                       lambda t,p:1/16*sqrt(105/(2*pi))*(cos(t)-cos(3*t))*e(2j*p))
ok &= chk("3+1",3, 1, lambda t,p:-1/8*sqrt(21/pi)*(5*cos(t)**2-1)*sin(t)*e(1j*p),
                       lambda t,p:-1/32*sqrt(21/pi)*(sin(t)+5*sin(3*t))*e(1j*p))
ok &= chk("30", 3, 0, lambda t,p:0.25*sqrt(7/pi)*(5*cos(t)**2-3)*cos(t),
                       lambda t,p:1/16*sqrt(7/pi)*(3*cos(t)+5*cos(3*t)))
ok &= chk("3-1",3,-1, lambda t,p:1/8*sqrt(21/pi)*(5*cos(t)**2-1)*sin(t)*e(-1j*p),
                       lambda t,p:1/32*sqrt(21/pi)*(sin(t)+5*sin(3*t))*e(-1j*p))
ok &= chk("3-2",3,-2, lambda t,p:0.25*sqrt(105/(2*pi))*cos(t)*sin(t)**2*e(-2j*p),
                       lambda t,p:1/16*sqrt(105/(2*pi))*(cos(t)-cos(3*t))*e(-2j*p))
ok &= chk("3-3",3,-3, lambda t,p:1/8*sqrt(35/pi)*sin(t)**3*e(-3j*p),
                       lambda t,p:1/32*sqrt(35/pi)*(3*sin(t)-sin(3*t))*e(-3j*p))
# l=4
ok &= chk("4+4",4, 4, lambda t,p:3/16*sqrt(35/(2*pi))*sin(t)**4*e(4j*p),
                       lambda t,p:3/128*sqrt(35/(2*pi))*(3-4*cos(2*t)+cos(4*t))*e(4j*p))
ok &= chk("4+3",4, 3, lambda t,p:-3/8*sqrt(35/pi)*sin(t)**3*cos(t)*e(3j*p),
                       lambda t,p:-3/64*sqrt(35/pi)*(2*sin(2*t)-sin(4*t))*e(3j*p))
ok &= chk("4+2",4, 2, lambda t,p:3/8*sqrt(5/(2*pi))*sin(t)**2*(7*cos(t)**2-1)*e(2j*p),
                       lambda t,p:3/64*sqrt(5/(2*pi))*(3+4*cos(2*t)-7*cos(4*t))*e(2j*p))
ok &= chk("4+1",4, 1, lambda t,p:-3/8*sqrt(5/pi)*sin(t)*(7*cos(t)**3-3*cos(t))*e(1j*p),
                       lambda t,p:-3/64*sqrt(5/pi)*(2*sin(2*t)+7*sin(4*t))*e(1j*p))
ok &= chk("40", 4, 0, lambda t,p:3/(16*sqrt(pi))*(35*cos(t)**4-30*cos(t)**2+3),
                       lambda t,p:3/(128*sqrt(pi))*(9+20*cos(2*t)+35*cos(4*t)))
ok &= chk("4-1",4,-1, lambda t,p:3/8*sqrt(5/pi)*sin(t)*(7*cos(t)**3-3*cos(t))*e(-1j*p),
                       lambda t,p:3/64*sqrt(5/pi)*(2*sin(2*t)+7*sin(4*t))*e(-1j*p))
ok &= chk("4-2",4,-2, lambda t,p:3/8*sqrt(5/(2*pi))*sin(t)**2*(7*cos(t)**2-1)*e(-2j*p),
                       lambda t,p:3/64*sqrt(5/(2*pi))*(3+4*cos(2*t)-7*cos(4*t))*e(-2j*p))
# 4-3: printed RHS has 'sin 4 theta' and 'e^{-i S varphi}' garbles -> sin4t, e^{-3i}
ok &= chk("4-3",4,-3, lambda t,p:3/8*sqrt(35/pi)*sin(t)**3*cos(t)*e(-3j*p),
                       lambda t,p:3/64*sqrt(35/pi)*(2*sin(2*t)-sin(4*t))*e(-3j*p))
ok &= chk("4-4",4,-4, lambda t,p:3/16*sqrt(35/(2*pi))*sin(t)**4*e(-4j*p),
                       lambda t,p:3/128*sqrt(35/(2*pi))*(3-4*cos(2*t)+cos(4*t))*e(-4j*p))
# l=5
ok &= chk("5+5",5, 5, lambda t,p:-3/32*sqrt(77/pi)*sin(t)**5*e(5j*p),
                       lambda t,p:-3/512*sqrt(77/pi)*(10*sin(t)-5*sin(3*t)+sin(5*t))*e(5j*p))
ok &= chk("5+4",5, 4, lambda t,p:3/16*sqrt(385/(2*pi))*sin(t)**4*cos(t)*e(4j*p),
                       lambda t,p:3/256*sqrt(385/(2*pi))*(2*cos(t)-3*cos(3*t)+cos(5*t))*e(4j*p))
# 5+9 -> 5+3
ok &= chk("5+3",5, 3, lambda t,p:-1/32*sqrt(385/pi)*sin(t)**3*(9*cos(t)**2-1)*e(3j*p),
                       lambda t,p:-1/512*sqrt(385/pi)*(6*sin(t)+13*sin(3*t)-9*sin(5*t))*e(3j*p))
ok &= chk("5+2",5, 2, lambda t,p:1/8*sqrt(1155/(2*pi))*sin(t)**2*(3*cos(t)**3-cos(t))*e(2j*p),
                       lambda t,p:1/128*sqrt(1155/(2*pi))*(2*cos(t)+cos(3*t)-3*cos(5*t))*e(2j*p))
# delta+1 -> 5+1 ; sin phi -> sin vartheta
ok &= chk("5+1",5, 1, lambda t,p:-1/16*sqrt(165/(2*pi))*sin(t)*(21*cos(t)**4-14*cos(t)**2+1)*e(1j*p),
                       lambda t,p:-1/256*sqrt(165/(2*pi))*(2*sin(t)+7*sin(3*t)+21*sin(5*t))*e(1j*p))
ok &= chk("50", 5, 0, lambda t,p:1/16*sqrt(11/pi)*(63*cos(t)**5-70*cos(t)**3+15*cos(t)),
                       lambda t,p:1/256*sqrt(11/pi)*(30*cos(t)+35*cos(3*t)+63*cos(5*t)))
ok &= chk("5-1",5,-1, lambda t,p:1/16*sqrt(165/(2*pi))*sin(t)*(21*cos(t)**4-14*cos(t)**2+1)*e(-1j*p),
                       lambda t,p:1/256*sqrt(165/(2*pi))*(2*sin(t)+7*sin(3*t)+21*sin(5*t))*e(-1j*p))
ok &= chk("5-3",5,-3, lambda t,p:1/32*sqrt(385/pi)*sin(t)**3*(9*cos(t)**2-1)*e(-3j*p),
                       lambda t,p:1/512*sqrt(385/pi)*(6*sin(t)+13*sin(3*t)-9*sin(5*t))*e(-3j*p))
ok &= chk("5-4",5,-4, lambda t,p:3/16*sqrt(385/(2*pi))*sin(t)**4*cos(t)*e(-4j*p),
                       lambda t,p:3/256*sqrt(385/(2*pi))*(2*cos(t)-3*cos(3*t)+cos(5*t))*e(-4j*p))
ok &= chk("5-5",5,-5, lambda t,p:3/32*sqrt(77/pi)*sin(t)**5*e(-5j*p),
                       lambda t,p:3/512*sqrt(77/pi)*(10*sin(t)-5*sin(3*t)+sin(5*t))*e(-5j*p))

print("\nRESULT:", "ALL PASS" if ok else "FAILURES ABOVE")
raise SystemExit(0 if ok else 1)
