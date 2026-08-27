#!/usr/bin/env python3
r"""Verify Sec 5.5 (behaviour of Y_lm under coordinate transformations) vs
mpmath.spherharm.  Rotations use Wigner D = e^{-i m alpha} d^l_{m m'}(beta)
e^{-i m' gamma} (VMK convention, sympy Rotation.d).
"""
import math, cmath
import mpmath as mp
import numpy as np
from sympy import Rational, N
from sympy.physics.quantum.spin import Rotation
mp.mp.dps = 25
pi = mp.pi
def Y(l,m,th,ph):
    return complex(mp.spherharm(l,m,th,ph)) if abs(m)<=l else 0j
def smalld(l,m,mp_,beta):
    return complex(N(Rotation.d(l,m,mp_,beta).doit(), 25))
def Dfun(l,m,mp_,al,be,ga):
    return cmath.exp(-1j*m*al)*smalld(l,m,mp_,be)*cmath.exp(-1j*mp_*ga)
def report(tag,w,tol=1e-11):
    okk = w<tol; print(f"  {tag:36s} {'PASS' if okk else 'FAIL'}  worst={w:.2e}"); return okk
ok = True

TH=[0.4,0.9,1.3,1.9,2.6]; PH=[0.3,1.8,4.0]
def sweep(tag, rel, ls=range(0,5)):
    w=0.0
    for l in ls:
        for m in range(-l,l+1):
            for th in TH:
                for ph in PH:
                    w=max(w, abs(rel(l,m,th,ph)))
    return report(tag,w)

print("Sec 5.5 special angle relations")
# 5.5.2 inversion:  Y(pi-th, ph+pi) = (-1)^l Y(th,ph)
ok &= sweep("5.5.2 inversion", lambda l,m,th,ph: Y(l,m,math.pi-th,ph+math.pi)-(-1)**l*Y(l,m,th,ph))
# 5.5.6 rot pi about x: Y(pi-th, 2pi-ph) = (-1)^l Y_{l,-m}
ok &= sweep("5.5.6 rot-pi x", lambda l,m,th,ph: Y(l,m,math.pi-th,2*math.pi-ph)-(-1)**l*Y(l,-m,th,ph))
# 5.5.7 rot pi about y: Y(pi-th, pi-ph) = (-1)^{l-m} Y_{l,-m}
ok &= sweep("5.5.7 rot-pi y", lambda l,m,th,ph: Y(l,m,math.pi-th,math.pi-ph)-(-1)**(l-m)*Y(l,-m,th,ph))
# 5.5.8 rot pi about z: Y(th, pi+ph) = (-1)^m Y
ok &= sweep("5.5.8 rot-pi z", lambda l,m,th,ph: Y(l,m,th,math.pi+ph)-(-1)**m*Y(l,m,th,ph))
# 5.5.9 rot about z by chi: Y(th, ph-chi) = e^{-i m chi} Y
CHI=0.7
ok &= sweep("5.5.9 rot-z chi", lambda l,m,th,ph: Y(l,m,th,ph-CHI)-cmath.exp(-1j*m*CHI)*Y(l,m,th,ph))
# 5.5.11 equatorial reflection: Y(pi-th, ph) = (-1)^{l+m} Y
ok &= sweep("5.5.11 equatorial reflection", lambda l,m,th,ph: Y(l,m,math.pi-th,ph)-(-1)**(l+m)*Y(l,m,th,ph))
# 5.5.12 meridian reflection: Y(th, 2phi0-ph) = e^{i 2 m phi0} (-1)^m Y_{l,-m}
P0=0.5
ok &= sweep("5.5.12 meridian reflection", lambda l,m,th,ph:
            Y(l,m,th,2*P0-ph)-cmath.exp(1j*2*m*P0)*(-1)**m*Y(l,-m,th,ph))

print("\nSec 5.5.1 general rotation via Wigner D")
# Y_{l m'}(th', ph') = sum_m Y_{l m}(th, ph) D^l_{m m'}(al,be,ga)
# (th',ph') = direction of the fixed point expressed in the rotated frame S'.
def rotmat(al,be,ga):  # active z-y-z? build passive: components in S' = Rz(ga)Ry(be)Rz(al) x_S ... test convention
    ca,sa=math.cos(al),math.sin(al); cb,sb=math.cos(be),math.sin(be); cg,sg=math.cos(ga),math.sin(ga)
    Rz=lambda a:np.array([[math.cos(a),-math.sin(a),0],[math.sin(a),math.cos(a),0],[0,0,1]])
    Ry=lambda b:np.array([[math.cos(b),0,math.sin(b)],[0,1,0],[-math.sin(b),0,math.cos(b)]])
    return Rz, Ry
def dir_to_vec(th,ph): return np.array([math.sin(th)*math.cos(ph), math.sin(th)*math.sin(ph), math.cos(th)])
def vec_to_dir(v):
    r=np.linalg.norm(v); th=math.acos(max(-1,min(1,v[2]/r))); ph=math.atan2(v[1],v[0])%(2*math.pi); return th,ph
def Rzyz(al,be,ga):
    Rz=lambda a:np.array([[math.cos(a),-math.sin(a),0],[math.sin(a),math.cos(a),0],[0,0,1]])
    Ry=lambda b:np.array([[math.cos(b),0,math.sin(b)],[0,1,0],[-math.sin(b),0,math.cos(b)]])
    return Rz(al)@Ry(be)@Rz(ga)
AL,BE,GA=0.6,1.1,2.3
w=0.0
for l in range(0,5):
    for mp_ in range(-l,l+1):
        for th in TH:
            for ph in PH:
                n=dir_to_vec(th,ph)
                # passive: point coords in S' = R^{-1} n  (S' axes = R applied to S)
                np_=Rzyz(AL,BE,GA).T@n
                thp,php=vec_to_dir(np_)
                lhs=Y(l,mp_,thp,php)
                rhs=sum(Y(l,m,th,ph)*Dfun(l,m,mp_,AL,BE,GA) for m in range(-l,l+1))
                w=max(w,abs(lhs-rhs))
ok &= report("5.5.1 rotation (passive R^-1)", w)

print("\nSec 5.5.10 infinitesimal rotation (linear in omega)")
# D Y_lm(n) = Y_lm(R_{-w}^{... } n); check d/dw at 0 matches the ladder-operator bracket
def axis_rot(nvec, axis, w):  # rotate nvec about unit axis by angle w (Rodrigues)
    axis=axis/np.linalg.norm(axis)
    return (nvec*math.cos(w) + np.cross(axis,nvec)*math.sin(w) + axis*np.dot(axis,nvec)*(1-math.cos(w)))
THETA_n, PHI_n = 0.8, 1.4
nax = dir_to_vec(THETA_n, PHI_n)
def bracket(l,m,th,ph):  # the RHS {...} in 5.5.10
    cT=math.cos(THETA_n); sT=math.sin(THETA_n)
    return (m*cT*Y(l,m,th,ph)
            + sT/2*(cmath.exp(-1j*PHI_n)*math.sqrt(l*(l+1)-m*(m+1))*Y(l,m+1,th,ph)
                    + cmath.exp(1j*PHI_n)*math.sqrt(l*(l+1)-m*(m-1))*Y(l,m-1,th,ph)))
w=0.0; h=1e-6
for l in range(1,5):
    for m in range(-l,l+1):
        for th in TH:
            for ph in PH:
                n=dir_to_vec(th,ph)
                # (D f)(n)=f(R^{-1}n); R = rot about nax by +w  => R^{-1} = rot by -w
                def f(wv):
                    nr=axis_rot(n, nax, -wv)
                    thr,phr=vec_to_dir(nr)
                    return Y(l,m,thr,phr)
                dref=(f(h)-f(-h))/(2*h)         # d/dw (D Y)|_0
                w=max(w, abs(dref - (-1j)*bracket(l,m,th,ph)))
ok &= report("5.5.10 infinitesimal (d/domega)", w, tol=1e-6)

print("\nSec 5.5.3 parallel translation (solid-harmonic addition theorem)")
# Y_{l'm'}(th',ph') = sum_{l=0}^{l'} (-1)^{l'+l} C(l',l) (a/r')^{l'} (r/a)^l
#   {Y_l(th,ph) (x) Y_{l'-l}(Theta,Phi)}_{l'm'},  r'=r-a.
# BOOK MISPRINT: printed C = [4pi (2l+1)(2l'-2l+1)/(2l'+1)]^{1/2} is wrong for
# l'>=2 (coincides only at l'=1).  Correct: [4pi (2l'+1)!/((2l+1)!(2l'-2l+1)!)]^{1/2}
# (derived from the regular-solid-harmonic addition theorem; verified below).
from sympy.physics.wigner import clebsch_gordan as _clg
_cg2={}
def CG(j1,m1,j2,m2,j3,m3):
    k=(j1,m1,j2,m2,j3,m3)
    if k not in _cg2: _cg2[k]=float(_clg(j1,j2,j3,m1,m2,m3))
    return _cg2[k]
def tens(l1,l2,L,M,t1,p1,t2,p2):
    s=0j
    for mu in range(-l1,l1+1):
        nu=M-mu
        if abs(nu)>l2: continue
        s+=CG(l1,mu,l2,nu,L,M)*Y(l1,mu,t1,p1)*Y(l2,nu,t2,p2)
    return s
def dvec(r,t,p): return np.array([r*math.sin(t)*math.cos(p), r*math.sin(t)*math.sin(p), r*math.cos(t)])
def vdir(v):
    r=np.linalg.norm(v); return r, math.acos(max(-1,min(1,v[2]/r))), math.atan2(v[1],v[0])%(2*math.pi)
def facn(n): return math.factorial(n)
w=0.0
for (r,th,ph,a,TH,PH) in [(2.0,0.7,1.1,0.5,1.9,4.2),(1.5,1.2,0.3,0.9,0.6,2.0)]:
    rp,thp,php=vdir(dvec(r,th,ph)-dvec(a,TH,PH))
    for lp in range(0,5):
        for mp_ in range(-lp,lp+1):
            s=sum((-1)**(lp+l)
                  *math.sqrt(4*math.pi*facn(2*lp+1)/(facn(2*l+1)*facn(2*lp-2*l+1)))
                  *(a/rp)**lp*(r/a)**l*tens(l,lp-l,lp,mp_,th,ph,TH,PH)
                  for l in range(lp+1))
            w=max(w, abs(Y(lp,mp_,thp,php)-s))
ok &= report("5.5.3 translation [corrected coeff]", w)

print("\nRESULT:", "ALL PASS" if ok else "FAILURES ABOVE")
raise SystemExit(0 if ok else 1)
