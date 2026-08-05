#!/usr/bin/env python3
r"""
Check of the Clebsch-Gordan "zero" selection rules of Section 8.10 of
Varshalovich, Moskalev & Khersonskii, "Quantum Theory of Angular Momentum".

Each rule (a)-(g) asserts that a Clebsch-Gordan coefficient vanishes whenever a
stated condition on its arguments holds (beyond the triangle rule).  For every
rule we enumerate all small, physically-valid argument sets that satisfy the
condition and confirm the coefficient is exactly zero (via sympy).  The count
of matching configurations is reported so the test is non-vacuous.

Three coefficients are OCR-damaged in the source; we test the reconstructed
form and flag the fix needed:
  (b) 3rd coeff  C_{c gamma, alpha, -alpha}^{a alpha}  -- the "alpha" momentum
                 should be "a":  C_{c gamma, a -alpha}^{a alpha}.
  (f)            C_{a a, b beta}^{c,c-1}  -- the m1="a" should be "alpha":
                 C_{a alpha, b beta}^{c,c-1}.
  (g)            \clebsch{a}{alpha}{a}{b}{b-1}{c}{gamma}  has 7 arguments; the
                 spurious 3rd token "{a}" should go: C_{a alpha, b,b-1}^{c gamma}.

Usage:  python3 check_8_10.py
"""
from sympy import Rational, S
from sympy.physics.wigner import clebsch_gordan as CG

HALF = Rational(1, 2)


def halfs(hi):
    """0, 1/2, 1, ... , hi."""
    return [Rational(i, 2) for i in range(0, int(2 * hi) + 1)]


def ints(hi):
    return list(range(0, int(hi) + 1))


def proj(j):
    return [-j + i for i in range(int(2 * j) + 1)]


def tri(a, b, c):
    return abs(a - b) <= c <= a + b and (a + b + c) == int(a + b + c)


def _int(x):
    return S(x).is_integer


def C(a, b, c, al, be, ga):
    """Clebsch-Gordan, 0 outside the physical domain."""
    if a < 0 or b < 0 or c < 0 or abs(al) > a or abs(be) > b or abs(ga) > c:
        return S.Zero
    if not tri(a, b, c) or al + be != ga:
        return S.Zero
    if not _int(a - al) or not _int(b - be) or not _int(c - ga):
        return S.Zero
    return CG(a, b, c, al, be, ga)


def report(tag, matched, nonzero, offnz=None):
    # a rule is confirmed if every on-condition coeff is zero (matched>0), and
    # -- for a meaningful (tight) rule -- some off-condition coeff is non-zero.
    ok = nonzero == 0 and matched > 0
    mark = "OK  " if ok else ("SKIP" if matched == 0 else "FAIL")
    tail = f", {nonzero} non-zero" + (f"; tight ({offnz} off-cond non-zero)" if offnz is not None else "")
    print(f"  [{mark}] {tag:50s} {matched} configs{tail}")
    return ok


def rule_a(J=8):
    # C_{a0,b0}^{c0} = 0 if a+b+c odd   (off: a+b+c even -> generally non-zero)
    m = z = off = 0
    for a in ints(J):
        for b in ints(J):
            for c in ints(J):
                if not tri(a, b, c):
                    continue
                cg = C(a, b, c, 0, 0, 0)
                if (a + b + c) % 2:
                    m += 1; z += (cg != 0)
                else:
                    off += (cg != 0)
    return report("(a)  C_{a0,b0}^{c0}=0, a+b+c odd", m, z, off)


def rule_b(hi=4):
    # 2a+c odd -> three (permuted) coefficients vanish
    m = [0, 0, 0]; z = [0, 0, 0]; off = [0, 0, 0]
    for a in halfs(hi):
        for c in halfs(hi):
            odd = int(2 * a + c) % 2
            for al in proj(a):
                ga = 2 * al
                if abs(ga) > c:
                    continue
                coeffs = []
                if tri(a, a, c):
                    coeffs.append((0, C(a, a, c, al, al, ga)))
                if tri(a, c, a):
                    coeffs.append((1, C(a, c, a, -al, ga, al)))
                    coeffs.append((2, C(c, a, a, ga, -al, al)))     # (b) 3rd, reconstructed
                for i, cg in coeffs:
                    if odd:
                        m[i] += 1; z[i] += (cg != 0)
                    else:
                        off[i] += (cg != 0)
    ok = report("(b1) C_{a al,a al}^{c ga}=0, 2a+c odd", m[0], z[0], off[0])
    ok &= report("(b2) C_{a -al,c ga}^{a al}=0", m[1], z[1], off[1])
    ok &= report("(b3) C_{c ga,a -al}^{a al}=0  [OCR: alpha->a]", m[2], z[2], off[2])
    return ok


def rule_c(J=17):
    # C_{a,J-3a, b,J-3b}^{c,-(J-3c)} = 0 if J odd   (off: J even)
    m = z = off = 0
    for a in ints(J):
        for b in ints(J):
            for c in ints(J):
                if not tri(a, b, c) or a + b + c > J:
                    continue
                Jv = a + b + c
                al, be = Jv - 3 * a, Jv - 3 * b
                if abs(al) > a or abs(be) > b:
                    continue
                cg = C(a, b, c, al, be, -(Jv - 3 * c))
                if Jv % 2:
                    m += 1; z += (cg != 0)
                else:
                    off += (cg != 0)
    return report("(c)  C_{a,J-3a,b,J-3b}^{c,-(J-3c)}=0, J odd", m, z, off)


def rule_d(hi=5):
    # C_{a al,b be}^{a+b-1} = 0 if b al = a be   (off: b al != a be)
    m = z = off = 0
    for a in halfs(hi):
        for b in halfs(hi):
            c = a + b - 1
            if c < abs(a - b) or c < 0:
                continue
            for al in proj(a):
                for be in proj(b):
                    if abs(al + be) > c:
                        continue
                    cg = C(a, b, c, al, be, al + be)
                    if b * al == a * be:
                        m += 1; z += (cg != 0)
                    else:
                        off += (cg != 0)
    return report("(d)  C_{a al,b be}^{a+b-1}=0, al/be=a/b", m, z, off)


def rule_e(hi=5):
    # C_{a al,b be}^{a-b+1} = 0 if b al = -(a+1) be   (off: otherwise)
    m = z = off = 0
    for a in halfs(hi):
        for b in halfs(hi):
            c = a - b + 1
            if c < abs(a - b) or c < 0:
                continue
            for al in proj(a):
                for be in proj(b):
                    if abs(al + be) > c:
                        continue
                    cg = C(a, b, c, al, be, al + be)
                    if b * al == -(a + 1) * be:
                        m += 1; z += (cg != 0)
                    else:
                        off += (cg != 0)
    return report("(e)  C_{a al,b be}^{a-b+1}=0, al/be=-(a+1)/b", m, z, off)


def rule_f(hi=5):
    # C_{a al,b be}^{c,c-1} = 0 if a(a+1)-b(b+1) = (al-be) c   [OCR: {a}{a}->{a}{al}]
    m = z = off = 0
    for a in halfs(hi):
        for b in halfs(hi):
            for c in halfs(hi):
                if c < 1 or not tri(a, b, c):
                    continue
                for al in proj(a):
                    be = (c - 1) - al                       # m3 = c-1
                    if abs(be) > b:
                        continue
                    cg = C(a, b, c, al, be, c - 1)
                    if a * (a + 1) - b * (b + 1) == (al - be) * c:
                        m += 1; z += (cg != 0)
                    else:
                        off += (cg != 0)
    return report("(f)  C_{a al,b be}^{c,c-1}=0  [OCR: m1 a->al]", m, z, off)


def rule_g(hi=5):
    # C_{a al,b,b-1}^{c ga} = 0 if a(a+1)-c(c+1) = -(al+ga) b   [OCR: drop extra {a}]
    m = z = off = 0
    for a in halfs(hi):
        for b in halfs(hi):
            if b < 1:
                continue
            be = b - 1
            for c in halfs(hi):
                if not tri(a, b, c):
                    continue
                for al in proj(a):
                    ga = al + be
                    if abs(ga) > c:
                        continue
                    cg = C(a, b, c, al, be, ga)
                    if a * (a + 1) - c * (c + 1) == -(al + ga) * b:
                        m += 1; z += (cg != 0)
                    else:
                        off += (cg != 0)
    return report("(g)  C_{a al,b,b-1}^{c ga}=0  [OCR: 7-arg -> 6]", m, z, off)


def run():
    print("Section 8.10 zero-selection-rule checks\n")
    all_ok = True
    all_ok &= rule_a()
    all_ok &= rule_b()
    all_ok &= rule_c()
    all_ok &= rule_d()
    all_ok &= rule_e()
    all_ok &= rule_f()
    all_ok &= rule_g()
    print("\nALL RULES HOLD" if all_ok else "\nSOME RULES FAILED -- see above")
    return all_ok


if __name__ == "__main__":
    raise SystemExit(0 if run() else 1)
