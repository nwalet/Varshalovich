#!/usr/bin/env python3
r"""
Check of the algebraic Clebsch-Gordan tables 8.1-8.3 of Varshalovich, Moskalev
& Khersonskii (b = 1/2, 1, 3/2).

Each table entry is a closed-form expression for  C_{a alpha, b beta}^{c gamma}
as a function of (c, gamma); the row fixes c relative to a (e.g. row "a+1/2"
means c = a+1/2, i.e. a = c-1/2) and the column fixes beta.  We substitute many
numeric (c, gamma), recover a = c+off and alpha = gamma-beta, and compare the
formula with sympy's exact clebsch_gordan.

Usage:  python3 check_cg_tables.py
"""
from sympy import Rational, sqrt, S, simplify
from sympy.physics.wigner import clebsch_gordan as CG

H = Rational(1, 2)


def gammas(c):
    return [-c + i for i in range(int(2 * c) + 1)]


def exact(a, b, c, al, be, ga):
    if a < 0 or abs(al) > a or abs(be) > b or abs(ga) > c:
        return None
    if c < abs(a - b) or c > a + b or (a + b + c) != int(a + b + c):
        return None
    if not S(a - al).is_integer or not S(b - be).is_integer or not S(c - ga).is_integer:
        return None
    return CG(a, b, c, al, be, ga)


# entries: (offset so that a = c + offset, beta, formula(c, gamma))
TABLE_1 = ("Table 8.1  b=1/2", H, [
    (-H, H,  lambda c, g: sqrt((c + g) / (2 * c))),
    (-H, -H, lambda c, g: sqrt((c - g) / (2 * c))),
    (H,  H,  lambda c, g: -sqrt((c - g + 1) / (2 * c + 2))),
    (H,  -H, lambda c, g: sqrt((c + g + 1) / (2 * c + 2))),
])

TABLE_2 = ("Table 8.2  b=1", 1, [
    (-1, 1,  lambda c, g: sqrt((c + g - 1) * (c + g) / ((2 * c - 1) * 2 * c))),
    (-1, 0,  lambda c, g: sqrt((c + g) * (c - g) / ((2 * c - 1) * c))),
    (-1, -1, lambda c, g: sqrt((c - g - 1) * (c - g) / ((2 * c - 1) * 2 * c))),
    (0,  1,  lambda c, g: -sqrt((c + g) * (c - g + 1) / (2 * c * (c + 1)))),
    (0,  0,  lambda c, g: g / sqrt(c * (c + 1))),
    (0,  -1, lambda c, g: sqrt((c + g + 1) * (c - g) / (2 * c * (c + 1)))),
    (1,  1,  lambda c, g: sqrt((c - g + 1) * (c - g + 2) / ((2 * c + 2) * (2 * c + 3)))),
    (1,  0,  lambda c, g: -sqrt((c + g + 1) * (c - g + 1) / ((c + 1) * (2 * c + 3)))),
    (1,  -1, lambda c, g: sqrt((c + g + 2) * (c + g + 1) / ((2 * c + 2) * (2 * c + 3)))),
])

TABLE_3 = ("Table 8.3  b=3/2", Rational(3, 2), [
    (-3 * H, 3 * H, lambda c, g: sqrt((c + g - 2) * (c + g - 1) * (c + g) / ((2 * c - 2) * (2 * c - 1) * 2 * c))),
    (-3 * H, H,     lambda c, g: sqrt(3 * (c + g - 1) * (c + g) * (c - g) / ((2 * c - 2) * (2 * c - 1) * 2 * c))),
    (-3 * H, -H,    lambda c, g: sqrt(3 * (c + g) * (c - g - 1) * (c - g) / ((2 * c - 2) * (2 * c - 1) * 2 * c))),
    (-H, 3 * H, lambda c, g: -sqrt(3 * (c + g - 1) * (c + g) * (c - g + 1) / ((2 * c - 1) * 2 * c * (2 * c + 2)))),
    (-H, H,     lambda c, g: -(c - 3 * g + 1) * sqrt((c + g) / ((2 * c - 1) * 2 * c * (2 * c + 2)))),
    (-H, -H,    lambda c, g: (c + 3 * g + 1) * sqrt((c - g) / ((2 * c - 1) * 2 * c * (2 * c + 2)))),
    (H, 3 * H, lambda c, g: sqrt(3 * (c + g) * (c - g + 1) * (c - g + 2) / (2 * c * (2 * c + 2) * (2 * c + 3)))),
    (H, H,     lambda c, g: -(c + 3 * g) * sqrt((c - g + 1) / (2 * c * (2 * c + 2) * (2 * c + 3)))),
    (H, -H,    lambda c, g: -(c - 3 * g) * sqrt((c + g + 1) / (2 * c * (2 * c + 2) * (2 * c + 3)))),
    (3 * H, 3 * H, lambda c, g: -sqrt((c - g + 1) * (c - g + 2) * (c - g + 3) / ((2 * c + 2) * (2 * c + 3) * (2 * c + 4)))),
    (3 * H, H,     lambda c, g: sqrt(3 * (c + g + 1) * (c - g + 1) * (c - g + 2) / ((2 * c + 2) * (2 * c + 3) * (2 * c + 4)))),
    (3 * H, -H,    lambda c, g: -sqrt(3 * (c + g + 1) * (c + g + 2) * (c - g + 1) / ((2 * c + 2) * (2 * c + 3) * (2 * c + 4)))),
])


def check(table):
    name, b, entries = table
    cvals = [Rational(i, 2) for i in range(2, 13)]         # c = 1 .. 6 (half steps)
    tested = bad = 0
    first = None
    for off, be, f in entries:
        for c in cvals:
            a = c + off
            for g in gammas(c):
                al = g - be
                ex = exact(a, b, c, al, be, g)
                if ex is None:
                    continue
                fv = f(c, g)
                if not fv.is_real:                          # sqrt of a negative -> skip
                    continue
                tested += 1
                if simplify(fv - ex) != 0:
                    bad += 1
                    if first is None:
                        first = (off, be, c, g, fv, ex)
    ok = bad == 0 and tested > 0
    print(f"  [{'OK  ' if ok else 'FAIL'}] {name:22s} {tested} entries checked, {bad} mismatched")
    if first:
        off, be, c, g, fv, ex = first
        print(f"        e.g. a=c{'+' if off >= 0 else ''}{off}, beta={be}, c={c}, gamma={g}: table={fv}  exact={ex}")
    return ok


def run():
    print("Algebraic Clebsch-Gordan tables 8.1, 8.2\n")
    ok = True
    for tb in (TABLE_1, TABLE_2):
        ok &= check(tb)
    print("\nALL TABLE ENTRIES CORRECT" if ok else "\nSOME ENTRIES MISMATCH -- see above")
    return ok


if __name__ == "__main__":
    raise SystemExit(0 if run() else 1)
