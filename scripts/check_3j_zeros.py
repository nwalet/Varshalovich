#!/usr/bin/env python3
r"""
Check of the table "Accidental zeroes of the 3j coefficients for J<=17" in
Chapter 8 of Varshalovich, Moskalev & Khersonskii.

Every 3jm symbol listed in the table is claimed to vanish.  For each entry we
verify:
  * it is a physical 3jm  (m1+m2+m3=0, |mi|<=ji, triangle, integer perimeter);
  * J = j1+j2+j3 matches the row it is listed under;
  * wigner_3j(...) == 0.
Any entry that is non-physical, mis-filed, or actually non-zero is flagged
(such an entry would signal an OCR error in its arguments).

Usage:  python3 check_3j_zeros.py
"""
from sympy import Rational, S
from sympy.physics.wigner import wigner_3j


def num(tok):
    if "/" in tok:
        n, d = tok.split("/")
        return Rational(int(n), int(d))
    return Rational(int(tok))


# One string per 3jm symbol: "j1 j2 j3 m1 m2 m3", grouped by the table's J row.
TABLE = {
    3: ["1 1 1 0 0 0"],
    5: ["1 2 2 0 0 0", "3/2 3/2 2 1/2 1/2 -1"],
    7: ["1 3 3 0 0 0", "2 2 3 0 0 0", "2 2 3 1 1 -2", "2 5/2 5/2 -1 1/2 1/2"],
    8: ["2 3 3 0 2 -2", "3/2 3 7/2 1/2 1 -3/2"],
    9: ["1 4 4 0 0 0", "2 3 4 0 0 0", "3 3 3 0 0 0", "3 3 3 1 1 -2",
        "2 7/2 7/2 -1 1/2 1/2", "5/2 3 7/2 3/2 0 -3/2",
        "5/2 5/2 4 1/2 1/2 -1", "5/2 5/2 4 3/2 3/2 -3"],
    11: ["1 5 5 0 0 0", "2 4 5 0 0 0", "2 4 5 1 2 -3", "3 3 5 0 0 0",
         "3 3 5 1 1 -2", "3 3 5 2 2 -4", "3 4 4 0 0 0", "3 4 4 -2 1 1",
         "2 9/2 9/2 -1 1/2 1/2", "3 7/2 9/2 2 1/2 -5/2",
         "7/2 7/2 4 1/2 1/2 -1", "7/2 7/2 4 3/2 3/2 -3",
         "3/2 9/2 5 1/2 3/2 -2", "5/2 4 9/2 1/2 3 -7/2"],
    13: ["1 6 6 0 0 0", "2 5 6 0 0 0", "3 4 6 0 0 0", "3 5 5 0 0 0",
         "3 5 5 -2 1 1", "4 4 5 0 0 0", "4 4 5 1 1 -2", "4 4 5 2 2 -4",
         "2 11/2 11/2 -1 1/2 1/2", "7/2 9/2 5 5/2 -1/2 -2",
         "7/2 4 11/2 5/2 1 -7/2", "4 9/2 9/2 -1 1/2 1/2",
         "4 9/2 9/2 -3 3/2 3/2", "7/2 7/2 6 1/2 1/2 -1",
         "7/2 7/2 6 3/2 3/2 -3", "7/2 7/2 6 5/2 5/2 -5"],
    14: ["3 5 6 1 -1 0", "3 5 6 2 1 -3", "3 5 6 1 4 -5", "4 4 6 0 2 -2",
         "4 5 5 1 3 -4", "4 5 5 2 1 -3", "3 9/2 13/2 1 3/2 -5/2",
         "5/2 5 13/2 1/2 1 -3/2", "5/2 5 13/2 3/2 3 -9/2",
         "4 9/2 11/2 0 7/2 -7/2", "3/2 6 13/2 1/2 2 -5/2",
         "5/2 11/2 6 3/2 1/2 -2"],
    15: ["1 7 7 0 0 0", "2 6 7 0 0 0", "2 6 7 1 3 -4", "3 5 7 0 0 0",
         "3 6 6 0 0 0", "3 6 6 0 5 -5", "3 6 6 -2 1 1", "4 5 6 0 0 0",
         "4 5 6 3 0 -3", "4 4 7 0 0 0", "4 4 7 1 1 -2", "4 4 7 2 2 -4",
         "4 4 7 3 3 -6", "5 5 5 0 0 0", "5 5 5 1 1 -2", "5 5 5 2 2 -4",
         "2 13/2 13/2 -1 1/2 1/2", "4 9/2 13/2 3 3/2 -9/2",
         "4 11/2 11/2 -1 1/2 1/2", "4 11/2 11/2 -3 3/2 3/2",
         "9/2 5 11/2 3/2 0 -3/2", "9/2 9/2 6 1/2 1/2 -1",
         "9/2 9/2 6 3/2 3/2 -3", "9/2 9/2 6 5/2 5/2 -5"],
    17: ["1 8 8 0 0 0", "2 7 8 0 0 0", "3 6 8 0 0 0", "3 6 8 2 4 -6",
         "3 7 7 0 0 0", "3 7 7 -2 1 1", "4 5 8 0 0 0", "4 6 7 0 0 0",
         "5 5 7 0 0 0", "5 5 7 1 1 -2", "5 5 7 2 2 -4", "5 5 7 3 3 -6",
         "5 6 6 0 0 0", "5 6 6 -2 1 1", "5 6 6 -4 2 2",
         "2 15/2 15/2 -1 1/2 1/2", "3 13/2 15/2 2 3/2 -7/2",
         "4 13/2 13/2 -1 1/2 1/2", "4 13/2 13/2 -3 3/2 3/2",
         "5 11/2 13/2 2 1/2 -5/2", "7/2 6 15/2 3/2 5 -13/2",
         "9/2 6 13/2 1/2 -5 9/2", "9/2 6 13/2 -7/2 1 5/2",
         "11/2 11/2 6 1/2 1/2 -1", "11/2 11/2 6 3/2 3/2 -3",
         "11/2 11/2 6 5/2 5/2 -5", "9/2 11/2 7 7/2 1/2 -4",
         "3/2 15/2 8 1/2 5/2 -3", "9/2 9/2 8 1/2 1/2 -1",
         "9/2 9/2 8 3/2 3/2 -3", "9/2 9/2 8 5/2 5/2 -5",
         "9/2 9/2 8 7/2 7/2 -7", "9/2 5 15/2 7/2 2 -11/2"],
}


def classify(entry, Jrow):
    j1, j2, j3, m1, m2, m3 = (num(t) for t in entry.split())
    problems = []
    if m1 + m2 + m3 != 0:
        problems.append("sum m != 0")
    if abs(m1) > j1 or abs(m2) > j2 or abs(m3) > j3:
        problems.append("|m|>j")
    if j3 < abs(j1 - j2) or j3 > j1 + j2 or not S(j1 + j2 + j3).is_integer:
        problems.append("triangle")
    if j1 + j2 + j3 != Jrow:
        problems.append(f"J={j1 + j2 + j3}!={Jrow}")
    if not problems:
        if wigner_3j(j1, j2, j3, m1, m2, m3) != 0:
            problems.append("NON-ZERO")
    return problems


def run():
    print("Accidental 3j zeros table (J <= 17)\n")
    total = flagged = 0
    for Jrow in sorted(TABLE):
        n = len(TABLE[Jrow]); bad = []
        for e in TABLE[Jrow]:
            total += 1
            probs = classify(e, Jrow)
            if probs:
                flagged += 1
                bad.append((e, probs))
        mark = "OK  " if not bad else "FAIL"
        print(f"  [{mark}] J={Jrow:<2d} {n} symbols" + ("" if not bad else ":"))
        for e, probs in bad:
            print(f"          ({e})  -> {', '.join(probs)}")
    print(f"\n{total} symbols checked, {flagged} flagged")
    print("ALL LISTED 3j SYMBOLS VANISH" if flagged == 0 else "SOME ENTRIES ARE PROBLEMATIC")
    return flagged == 0


if __name__ == "__main__":
    raise SystemExit(0 if run() else 1)
