#!/usr/bin/env python3
r"""
Check for Section 10.4 (symmetry properties of the 9j symbol), Chapter 10,
Varshalovich, Moskalev & Khersonskii.

The big table eq 10.4.6 lists 72 (= 3!x3!x2) 9j symbols all equal to the
reference {a b c / d e f / g h j} up to the sign  eps = 1 (even permutation)
or (-1)^R (odd), R = sum of all nine arguments.  We extract every \ninej call
(and whether it carries an \varepsilon prefactor) straight from Chap10.tex and
verify each equals eps * (reference 9j) for several numeric assignments.  A
mismatch flags an OCR error in the printed table.

Usage:  python3 check_10_4.py
"""
import re
from sympy import Rational, S
from sympy.physics.wigner import wigner_9j

H = Rational(1, 2)


def tri(a, b, c):
    return abs(a - b) <= c <= a + b and (a + b + c) == int(a + b + c)


def valid9(v):
    a, b, c, d, e, f, g, h, j = v
    return all([tri(a, b, c), tri(d, e, f), tri(g, h, j),
                tri(a, d, g), tri(b, e, h), tri(c, f, j)])


def w9(v):
    return wigner_9j(*v) if valid9(v) else S.Zero


def extract_table():
    """Pull the eq 10.4.6 block out of Chap10.tex and return list of
    (eps_flag, [9 letters]) for every \\ninej in it."""
    text = open('../Chap10.tex').read()
    # region: from the paragraph introducing 729/72 symbols up to the r-symbol subsection
    start = text.index('are arranged horizontally')
    end = text.index('Symmetries of the $r$ Symbol')
    region = text[start:end]
    out = []
    for m in re.finditer(r'(\\varepsilon\s*)?\\ninej((?:\{[^{}]*\}){9})', region):
        eps = m.group(1) is not None
        args = re.findall(r'\{([^{}]*)\}', m.group(2))
        out.append((eps, [a.strip() for a in args]))
    return out


ASSIGN = [
    {'a': 1, 'b': 2, 'c': 2, 'd': 2, 'e': 1, 'f': 2, 'g': 2, 'h': 2, 'j': 1},
    {'a': 1, 'b': H, 'c': H, 'd': H, 'e': 1, 'f': H, 'g': H, 'h': H, 'j': 1},
    {'a': 2, 'b': 1, 'c': 2, 'd': 1, 'e': 2, 'f': 2, 'g': 2, 'h': 1, 'j': 2},
    {'a': Rational(3, 2), 'b': 1, 'c': H, 'd': 1, 'e': 1, 'f': 1, 'g': H, 'h': 1, 'j': Rational(3, 2)},
]


def run():
    print("Section 10.4 : eq 10.4.6 symmetry table\n")
    table = extract_table()
    print(f"  extracted {len(table)} \\ninej entries from the eq 10.4.6 block")
    ok = True
    bad_list = []
    for sub in ASSIGN:
        ref_v = [sub[k] for k in 'abcdefghj']
        R = sum(ref_v)
        ref = w9(ref_v)
        for idx, (eps, args) in enumerate(table):
            try:
                v = [sub[a] for a in args]
            except KeyError:
                bad_list.append((idx, args, 'non-letter arg'))
                ok = False
                continue
            eps_val = (-1) ** R if eps else 1
            if w9(v) != eps_val * ref:
                ok = False
                if (idx, tuple(args)) not in [(b[0], tuple(b[1])) for b in bad_list]:
                    bad_list.append((idx, args, f"eps={'(-1)^R' if eps else '1'}"))
    print(f"  [{'OK  ' if ok else 'FAIL'}] all {len(table)} entries equal eps * reference")
    if bad_list:
        print("  mismatches:")
        for idx, args, note in bad_list[:20]:
            print(f"    #{idx}: {{{' '.join(args)}}}  ({note})")
    print("\nALL 10.4.6 ENTRIES CORRECT" if ok else "\nSOME 10.4.6 ENTRIES MISMATCH")
    return ok


if __name__ == "__main__":
    raise SystemExit(0 if run() else 1)
