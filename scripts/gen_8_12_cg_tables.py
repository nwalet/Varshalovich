#!/usr/bin/env python3
r"""Generate any algebraic Clebsch-Gordan entry of the VMK Sec. 8.12 tables.

Each table fixes a small b and tabulates, with a and alpha left SYMBOLIC,

    <a alpha, b beta | c gamma>,   c = a + k,   gamma = alpha + beta,

as a function of the book's free variables c and gamma.  Rows are indexed by
k = c - a  (k = b, b-1, ..., -b), columns by beta (beta = b, b-1, ..., -b).

Method
------
Racah's formula has only min(b-k, b+beta) - max(0, beta-k) + 1 <= 2b+1 terms
once b, beta and k are fixed, so the whole entry is a short symbolic sum.

The one trick that makes it robust: change variables to

    p = a + alpha,   q = a - alpha

which are non-negative INTEGERS even when a and alpha are half-integral.
Every factorial argument is then an explicit integer offset from p, q or p+q,
so each (x+m)!/x! collapses to a rising/falling factorial and the factorials
disappear entirely up front.  What remains is pure rational-function algebra:
no combsimp, no gammasimp, no simplify() in the hot path.  The p!q! generated
by the Racah sum cancels exactly against the one under the square root,
leaving  CG = T(p,q) * sqrt(W(p,q))  with T and W explicit rational functions.

Book form
---------
CG = T*sqrt(W) is not a unique split -- any factor may be moved inside the
radical as its square.  The book absorbs every factor that is sign-definite on
the physical domain |gamma| <= c (c+gamma, c-gamma+1, ...) and leaves outside
only the sign-indefinite ones (gamma, c+3*gamma+1, ...).  decompose() applies
that rule, sampling each factor over the physical wedge to classify it and
fixing the overall sign at a physical point.

Usage
-----
    python3 gen_8_12_cg_tables.py --b 1 --beta=0 --k=0
    python3 gen_8_12_cg_tables.py --b 3/2                 # whole table
    python3 gen_8_12_cg_tables.py --b 1 --beta=0 --k=0 --vars a --format sympy

Note the '=' in --beta=-1/2: argparse would otherwise read a leading '-' as a
flag.  b may be any non-negative half-integer; the book's b <= 5 is not a limit.
"""
import sys

from sympy import (Rational, symbols, RisingFactorial, FallingFactorial, factorial,
                   cancel, factor, factor_list, sqrt, latex, fraction, Integer,
                   lambdify, nsimplify)

p, q = symbols('p q')                    # p = a+alpha, q = a-alpha
c, g = symbols('c gamma')                # the book's free variables
a, al = symbols('a alpha')               # optional alternative variables

__all__ = ['cg_algebraic', 'decompose', 'entry', 'render', 'p', 'q', 'c', 'g']


# --------------------------------------------------------------------------
# the Racah engine
# --------------------------------------------------------------------------
def _fr(x, m):
    """(x+m)! / x!  as an explicit rational function; m an explicit integer."""
    m = int(m)
    return RisingFactorial(x + 1, m) if m >= 0 else 1 / FallingFactorial(x, -m)


def _nfact(n):
    """n! for an explicit integer n, or None if n < 0 (term vanishes)."""
    n = int(n)
    return None if n < 0 else factorial(n)


def cg_algebraic(b, beta, k):
    """Return (T, W) with  CG = T*sqrt(W),  as rational functions of p and q.

    Returns None when the entry is identically zero.
    """
    b, beta, k = Rational(b), Rational(beta), Rational(k)
    if abs(k) > b or abs(beta) > b or (b - k) % 1 or (b - beta) % 1:
        return None
    s = p + q                                            # = 2a

    # Racah sum, with the factorials of p and q divided out as 1/(p! q!).
    T = Integer(0)
    z, zhi = max(Integer(0), beta - k), min(b - k, b + beta)
    while z <= zhi:
        denom, ok = Integer(1), True
        for n in (z, b - k - z, b + beta - z, k - beta + z):
            f = _nfact(n)
            if f is None:
                ok = False
                break
            denom *= f
        if ok:
            #  1/(q-z)!     = FallingFactorial(q, z) / q!
            #  1/(p+k-b+z)! = 1 / (_fr(p, k-b+z) * p!)
            T += Integer(-1)**int(z) * FallingFactorial(q, z) / (denom * _fr(p, k - b + z))
        z += 1
    T = cancel(T)
    if T == 0:
        return None

    W = ((s + 2*k + 1)                                   # 2c+1
         * _nfact(b - k) * _nfact(b + k)                 # (a+b-c)! (-a+b+c)!
         / _fr(s - b + k, 2*b + 1)                       # (a-b+c)! / (a+b+c+1)!
         * _nfact(b + beta) * _nfact(b - beta)
         * _fr(p, k + beta)                              # (c+gamma)! / p!
         * _fr(q, k - beta))                             # (c-gamma)! / q!
    return T, cancel(W)


def _to_book_vars(expr, beta, k):
    """p,q -> c,gamma  via  p = c-k+gamma-beta,  q = c-k-gamma+beta."""
    return cancel(expr.subs({p: c - k + g - beta, q: c - k - g + beta}))


def _definite(f):
    """+1 if f >= 0 on the physical domain |gamma| <= c, -1 if f <= 0, else 0."""
    fn = lambdify((c, g), f, 'math')
    pos = neg = False
    for cv in range(6, 16):
        for gv in range(-cv, cv + 1):
            try:
                v = fn(float(cv), float(gv))
            except (ValueError, ZeroDivisionError):
                continue
            if v > 1e-12:
                pos = True
            elif v < -1e-12:
                neg = True
    if pos and neg:
        return 0
    return -1 if neg else 1


# --------------------------------------------------------------------------
# book-form decomposition -- the single source of truth for render() and the
# verifier in check_8_12.py, so the checked path is the shipped path
# --------------------------------------------------------------------------
def decompose(b, beta, k):
    """Return (sign, outside, inside) with  CG = sign * outside * sqrt(inside).

    Sign-definite factors are absorbed into the radical, matching the book's
    +/-[...]^{1/2} presentation.  Returns None for an identically-zero entry.
    """
    res = cg_algebraic(b, beta, k)
    if res is None:
        return None
    beta, k = Rational(beta), Rational(k)
    T, W = (_to_book_vars(e, beta, k) for e in res)

    const, facs = factor_list(T)
    inside, outside, sign = W * const**2, Integer(1), (-1 if const < 0 else 1)
    for f, mult in facs:
        d = _definite(f)
        if d == 0:
            outside *= f**mult
        else:
            inside *= (f**2)**mult
            if d < 0 and mult % 2:
                sign = -sign
    return sign, factor(outside), factor(cancel(inside))


def entry(b, beta, k, vars='c'):
    """The entry as an exact sympy expression.

    vars='c' -> in terms of c and gamma (the book's form)
    vars='a' -> in terms of a and alpha  (c = a+k, gamma = alpha+beta)
    """
    dec = decompose(b, beta, k)
    if dec is None:
        return Integer(0)
    sign, outside, inside = dec
    expr = sign * outside * sqrt(inside)
    if vars == 'a':
        expr = expr.subs({c: a + Rational(k), g: al + Rational(beta)}, simultaneous=True)
    return expr


def render(b, beta, k):
    """LaTeX for one entry, in the book's  +/-[ ... ]^{1/2}  style."""
    dec = decompose(b, beta, k)
    if dec is None:
        return "0"
    sign, outside, inside = dec

    num, den = fraction(inside)
    num, den = factor(num), factor(den)
    pre = "-" if sign < 0 else ""
    out = "" if outside == 1 else latex(outside)
    if outside != 1 and outside.is_Add:              # a multi-term factor needs brackets
        out = r"\left(%s\right)" % out

    if num == 1 and den != 1:                        # the  gamma/[c(c+1)]^{1/2}  shape
        return r"%s\frac{%s}{\left[%s\right]^{1/2}}" % (pre, out or "1", latex(den))
    body = (r"\left[\frac{%s}{%s}\right]^{1/2}" % (latex(num), latex(den))
            if den != 1 else r"\left[%s\right]^{1/2}" % latex(num))
    return f"{pre}{out}{body}" if out == "" else f"{pre}{out}\\,{body}"


# --------------------------------------------------------------------------
def _fmt(x):
    """'3/2', '-1', '0' ..."""
    x = Rational(x)
    return str(x)


def _row_label(k):
    k = Rational(k)
    if k == 0:
        return "a"
    return f"a{'+' if k > 0 else '-'}{_fmt(abs(k))}"


def _one(b, beta, k, vars, fmt):
    beta = Rational(beta)
    gam = "alpha" if beta == 0 else f"alpha{'+' if beta > 0 else '-'}{_fmt(abs(beta))}"
    print(f"<a alpha, {_fmt(b)} {_fmt(beta)} | c={_row_label(k)}, gamma={gam}>")
    if fmt in ('sympy', 'both'):
        print("  expr :", entry(b, beta, k, vars))
    if fmt in ('latex', 'both'):
        print("  latex:", render(b, beta, k))


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(
        description="Generate algebraic Clebsch-Gordan entries (VMK Sec. 8.12).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Use '=' for negative values:  --beta=-1/2 --k=-1\n"
               "b may be any non-negative half-integer (the book stops at 5).")
    ap.add_argument('--b', required=True, help="the fixed momentum b, e.g. 1 or 3/2")
    ap.add_argument('--beta', help="column: beta (omit with --k to dump the whole table)")
    ap.add_argument('--k', help="row: k = c - a, e.g. 0 or -1/2")
    ap.add_argument('--vars', choices=['c', 'a'], default='c',
                    help="express in c,gamma (default, the book's form) or a,alpha")
    ap.add_argument('--format', choices=['latex', 'sympy', 'both'], default='both')
    args = ap.parse_args(argv)

    b = Rational(args.b)
    if b < 0 or (2*b) % 1:
        ap.error("b must be a non-negative half-integer")

    if (args.beta is None) != (args.k is None):
        ap.error("give both --beta and --k, or neither (to dump the table)")

    if args.beta is not None:
        bb, kk = Rational(args.beta), Rational(args.k)
        if abs(bb) > b or abs(kk) > b or (b - bb) % 1 or (b - kk) % 1:
            ap.error(f"beta and k must run over {_fmt(b)}, {_fmt(b-1)}, ..., {_fmt(-b)}")
        _one(b, bb, kk, args.vars, args.format)
        return 0

    vals = [b - i for i in range(int(2*b) + 1)]
    print(f"Table for b = {_fmt(b)}   (rows c = a+k, columns beta)\n")
    for kk in vals:
        for bb in vals:
            print(f"c = {_row_label(kk):<8} beta = {_fmt(bb):<5} "
                  f"{render(b, bb, kk)}")
        print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
