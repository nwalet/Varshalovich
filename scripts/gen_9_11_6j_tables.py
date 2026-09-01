#!/usr/bin/env python3
r"""Generate any algebraic 6j entry of the VMK Sec. 9.11 tables (Tables 9.1-9.8).

Each table fixes a small d and tabulates, with a, b, c left SYMBOLIC,

    {a b c; d e f},    f = b + m,    e = c + n,

as a function of the book's notation

    s = a+b+c,    X = -a(a+1) + b(b+1) + c(c+1).

Rows are indexed by m (f = b+m), columns by n (e = c+n); both run over
d, d-1, ..., -d.  Tables 9.1-9.8 are d = 1/2, 1, 3/2, 2, 5/2, 3, 7/2, 4.

Method
------
Work in the triangle deficits

    u = s - 2a = -a+b+c,   v = s - 2b = a-b+c,   w = s - 2c = a+b-c,

all non-negative INTEGERS even for half-integral a, b, c, with s = u+v+w.
Putting z = s + t in Racah's single sum turns every denominator factorial into
an explicit integer offset from t, u, v or w, e.g.

    z-a-b-c = t          a+b+d+e-z = d+n-t
    z-a-e-f = t-n-m      b+c+e+f-z = u+n+m-t
    z-d-b-f = v+t-d-m    a+c+d+f-z = d+m-t
    z-d-e-c = w+t-d-n    (z+1)!     = (s+t+1)!

and t runs over the explicit finite range max(0, n+m) <= t <= min(d+n, d+m),
so the sum has at most 2d+1 terms.  Each (x+j)!/x! then collapses to a
rising/falling factorial: no combsimp, no gammasimp, no simplify().

The symbolic factorials cancel exactly.  Delta(abc) supplies
sqrt(u!v!w!/(s+1)!), the sum supplies (s+1)!/(u!v!w!), and Delta(aef) supplies
the inverse, leaving

    {6j} = (-1)^s * T(u,v,w) * sqrt(W(u,v,w))

with T and W explicit rational functions.  The global (-1)^s appears because
(-1)^z = (-1)^s (-1)^t.

Book form
---------
T*sqrt(W) is not a unique split; as in Sec. 8.12 every sign-definite factor is
absorbed into the radical and only sign-indefinite factors stay outside.  The
overall sign is folded into the phase, giving (-1)^s or (-1)^(s+1).  Linear
factors are then printed in the book's own notation:

    a+b+c+k -> s+k       -a+b+c+k -> s-2a+k     2b+k -> 2b+k
    a-b+c+k -> s-2b+k     a+b-c+k -> s-2c+k     2c+k -> 2c+k

and the factor left outside the radical is rewritten as a polynomial in X over
Q[b,c] when that is possible (it is what the book does from Table 9.2 on).

Usage
-----
    python3 gen_9_11_6j_tables.py --d 1 --m=0 --n=0
    python3 gen_9_11_6j_tables.py --d 1/2                 # whole table
    python3 gen_9_11_6j_tables.py --d 1 --m=0 --n=0 --format sympy

Note the '=' in --m=-1: argparse would otherwise read a leading '-' as a flag.
d may be any positive half-integer; the book's d <= 4 is not a limit.
"""
import sys

from sympy import (Rational, symbols, RisingFactorial, FallingFactorial, factorial,
                   cancel, factor, factor_list, sqrt, latex, fraction, Integer,
                   lambdify, Poly, div, expand, simplify)

u, v, w = symbols('u v w')               # triangle deficits s-2a, s-2b, s-2c
a, b, c = symbols('a b c')               # the book's momenta
Xsym = symbols('X')

# X in terms of a, b, c
XABC = -a*(a + 1) + b*(b + 1) + c*(c + 1)

__all__ = ['sixj_algebraic', 'decompose', 'entry', 'render', 'u', 'v', 'w', 'a', 'b', 'c']


# --------------------------------------------------------------------------
# the Racah engine
# --------------------------------------------------------------------------
def _fr(x, j):
    """(x+j)! / x!  as an explicit rational function; j an explicit integer."""
    j = int(j)
    return RisingFactorial(x + 1, j) if j >= 0 else 1 / FallingFactorial(x, -j)


def _nfact(n):
    """n! for an explicit integer n, or None if n < 0 (term vanishes)."""
    n = int(n)
    return None if n < 0 else factorial(n)


def sixj_algebraic(d, m, n):
    """Return (T, W) with  {6j} = (-1)^s * T * sqrt(W),  rational in u, v, w.

    Returns None when the entry is identically zero.
    """
    d, m, n = Rational(d), Rational(m), Rational(n)
    if d < 0 or abs(m) > d or abs(n) > d or (d - m) % 1 or (d - n) % 1:
        return None
    s = u + v + w

    # ---- Racah sum, with u!, v!, w! and (s+1)! divided out ----
    T = Integer(0)
    t, thi = max(Integer(0), n + m), min(d + n, d + m)
    while t <= thi:
        denom, ok = Integer(1), True
        for j in (t, t - n - m, d + n - t, d + m - t):
            f = _nfact(j)
            if f is None:
                ok = False
                break
            denom *= f
        if ok:
            #  1/(v+t-d-m)! = 1/(v! * _fr(v, t-d-m)),  likewise w and u
            T += (Integer(-1)**int(t) * _fr(s + 1, t)
                  / (denom * _fr(v, t - d - m) * _fr(w, t - d - n) * _fr(u, n + m - t)))
        t += 1
    T = cancel(T)
    if T == 0:
        return None

    # ---- the square-rooted factor; all symbolic factorials have cancelled ----
    W = (_fr(v, n - m) * _fr(w, m - n) * _fr(u, n + m) / _fr(s + 1, n + m)      # Delta(abc)Delta(aef)
         * _nfact(d - m) * _nfact(d + m)                                        # Delta(dbf)
         * _fr(u + w, m - d) / _fr(u + w, d + m + 1)
         * _nfact(d + n) * _nfact(d - n)                                        # Delta(dec)
         * _fr(u + v, n - d) / _fr(u + v, d + n + 1))
    return T, cancel(W)


def _to_abc(expr):
    """u,v,w -> a,b,c."""
    return cancel(expr.subs({u: -a + b + c, v: a - b + c, w: a + b - c}))


def _definite(f):
    """+1 if f >= 0 on the physical domain, -1 if f <= 0, 0 if sign-indefinite.

    Sampled over a wide, deliberately asymmetric spread of the deficit lattice:
    factors like X only change sign when one momentum dominates the other two,
    so a narrow symmetric window silently misclassifies them as definite.
    """
    fn = lambdify((a, b, c), f, 'math')
    pos = neg = False
    for uu in (1, 2, 3, 5, 9, 17, 33, 65):
        for vv in (1, 2, 3, 5, 9, 17, 33, 65):
            for ww in (1, 2, 3, 5, 9, 17, 33, 65):
                av, bv, cv = (vv + ww) / 2, (uu + ww) / 2, (uu + vv) / 2
                try:
                    val = fn(av, bv, cv)
                except (ValueError, ZeroDivisionError):
                    continue
                if val > 1e-12:
                    pos = True
                elif val < -1e-12:
                    neg = True
    if pos and neg:
        return 0
    return -1 if neg else 1


# --------------------------------------------------------------------------
# book-form decomposition -- shared by render() and check_9_11.py so that the
# verified path is the shipped path
# --------------------------------------------------------------------------
def decompose(d, m, n):
    """Return (parity, outside, inside) with

        {6j} = (-1)**(s + parity) * outside * sqrt(inside),   s = a+b+c

    where sign-definite factors have been absorbed into `inside`.
    Returns None for an identically-zero entry.
    """
    res = sixj_algebraic(d, m, n)
    if res is None:
        return None
    T, W = (_to_abc(e) for e in res)

    const, facs = factor_list(T)
    inside, outside, parity = W * const**2, Integer(1), (1 if const < 0 else 0)
    for f, mult in facs:
        dsign = _definite(f)
        if dsign == 0:
            outside *= f**mult
        else:
            inside *= (f**2)**mult
            if dsign < 0 and mult % 2:
                parity ^= 1
    outside = factor(outside)
    if outside != 1 and _leading_negative(outside):       # keep X, not -X
        outside, parity = factor(-outside), parity ^ 1
    return parity, outside, factor(cancel(inside))


def _leading_negative(expr):
    """True if `expr` is naturally written negated -- judged by the leading
    coefficient of its X-polynomial when it has one, else by its content."""
    xs = _as_x_poly(expr)
    lead = xs[-1] if (xs is not None and len(xs) > 1) else expr
    if lead.is_number:
        return lead < 0
    return factor_list(lead)[0] < 0


def entry(d, m, n):
    """The entry as an exact sympy expression in a, b, c."""
    dec = decompose(d, m, n)
    if dec is None:
        return Integer(0)
    parity, outside, inside = dec
    return Integer(-1)**(a + b + c + parity) * outside * sqrt(inside)


# --------------------------------------------------------------------------
# rendering in the book's notation
# --------------------------------------------------------------------------
_LINEAR_NAMES = {                      # (coeff_a, coeff_b, coeff_c) -> name
    (1, 1, 1):  "s",
    (-1, 1, 1): "s-2a",
    (1, -1, 1): "s-2b",
    (1, 1, -1): "s-2c",
    (1, 0, 0):  "a",   (2, 0, 0): "2a",
    (0, 1, 0):  "b",   (0, 2, 0): "2b",
    (0, 0, 1):  "c",   (0, 0, 2): "2c",
}


def _name_linear(f):
    """The book's name for a linear factor (s+2, s-2a+1, 2b+1, ...), or None."""
    pf = Poly(f, a, b, c)
    if pf.total_degree() > 1:
        return None
    key = tuple(int(pf.coeff_monomial(x)) for x in (a, b, c))
    name = _LINEAR_NAMES.get(key)
    if name is None:
        return None
    k = int(pf.coeff_monomial(1))
    return name if k == 0 else f"{name}{'+' if k > 0 else '-'}{abs(k)}"


def _factors(expr):
    """factor_list(), with each factor sign-flipped when that makes it match a
    book name (s-2a+1 rather than -(a-b-c-1)); the sign moves into the const."""
    const, facs = factor_list(expr)
    out = []
    for f, mult in facs:
        if _name_linear(f) is None and _name_linear(-f) is not None:
            f, const = -f, const * Integer(-1)**mult
        out.append((f, mult))
    return const, out


def _split_square(r):
    """r = g**2 * rest, with g the largest rational square root of r."""
    from sympy import factorint
    r = Rational(r)
    g = Rational(1)
    for prime, e in factorint(r.p).items():
        g *= Rational(prime)**(e // 2)
    for prime, e in factorint(r.q).items():
        g /= Rational(prime)**(e // 2)
    return g, r / g**2


def _render_facs(const, facs):
    """A product of named factors, e.g. '(s+2)(s-2a+1)'."""
    parts = []
    if const != 1:
        parts.append(latex(const))
    for f, mult in facs:
        nm = _name_linear(f)
        txt = nm if nm else latex(f)
        parts.extend([txt if (nm and nm.isalnum()) else f"({txt})"] * int(mult))
    return "".join(parts) if parts else "1"


def _fmt_product(expr):
    """Factor an expression and print it with the book's names."""
    const, facs = _factors(expr)
    return _render_facs(const, facs)


def _as_x_poly(P):
    """P(a,b,c) as coefficients in X over Q[b,c], or None.

    X is quadratic in a, so divide repeatedly by X and require each remainder
    to be free of a.
    """
    coeffs, R = [], expand(P)
    while True:
        if not R.has(a):
            coeffs.append(R)
            return coeffs
        Q, rem = div(R, XABC, a)
        rem = expand(rem)
        if rem.has(a):
            return None
        coeffs.append(rem)
        R = expand(Q)


def _fmt_outside(outside):
    """The factor outside the radical, in X-notation when that is possible."""
    if outside == 1:
        return ""
    xs = _as_x_poly(outside)
    if xs is not None and len(xs) > 1:                 # genuinely involves X
        terms = []
        for i, co in reversed(list(enumerate(xs))):
            if co == 0:
                continue
            xp = "" if i == 0 else ("X" if i == 1 else f"X^{{{i}}}")
            if co in (1, -1) and xp:
                terms.append(("-" if co == -1 else "+", xp))
            else:
                neg = (co < 0) if co.is_number else (factor_list(co)[0] < 0)
                cl = _fmt_product(-co if neg else co)
                terms.append(("-" if neg else "+", f"{cl}{xp}" if xp else cl))
        out = ""
        for i, (sgn, txt) in enumerate(terms):
            out += (txt if sgn == "+" else f"-{txt}") if i == 0 else f" {sgn} {txt}"
        return out if len(terms) == 1 else r"\left\{%s\right\}" % out
    if outside.is_Add:
        return r"\left(%s\right)" % latex(outside)
    return _fmt_product(outside)


def render(d, m, n):
    """LaTeX for one entry, in the book's  (-1)^s (1/2)[ ... ]^{1/2}  style."""
    dec = decompose(d, m, n)
    if dec is None:
        return "0"
    parity, outside, inside = dec
    phase = r"(-1)^{s+1}" if parity else r"(-1)^{s}"

    num, den = fraction(inside)
    cn, fn = _factors(num)
    cd, fd = _factors(den)
    g, rest = _split_square(Rational(cn) / Rational(cd))   # pull the square part out
    pre = "" if g == 1 else latex(g)
    numstr = _render_facs(Rational(rest.p), fn)
    denstr = _render_facs(Rational(rest.q), fd)

    out = _fmt_outside(outside)
    if numstr == "1" and denstr != "1":                   # the  X/[...]^{1/2}  shape
        return r"%s%s\frac{%s}{\left[%s\right]^{1/2}}" % (phase, pre, out or "1", denstr)
    body = (r"\left[\frac{%s}{%s}\right]^{1/2}" % (numstr, denstr)
            if denstr != "1" else r"\left[%s\right]^{1/2}" % numstr)
    return f"{phase}{pre}{out}{body}"


def _fmt(x):
    return str(Rational(x))


def _lbl(base, k):
    k = Rational(k)
    if k == 0:
        return base
    return f"{base}{'+' if k > 0 else '-'}{_fmt(abs(k))}"


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(
        description="Generate algebraic 6j entries (VMK Sec. 9.11, Tables 9.1-9.8).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Use '=' for negative values:  --m=-1 --n=-1/2\n"
               "d may be any positive half-integer (the book stops at 4).")
    ap.add_argument('--d', required=True, help="the fixed momentum d, e.g. 1 or 3/2")
    ap.add_argument('--m', help="row: f = b + m")
    ap.add_argument('--n', help="column: e = c + n")
    ap.add_argument('--format', choices=['latex', 'sympy', 'both'], default='both')
    args = ap.parse_args(argv)

    d = Rational(args.d)
    if d <= 0 or (2*d) % 1:
        ap.error("d must be a positive half-integer")
    if (args.m is None) != (args.n is None):
        ap.error("give both --m and --n, or neither (to dump the table)")

    if args.m is not None:
        mm, nn = Rational(args.m), Rational(args.n)
        if abs(mm) > d or abs(nn) > d or (d - mm) % 1 or (d - nn) % 1:
            ap.error(f"m and n must run over {_fmt(d)}, {_fmt(d-1)}, ..., {_fmt(-d)}")
        print(f"{{a b c; {_fmt(d)} {_lbl('c', nn)} {_lbl('b', mm)}}}"
              f"   (e = {_lbl('c', nn)}, f = {_lbl('b', mm)})")
        if args.format in ('sympy', 'both'):
            print("  expr :", entry(d, mm, nn))
        if args.format in ('latex', 'both'):
            print("  latex:", render(d, mm, nn))
        return 0

    vals = [d - i for i in range(int(2*d) + 1)]
    print(f"Table for d = {_fmt(d)}   (rows f = b+m, columns e = c+n)\n")
    for nn in vals:
        for mm in vals:
            print(f"e = {_lbl('c', nn):<8} f = {_lbl('b', mm):<8} {render(d, mm, nn)}")
        print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
