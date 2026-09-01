#!/usr/bin/env python3
r"""Generate any algebraic 9j entry of the VMK Sec. 10.11 tables (10.1-10.12).

The tables give, with a, b, c left SYMBOLIC,

    { a+lambda  b+mu  c+nu }
    {    a        b     c  }
    {  alpha    beta  gamma}

for 0 <= alpha <= 3, 0 <= beta <= 2, gamma = 0 or 1 (and alpha >= beta; the
rest follow from the symmetries), with -alpha <= lambda <= alpha etc., in the
book's notation

    S = a+b+c,    Z = -c(c+1) + a(a+1) + b(b+1).

Method
------
Unlike the 6j of Sec. 9.11 this is not a single Racah sum.  It is expanded as

    9j = sum_x (-1)^(2x) (2x+1) {j1 j4 j7; j8 j9 x}{j2 j5 j8; j4 x j6}
                                {j3 j6 j9; x j1 j2}

and the triad (beta, a, x) forces x = a + xi with xi = -beta ... beta, an
explicit sum of at most 2*beta+1 terms.  Each of the three 6j symbols is then
carried by column permutations plus an upper/lower swap in two columns into
the Sec. 9.11 template {A B C; d C+n B+m} with d small:

    {a+lam, a, alpha; beta, gamma, a+xi} = {a+xi, beta, a; alpha, a+lam, gamma}
    {b+mu, b, beta; a, a+xi, c}          = {c, a, b; beta, b+mu, a+xi}
    {c+nu, c, gamma; a+xi, a+lam, b+mu}  = {b+mu, a+xi, c; gamma, c+nu, a+lam}

so the verified Sec. 9.11 engine supplies all three; only its deficits u,v,w
need remapping onto the global U = -a+b+c, V = a-b+c, W = a+b-c.

The phases collapse completely.  Summing the three (-1)^(A+B+C) with (-1)^(2x)
gives exponent 4a + 4xi + 2(a+b+c) + beta + mu, and 2a, 2xi and a+b+c are all
integers, so every term dies except an overall explicit (-1)^(beta+mu).

Each xi term is (rational)*sqrt(radicand).  They share a common radical: the
ratio of any two radicands is a perfect square rational function (checked, not
assumed), so it is pulled out and the entry collapses to

    9j = (-1)^(beta+mu) * (rational) * sqrt(common radicand).

Sign-definite factors are then absorbed into the radical exactly as in
Secs. 8.12 and 9.11, and what is left outside is rewritten as a polynomial in
Z over Q[a,b] when possible -- the {3Z - 2ab} forms the book uses.

Usage
-----
    python3 gen_10_11_9j_tables.py --alpha 1/2 --beta 1/2 --gamma 0 \
                                   --lam=1/2 --mu=1/2 --nu=0
    python3 gen_10_11_9j_tables.py --alpha 1/2 --beta 1/2 --gamma 0    # table
Use '=' for negative values (--lam=-1/2); argparse reads a bare '-' as a flag.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from sympy import (Rational, symbols, cancel, factor, factor_list, fraction, sqrt,
                   latex, Integer, lambdify, nsimplify, Poly, div, expand)

from gen_9_11_6j_tables import sixj_algebraic, u as _u9, v as _v9, w as _w9

a, b, c = symbols('a b c')
U, V, W = -a + b + c, a - b + c, a + b - c          # global triangle deficits
ZABC = -c*(c + 1) + a*(a + 1) + b*(b + 1)           # the book's Z

__all__ = ['ninej_terms', 'combine', 'decompose', 'entry', 'render', 'a', 'b', 'c']

_PT = (7.5, 6.5, 5.0)                               # a physical (a,b,c) sample point


# --------------------------------------------------------------------------
# the engine
# --------------------------------------------------------------------------
def _sixj(d, m, n, uu, vv, ww):
    """Sec. 9.11 engine mapped onto given deficits: {6j} = (-1)^(A+B+C) T sqrt(Wr)."""
    res = sixj_algebraic(d, m, n)
    if res is None:
        return None
    sub = {_u9: uu, _v9: vv, _w9: ww}
    return (cancel(res[0].subs(sub, simultaneous=True)),
            cancel(res[1].subs(sub, simultaneous=True)))


def ninej_terms(al, be, ga, lam, mu, nu):
    """[(coef, radicand)] with  9j = (-1)**(be+mu) * sum coef*sqrt(radicand)."""
    al, be, ga = Rational(al), Rational(be), Rational(ga)
    lam, mu, nu = Rational(lam), Rational(mu), Rational(nu)
    out, xi = [], -be
    while xi <= be:
        r1 = _sixj(al, ga - be, lam,     be - xi,     V + W + xi - be, xi + be)
        r2 = _sixj(be, xi,      mu,      W,           U,               V)
        r3 = _sixj(ga, lam - xi, nu,     V + xi - mu, U + mu - xi,     W + mu + xi)
        if r1 and r2 and r3:
            coef = cancel((V + W + 2*xi + 1) * r1[0] * r2[0] * r3[0])
            if coef != 0:
                out.append((coef, cancel(r1[1] * r2[1] * r3[1])))
        xi += 1
    return out


def _sqrt_poly(p):
    """(g, leftover, const) with p == const * g**2 * leftover."""
    const, facs = factor_list(p)
    g, rest = Integer(1), Integer(1)
    for fac, mult in facs:
        g *= fac**(mult // 2)
        if mult % 2:
            rest *= fac
    return g, rest, Rational(const)


def _sqrt_rational(r):
    """g with g**2 == r, or None when r is not a perfect square rational function."""
    num, den = fraction(cancel(r))
    gn, rn, cn = _sqrt_poly(num)
    gd, rd, cd = _sqrt_poly(den)
    if rn != 1 or rd != 1:
        return None
    cs = nsimplify(sqrt(abs(cn / cd)))
    if not cs.is_Rational:
        return None
    return cancel(gn / gd) * cs


def combine(terms):
    """sum coef*sqrt(rad)  ->  (coef_total, common_radicand).

    Raises ValueError if the radicands do not share a common radical, which
    would mean the entry is not of the form (rational)*sqrt(rational).
    """
    if not terms:
        return Integer(0), Integer(1)
    rad0 = terms[0][1]
    tot = Integer(0)
    for coef, rad in terms:
        g = _sqrt_rational(cancel(rad / rad0))
        if g is None:
            raise ValueError(f"radicand ratio is not a perfect square: {factor(cancel(rad/rad0))}")
        if float(lambdify((a, b, c), g, 'math')(*_PT)) < 0:   # positive branch
            g = -g
        tot += coef * g
    return cancel(tot), rad0


def _definite(f):
    """+1 if f >= 0 on the physical domain, -1 if f <= 0, 0 if sign-indefinite."""
    fn = lambdify((a, b, c), f, 'math')
    pos = neg = False
    for uu in (1, 2, 3, 5, 9, 17, 33, 65):
        for vv in (1, 2, 3, 5, 9, 17, 33, 65):
            for ww in (1, 2, 3, 5, 9, 17, 33, 65):
                av, bv, cv = (vv+ww)/2, (uu+ww)/2, (uu+vv)/2
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


def _as_z_poly(P):
    """P(a,b,c) as coefficients in Z over Q[a,b], or None (Z is quadratic in c)."""
    coeffs, R = [], expand(P)
    while True:
        if not R.has(c):
            coeffs.append(R)
            return coeffs
        Q, rem = div(R, ZABC, c)
        rem = expand(rem)
        if rem.has(c):
            return None
        coeffs.append(rem)
        R = expand(Q)


def _leading_negative(expr):
    zs = _as_z_poly(expr)
    lead = zs[-1] if (zs is not None and len(zs) > 1) else expr
    return lead < 0 if lead.is_number else factor_list(lead)[0] < 0


def decompose(al, be, ga, lam, mu, nu):
    """(sign, outside, inside) with  9j = sign * outside * sqrt(inside).

    Sign-definite factors are absorbed into the radical; the (-1)**(be+mu)
    phase is folded into `sign`.  Returns None for an identically-zero entry.
    """
    coef, rad = combine(ninej_terms(al, be, ga, lam, mu, nu))
    if coef == 0:
        return None
    sign = (-1)**int(Rational(be) + Rational(mu))

    const, facs = factor_list(coef)
    inside, outside = rad * const**2, Integer(1)
    if const < 0:
        sign = -sign
    for f, mult in facs:
        d = _definite(f)
        if d == 0:
            outside *= f**mult
        else:
            inside *= (f**2)**mult
            if d < 0 and mult % 2:
                sign = -sign
    outside = factor(outside)
    if outside != 1 and _leading_negative(outside):          # keep Z, not -Z
        outside, sign = factor(-outside), -sign
    return sign, outside, factor(cancel(inside))


def entry(al, be, ga, lam, mu, nu):
    """The entry as an exact sympy expression in a, b, c."""
    dec = decompose(al, be, ga, lam, mu, nu)
    if dec is None:
        return Integer(0)
    sign, outside, inside = dec
    return sign * outside * sqrt(inside)


# --------------------------------------------------------------------------
# rendering in the book's notation
# --------------------------------------------------------------------------
_NAMES = {                             # (coeff_a, coeff_b, coeff_c) -> name
    (1, 1, 1):  "S",
    (-1, 1, 1): "S-2a",
    (1, -1, 1): "S-2b",
    (1, 1, -1): "S-2c",
    (1, 0, 0):  "a",   (2, 0, 0): "2a",
    (0, 1, 0):  "b",   (0, 2, 0): "2b",
    (0, 0, 1):  "c",   (0, 0, 2): "2c",
}


def _name_linear(f):
    """The book's name for a linear factor (S+2, S-2c+1, 2a+1, ...), or None."""
    pf = Poly(f, a, b, c)
    if pf.total_degree() > 1:
        return None
    name = _NAMES.get(tuple(int(pf.coeff_monomial(x)) for x in (a, b, c)))
    if name is None:
        return None
    k = int(pf.coeff_monomial(1))
    return name if k == 0 else f"{name}{'+' if k > 0 else '-'}{abs(k)}"


def _factors(expr):
    """factor_list(), sign-flipping factors so they match a book name."""
    const, facs = factor_list(expr)
    out = []
    for f, mult in facs:
        if _name_linear(f) is None and _name_linear(-f) is not None:
            f, const = -f, const * Integer(-1)**mult
        out.append((f, mult))
    return const, out


def _render_facs(const, facs):
    parts = []
    if const != 1:
        parts.append(latex(const))
    for f, mult in facs:
        nm = _name_linear(f)
        txt = nm if nm else latex(f)
        parts.extend([txt if (nm and nm.isalnum()) else f"({txt})"] * int(mult))
    return "".join(parts) if parts else "1"


def _fmt_product(expr):
    return _render_facs(*_factors(expr))


def _fmt_outside(outside):
    """The factor outside the radical, in Z-notation when that is possible."""
    if outside == 1:
        return ""
    zs = _as_z_poly(outside)
    if zs is not None and len(zs) > 1:                  # genuinely involves Z
        terms = []
        for i, co in reversed(list(enumerate(zs))):
            if co == 0:
                continue
            zp = "" if i == 0 else ("Z" if i == 1 else f"Z^{{{i}}}")
            if co in (1, -1) and zp:
                terms.append(("-" if co == -1 else "+", zp))
            else:
                neg = (co < 0) if co.is_number else (factor_list(co)[0] < 0)
                cl = _fmt_product(-co if neg else co)
                terms.append(("-" if neg else "+", f"{cl}{zp}" if zp else cl))
        out = ""
        for i, (sgn, txt) in enumerate(terms):
            out += (txt if sgn == "+" else f"-{txt}") if i == 0 else f" {sgn} {txt}"
        return out if len(terms) == 1 else r"\left\{%s\right\}" % out
    if outside.is_Add:
        return r"\left(%s\right)" % latex(outside)
    return _fmt_product(outside)


def render(al, be, ga, lam, mu, nu):
    """LaTeX for one entry, in the book's  +/-{3Z-2ab}[ ... ]^{1/2}  style."""
    dec = decompose(al, be, ga, lam, mu, nu)
    if dec is None:
        return "0"
    sign, outside, inside = dec

    num, den = fraction(inside)
    cn, fn = _factors(num)
    cd, fd = _factors(den)
    from gen_9_11_6j_tables import _split_square
    g, rest = _split_square(Rational(cn) / Rational(cd))
    pre = ("-" if sign < 0 else "") + ("" if g == 1 else latex(g))
    numstr = _render_facs(Rational(rest.p), fn)
    denstr = _render_facs(Rational(rest.q), fd)
    out = _fmt_outside(outside)

    if numstr == "1" and denstr != "1":
        return r"%s\frac{%s}{\left[%s\right]^{1/2}}" % (pre, out or "1", denstr)
    body = (r"\left[\frac{%s}{%s}\right]^{1/2}" % (numstr, denstr)
            if denstr != "1" else r"\left[%s\right]^{1/2}" % numstr)
    return f"{pre}{out}{body}"


# --------------------------------------------------------------------------
def _fmt(x):
    return str(Rational(x))


def _lbl(base, k):
    k = Rational(k)
    return base if k == 0 else f"{base}{'+' if k > 0 else '-'}{_fmt(abs(k))}"


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(
        description="Generate algebraic 9j entries (VMK Sec. 10.11, Tables 10.1-10.12).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Use '=' for negative values:  --lam=-1/2 --nu=-1")
    ap.add_argument('--alpha', required=True)
    ap.add_argument('--beta', required=True)
    ap.add_argument('--gamma', required=True)
    ap.add_argument('--lam')
    ap.add_argument('--mu')
    ap.add_argument('--nu')
    ap.add_argument('--format', choices=['latex', 'sympy', 'both'], default='both')
    args = ap.parse_args(argv)

    al, be, ga = (Rational(x) for x in (args.alpha, args.beta, args.gamma))
    given = [args.lam, args.mu, args.nu]
    if any(x is None for x in given) and any(x is not None for x in given):
        ap.error("give all of --lam --mu --nu, or none (to dump the table)")

    if args.lam is not None:
        lam, mu, nu = (Rational(x) for x in given)
        for nm, val, lim in (("lam", lam, al), ("mu", mu, be), ("nu", nu, ga)):
            if abs(val) > lim or (lim - val) % 1:
                ap.error(f"{nm} must run over {_fmt(lim)}, ..., {_fmt(-lim)}")
        print(f"{{{_lbl('a',lam)} {_lbl('b',mu)} {_lbl('c',nu)}; a b c; "
              f"{_fmt(al)} {_fmt(be)} {_fmt(ga)}}}")
        if args.format in ('sympy', 'both'):
            print("  expr :", entry(al, be, ga, lam, mu, nu))
        if args.format in ('latex', 'both'):
            print("  latex:", render(al, be, ga, lam, mu, nu))
        return 0

    print(f"Table for alpha={_fmt(al)}, beta={_fmt(be)}, gamma={_fmt(ga)}\n")
    for lam in [al - i for i in range(int(2*al) + 1)]:
        for mu in [be - i for i in range(int(2*be) + 1)]:
            for nu in [ga - i for i in range(int(2*ga) + 1)]:
                print(f"lam={_fmt(lam):<5} mu={_fmt(mu):<5} nu={_fmt(nu):<5} "
                      f"{render(al, be, ga, lam, mu, nu)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
