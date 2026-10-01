#!/usr/bin/env python3
"""
Usage:
    gen_kernel_hydro.py [--check] [path/to/kernel_hydro.h]

Regenerates the kernel tables of src/kernel_hydro.h, i.e. the region between
the "/* clang-format off */" and "/* clang-format on */" guards, from the exact
definitions of the kernels. Everything outside that region is left untouched.
Without a path the header is looked up relative to this script. With --check
nothing is written and the exit status is 1 if the header is not up to date.

The kernels are those of table 1 of Dehnen & Aly, MNRAS 425, 1068 (2012):
W(r, h) = C / H^d f(r / H) with H = gamma h the compact support and
gamma = H / h = 1 / (2 sigma) fixed by the second moment sigma^2 of the kernel
(h = 2 sigma), so that a given resolution_eta gives the same number of
neighbours for all kernels. gamma is computed exactly, written as a double
literal and cast to float by the compiler; all tables are built for that
FLOAT value of gamma, so that the support is identical in float and double.

For each kernel and dimension, [0, gamma) is split into kernel_poly_ivals
uniform sub-intervals and the kernel is stored on each as a polynomial in
t = u - origin, with the normalisation folded in (kernel_poly_coeffs, highest
degree first; a last all-zero row is selected for u >= gamma). The origin is
one of the two ends of the sub-interval, the one giving the smaller
condition number of Horner's scheme for W and dW/du: this is u = 0 for the
first sub-interval (W(0) is then a table entry and dW/du has no cancellation
at the centre) and the right end for all others, so the last origin is the
support edge and W(gamma) = 0 exactly. Both ends being floats, t is exact in
float arithmetic (Sterbenz). For the splines the sub-intervals coincide with
or subdivide the natural branches. See the header for the choice of the
number of sub-intervals.

Every constant is computed exactly (sympy) and rounded ONCE to float or to
double (mpmath, round to nearest even). kernel_poly_ivals_over_gamma is
rounded down until no u < gamma maps to the zero row.

Requires numpy, sympy and mpmath.

This file is part of SWIFT.
Copyright (c) 2026 Matthieu Schaller (schaller@strw.leidenuniv.nl)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import os
import sys

import mpmath as mp
import numpy as np
import sympy as sp

mp.mp.dps = 60
q, u, t = sp.symbols("q u t")
R = sp.Rational
PI = sp.pi

GUARD_OPEN = "/* clang-format off */\n"
GUARD_CLOSE = "/* clang-format on */\n"
RULE = "/* " + "-" * 73 + " */"

# Volume element of the unit sphere surface, per dimension
VOL = {1: 2, 2: 2 * PI, 3: 4 * PI}


def spline_pieces(terms, breaks):
    """Polynomials of a B-spline sum_j c_j (a_j - q)_+^n on each interval
    between consecutive breaks (the term j is active where q < a_j)."""
    pieces = []
    for i in range(len(breaks) - 1):
        hi = breaks[i + 1]
        pieces.append(sp.expand(sum(c * (a - q) ** n for c, a, n in terms if a >= hi)))
    return pieces


def wendland(d_low, d_high):
    """Wendland kernel: (1 - q)^l p(q); the 1D form differs from 2D/3D."""

    def pieces(d):
        return [sp.expand(d_low if d == 1 else d_high)]

    return pieces


# Kernel definitions: shape f(q) on [0, 1] as polynomial pieces between the
# breaks, normalisation C per dimension (exact, and SWIFT's spelling of it),
# dimensions in which the kernel is defined, number of sub-intervals.
KERNELS = [
    dict(
        key="CUBIC_SPLINE_KERNEL",
        name="Cubic spline (M4)",
        breaks=[0, R(1, 2), 1],
        pieces=lambda d: spline_pieces([(1, 1, 3), (-4, R(1, 2), 3)], [0, R(1, 2), 1]),
        C={3: 16 / PI, 2: R(80, 7) / PI, 1: R(8, 3)},
        C_expr={3: "16. * M_1_PI", 2: "80. * M_1_PI / 7.", 1: "8. / 3."},
        dims=(3, 2, 1),
        nsub=2,
    ),
    dict(
        key="QUARTIC_SPLINE_KERNEL",
        name="Quartic spline (M5)",
        breaks=[0, R(1, 5), R(3, 5), 1],
        pieces=lambda d: spline_pieces(
            [(1, 1, 4), (-5, R(3, 5), 4), (10, R(1, 5), 4)],
            [0, R(1, 5), R(3, 5), 1],
        ),
        C={3: R(15625, 512) / PI, 2: R(46875, 2398) / PI, 1: R(3125, 768)},
        C_expr={
            3: "15625. * M_1_PI / 512.",
            2: "46875. * M_1_PI / 2398.",
            1: "3125. / 768.",
        },
        dims=(3, 2, 1),
        nsub=5,
    ),
    dict(
        key="QUINTIC_SPLINE_KERNEL",
        name="Quintic spline (M6)",
        breaks=[0, R(1, 3), R(2, 3), 1],
        pieces=lambda d: spline_pieces(
            [(1, 1, 5), (-6, R(2, 3), 5), (15, R(1, 3), 5)],
            [0, R(1, 3), R(2, 3), 1],
        ),
        C={3: R(2187, 40) / PI, 2: R(15309, 478) / PI, 1: R(243, 40)},
        C_expr={
            3: "2187. * M_1_PI / 40.",
            2: "15309. * M_1_PI / 478.",
            1: "243. / 40.",
        },
        dims=(3, 2, 1),
        nsub=3,
    ),
    dict(
        key="WENDLAND_C2_KERNEL",
        name="Wendland C2",
        breaks=[0, 1],
        pieces=wendland((1 - q) ** 3 * (1 + 3 * q), (1 - q) ** 4 * (1 + 4 * q)),
        C={3: R(21, 2) / PI, 2: 7 / PI, 1: R(5, 4)},
        C_expr={3: "21. * M_1_PI / 2.", 2: "7. * M_1_PI", 1: "5. / 4."},
        dims=(3, 2, 1),
        nsub=4,
    ),
    dict(
        key="WENDLAND_C4_KERNEL",
        name="Wendland C4",
        breaks=[0, 1],
        pieces=wendland(None, (1 - q) ** 6 * (1 + 6 * q + R(35, 3) * q**2)),
        C={3: R(495, 32) / PI, 2: 9 / PI},
        C_expr={3: "495. * M_1_PI / 32.", 2: "9. * M_1_PI"},
        dims=(3, 2),
        nsub=4,
    ),
    dict(
        key="WENDLAND_C6_KERNEL",
        name="Wendland C6",
        breaks=[0, 1],
        pieces=wendland(None, (1 - q) ** 8 * (1 + 8 * q + 25 * q**2 + 32 * q**3)),
        C={3: R(1365, 64) / PI, 2: R(78, 7) / PI},
        C_expr={3: "1365. * M_1_PI / 64.", 2: "78. * M_1_PI / 7."},
        dims=(3, 2),
        nsub=4,
    ),
]


# --- Rounding and C literals -------------------------------------------------


def to_mpf(x):
    """Exact sympy expression (or number) to a 60-digit mpf."""
    return mp.mpf(sp.N(x, 60)) if isinstance(x, sp.Basic) else mp.mpf(x)


def f32(x):
    """Round once (nearest even) to float."""
    with mp.workprec(24):
        y = +to_mpf(x)
    return np.float32(float(y))


def f64(x):
    """Round once (nearest even) to double."""
    with mp.workprec(53):
        y = +to_mpf(x)
    return float(y)


def exact(x):
    """Exact rational value of a float or double."""
    return sp.Rational(*float(x).as_integer_ratio())


def flit(x):
    """Float literal that round-trips (e.g. 0.418429196f)."""
    s = "%.9g" % float(x)
    if "e" not in s and "." not in s:
        s += "."
    return s + "f"


def dlit(x):
    """Double literal that round-trips (e.g. 0.4184291920938659)."""
    s = repr(float(x))
    if "e" not in s and "." not in s:
        s += "."
    return s


# --- Tables -------------------------------------------------------------------


def condition_number(coeffs, tt):
    """max over the points tt of sum |c_k| |t|^k / |sum c_k t^k| for the
    polynomial and for its derivative (mpf arithmetic)."""
    c = [to_mpf(x) for x in coeffs]
    dc = [k * x for k, x in enumerate(c)][1:]
    worst = mp.mpf(0)
    for x in tt:
        for cc in (c, dc):
            val = abs(sum(ck * x**k for k, ck in enumerate(cc)))
            mag = sum(abs(ck) * abs(x) ** k for k, ck in enumerate(cc))
            worst = max(worst, mag / val if val else mp.inf)
    return worst


def build(kernel, d):
    """All numbers of the tables of one kernel in dimension d."""
    pieces = kernel["pieces"](d)
    breaks = kernel["breaks"]
    C = kernel["C"][d]
    nsub = kernel["nsub"]

    # gamma = H / h = 1 / (2 sigma), with sigma^2 = <r^2> / d for H = 1
    m2 = sum(
        sp.integrate(
            C * p * q**2 * VOL[d] * q ** (d - 1), (q, breaks[i], breaks[i + 1])
        )
        for i, p in enumerate(pieces)
    )
    gamma2 = sp.simplify(d / (4 * m2))
    assert gamma2.is_Rational, gamma2
    gamma_d = f64(sp.sqrt(gamma2))  # the double literal of the header
    gamma_f = np.float32(gamma_d)  # (float)(double literal), as the compiler
    G = exact(gamma_f)  # the tables are built for this exact value

    degree = max(sp.Poly(p, q).degree() for p in pieces)

    def piece_at(x):
        return pieces[max(i for i in range(len(breaks) - 1) if breaks[i] <= x)]

    origins, rows_f, rows_d = [], [], []
    for j in range(nsub):
        a, b = G * R(j, nsub), G * R(j + 1, nsub)
        W_u = sp.expand(C / G**d * piece_at((a + b) / (2 * G)).subs(q, u / G))
        best = None
        for end in (a, b):
            origin = np.float32(float(end))  # a and b are exactly floats
            P = exact(origin)
            poly = sp.Poly(sp.expand(W_u.subs(u, P + t)), t)
            coeffs = [poly.coeff_monomial(t**k) for k in range(degree + 1)]
            tt = [to_mpf(a + (b - a) * R(s, 60) - P) for s in range(1, 60)]
            cond = condition_number(coeffs, tt)
            if best is None or cond < best[0]:
                best = (cond, origin, coeffs)
        _, origin, coeffs = best
        origins.append(origin)
        rows_f.append([f32(c) for c in reversed(coeffs)])
        rows_d.append([f64(c) for c in reversed(coeffs)])
    origins.append(gamma_f)
    rows_f.append([np.float32(0)] * (degree + 1))
    rows_d.append([0.0] * (degree + 1))

    # Sanity: the first origin is 0 so that W(0) is the last entry of row 0
    assert origins[0] == 0, "first sub-interval not expanded about u = 0"
    root = f32(C / G**d * pieces[0].subs(q, 0))
    assert rows_f[0][-1] == root

    # nsub / gamma rounded down until every u < gamma maps to a row < nsub
    S = np.float32(float(nsub / G))
    below = np.nextafter(gamma_f, np.float32(0))
    while np.float32(below * S) >= nsub:
        S = np.nextafter(S, np.float32(0))
    S_d = float(nsub / G)
    while np.nextafter(float(gamma_f), 0.0) * S_d >= nsub:
        S_d = np.nextafter(S_d, 0.0)

    return dict(
        gamma2=gamma2,
        gamma_d=gamma_d,
        degree=degree,
        nsub=nsub,
        S=S,
        S_d=S_d,
        root=root,
        origins=origins,
        rows_f=rows_f,
        rows_d=rows_d,
    )


def dimension_block(kernel, d, first):
    """Lines of the #if HYDRO_DIMENSION_dD block of one kernel."""
    r = build(kernel, d)
    out = [("#if" if first else "#elif") + f" defined(HYDRO_DIMENSION_{d}D)"]
    out.append(
        f"#define kernel_gamma ((float)({dlit(r['gamma_d'])})) /* sqrt({r['gamma2']}) */"
    )
    out.append(f"#define kernel_constant ((float)({kernel['C_expr'][d]}))")
    out.append(f"#define kernel_poly_degree {r['degree']}")
    out.append(f"#define kernel_poly_ivals {r['nsub']}")
    out.append(f"#define kernel_poly_ivals_over_gamma {flit(r['S'])}")
    out.append(f"#define kernel_poly_ivals_over_gamma_d {dlit(r['S_d'])}")
    out.append(f"#define kernel_poly_root {flit(r['root'])}")
    out.append(
        "static const float kernel_poly_origin[kernel_poly_ivals + 1] = {\n    "
        + ", ".join(flit(x) for x in r["origins"])
        + "};"
    )
    out.append(
        "static const double kernel_poly_origin_d[kernel_poly_ivals + 1] = {\n    "
        + ", ".join(dlit(x) for x in r["origins"])
        + "};"
    )
    out.append(
        "static const float\n"
        "    kernel_poly_coeffs[(kernel_poly_degree + 1) * (kernel_poly_ivals + 1)] = {"
    )
    for row in r["rows_f"]:
        out.append("        " + ", ".join(flit(c) for c in row) + ",")
    out.append("};")
    out.append(
        "static const double\n"
        "    kernel_poly_coeffs_d[(kernel_poly_degree + 1) * (kernel_poly_ivals + 1)] = {"
    )
    for row in r["rows_d"]:
        out.append("        " + ", ".join(dlit(c) for c in row) + ",")
    out.append("};")
    return out


def generate():
    """The text between the clang-format guards."""
    body = []
    for i, kernel in enumerate(KERNELS):
        body.append(RULE)
        body.append(("#if" if i == 0 else "#elif") + f" defined({kernel['key']})\n")
        body.append(f'#define kernel_name "{kernel["name"]}"')
        for j, d in enumerate(kernel["dims"]):
            body += dimension_block(kernel, d, j == 0)
        if 1 not in kernel["dims"]:
            body.append(
                "#elif defined(HYDRO_DIMENSION_1D)\n"
                f'#error "{kernel["name"]} kernel not defined in 1D."'
            )
        body.append("#endif\n")
    body.append(RULE)
    body.append(
        '#else\n\n#error "A kernel function must be chosen at configure time !!"\n\n'
        + RULE
        + "\n#endif"
    )
    return "\n".join(body) + "\n"


def main(argv):
    check = "--check" in argv
    args = [a for a in argv if a != "--check"]
    if len(args) > 1:
        sys.exit(__doc__.split("\n\n")[0])
    path = (
        args[0]
        if args
        else os.path.join(os.path.dirname(__file__), "..", "src", "kernel_hydro.h")
    )

    header = open(path).read()
    start = header.index(GUARD_OPEN) + len(GUARD_OPEN)
    end = header.index(GUARD_CLOSE)
    new = header[:start] + generate() + header[end:]

    if new == header:
        print(f"{path}: tables are up to date", file=sys.stderr)
        return 0
    if check:
        print(f"{path}: tables are NOT up to date", file=sys.stderr)
        return 1
    open(path, "w").write(new)
    print(f"{path}: tables regenerated", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
