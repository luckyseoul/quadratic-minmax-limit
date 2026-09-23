"""Exact algebra regressions for the finite interior-gap proof, not a census."""

import hashlib
import json
from pathlib import Path

import sympy as sp


checks = []


def identity(name, actual, expected):
    assert sp.simplify(actual - expected) == 0, name
    checks.append(name)


b, m, k = sp.symbols("b m k", integer=True, positive=True)
H, L, p, alpha, n = sp.symbols("H L p alpha n", positive=True)
e, f, eta = sp.symbols("e f eta", real=True)
identity("Rademacher fourth-moment coefficients", b + 3 * b * (b - 1), 3 * b**2 - 2 * b)
identity("even half-split", 8 * m**3, (2 * m) ** 3)
identity("odd half-split remainder", 8 * (m + 1) ** 2 * m - (2 * m + 1) ** 3, 4 * m**2 + 2 * m - 1)
identity("odd remainder positive decomposition", 4 * m**2 + 2 * m - 1, 4 * m * (m - 1) + 6 * m - 1)
t = L / (2 * H)
identity("admissible quadratic optimizer value", t * L - t**2 * H, L**2 / (4 * H))
identity("finite product constant", sp.Rational(1, 8) / 12, sp.Rational(1, 96))

# A generic zero-diagonal quadratic: all coefficients and coordinates symbolic.
a01, a02, a12 = sp.symbols("a01 a02 a12", real=True)
z = sp.symbols("z0:3", real=True)
v = sp.symbols("v0:3", real=True)
epsilon = sp.symbols("epsilon", real=True)


def q(x):
    return a01 * x[0] * x[1] + a02 * x[0] * x[2] + a12 * x[1] * x[2]


identity(
    "symmetric quadratic perturbation",
    (q([z[i] + epsilon * v[i] for i in range(3)]) + q([z[i] - epsilon * v[i] for i in range(3)])) / 2,
    q(z) + epsilon**2 * q(v),
)
identity("layer-cake integral", sp.integrate(epsilon ** sp.Rational(-2, 3), (epsilon, 0, 1)), 3)
identity("cubed softness constant", 3**3 * 96, 2592)
F = alpha * n ** sp.Rational(3, 2)
identity("normalized fixed-update deficit", (2 * p) ** 2 * k**3 / (96 * F * sp.sqrt(n)), p**2 * k**3 / (24 * alpha * n**2))
identity(
    "phase-averaged update rearrangement",
    alpha - p**2 * eta / (24 * alpha) - ((1 - p) ** 2 * e + p * (1 - p) * f - p**2 * alpha),
    (1 + p**2) * alpha - (1 - p) ** 2 * e - p * (1 - p) * f - p**2 * eta / (24 * alpha),
)

# One retained C5 boundary point. Phi(C5)=4 is prior evidence, not recomputed.
# It checks why an isolated zero-field coordinate must remain admissible.
C5 = [[0 if i == j else (1 if (i - j) % 5 in (1, 4) else -1) for j in range(5)] for i in range(5)]
x = [-1, -1, 1, 1, 1]
soft = [-1, -1, 0, 1, 1]


def energy(matrix, point):
    return sum(matrix[i][j] * point[i] * point[j] for i in range(len(point)) for j in range(i + 1, len(point)))


identity("retained C5 zero field", sum(C5[2][j] * x[j] for j in range(5)), 0)
identity("retained C5 fractional face value", energy(C5, soft), energy(C5, x))
identity("retained C5 value normalization", energy(C5, soft), 4)
assert sum(1 - abs(value) for value in soft) == 1
checks.append("single fractional coordinate allowed")

print(json.dumps({
    "classification": "EXACT_ALGEBRA_REGRESSION_NOT_ALL_ORDERS_PROOF",
    "sympy_version": sp.__version__,
    "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "passed": len(checks),
    "checks": checks,
    "original_problem_status": "OPEN",
}, indent=2))
