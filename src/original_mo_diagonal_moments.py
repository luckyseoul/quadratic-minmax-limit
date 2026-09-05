"""Exact mixed moments for the original-MO cross/diagonal decomposition.

Scalars are integers, Fractions, or ``(real, imaginary)`` pairs of those.
No floating-point or built-in complex arithmetic is used.  The identities
hold for fixed blocks and uniform independent fourth roots of unity; they
do not supply the uniform phase bound needed by the multiplier-two route.
"""

from __future__ import annotations

from collections.abc import Sequence
from fractions import Fraction
from math import comb

Rational = int | Fraction
ExactScalar = Rational | tuple[Rational, Rational]
Gaussian = tuple[Fraction, Fraction]
Matrix = list[list[Gaussian]]

original_mo_limit_closed = False
multiplier_two_closed = False

_ZERO = (Fraction(0), Fraction(0))
_PHASES = {(1, 0), (0, 1), (-1, 0), (0, -1)}


def _gaussian(value: ExactScalar) -> Gaussian:
    if isinstance(value, (int, Fraction)) and not isinstance(value, bool):
        return Fraction(value), Fraction(0)
    if (
        isinstance(value, tuple)
        and len(value) == 2
        and all(isinstance(part, (int, Fraction)) and not isinstance(part, bool)
                for part in value)
    ):
        return Fraction(value[0]), Fraction(value[1])
    raise TypeError("use integers, Fractions, or exact (real, imaginary) pairs")


def _matrix(values: Sequence[Sequence[ExactScalar]]) -> Matrix:
    if not values or not values[0]:
        raise ValueError("matrices must have positive dimensions")
    width = len(values[0])
    if any(len(row) != width for row in values):
        raise ValueError("matrix rows must have equal lengths")
    return [[_gaussian(value) for value in row] for row in values]


def _add(a: Gaussian, b: Gaussian) -> Gaussian:
    return a[0] + b[0], a[1] + b[1]


def _mul(a: Gaussian, b: Gaussian) -> Gaussian:
    return a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]


def _conj(a: Gaussian) -> Gaussian:
    return a[0], -a[1]


def _adjoint(a: Matrix) -> Matrix:
    return [[_conj(a[i][j]) for i in range(len(a))] for j in range(len(a[0]))]


def _matmul(a: Matrix, b: Matrix) -> Matrix:
    if len(a[0]) != len(b):
        raise ValueError("incompatible matrix dimensions")
    out = [[_ZERO for _ in b[0]] for _ in a]
    for i in range(len(a)):
        for j in range(len(b[0])):
            for k in range(len(b)):
                out[i][j] = _add(out[i][j], _mul(a[i][k], b[k][j]))
    return out


def _trace(a: Matrix) -> Gaussian:
    total = _ZERO
    for i in range(len(a)):
        total = _add(total, a[i][i])
    return total


def _off_diagonal(a: Matrix) -> Matrix:
    return [[_ZERO if i == j else value for j, value in enumerate(row)]
            for i, row in enumerate(a)]


def _real(a: Gaussian) -> Fraction:
    if a[1] != 0:
        raise ArithmeticError("a theoretically real mixed moment was not real")
    return a[0]


def _zero_diagonal_hermitian(a: Matrix) -> None:
    if len(a) != len(a[0]):
        raise ValueError("diagonal-energy blocks must be square")
    if any(a[i][i] != _ZERO for i in range(len(a))):
        raise ValueError("diagonal-energy blocks must have zero diagonal")
    if a != _adjoint(a):
        raise ValueError("diagonal-energy blocks must be Hermitian")


def _energy(a: Matrix, v: Sequence[Gaussian]) -> Fraction:
    total = _ZERO
    for i in range(len(v)):
        for j in range(len(v)):
            total = _add(total, _mul(_mul(_conj(v[i]), a[i][j]), v[j]))
    return _real(total)


def lee_weight_two_coefficient(r: int) -> Fraction:
    """Return a_r, the balanced Lee-weight-two coefficient of h_r.

    Here h_r(w)=rho((1+i) sum_j w_j), rho(t)=max(|Re t|, |Im t|),
    and balanced means that the character exponents sum to zero modulo 4.
    """

    if not isinstance(r, int) or r < 2:
        raise ValueError("r must be an integer at least two")
    return Fraction(comb(2 * r - 2, r - 1), 2 ** (2 * r - 2))


def lee_weight_four_coefficient(r: int) -> Fraction:
    """Return -a_r/(2r-3), including repeated-index Lee-weight-four terms.

    All balanced characters of Lee weight four have this coefficient;
    the unbalanced characters of that weight have coefficient zero.
    """

    return -lee_weight_two_coefficient(r) / (2 * r - 3)


def row_energy_covariance(
    q: Sequence[Sequence[ExactScalar]],
    b: Sequence[Sequence[ExactScalar]],
) -> Fraction:
    """Return Cov(F_Q(w), (w* B w)^2) exactly (mixed-moment identity D8).

    Q is ell by r, r>=2, with fourth-root entries; B is zero-diagonal
    Hermitian of order r.  For uniform w in mu_4^r, define
    F_Q(w)=sum_a rho((1+i) sum_j Q_aj w_j).  With T=Q*Q,
    v_a=conj(Q_a), a=a_r and d=a/(2r-3), the returned expression is

      2(a+d) tr[B^2(T-ell I)]
        - d (sum_a (v_a* B v_a)^2 - ell tr B^2).

    The row vectors v_a are conjugated: replacing them by Q_a generally
    changes the actual diagonal-energy contribution.
    """

    qm, bm = _matrix(q), _matrix(b)
    ell, r = len(qm), len(qm[0])
    a = lee_weight_two_coefficient(r)
    if any(value not in _PHASES for row in qm for value in row):
        raise ValueError("Q entries must be fourth roots of unity")
    _zero_diagonal_hermitian(bm)
    if len(bm) != r:
        raise ValueError("B order must equal the number of Q columns")
    d = -lee_weight_four_coefficient(r)
    b2 = _matmul(bm, bm)
    trace_b2 = _real(_trace(b2))
    # diag(Q*Q)=ell I because every entry of Q has unit modulus.
    gram_defect = _off_diagonal(_matmul(_adjoint(qm), qm))
    gram_term = _real(_trace(_matmul(b2, gram_defect)))
    row_squares = sum(_energy(bm, [_conj(x) for x in row]) ** 2 for row in qm)
    return 2 * (a + d) * gram_term - d * (row_squares - ell * trace_b2)


def joint_diagonal_covariance(
    left: Sequence[Sequence[ExactScalar]],
    cross: Sequence[Sequence[ExactScalar]],
    right: Sequence[Sequence[ExactScalar]],
) -> Fraction:
    """Return Cov(|z* G w|^2, (z* L z + w* B w)^2) exactly.

    L and B are zero-diagonal Hermitian blocks; G is any compatible
    rectangular Gaussian-rational matrix.  z and w are independent
    uniform fourth-root vectors.  Writing off(K)=K-diag(K), return

      2 tr[L^2 off(GG*)] + 2 tr[B^2 off(G*G)] + 2 tr[L G B G*].

    For cross entries of squared modulus two, the off-diagonal Gram
    terms equal GG*-2r I and G*G-2ell I, giving identity D17.
    """

    lm, gm, bm = _matrix(left), _matrix(cross), _matrix(right)
    _zero_diagonal_hermitian(lm)
    _zero_diagonal_hermitian(bm)
    if len(gm) != len(lm) or len(gm[0]) != len(bm):
        raise ValueError("cross dimensions must match the two diagonal blocks")
    ga = _adjoint(gm)
    left_term = _trace(_matmul(_matmul(lm, lm), _off_diagonal(_matmul(gm, ga))))
    right_term = _trace(_matmul(_matmul(bm, bm), _off_diagonal(_matmul(ga, gm))))
    mixed_term = _trace(_matmul(_matmul(_matmul(lm, gm), bm), ga))
    return 2 * _real(_add(_add(left_term, right_term), mixed_term))
