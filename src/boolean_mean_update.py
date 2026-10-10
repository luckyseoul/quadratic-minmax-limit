"""Exact conditional moments of an independent Boolean correction layer.

All probabilities and answers are rational. Pair fields share their masks:
their dependence and endpoint tie rules are retained, not Gaussianized.
See evidence/NOTE_2026-09-23_CONDITIONAL_BOOLEAN_UPDATE.md for the proof.
This does not close an arbitrary outer mixture or an iterated update law.
"""

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from operator import index


def _rational(value):
    if isinstance(value, Fraction):
        return value
    # Deliberately reject floats, even when their binary value is exact.
    try:
        return Fraction(index(value))
    except TypeError as exc:
        raise TypeError("use integers or Fraction for exact inputs") from exc


def _probability(value):
    value = _rational(value)
    if not 0 <= value <= 1:
        raise ValueError("probabilities must lie in [0, 1]")
    return value


def _validated(matrix, spins):
    a = tuple(tuple(index(v) for v in row) for row in matrix)
    x = tuple(index(v) for v in spins)
    n = len(a)
    if n < 2 or len(x) != n or any(len(row) != n for row in a):
        raise ValueError("need an n by n signing and n spins, n >= 2")
    if any(v not in (-1, 1) for v in x):
        raise ValueError("spins must be +/-1")
    for i in range(n):
        if a[i][i] != 0:
            raise ValueError("the diagonal must be zero")
        for j in range(i):
            if a[i][j] not in (-1, 1) or a[i][j] != a[j][i]:
                raise ValueError("need a symmetric complete signing")
    return a, x


def _probabilities(values, n):
    if isinstance(values, Fraction) or hasattr(values, "__index__"):
        return (_probability(values),) * n
    result = tuple(_probability(v) for v in values)
    if len(result) != n:
        raise ValueError("need one flip probability per spin")
    return result


def mean_update_probabilities(matrix, spins, probability):
    """Flip only strictly unstable spins; a zero field keeps the old sign."""
    a, x = _validated(matrix, spins)
    p = _probability(probability)
    return tuple(p if x[i] * sum(v * s for v, s in zip(row, x)) < 0
                 else Fraction(0) for i, row in enumerate(a))


@lru_cache(maxsize=512)
def _difference_distribution(positive, negative):
    """Integer masses for sum Bernoulli(positive) - sum Bernoulli(negative).

    Arguments are sorted tuples of Fraction probabilities. Shared masks
    must not occur in both groups. Caches are bounded across calls.
    """
    mass, denominator = {0: 1}, 1
    for sign, group in ((1, positive), (-1, negative)):
        for p in group:
            u, v = p.numerator, p.denominator
            if not u:
                continue
            following = {}
            for k, weight in mass.items():
                if v != u:
                    following[k] = following.get(k, 0) + weight * (v - u)
                following[k + sign] = following.get(k + sign, 0) + weight * u
            mass = following
            denominator *= v
    return tuple(sorted(mass.items())), denominator


def _distribution(positive, negative):
    return _difference_distribution(tuple(sorted(positive)), tuple(sorted(negative)))


def _bernoulli(p):
    return tuple((b, w) for b, w in ((0, p.denominator - p.numerator),
                                     (1, p.numerator)) if w)


def _target(field, old_spin):
    return 1 if field > 0 else -1 if field < 0 else old_spin


def _energy(a, x):
    return sum(a[i][j] * x[i] * x[j]
               for i in range(len(a)) for j in range(i))


def instability_cubic_lower(mean_count):
    """Convex interpolation of g(0)=g(1)=0, g(k)=k^3 for k >= 2."""
    q = _rational(mean_count)
    if q < 0:
        raise ValueError("a mean count cannot be negative")
    k = q.numerator // q.denominator
    left = k ** 3 if k >= 2 else 0
    right = (k + 1) ** 3 if k + 1 >= 2 else 0
    return left + (q - k) * (right - left)


@dataclass(frozen=True)
class NormCertificate:
    """Phi satisfies quadratic*Phi^2 - linear*Phi - constant >= 0."""

    quadratic: Fraction
    linear: Fraction
    constant: Fraction

    def polynomial(self, value):
        f = _rational(value)
        return self.quadratic * f * f - self.linear * f - self.constant

    def root_interval(self, digits=12):
        """Rigorous rational enclosure of the nonnegative root, no float sqrt."""
        digits = index(digits)
        if digits < 0 or self.quadratic <= 0 or self.constant < 0:
            raise ValueError("need digits >= 0, quadratic > 0, constant >= 0")
        if not self.constant:
            root = max(Fraction(0), self.linear / self.quadratic)
            return root, root
        scale = 10 ** digits
        bound = (abs(self.linear) + self.constant) / self.quadratic + 1
        lo, hi = 0, (bound.numerator // bound.denominator + 1) * scale
        while hi - lo > 1:
            mid = (lo + hi) // 2
            value = self.polynomial(Fraction(mid, scale))
            if not value:
                root = Fraction(mid, scale)
                return root, root
            if value < 0:
                lo = mid
            else:
                hi = mid
        return Fraction(lo, scale), Fraction(hi, scale)


@dataclass(frozen=True)
class FieldMoments:
    law: tuple
    mean: Fraction
    second: Fraction
    absolute: Fraction
    positive: Fraction
    negative: Fraction
    zero: Fraction
    target_mean: Fraction
    unstable: Fraction


@dataclass(frozen=True)
class ConditionalMoments:
    flip_probabilities: tuple
    spin_mean: tuple
    fields: tuple
    target_second: tuple
    energy: Fraction
    target_energy: Fraction
    field_l1: Fraction
    unstable_count: Fraction

    def second_mean_energy(self, weight):
        t = _probability(weight)
        return ((1 - t) ** 2 * self.energy + t * (1 - t) * self.field_l1
                + t ** 2 * self.target_energy)

    def norm_certificate(self, weight, *, coarse_target=False):
        """Include the finite interior gap; weight must be at most 1/2.

        coarse_target=True restores the old E Q(Y) >= -Phi surrogate,
        for an exact like-for-like comparison without parameter tuning.
        """
        t = _probability(weight)
        if t > Fraction(1, 2):
            raise ValueError("the interior-gap certificate needs weight <= 1/2")
        d = t ** 2 * instability_cubic_lower(self.unstable_count) / 24
        if coarse_target:
            c = (1 - t) ** 2 * self.energy + t * (1 - t) * self.field_l1
            return NormCertificate(1 + t ** 2, c, d)
        return NormCertificate(Fraction(1), self.second_mean_energy(t), d)


def independent_flip_moments(matrix, spins, probabilities):
    """Exact moments after independent flips from one fixed Boolean source.

    For B_j ~ Bernoulli(p_j), X'_j=x_j(1-2B_j) and Y'_i is the sign
    of (A X')_i, keeping X'_i at a zero field. Return row laws and the
    complete matrix E[Y' Y'^T], including dependent-field tie effects.
    General rational probabilities cost O(n^4) integer arithmetic
    operations; integer bit sizes depend on the supplied denominators.
    """
    a, x = _validated(matrix, spins)
    n = len(a)
    p = _probabilities(probabilities, n)
    m = tuple(xi * (1 - 2 * pi) for xi, pi in zip(x, p))
    h = tuple(sum(a[i][j] * x[j] for j in range(n)) for i in range(n))
    rows = []
    for i in range(n):
        groups = {1: [], -1: []}
        for j in range(n):
            if j != i:
                groups[a[i][j] * x[j]].append(p[j])
        mass, den = _distribution(groups[1], groups[-1])
        law = tuple(sorted((h[i] - 2 * k, Fraction(w, den)) for k, w in mass))
        positive = sum((w for f, w in law if f > 0), Fraction(0))
        negative = sum((w for f, w in law if f < 0), Fraction(0))
        zero = sum((w for f, w in law if f == 0), Fraction(0))
        rows.append(FieldMoments(
            law, sum(f * w for f, w in law), sum(f * f * w for f, w in law),
            sum(abs(f) * w for f, w in law), positive, negative, zero,
            positive - negative + m[i] * zero,
            ((1 + m[i]) * negative + (1 - m[i]) * positive) / 2))

    second = [[Fraction(int(i == j)) for j in range(n)] for i in range(n)]
    for i in range(n):
        for j in range(i):
            groups = {(s, t): [] for s in (-1, 1) for t in (-1, 1)}
            for k in range(n):
                if k != i and k != j:
                    groups[a[i][k] * x[k], a[j][k] * x[k]].append(p[k])
            u_mass, u_den = _distribution(groups[1, 1], groups[-1, -1])
            v_mass, v_den = _distribution(groups[1, -1], groups[-1, 1])
            numerator = 0
            for bi, wi in _bernoulli(p[i]):
                for bj, wj in _bernoulli(p[j]):
                    endpoint_weight = wi * wj
                    for u, wu in u_mass:
                        for v, wv in v_mass:
                            fi = h[i] - 2 * (u + v + a[i][j] * x[j] * bj)
                            fj = h[j] - 2 * (u - v + a[i][j] * x[i] * bi)
                            yi = _target(fi, x[i] * (1 - 2 * bi))
                            yj = _target(fj, x[j] * (1 - 2 * bj))
                            numerator += endpoint_weight * wu * wv * yi * yj
            denominator = u_den * v_den * p[i].denominator * p[j].denominator
            second[i][j] = second[j][i] = Fraction(numerator, denominator)

    return ConditionalMoments(
        p, m, tuple(rows), tuple(tuple(row) for row in second), _energy(a, m),
        sum(a[i][j] * second[i][j] for i in range(n) for j in range(i)),
        sum(row.absolute for row in rows), sum(row.unstable for row in rows))
