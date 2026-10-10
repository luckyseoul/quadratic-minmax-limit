"""Direct mask-oracle checks on fixed fixtures, not a signing census."""

from fractions import Fraction as F
from itertools import product
import unittest

from src.boolean_mean_update import (
    NormCertificate, independent_flip_moments, instability_cubic_lower,
    mean_update_probabilities,
)


def cycle_five():
    return tuple(tuple(0 if i == j else 1 if (i - j) % 5 in (1, 4) else -1
                       for j in range(5)) for i in range(5))


MIXED_SIX = (
    (0, 1, 1, 1, -1, -1), (1, 0, 1, -1, 1, -1),
    (1, 1, 0, 1, -1, 1), (1, -1, 1, 0, 1, -1),
    (-1, 1, -1, 1, 0, 1), (-1, -1, 1, -1, 1, 0),
)


def energy(a, x):
    return sum(a[i][j] * x[i] * x[j]
               for i in range(len(x)) for j in range(i))


def mask_oracle(a, x, p, t):
    """Independent implementation: literally visit every mask."""
    n = len(x)
    rows = [dict() for _ in x]
    mean = [F(0) for _ in x]
    unstable = [F(0) for _ in x]
    second = [[F(0) for _ in x] for _ in x]
    total = e = ey = l1 = k3 = z_energy = F(0)
    for mask in product((0, 1), repeat=n):
        w = F(1)
        for bit, prob in zip(mask, p):
            w *= prob if bit else 1 - prob
        if not w:
            continue
        xp = tuple(s * (-1 if bit else 1) for s, bit in zip(x, mask))
        fields = tuple(sum(v * s for v, s in zip(row, xp)) for row in a)
        y = tuple(1 if h > 0 else -1 if h < 0 else s
                  for h, s in zip(fields, xp))
        k = sum(s * h < 0 for s, h in zip(xp, fields))
        z = tuple((1 - t) * s + t * r for s, r in zip(xp, y))
        total += w
        e += w * energy(a, xp)
        ey += w * energy(a, y)
        l1 += w * sum(map(abs, fields))
        k3 += w * (k ** 3 if k >= 2 else 0)
        z_energy += w * energy(a, z)
        for i in range(n):
            rows[i][fields[i]] = rows[i].get(fields[i], F(0)) + w
            mean[i] += w * y[i]
            unstable[i] += w * (xp[i] * fields[i] < 0)
            for j in range(n):
                second[i][j] += w * y[i] * y[j]
    return dict(total=total, rows=rows, mean=mean, unstable=unstable,
                second=second, e=e, ey=ey, l1=l1, k3=k3, z_energy=z_energy)


class ConditionalUpdateTests(unittest.TestCase):
    def check_fixture(self, a, x, p):
        t = F(1, 4)
        actual = independent_flip_moments(a, x, p)
        expected = mask_oracle(a, x, p, t)
        self.assertEqual(expected['total'], 1)
        self.assertEqual(actual.energy, expected['e'])
        self.assertEqual(actual.target_energy, expected['ey'])
        self.assertEqual(actual.field_l1, expected['l1'])
        self.assertEqual(actual.second_mean_energy(t), expected['z_energy'])
        self.assertEqual(actual.unstable_count, sum(expected['unstable']))
        self.assertEqual(actual.target_second, tuple(map(tuple, expected['second'])))
        for i, row in enumerate(actual.fields):
            law = expected['rows'][i]
            self.assertEqual(dict(row.law), law)
            self.assertEqual(row.mean, sum(h * w for h, w in law.items()))
            self.assertEqual(row.second, sum(h * h * w for h, w in law.items()))
            self.assertEqual(row.absolute, sum(abs(h) * w for h, w in law.items()))
            self.assertEqual(row.positive, sum(w for h, w in law.items() if h > 0))
            self.assertEqual(row.negative, sum(w for h, w in law.items() if h < 0))
            self.assertEqual(row.zero, law.get(0, 0))
            self.assertEqual(row.target_mean, expected['mean'][i])
            self.assertEqual(row.unstable, expected['unstable'][i])
        self.assertLessEqual(instability_cubic_lower(actual.unstable_count), expected['k3'])
        # A small-fixture norm oracle only; no large input is rescored here.
        phi = max(abs(energy(a, s)) for s in product((-1, 1), repeat=len(x)))
        for weight in (F(0), F(1, 4), F(1, 2)):
            new = actual.norm_certificate(weight)
            old = actual.norm_certificate(weight, coarse_target=True)
            for cert in (new, old):
                self.assertGreaterEqual(cert.polynomial(phi), 0)
                lo, hi = cert.root_interval(10)
                self.assertLessEqual(lo, phi)
                self.assertLessEqual(hi - lo, F(1, 10 ** 10))
                self.assertLessEqual(cert.polynomial(lo), 0)
                self.assertGreaterEqual(cert.polynomial(hi), 0)
            # Retain elementary moment lower bounds when comparing certificates.
            new_lo = max(new.root_interval(10)[0], abs(actual.target_energy))
            self.assertLessEqual(old.root_interval(10)[0], new_lo + F(1, 10 ** 10))

    def test_two_vertices_endpoints(self):
        for p in ((F(0), F(1)), (F(1, 3), F(2, 5)), (F(1), F(1))):
            with self.subTest(p=p):
                self.check_fixture(((0, -1), (-1, 0)), (1, -1), p)

    def test_odd_order_ties_and_nonuniform_probabilities(self):
        self.check_fixture(cycle_five(), (1, 1, 1, 1, 1),
                           (F(1, 2), F(1, 3), F(2, 3), F(1, 5), F(3, 4)))

    def test_deterministic_ties_keep_old_spin(self):
        a = cycle_five()
        result = independent_flip_moments(a, (1,) * 5, 0)
        self.assertEqual(tuple(row.zero for row in result.fields), (1,) * 5)
        self.assertEqual(tuple(row.target_mean for row in result.fields), (1,) * 5)
        self.assertEqual(result.unstable_count, 0)
        self.check_fixture(a, (-1, -1, 1, 1, 1), (F(0),) * 5)

    def test_all_unstable_uniform_probabilities(self):
        a = tuple(tuple(0 if i == j else -1 for j in range(6)) for i in range(6))
        for p in (F(0), F(1), F(1, 2), F(1, 10)):
            with self.subTest(p=p):
                self.assertEqual(mean_update_probabilities(a, (1,) * 6, p), (p,) * 6)
                self.check_fixture(a, (1,) * 6, (p,) * 6)

    def test_four_pair_classes_and_random_endpoints(self):
        self.assertEqual({(MIXED_SIX[0][k], MIXED_SIX[1][k]) for k in range(2, 6)},
                         {(1, 1), (1, -1), (-1, 1), (-1, -1)})
        self.check_fixture(MIXED_SIX, (1,) * 6,
                           (F(1, 3), F(2, 5), F(1, 2), F(1, 4), F(3, 5), F(2, 3)))

    def test_zero_and_one_in_nonuniform_vector(self):
        self.check_fixture(MIXED_SIX, (1, -1, 1, -1, -1, 1),
                           (F(0), F(1), F(1, 3), F(2, 3), F(1, 2), F(1, 5)))

    def test_actual_mean_update_policy(self):
        a, x = cycle_five(), (1, -1, 1, -1, 1)
        p = mean_update_probabilities(a, x, F(1, 4))
        self.assertIn(F(0), p)
        self.assertIn(F(1, 4), p)
        self.check_fixture(a, x, p)

    def test_uniform_scalar_matches_vector(self):
        x = (1, -1, -1, 1, 1, -1)
        self.assertEqual(independent_flip_moments(MIXED_SIX, x, F(1, 3)),
                         independent_flip_moments(MIXED_SIX, x, (F(1, 3),) * 6))

    def test_switching_covariance(self):
        a, x = cycle_five(), (1, -1, 1, -1, 1)
        d, p = (-1, 1, -1, 1, 1), (F(1, 3), F(1, 5), F(1), F(0), F(2, 3))
        switched = tuple(tuple(d[i] * d[j] * a[i][j] for j in range(5)) for i in range(5))
        first = independent_flip_moments(a, x, p)
        second = independent_flip_moments(switched, tuple(s * r for s, r in zip(x, d)), p)
        self.assertEqual(first.energy, second.energy)
        self.assertEqual(first.target_energy, second.target_energy)
        self.assertEqual(first.field_l1, second.field_l1)
        self.assertEqual(first.unstable_count, second.unstable_count)
        for i in range(5):
            self.assertEqual(dict(second.fields[i].law),
                             {d[i] * h: w for h, w in first.fields[i].law})
            for j in range(5):
                self.assertEqual(second.target_second[i][j], d[i] * d[j] * first.target_second[i][j])

    def test_cubic_interpolation_boundaries(self):
        for q, expected in ((F(0), 0), (F(1, 2), 0), (F(1), 0),
                            (F(3, 2), 4), (F(2), 8), (F(5, 2), F(35, 2)), (F(3), 27)):
            self.assertEqual(instability_cubic_lower(q), expected)

    def test_root_brackets(self):
        for cert in (NormCertificate(F(1), F(2), F(0)),
                     NormCertificate(F(1), F(-2), F(0)),
                     NormCertificate(F(1), F(-5), F(17, 3)),
                     NormCertificate(F(17, 16), F(7, 2), F(11, 7)),
                     NormCertificate(F(1, 100), F(0), F(1, 100))):
            lo, hi = cert.root_interval(8)
            self.assertLessEqual(cert.polynomial(lo), 0)
            self.assertGreaterEqual(cert.polynomial(hi), 0)
            self.assertLessEqual(hi - lo, F(1, 10 ** 8))

    def test_invalid_inputs(self):
        a, x = cycle_five(), (1,) * 5
        for bad in (F(-1, 2), F(3, 2), (F(1, 2),) * 4):
            with self.assertRaises(ValueError):
                independent_flip_moments(a, x, bad)
        with self.assertRaises(TypeError):
            independent_flip_moments(a, x, [0.25] * 5)
        for bad_a, bad_x in ((((0,),), (1,)), (a, (1,) * 4),
                             (a, (0,) * 5), (((0, 1), (-1, 0)), (1, 1)),
                             (((1, 1), (1, 0)), (1, 1)),
                             (((0, 2), (2, 0)), (1, 1))):
            with self.assertRaises(ValueError):
                independent_flip_moments(bad_a, bad_x, 0)
        result = independent_flip_moments(a, x, 0)
        with self.assertRaises(ValueError):
            result.norm_certificate(F(3, 4))
        with self.assertRaises(ValueError):
            instability_cubic_lower(F(-1))


if __name__ == '__main__':
    unittest.main(verbosity=2)
