import math


def bound(t: float):
    kappa = 2.0 / math.pi
    z = t * t / (1.0 + t * t)
    a = 1.0 - kappa * math.asin(z) + kappa * z
    e = kappa * t / (1.0 + t * t)
    f = math.sqrt(kappa * a)
    disc = math.sqrt(f * f + (f - 2.0 * e) ** 2)
    p = (f - 2.0 * e) / (f + disc)
    B = e - f / 2.0 + disc / 2.0
    return B, p, e, f


def test_tilted_paired_field_improves_t_equals_one():
    B1, p1, _, _ = bound(1.0)
    Bt, pt, _, _ = bound(993.0 / 1000.0)
    assert 0.0 < p1 < 1.0
    assert 0.0 < pt < 1.0
    assert B1 > 0.3258474
    assert Bt > B1 + 5.6e-6
    assert Bt > 0.32585303


def test_envelope_derivative_at_one_points_toward_smaller_t():
    kappa = 2.0 / math.pi
    _, p1, _, f1 = bound(1.0)
    z = 0.5
    zprime = 0.5
    aprime = kappa * (1.0 - 1.0 / math.sqrt(1.0 - z * z)) * zprime
    fprime = kappa * aprime / (2.0 * f1)
    Bprime = p1 * (1.0 - p1) * fprime / (1.0 + p1 * p1)
    assert Bprime < 0.0
