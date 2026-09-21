"""Exact local algebra only; not an all-orders theorem verifier or census."""

import json

import sympy as s


def main():
    a, b, q = s.symbols("a b q", positive=True)
    spins = [(x, y) for x in (-1, 1) for y in (-1, 1)]
    weights = {(x, y): a**x * b**y * q ** (x * y) for x, y in spins}
    normalizer = sum(weights.values())

    def mean(values):
        return s.factor(sum(weights[p] * values[p] for p in spins) / normalizer)

    def tanh_exp(z):
        return (z**2 - 1) / (z**2 + 1)

    m1 = {(x, y): tanh_exp(a * q**y) for x, y in spins}
    m2 = {(x, y): tanh_exp(b * q**x) for x, y in spins}
    r1 = {(x, y): x - m1[x, y] for x, y in spins}
    r2 = {(x, y): y - m2[x, y] for x, y in spins}
    d1 = {p: 1 - m1[p] ** 2 for p in spins}
    d2 = {p: 1 - m2[p] ** 2 for p in spins}
    delta2 = (tanh_exp(b * q) - tanh_exp(b / q)) / 2
    passed = []

    def check(name, expression):
        residual = s.factor(expression)
        if residual != 0:
            raise AssertionError((name, residual))
        passed.append(name)

    check("score_1_centered", mean(r1))
    check("score_2_centered", mean(r2))
    check("off_diagonal_spin_score", mean({(x, y): y * r1[x, y] for x, y in spins}))
    check("diagonal_spin_score", mean({(x, y): x * r1[x, y] for x, y in spins}) - mean(d1))
    check("diagonal_score_covariance", mean({p: r1[p] ** 2 for p in spins}) - mean(d1))
    check("exact_pair_score_identity", mean({p: r1[p] * r2[p] for p in spins}) + delta2 * mean(d1))
    delta1 = (tanh_exp(a * q) - tanh_exp(a / q)) / 2
    check("pair_identity_symmetric", delta2 * mean(d1) - delta1 * mean(d2))

    z = s.symbols("z", real=True)
    # d/dh sech(h)^2 with z=tanh(h); the proof uses |z|<=1.
    h = s.symbols("h", real=True)
    check("sech_derivative", s.diff(1 - s.tanh(h) ** 2, h) + 2 * s.tanh(h) * (1 - s.tanh(h) ** 2))
    for endpoint in (-1, 1):
        check(f"lipschitz_integral_{endpoint}", s.integrate(1 - endpoint * z, (z, -1, 1)) - 2)

    u = s.symbols("u", positive=True)
    rb, rd = s.symbols("rbar rdiff", real=True)
    ch, sh = (u + 1 / u) / 2, (u - 1 / u) / 2
    check("phase_flip_product", (ch - sh * (rb + rd)) * (ch - sh * (rb - rd)) - ((ch - sh * rb) ** 2 - sh**2 * rd**2))
    check("edge_threshold", ch - 1 - sh * (u - 1) / (u + 1))
    print(json.dumps({"classification": "exact symbolic local algebra regression", "sympy_version": s.__version__, "checks": passed, "count": len(passed), "status": "PASS", "limitations": "No matrix census; no formal verification of PSD, probabilistic, all-orders, or convergence claims."}, indent=2))


if __name__ == "__main__":
    main()
