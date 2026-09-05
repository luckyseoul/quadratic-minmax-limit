"""Independent exact checks; finite checks supplement the analytic proof."""
from fractions import Fraction
from itertools import combinations, product

import pytest

from original_mo_diagonal_compatibility import (
    SCOPE, additive_interval_certificate, conjugate_pair_parameters,
    relaxed_payment_certificate,
)


def test_conjugate_pair_identity_and_coherent_specialization():
    for a, gp, gm, d in product(range(-4, 5), range(4), range(4), range(-6, 7)):
        c, d0 = conjugate_pair_parameters(a, gp, gm)
        assert max(abs(a + d) + 2 * gp, abs(a - d) + 2 * gm) == c + abs(d - d0)
        if gp == gm:
            assert d0 == 0
            assert c + abs(d - d0) == abs(a) + abs(d) + 2 * gp


@pytest.mark.parametrize("cross_sign", [1, -1])
def test_actual_coherent_matrix_cannot_be_rescued_by_diagonal_skew(cross_sign):
    phases = list(product((1, 1j, -1, -1j), repeat=2))
    c = ((1, 1), (-1, 1))

    def bilinear(z, matrix, w):
        return sum(complex(z[i]).conjugate() * matrix[i][j] * w[j]
                   for i in range(2) for j in range(2))

    def rho(value):
        return max(abs(value.real), abs(value.imag))

    g = [[(1 + cross_sign * 1j) * x for x in row] for row in c]
    for left_skew, right_skew in product((-1, 1), repeat=2):
        left = ((0, 1 + 1j * left_skew), (1 - 1j * left_skew, 0))
        right = ((0, -1 + 1j * right_skew), (-1 - 1j * right_skew, 0))
        direct_norm = baseline = pair_max = 0
        for z, w in product(phases, repeat=2):
            actual_diagonal = bilinear(z, left, z) + bilinear(w, right, w)
            assert actual_diagonal.imag == 0
            a = 2 * (complex(z[0]).conjugate() * z[1]).real
            a -= 2 * (complex(w[0]).conjugate() * w[1]).real
            d = actual_diagonal.real - a
            s, cross = bilinear(z, c, w), bilinear(z, g, w)
            coherent_cross = abs(s.real) + abs(s.imag)
            assert rho(cross) == coherent_cross
            direct_norm = max(direct_norm, abs(actual_diagonal.real + 2 * cross.real))
            baseline = max(baseline, abs(a) + 2 * coherent_cross)
            pair_max = max(pair_max, abs(a) + abs(d) + 2 * coherent_cross)
        assert direct_norm == pair_max
        assert direct_norm >= baseline


def test_physical_global_minimizer_separates_pair_cycle_and_sign_costs():
    # This is one algebraic fixture, not an orientation or signing census.
    phases = (1, -1, 1j, -1j)
    cross_real = ((1, 1), (1, -1))
    cross_skew = ((-1, 1), (1, 1))

    def scalar(zeta, matrix, omega):
        return (matrix[0][0] + matrix[0][1] * omega
                + complex(zeta).conjugate() * (matrix[1][0] + matrix[1][1] * omega))

    def rho(value):
        return int(max(abs(value.real), abs(value.imag)))

    a, gp, gm = [], [], []
    for zeta in phases:
        a.append([])
        gp.append([])
        gm.append([])
        for omega in phases:
            real, skew = scalar(zeta, cross_real, omega), scalar(zeta, cross_skew, omega)
            a[-1].append(int(2 * complex(zeta).real + 2 * complex(omega).real))
            gp[-1].append(rho(real + 1j * skew))
            gm[-1].append(rho(real - 1j * skew))
    parameters = [[conjugate_pair_parameters(a[i][j], gp[i][j], gm[i][j])
                   for j in range(4)] for i in range(4)]
    assert max(c for row in parameters for c, _ in row) == 8
    assert [[parameters[i][j] for j in (2, 3)] for i in (0, 1)] == [
        [(8, 2), (8, -2)], [(8, -2), (8, 2)],
    ]
    assert not relaxed_payment_certificate(a, gp, gm, Fraction(99, 10))["feasible"]
    assert relaxed_payment_certificate(a, gp, gm, 10)["feasible"]
    assert max(c + abs(d0) for row in parameters for c, d0 in row) == 10

    # Independent direct Hermitian energies, including relative block phase.
    for left_sign, right_sign in product((-1, 1), repeat=2):
        direct_norm = 0
        for zeta, omega, rotation in product(phases, repeat=3):
            diagonal = (2 * complex(zeta).real + 2 * complex(omega).real
                        - 2 * left_sign * complex(zeta).imag
                        - 2 * right_sign * complex(omega).imag)
            cross = scalar(zeta, cross_real, omega) + 1j * scalar(zeta, cross_skew, omega)
            direct_norm = max(direct_norm, abs(diagonal + 2 * (rotation * cross).real))
        assert direct_norm == 12

    # The real signing has precisely one negative edge and edge-energy norm 4.
    energies = [x[0] * x[1] + x[2] * x[3]
                + sum(x[i] * cross_real[i][j] * x[2 + j]
                      for i in range(2) for j in range(2))
                for x in product((-1, 1), repeat=4)]
    assert max(map(abs, energies)) == 4
    assert sum(q * q for q in energies) == 16 * 6
    # Even edge energies and second moment six prove every order-four norm>=4.
    assert 10**2 < (2 * 4)**2 * 2 < 12**2


def test_feasible_certificate_replays_against_input():
    u, v = [Fraction(1, 3), -2, 5], [0, Fraction(-2, 7), 4]
    lo = [[x + y - Fraction(1, 11) for y in v] for x in u]
    hi = [[x + y + Fraction(1, 13) for y in v] for x in u]
    out = additive_interval_certificate(lo, hi)
    assert out["feasible"]
    assert all(lo[i][j] <= out["u"][i] + out["v"][j] <= hi[i][j]
               for i in range(3) for j in range(3))


def test_four_cycles_do_not_replace_all_cycles():
    # Generic intervals, NOT claimed to arise from an actual signing.
    lo = [[0, 1, -2], [-2, 0, 1], [1, -2, 0]]
    hi = [[0, 2, 2], [2, 0, 2], [2, 2, 0]]
    for rows in combinations(range(3), 2):
        for cols in combinations(range(3), 2):
            assert additive_interval_certificate(
                [[lo[i][j] for j in cols] for i in rows],
                [[hi[i][j] for j in cols] for i in rows],
            )["feasible"]
    out = additive_interval_certificate(lo, hi)
    assert not out["feasible"]
    cycle = out["negative_cycle"]
    assert len(cycle) == 6
    for k, (source, target, weight, side, i, j) in enumerate(cycle):
        assert target == cycle[(k + 1) % len(cycle)][0]
        if side == "upper":
            assert (source, target, weight) == (3 + j, i, hi[i][j])
        else:
            assert (source, target, weight) == (i, 3 + j, -lo[i][j])
    assert sum(arc[2] for arc in cycle) == out["weight"] < 0


def test_empty_edge_and_exact_threshold():
    out = additive_interval_certificate([[2]], [[1]])
    assert out["feasible"] is False
    assert len(out["negative_cycle"]) == 2
    assert out["weight"] == -1
    # Pairwise optimal payment is d0=-1 and its cost is 8.
    assert conjugate_pair_parameters(3, 3, 2) == (8, -1)
    assert relaxed_payment_certificate([[3]], [[3]], [[2]], 8)["feasible"]
    assert not relaxed_payment_certificate([[3]], [[3]], [[2]],
                                           Fraction(79, 10))["feasible"]


@pytest.mark.parametrize("lower,upper", [([], []), ([[]], [[]]),
                                          ([[0]], [[0, 1]]), ([[0], [1]], [[0]])])
def test_bad_shapes_rejected(lower, upper):
    with pytest.raises(ValueError):
        additive_interval_certificate(lower, upper)


@pytest.mark.parametrize("bad", [True, 0.5, float("nan"), "1"])
def test_inexact_inputs_rejected(bad):
    with pytest.raises(TypeError):
        additive_interval_certificate([[bad]], [[1]])
    with pytest.raises(TypeError):
        additive_interval_certificate([[0]], [[bad]])
    with pytest.raises(TypeError):
        relaxed_payment_certificate([[0]], [[1]], [[1]], bad)


def test_invalid_cross_norm_and_payment_shapes_rejected():
    with pytest.raises(ValueError):
        conjugate_pair_parameters(0, -1, 1)
    with pytest.raises(ValueError):
        conjugate_pair_parameters(0, 1, -1)
    with pytest.raises(ValueError):
        relaxed_payment_certificate([[0]], [[1, 1]], [[1]], 3)


def test_scope_does_not_promote_relaxation_to_a_signing():
    assert SCOPE["all_cycle_relaxed_criterion"]
    assert not SCOPE["skew_sign_realization_proved"]
    assert not SCOPE["multiplier_two_closed"]
    assert not SCOPE["original_mo_limit_closed"]


def test_new_proofs_are_linked_without_changing_the_live_gate():
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    notes = ("NOTE_2026-09-05_DIAGONAL_PAYMENT_COMPATIBILITY.md",
             "NOTE_2026-09-05_ACTUAL_DIAGONAL_MIXED_MOMENTS.md")
    for note in notes:
        assert (root / "evidence" / note).is_file()
        for doc in ("AGENTS.md", "STATUS.md", "HANDOFF.md", "README.md",
                    "solution.md", "evidence/PROPOSITION_DEDUP_AUDIT_2026-08-30.md"):
            assert note in (root / doc).read_text(encoding="utf-8")
    status = (root / "STATUS.md").read_text(encoding="utf-8")
    assert "Multiplier two and the original MO limit remain **OPEN**." in status
