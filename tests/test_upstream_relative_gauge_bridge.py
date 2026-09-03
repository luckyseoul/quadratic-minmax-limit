from itertools import product


def qform(a, x):
    return sum(
        a[i][j] * x[i] * x[j]
        for i in range(len(a))
        for j in range(i + 1, len(a))
    )


def matmul_diag_left_right(alpha, c, beta):
    return [
        [alpha[i] * c[i][j] * beta[j] for j in range(len(beta))]
        for i in range(len(alpha))
    ]


def bilinear(c, x, y):
    return sum(
        c[i][j] * x[i] * y[j]
        for i in range(len(x))
        for j in range(len(y))
    )


def block_energy(a, b, c, alpha, beta, tau, x, y):
    xa = tuple(alpha[i] * x[i] for i in range(len(x)))
    yb = tuple(beta[j] * y[j] for j in range(len(y)))
    return qform(a, xa) + tau * qform(b, yb) + bilinear(c, x, y)


def block_max(a, b, c, alpha, beta, tau):
    return max(
        abs(block_energy(a, b, c, alpha, beta, tau, x, y))
        for x in product((-1, 1), repeat=len(a))
        for y in product((-1, 1), repeat=len(b))
    )


def diamond(a, b, d, tau):
    return max(
        abs(qform(a, x) + tau * qform(b, y)) + abs(bilinear(d, x, y))
        for x in product((-1, 1), repeat=len(a))
        for y in product((-1, 1), repeat=len(b))
    )


def sign_matrices(n, k):
    for entries in product((-1, 1), repeat=n * k):
        yield [list(entries[i * k : (i + 1) * k]) for i in range(n)]


def character_half_statistics(occupancy, character):
    size = len(occupancy)
    mean = sum(occupancy) / size
    coefficient = sum(b * chi for b, chi in zip(occupancy, character)) / size
    sign = 1 if coefficient >= 0 else -1
    opposite_mass = sum(b for b, chi in zip(occupancy, character) if chi == -sign)
    opposite_zeros = sum(
        1
        for b, chi in zip(occupancy, character)
        if chi == -sign and b == 0
    )
    return mean, coefficient, opposite_mass, opposite_zeros


def test_gauge_absorption_identity():
    a = [[0, 1, -1], [1, 0, 1], [-1, 1, 0]]
    b = [[0, -1], [-1, 0]]
    c = [[1, -1], [-1, -1], [1, 1]]

    for alpha in product((-1, 1), repeat=3):
        for beta in product((-1, 1), repeat=2):
            d = matmul_diag_left_right(alpha, c, beta)
            for tau in (-1, 1):
                assert block_max(a, b, c, alpha, beta, tau) == diamond(a, b, d, tau)


def test_free_seed_makes_gauge_orbit_redundant():
    a = [[0, 1], [1, 0]]
    b = [[0, -1], [-1, 0]]

    left = min(
        block_max(a, b, c, alpha, beta, tau)
        for c in sign_matrices(2, 2)
        for alpha in product((-1, 1), repeat=2)
        for beta in product((-1, 1), repeat=2)
        for tau in (-1, 1)
    )
    right = min(
        diamond(a, b, d, tau)
        for d in sign_matrices(2, 2)
        for tau in (-1, 1)
    )
    assert left == right


def test_opposite_diagonal_specialization():
    a = [[0, 1, -1], [1, 0, 1], [-1, 1, 0]]
    c = [[1, -1, 1], [1, -1, -1], [-1, 1, 1]]
    ones = (1, 1, 1)

    assert block_max(a, a, c, ones, ones, -1) == diamond(a, a, c, -1)


def test_one_character_abundance_identity():
    occupancy = [0, 0, 0, 2, 4, 1, 3, 2]
    character = [1, 1, 1, 1, -1, -1, -1, -1]
    mean, coefficient, opposite_mass, opposite_zeros = character_half_statistics(
        occupancy, character
    )
    delta = mean - abs(coefficient)

    assert opposite_mass == len(occupancy) * delta / 2
    assert opposite_zeros >= len(occupancy) * (1 - delta) / 2
