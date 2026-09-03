from __future__ import annotations

import itertools


A = (
    (0, 1, -1, 1),
    (1, 0, 1, -1),
    (-1, 1, 0, 1),
    (1, -1, 1, 0),
)


def _quadratic(y: tuple[int, ...]) -> int:
    return sum(
        A[i][j] * y[i] * y[j]
        for i, j in itertools.combinations(range(len(A)), 2)
    )


def test_prefix_identity_and_random_order_variance() -> None:
    """Exact factor check for equations (5), (12), and (16) of the note."""
    n = len(A)
    states = tuple(itertools.product((-1, 1), repeat=n))
    orders = tuple(itertools.permutations(range(n)))

    # Prefix jump and total-variation normalization, for every state and order.
    for y in states:
        for order in orders:
            position = {vertex: index for index, vertex in enumerate(order)}
            R = tuple(
                tuple(
                    0
                    if i == j
                    else A[i][j] * (1 if position[i] < position[j] else -1)
                    for j in range(n)
                )
                for i in range(n)
            )
            local = tuple(
                y[i] * sum(R[i][j] * y[j] for j in range(n))
                for i in range(n)
            )

            flipped = list(y)
            prefix_energies = [_quadratic(tuple(flipped))]
            for vertex in order:
                flipped[vertex] *= -1
                prefix_energies.append(_quadratic(tuple(flipped)))

            jumps = tuple(
                prefix_energies[k - 1] - prefix_energies[k]
                for k in range(1, n + 1)
            )
            assert jumps == tuple(2 * local[vertex] for vertex in order)

            direct_max = max(
                abs(
                    sum(
                        x[i] * R[i][j] * y[j]
                        for i in range(n)
                        for j in range(n)
                    )
                )
                for x in states
            )
            assert direct_max == sum(abs(value) for value in local)
            assert 2 * direct_max == sum(abs(jump) for jump in jumps)

    # Fixed-cut permutation variance and its Boolean-state average.
    for mask in range(1, (1 << n) - 1):
        U = tuple(i for i in range(n) if mask & (1 << i))
        V = tuple(i for i in range(n) if not mask & (1 << i))
        h = len(U) * len(V)
        sum_D = 0
        sum_Z_squared = 0

        for y in states:
            W = tuple(
                tuple(A[i][j] * y[i] * y[j] for j in V)
                for i in U
            )
            row_sums = tuple(sum(row) for row in W)
            column_sums = tuple(
                sum(W[a][b] for a in range(len(U))) for b in range(len(V))
            )
            D = sum(value * value for value in row_sums + column_sums)

            values = []
            for order in orders:
                position = {vertex: index for index, vertex in enumerate(order)}
                # Z=G-F: an edge contributes positively when V precedes U.
                Z = sum(
                    W[a][b]
                    * (1 if position[U[a]] > position[V[b]] else -1)
                    for a in range(len(U))
                    for b in range(len(V))
                )
                values.append(Z)

            assert sum(values) == 0
            assert 3 * sum(value * value for value in values) == len(orders) * (
                h + D
            )
            sum_D += D
            sum_Z_squared += sum(value * value for value in values)

        assert sum_D == len(states) * 2 * h
        assert sum_Z_squared == len(states) * len(orders) * h
