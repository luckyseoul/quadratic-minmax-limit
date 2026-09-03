import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
NOTE = ROOT / "evidence" / "NOTE_2026-09-02_BANASZCZYK_WEIGHTED_ANCHOR_ROUNDING.md"


def _skew_from_upper(n, upper):
    matrix = [[0] * n for _ in range(n)]
    for (i, j), value in upper.items():
        matrix[i][j] = value
        matrix[j][i] = -value
    return matrix


def _matvec(matrix, vector):
    return [sum(row[j] * vector[j] for j in range(len(vector))) for row in matrix]


def test_edge_block_embedding_is_exact():
    n = 6
    anchors = [
        [1, 1, -1, 1, -1, -1],
        [-1, 1, 1, -1, 1, -1],
        [1, -1, 1, -1, -1, 1],
    ]
    weights = [0.75, 1.25, 2.0]
    upper = {
        (i, j): (1 if (3 * i + 5 * j) % 4 < 2 else -1)
        for i in range(n)
        for j in range(i + 1, n)
    }
    matrix = _skew_from_upper(n, upper)

    embedded = [[0.0] * n for _ in anchors]
    for (i, j), sign in upper.items():
        for a, anchor in enumerate(anchors):
            embedded[a][i] += sign * weights[a] * anchor[j]
            embedded[a][j] -= sign * weights[a] * anchor[i]

    for a, anchor in enumerate(anchors):
        expected = [weights[a] * value for value in _matvec(matrix, anchor)]
        assert embedded[a] == expected


def test_uniform_budget_is_the_weighted_capacity_endpoint():
    for n, k in [(20, 1), (100, 7), (1000, 31)]:
        rho = n * math.sqrt(2 / math.pi) + math.sqrt(2 * n * math.log(2 * k))
        budget = 5 * math.sqrt(2 * k) * rho
        capacity = 50 * k * rho**2 / budget**2
        assert math.isclose(capacity, 1.0, rel_tol=1e-12, abs_tol=1e-12)
        expanded = 10 * n * math.sqrt(k / math.pi) + 10 * math.sqrt(
            k * n * math.log(2 * k)
        )
        assert math.isclose(budget, expanded, rel_tol=1e-12, abs_tol=1e-12)


def test_linear_anchor_capacity_constant():
    constant = (3 - 2 * math.sqrt(2)) / (25 * math.pi)
    assert math.isclose(constant, 0.0021845336957706, rel_tol=1e-13)


def test_anchor_incident_energy_budget_closes_exact_diamond():
    # This checks the worst-endpoint algebra in (10)--(11), independently of
    # any asymptotic replacement for M.
    M = 137.0
    for qx in [-M, -73.0, 0.0, 41.0, M]:
        budget = (2 * math.sqrt(2) - 1) * M - abs(qx)
        assert budget > 0
        for qy in [-M, -51.0, 0.0, 89.0, M]:
            assert abs(qx - qy) + budget <= 2 * math.sqrt(2) * M + 1e-12


def test_note_preserves_open_status_and_integrality_scope():
    text = NOTE.read_text()
    assert "50\\sum_{a=1}^k {\\rho_a^2\\over B_a^2}\\le1" in text
    assert "every upper-triangular\nentry of `R` is `+1` or `-1`" in text
    assert "does not orient every constraint" in text
    assert "still-open\ncover condition" in text
    assert "does not bound `Gamma(R)`" in text
    assert "Define the tournament `S=A circ R`" in text
    assert "where `0<delta<c_*` is arbitrary" in text
    assert "only deep external input" in text
    assert "MathOverflow limit remain open" in text
