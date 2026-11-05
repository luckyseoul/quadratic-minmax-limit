"""Constructive upper bounds on F(n) by local search on edge signs."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np


def _score_from_matrix(a: np.ndarray, x: np.ndarray) -> int:
    # a is symmetric ±1 with zero diagonal; score = (1/2) x^T a x
    return int((x @ (a @ x)) // 2)


def maximize_abs_score(
    a: np.ndarray,
    rng: np.random.Generator,
    n_restarts: int = 32,
    max_passes: int = 40,
) -> int:
    """
    Approximate max_x |x^T A x|/2 by multi-start greedy coordinate ascent.

    For each restart, flip coordinates that increase |score| until a local peak.
    """
    n = a.shape[0]
    best = 0
    for _ in range(n_restarts):
        x = rng.choice(np.array([-1, 1], dtype=np.int8), size=n)
        # Keep a representative with x[-1] = +1
        if x[-1] == -1:
            x = -x
        cur = abs(_score_from_matrix(a, x))
        for _pass in range(max_passes):
            improved = False
            order = rng.permutation(n - 1)  # last coord fixed by switching
            for i in order:
                x[i] = -x[i]
                val = abs(_score_from_matrix(a, x))
                if val > cur:
                    cur = val
                    improved = True
                else:
                    x[i] = -x[i]
            if not improved:
                break
        best = max(best, cur)
    return int(best)


def random_seidel_matrix(n: int, rng: np.random.Generator) -> np.ndarray:
    a = np.zeros((n, n), dtype=np.int8)
    for i in range(n):
        for j in range(i + 1, n):
            s = rng.choice([-1, 1])
            a[i, j] = s
            a[j, i] = s
    return a


def flip_edge(a: np.ndarray, i: int, j: int) -> None:
    a[i, j] *= -1
    a[j, i] *= -1


@dataclass
class ConstructiveResult:
    n: int
    upper_f: int
    normalized: float
    method: str
    seed: int
    evals: int


def local_search_upper(
    n: int,
    seed: int = 0,
    n_outer: int = 200,
    n_restarts: int = 16,
    temperature0: float = 2.0,
) -> ConstructiveResult:
    """
    Simulated-annealing style search for a signing with small max_|score|.

    Returns an *upper bound* on F(n): any concrete a gives
    F(n) <= max_x |score(a,x)|, and the inner max is itself approximated
    from below, so the reported value can be slightly optimistic. We
    mitigate by many restarts and a final dense re-evaluation.
    """
    rng = np.random.default_rng(seed)
    a = random_seidel_matrix(n, rng)
    cur = maximize_abs_score(a, rng, n_restarts=n_restarts)
    best_a = a.copy()
    best = cur
    evals = 1

    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]

    for t in range(1, n_outer + 1):
        temperature = temperature0 * (1.0 - t / (n_outer + 1))
        i, j = pairs[int(rng.integers(len(pairs)))]
        flip_edge(a, i, j)
        val = maximize_abs_score(a, rng, n_restarts=max(4, n_restarts // 2))
        evals += 1
        delta = val - cur
        if delta <= 0 or rng.random() < math.exp(-delta / max(temperature, 1e-9)):
            cur = val
            if val < best:
                best = val
                best_a = a.copy()
        else:
            flip_edge(a, i, j)  # revert

    # Final denser evaluation of the best matrix found
    final = maximize_abs_score(best_a, rng, n_restarts=max(64, 4 * n_restarts), max_passes=80)
    evals += 1
    upper = max(best, final)

    return ConstructiveResult(
        n=n,
        upper_f=int(upper),
        normalized=upper / (n ** 1.5),
        method="sa_edge_flips+greedy_x",
        seed=seed,
        evals=evals,
    )


def conference_like_upper(n: int, seed: int = 0) -> ConstructiveResult | None:
    """
    Try a deterministic ±1 Seidel pattern from a cyclic tournament / Paley-like
    quadratic residue construction when n ≡ 3 (mod 4) is prime, else skip.
    """
    if n % 4 != 3 or not _is_prime(n):
        return None
    a = np.zeros((n, n), dtype=np.int8)
    residues = {pow(k, 2, n) for k in range(1, n)}
    for i in range(n):
        for j in range(i + 1, n):
            d = (j - i) % n
            s = 1 if d in residues else -1
            a[i, j] = s
            a[j, i] = s
    rng = np.random.default_rng(seed)
    upper = maximize_abs_score(a, rng, n_restarts=128, max_passes=100)
    return ConstructiveResult(
        n=n,
        upper_f=int(upper),
        normalized=upper / (n ** 1.5),
        method="paley_seidel",
        seed=seed,
        evals=1,
    )


def _is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    r = int(math.isqrt(n))
    for d in range(3, r + 1, 2):
        if n % d == 0:
            return False
    return True
