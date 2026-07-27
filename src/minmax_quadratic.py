"""
Min-max quadratic form of ±1 coefficients.

m_n = min_{a_ij=±1} max_{x=±1} |sum_{i<j} a_ij x_i x_j|

Equivalently m_n = (1/2) min_A max_x |x^T A x| over symmetric zero-diagonal ±1 matrices.
"""
from __future__ import annotations

import numpy as np
from typing import Optional, Tuple, List


def form_Q(A: np.ndarray, x: np.ndarray) -> float:
    """Q = sum_{i<j} A_ij x_i x_j = (1/2) x^T A x for zero-diagonal symmetric A."""
    return 0.5 * float(x @ A @ x)


def phi(A: np.ndarray, x_batch: Optional[np.ndarray] = None) -> float:
    """max_x |Q(A,x)|. If x_batch given, max over those rows; else brute force for n<=20."""
    n = A.shape[0]
    if x_batch is not None:
        # x_batch: (m, n) of ±1
        # Q = 0.5 * sum_{ab} x_a A_ab x_b
        vals = 0.5 * np.einsum('ia,ab,ib->i', x_batch, A, x_batch)
        return float(np.max(np.abs(vals)))
    if n > 20:
        raise ValueError("brute force phi only for n<=20; pass x_batch or use phi_local")
    npat = 1 << (n - 1)
    best = 0.0
    for xb in range(npat):
        x = np.ones(n, dtype=np.float64)
        for i in range(n - 1):
            x[i + 1] = 1.0 if (xb >> i) & 1 else -1.0
        best = max(best, abs(form_Q(A, x)))
    return best


def phi_local(A: np.ndarray, restarts: int = 40, rng: Optional[np.random.Generator] = None) -> float:
    """Local search lower bound on phi(A) = max |Q|."""
    if rng is None:
        rng = np.random.default_rng(0)
    n = A.shape[0]
    best = 0.0
    for _ in range(restarts):
        x = rng.choice([-1.0, 1.0], size=n)
        Ax = A @ x
        qfull = float(x @ Ax)
        cur = abs(qfull) / 2.0
        improved = True
        while improved:
            improved = False
            for i in range(n):
                delta = -4.0 * x[i] * Ax[i]
                new_abs = abs(qfull + delta) / 2.0
                if new_abs > cur + 1e-12:
                    xi = x[i]
                    x[i] = -xi
                    Ax -= 2.0 * xi * A[:, i]
                    qfull += delta
                    cur = new_abs
                    improved = True
        best = max(best, cur)
    return best


def greedy_Q(A: np.ndarray, order: Optional[np.ndarray] = None) -> float:
    """
    Online greedy signing: produces some x with
    |Q| = sum_t |sum_{u<t} A[order_t, order_u] x[order_u]|.
    This is a constructive lower bound on phi(A).
    """
    n = A.shape[0]
    if order is None:
        order = np.arange(n)
    x = np.zeros(n, dtype=np.float64)
    total = 0.0
    for t, i in enumerate(order):
        s = 0.0
        for u in range(t):
            j = order[u]
            s += A[i, j] * x[j]
        x[i] = 1.0 if s >= 0 else -1.0
        total += abs(s)
    return total


def exact_m(n: int) -> int:
    """
    Exact m_n by exhaustive search over matrices (with vertex-folding a[0,j]=+1)
    and exhaustive x (with x[0]=+1). Feasible for n <= 9 (n=9 is 2^28; use
    exact_m9_parallel.py for sharded 80-worker run).
    """
    if n < 2:
        raise ValueError("n>=2")
    if n > 9:
        raise ValueError("exact_m only supported for n<=9 in this implementation")
    free = [(i, j) for i in range(1, n) for j in range(i + 1, n)]
    nf = len(free)
    npat = 1 << (n - 1)
    X = np.ones((npat, n), dtype=np.int16)
    for xb in range(npat):
        for i in range(n - 1):
            X[xb, i + 1] = 1 if (xb >> i) & 1 else -1
    XXf = np.stack([X[:, 0] * X[:, j] for j in range(1, n)], axis=1).astype(np.int32)
    if nf:
        XXu = np.stack([X[:, i] * X[:, j] for i, j in free], axis=1).astype(np.int32)
    else:
        XXu = np.zeros((npat, 0), dtype=np.int32)
    Q = XXf.sum(axis=1).astype(np.int32)
    if nf:
        Q = Q + XXu @ np.ones(nf, dtype=np.int32)
    best = int(np.max(np.abs(Q)))
    signs = np.ones(nf, dtype=np.int32)
    gray = 0
    total = 1 << nf
    for i in range(1, total):
        ng = i ^ (i >> 1)
        bit = (gray ^ ng).bit_length() - 1
        gray = ng
        os = int(signs[bit])
        signs[bit] = -os
        Q -= 2 * os * XXu[:, bit]
        mx = int(np.max(np.abs(Q)))
        if mx < best:
            best = mx
    return best


def random_sym_pm1(n: int, rng: Optional[np.random.Generator] = None) -> np.ndarray:
    if rng is None:
        rng = np.random.default_rng()
    U = rng.choice([-1.0, 1.0], size=(n, n))
    A = np.triu(U, 1)
    return A + A.T


def _is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def paley_conference_matrix(q: int) -> np.ndarray:
    """
    Symmetric conference matrix of order n=q+1 via Paley construction.
    Requires q prime, q ≡ 1 (mod 4).
    Satisfies C C^T = (n-1) I, zero diagonal, off-diagonal ±1.

    For prime-power fields F_{p^2} (orders n=p^2+1), use
    ``paley_conference_prime_power(p)`` instead.
    """
    if q % 4 != 1:
        raise ValueError("Paley conference requires q ≡ 1 (mod 4)")
    if not _is_prime(q):
        raise ValueError("q must be prime for this implementation; use paley_conference_prime_power for p^2")
    n = q + 1
    C = np.zeros((n, n), dtype=np.float64)
    C[0, 1:] = 1.0
    C[1:, 0] = 1.0
    for i in range(q):
        for j in range(i + 1, q):
            # Legendre(j-i mod q)
            val = pow((j - i) % q, (q - 1) // 2, q)
            s = 1.0 if val == 1 else -1.0
            C[i + 1, j + 1] = s
            C[j + 1, i + 1] = s
    return C


def paley_conference_prime_power(p: int) -> np.ndarray:
    """
    Paley conference matrix of order n = p^2 + 1 over F_{p^2}.

    Requires odd prime p. Encodes F_p(ω) elements as a + b*p with
    ω^2 = ia*ω + ib for a fixed irreducible quadratic. Quadratic character
    is the field multiplicative character of order 2.
    """
    if p < 3 or p % 2 == 0 or not _is_prime(p):
        raise ValueError("p must be an odd prime")
    q = p * p

    def is_irr(a: int, b: int) -> bool:
        return all((x * x - a * x - b) % p != 0 for x in range(p))

    ia = ib = None
    for a in range(p):
        for b in range(p):
            if is_irr(a, b):
                ia, ib = a, b
                break
        if ia is not None:
            break
    assert ia is not None

    def mul(u: int, v: int) -> int:
        c0, c1 = u % p, u // p
        d0, d1 = v % p, v // p
        e0 = (c0 * d0 + c1 * d1 * ib) % p
        e1 = (c0 * d1 + c1 * d0 + c1 * d1 * ia) % p
        return e0 + e1 * p

    def powm(u: int, e: int) -> int:
        r, base = 1, u
        while e:
            if e & 1:
                r = mul(r, base)
            base = mul(base, base)
            e >>= 1
        return r

    def chi(x: int) -> int:
        if x == 0:
            return 0
        return 1 if powm(x, (q - 1) // 2) == 1 else -1

    n = q + 1
    C = np.zeros((n, n), dtype=np.float64)
    C[0, 1:] = 1.0
    C[1:, 0] = 1.0
    for i in range(q):
        for j in range(i + 1, q):
            e0 = (i % p - j % p) % p
            e1 = (i // p - j // p) % p
            s = float(chi(e0 + e1 * p))
            C[i + 1, j + 1] = s
            C[j + 1, i + 1] = s
    return C


def halfspace_boolean_vector(p: int) -> np.ndarray:
    """
    Boolean eigenvector for Paley order n=p^2+1 (ρ=1 theorem).

    With F_p-linear form L(a+bω)=b and S={0,1,...,⌊p/2⌋}, set
    x_∞=1 and x_u = +1 iff L(u)∈S. Then C x = p x for the Paley
    conference matrix over F_{p^2} (see solution.md / HANDOFF.md).
    Encoding of field elements matches ``paley_conference_prime_power``.
    """
    if p < 3 or p % 2 == 0 or not _is_prime(p):
        raise ValueError("p must be an odd prime")
    q = p * p
    n = q + 1
    S = set(range((p + 1) // 2))  # 0 .. floor(p/2)
    x = np.ones(n, dtype=np.float64)
    for i in range(q):
        b = i // p  # L(a + b ω) with encoding a + b*p
        x[i + 1] = 1.0 if b in S else -1.0
    return x


def rho_from_boolean_evec(C: np.ndarray, x: np.ndarray) -> float:
    """
    If C x = λ x with x ∈ {±1}^n and C conference of order n, then
    |x^T C x| = |λ| n and ρ(C) = |x^T C x| / (n √(n-1)) = |λ| / √(n-1).
    For Paley p^2+1 one has λ=p=√(n-1), hence ρ=1.
    """
    n = C.shape[0]
    # residual of eigenvector equation
    Cx = C @ x
    # best scalar λ in least squares: λ = (x·Cx) / (x·x) = (x·Cx)/n
    lam = float(x @ Cx) / n
    residual = float(np.max(np.abs(Cx - lam * x)))
    if residual > 1e-6:
        # not an eigenvector; fall back to achieved cube/sphere ratio
        return abs(float(x @ Cx)) / (n * op_norm_sym(C))
    return abs(lam) / float(np.sqrt(n - 1))


def is_conference_matrix(C: np.ndarray, tol: float = 1e-8) -> bool:
    n = C.shape[0]
    # zero diagonal
    if not np.allclose(np.diag(C), 0, atol=tol):
        return False
    # off-diagonal ±1
    off = C[~np.eye(n, dtype=bool)]
    if not np.allclose(np.abs(off), 1.0, atol=tol):
        return False
    # C C^T = (n-1) I
    return np.allclose(C @ C.T, (n - 1) * np.eye(n), atol=tol * n)


def op_norm_sym(A: np.ndarray) -> float:
    ev = np.linalg.eigvalsh(A)
    return float(max(abs(ev[0]), abs(ev[-1])))


def spherical_half_bound(n: int) -> float:
    """(1/2) * n * sqrt(n-1): upper bound on m_n when a conference matrix of order n exists."""
    return 0.5 * n * np.sqrt(n - 1)


def alpha(m: float, n: int) -> float:
    return m / (n ** 1.5)


# --- Theoretical constants ---
# DMP/BH style lower bound constant claimed by Ivanisvili (comments on MO 413935)
DMP_LIminF_CONSTANT = 2 ** (-2.5)  # 2^{-5/2}
# Dual-Gaussian arcsine floor (Prop 5.2): m_n >= n sqrt(n-1)/pi for every Seidel A
DUAL_GAUSS_LIMINF = 1.0 / np.pi


def dmp_lower_bound(n: int) -> float:
    """Proven-shape lower bound m_n >= c n^{3/2} with c = 2^{-5/2} (author's calculus)."""
    return DMP_LIminF_CONSTANT * (n ** 1.5)


def dual_gaussian_lower_bound(n: int) -> float:
    """
    Prop 5.2: for every Seidel matrix Phi(A) >= n * sqrt(n-1) / pi.
    Hence m_n >= n * sqrt(n-1) / pi and liminf alpha_n >= 1/pi.
    """
    if n < 2:
        raise ValueError("n>=2")
    return n * float(np.sqrt(n - 1)) / float(np.pi)


def edge_hamming(A: np.ndarray, B: np.ndarray) -> int:
    """Number of undirected disagreeing edges between two Seidel matrices."""
    if A.shape != B.shape:
        raise ValueError("shape mismatch")
    # each undirected edge counted twice in the symmetric difference
    return int(np.sum(A != B) // 2)


def phi_edge_lipschitz_lower(phi_ref: float, k: int) -> float:
    """
    Prop 15.20b (edge-counting Lipschitz): if A and C differ in k edges then
    Phi(A) >= Phi(C) - 2k.

    Sharper than Frobenius Lipschitz (Prop 15.20: gap <= n sqrt(k)) whenever
    2k < n sqrt(k), i.e. k < n^2/4 (the sparse-disagreement regime).
    Consequence for E(1): if a Phi-minimiser is within k = o(n^{3/2}) edges of a
    rho=1 conference matrix (after switching), then m_n = Phi(C) - o(n^{3/2}).
    """
    if k < 0:
        raise ValueError("k>=0")
    return float(phi_ref) - 2.0 * float(k)


def phi_degree_lipschitz_lower(phi_ref: float, max_degree: int, n: int) -> float:
    """
    Prop 15.20c (degree Lipschitz): if the disagreement graph has max degree D
    then Phi(A) >= Phi(C) - D*n.

    Proof: E=A-C has row l1-norm <= 2D, hence ||E||_op <= 2D, so
    |x^T E x| <= 2 D n and |Q_A - Q_C| <= D n.
    For matchings D=1 this recovers gap <= n (same order as edge lip with k=n/2).
    """
    if max_degree < 0 or n < 2:
        raise ValueError("bad args")
    return float(phi_ref) - float(max_degree) * float(n)


def random_method_upper_bound(n: int) -> float:
    """
    Crude union-bound upper bound: for random A, with positive probability
    max_x |Q| <= sqrt(2 n^2 log(2^n)) = n^{3/2} sqrt(2 log 2), roughly.
    More carefully: each Q is Rademacher sum of binom(n,2) terms, subgauss with
    proxy binom(n,2); union bound over 2^{n-1} orbits gives
    m_n <= sqrt( n(n-1) log 2 ) * something — we use the standard
    m_n <= sqrt(log 2) * n^{3/2} as a convenient closed form that holds for n large
    (and we only use it as an O(n^{3/2}) envelope in tests).
    """
    return np.sqrt(np.log(2.0)) * (n ** 1.5)
