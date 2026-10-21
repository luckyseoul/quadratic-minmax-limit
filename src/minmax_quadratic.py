"""
Min-max quadratic form of ±1 coefficients.

m_n = min_{a_ij=±1} max_{x=±1} |sum_{i<j} a_ij x_i x_j|

Equivalently m_n = (1/2) min_A max_x |x^T A x| over symmetric zero-diagonal ±1 matrices.
"""
from __future__ import annotations

import math

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


def phi_mitm(A: np.ndarray) -> float:
    """
    Exact Phi(A)=max|Q| via meet-in-the-middle for even n.

    Fix x[0]=+1; split the remaining n-1 coordinates into two halves of size
    floor((n-1)/2) and ceil((n-1)/2). Enumerate both halves (2^{ceil((n-1)/2)}
    states) and combine. Feasible through n=26 (2^13 states per half); memory
    O(2^{n/2}).
    """
    n = A.shape[0]
    if n % 2 != 0:
        raise ValueError("phi_mitm requires even n")
    if n > 28:
        raise ValueError("phi_mitm only supported for n<=28")
    free = list(range(1, n))
    L = free[: len(free) // 2]
    R = free[len(free) // 2 :]
    L0 = [0] + L
    nL, nR = len(L), len(R)
    mR = 1 << nR
    bitsR = ((np.arange(mR)[:, None] >> np.arange(nR)[None, :]) & 1) * 2.0 - 1.0
    cR = 0.5 * np.einsum("ia,ab,ib->i", bitsR, A[np.ix_(R, R)], bitsR)
    mL = 1 << nL
    bitsL = ((np.arange(mL)[:, None] >> np.arange(nL)[None, :]) & 1) * 2.0 - 1.0
    xL0 = np.ones((mL, 1 + nL), dtype=np.float64)
    xL0[:, 1:] = bitsL
    cL = 0.5 * np.einsum("ia,ab,ib->i", xL0, A[np.ix_(L0, L0)], xL0)
    V = xL0 @ A[np.ix_(L0, R)]
    q = cL[:, None] + cR[None, :] + V @ bitsR.T
    return float(np.max(np.abs(q)))


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


def root_scale(height: float, order: int) -> float:
    """m / sqrt(n). Fekete would need this to be superadditive or subadditive."""
    if order < 1 or height < 0.0:
        raise ValueError("order>=1 and height>=0 required")
    return height / math.sqrt(order)


def fourth_moment_deletion_floor(deleted: int) -> float:
    """Lower bound c4*sqrt(binom(t,2)) on m_n - m_{n-t} for t>=2."""
    if deleted < 2:
        raise ValueError("deleted>=2 required")
    coefficient = (math.sqrt(18.0) - math.sqrt(14.0)) / 2.0
    return coefficient * math.sqrt(deleted * (deleted - 1) / 2.0)


def half_vertex_deletion_floor(deleted: int) -> int:
    """Host-free deletion floor floor(t/2) on m_n - m_{n-t}."""
    if deleted < 2:
        raise ValueError("deleted>=2 required")
    return deleted // 2


def hypercontractive_deletion_floor(deleted: int) -> float:
    """Host-free deletion floor sqrt(binom(t,2))/(2 e^2)."""
    if deleted < 2:
        raise ValueError("deleted>=2 required")
    return math.sqrt(deleted * (deleted - 1) / 2.0) / (2.0 * math.e * math.e)


def fractional_half_cut_gap(height: float) -> float:
    """Distance from the convex half-cut value to the integral target sqrt(2)*height."""
    if height < 0.0:
        raise ValueError("height>=0 required")
    return (math.sqrt(2.0) - 1.0) * height


def half_cut_witness_slack() -> float:
    """Additive drop forced by one anti-aligned near-maximizer."""
    return 2.0


def neutral_mass_increment(level: float, order: int, deleted: int) -> float:
    """m-increment that holds alpha fixed across a deletion of this size."""
    if order <= deleted or deleted < 1 or level < 0.0:
        raise ValueError("order>deleted>=1 and level>=0 required")
    return level * (order ** 1.5 - (order - deleted) ** 1.5)


def partition_extension_ceiling(order: int, lam: float) -> float:
    """Full-cube form of the partition extension bound.

    The lemma charges n log cosh(lambda) + log(2 S). With S = 2^n this is
    an explicit upper bound on the one-vertex increment the lemma can certify.
    """
    if order < 2 or lam <= 0.0:
        raise ValueError("order>=2 and lambda>0 required")
    return (order * math.log(math.cosh(lam)) + (order + 1) * math.log(2.0)) / lam


def partition_extension_full_cube_min(order: int) -> float:
    """Minimum of the full-cube partition ceiling over lambda in (0, 4]."""
    if order < 2:
        raise ValueError("order>=2 required")
    best = None
    for step in range(1, 4001):
        value = partition_extension_ceiling(order, step / 1000.0)
        if best is None or value < best:
            best = value
    return float(best)


def odd_order_forced_increment() -> float:
    """Smallest positive Boolean dot product on an odd order.

    A dot product of two ±1 vectors of odd length is a nonzero odd integer,
    so a flat step m_{n+1}=m_n is impossible when n is odd.
    """
    return 1.0


def neutral_one_step(level: float, order: int) -> float:
    """m_n*((1+1/n)^{3/2}-1), the increment that holds alpha fixed."""
    if order < 1 or level < 0.0:
        raise ValueError("order>=1 and level>=0 required")
    return level * (order ** 1.5) * ((1.0 + 1.0 / order) ** 1.5 - 1.0)


def deletion_shortfall_density(order: int, level: float) -> float:
    """Normalized gap between the neutral two-step increment and the floor."""
    if order < 4 or level < 0.0:
        raise ValueError("order>=4 and level>=0 required")
    gap = neutral_mass_increment(level, order, 2) - fourth_moment_deletion_floor(2)
    return gap / (order ** 1.5)


def doubling_image(height: float) -> float:
    """2*sqrt(2)*height, the zero-error double of one block norm."""
    if height < 0.0:
        raise ValueError("height>=0 required")
    return 2.0 * math.sqrt(2.0) * height


def relaxation_ceiling_exceeds_doubling(pair_cost: int, height: int) -> bool:
    """Whether (3/2)*pair_cost sits strictly above 2*sqrt(2)*height.

    Squares remove the root: 9 pair_cost^2 > 32 height^2.
    """
    if pair_cost < 0 or height < 0:
        raise ValueError("pair_cost>=0 and height>=0 required")
    return 9 * pair_cost * pair_cost > 32 * height * height


def bivector_elliptope_score(order: int) -> float:
    """Critical elliptope ratio bound pi^2 / (4 (sqrt(2)-1)^2 (n-1))."""
    if order < 2:
        raise ValueError("order>=2 required")
    gap = math.sqrt(2.0) - 1.0
    return (math.pi * math.pi) / (4.0 * gap * gap * (order - 1))


def bivector_critical_excess(order: int) -> float:
    """(sqrt(2)-1) times the Gaussian floor, divided by n^{3/2}."""
    if order < 2:
        raise ValueError("order>=2 required")
    return (math.sqrt(2.0) - 1.0) * dual_gaussian_lower_bound(order) / (order ** 1.5)


def degree_four_width_dominates(order: int) -> bool:
    """Rational certificate that critical widths beat degree-four row support.

    Uses sqrt(2)<5657/4000 and pi^2<10. The resulting comparison
    1372*(n-1) > 60000 holds for every n>=45 and fails for n<=44.
    """
    if order < 2:
        raise ValueError("order>=2 required")
    return 1372 * (order - 1) > 60000


def endpoint_annihilator_signing(level: int) -> np.ndarray:
    """Complete signing A from the replicated near-annihilator family.

    level ell>=1 gives order n=2*4^ell. The Boolean norm is n*d/2 with
    d=2*sqrt(4^ell)+1, which is above the conference cap.
    """
    if level < 1:
        raise ValueError("level>=1 required")
    block = 2.0 * np.eye(4) - np.ones((4, 4))
    hadamard = block
    for _ in range(level - 1):
        hadamard = np.kron(hadamard, block)
    width = 4 ** level
    paired = np.kron(hadamard, np.ones((2, 2)))
    return paired - np.eye(2 * width)


def endpoint_annihilator_height(level: int) -> float:
    """Boolean norm n*d/2 of the replicated near-annihilator."""
    if level < 1:
        raise ValueError("level>=1 required")
    width = 4 ** level
    order = 2 * width
    singular = 2.0 * math.sqrt(width) + 1.0
    return order * singular / 2.0


def boolean_operator_product(order: int, height: float) -> float:
    """Prop 15.6 product rho(A) * ||A||_op = 2 * Phi(A) / n.

    The product is fixed by the Boolean height. It does not take the
    operator norm as an input.
    """
    if order < 2:
        raise ValueError("order>=2 required")
    if height < 0:
        raise ValueError("height>=0 required")
    return 2.0 * height / order


def character_block_indices(q: int) -> tuple[np.ndarray, np.ndarray]:
    """Infinity-squares split for a prime q ≡ 1 (mod 4).

    Index 0 is infinity and index a+1 is the field element a. The second
    block is 0 followed by nu times the same squares.
    """
    if q < 5 or q % 4 != 1 or not _is_prime(q):
        raise ValueError("q>=5 prime and q≡1 mod 4 required")
    squares = sorted({(r * r) % q for r in range(1, q)})
    square_set = set(squares)
    nu = next(element for element in range(1, q) if element not in square_set)
    residues = np.asarray(squares, dtype=int)
    squares_block = np.concatenate(([0], residues + 1))
    nonsquares_block = np.concatenate(([1], (nu * residues) % q + 1))
    return squares_block, nonsquares_block


def character_block_gram_gap(q: int) -> float:
    """Max entry of A^2 + B B^T - q I on the character blocks.

    C^2 = q I and the two-block form force the gap to vanish.
    """
    conference = paley_conference_matrix(q)
    squares_block, nonsquares_block = character_block_indices(q)
    block = conference[np.ix_(squares_block, squares_block)]
    cross = conference[np.ix_(squares_block, nonsquares_block)]
    gram = block @ block + cross @ cross.T
    eye = q * np.eye(block.shape[0])
    return float(np.max(np.abs(gram - eye)))


def character_block_height_cap(q: int) -> float:
    """Spectral height cap (1/2)*k*sqrt(q) on the character block.

    k=(q+1)/2. Positive-semidefiniteness of B B^T and A^2+B B^T=q I
    give ||A||_op <= sqrt(q), hence Phi(A) <= (1/2)*k*sqrt(q).
    """
    if q < 5 or q % 4 != 1:
        raise ValueError("q>=5 and q≡1 mod 4 required")
    block_order = (q + 1) // 2
    return 0.5 * block_order * math.sqrt(q)


def character_block_cap_ratio(q: int) -> float:
    """Block cap divided by the spherical cap on the same block order."""
    block_order = (q + 1) // 2
    return character_block_height_cap(q) / spherical_half_bound(block_order)


def chsh_scalar_bilinear() -> float:
    """Maximum of |x^T H y| over Boolean vectors for the CHSH matrix.

    H is [[1, 1], [1, -1]]. The tensorization audit uses this value as
    the scalar bilinear norm.
    """
    matrix = np.array([[1.0, 1.0], [1.0, -1.0]])
    best = 0.0
    for left in (-1.0, 1.0):
        for right in (-1.0, 1.0):
            row = np.array([left, right])
            for first in (-1.0, 1.0):
                for second in (-1.0, 1.0):
                    column = np.array([first, second])
                    best = max(best, abs(float(row @ matrix @ column)))
    return best


def chsh_vector_value() -> float:
    """CHSH vector value on the same matrix.

    The first pair is the standard basis. The second pair is that basis
    rotated by a quarter turn. The pairing sum is H_ij times the inner
    product of the i-th and j-th vectors.
    """
    matrix = np.array([[1.0, 1.0], [1.0, -1.0]])
    left = np.array([[1.0, 0.0], [0.0, 1.0]])
    scale = math.sqrt(2.0)
    right = np.array([[1.0 / scale, 1.0 / scale], [1.0 / scale, -1.0 / scale]])
    total = 0.0
    for i in range(2):
        for j in range(2):
            total += float(matrix[i, j]) * float(left[i] @ right[j])
    return total


def random_padding_increment(order: int, added: int) -> float:
    """Hoeffding union-bound increment for filling ``added`` new vertices.

    E = added*(2*order+added-1)/2 new edges and N = order+added. The
    existence bound is m_N <= m_order + sqrt(2*E*(N*log(2)+1)).
    """
    if order < 2 or added < 1:
        raise ValueError("order>=2 and added>=1 required")
    edges = added * (2 * order + added - 1) / 2.0
    total = order + added
    return math.sqrt(2.0 * edges * (total * math.log(2.0) + 1.0))


def lexicographic_blowup(seed: np.ndarray, block: int) -> np.ndarray:
    """Lexicographic blow-up with constant cross blocks and positive internals.

    Each vertex of the seed becomes ``block`` vertices. Distinct seed
    vertices contribute the rank-one block seed_ij * J. Each diagonal
    block is the Seidel matrix with off-diagonal entries +1.
    """
    host = np.asarray(seed, dtype=np.float64)
    seed_order = int(host.shape[0])
    if block < 2 or host.shape != (seed_order, seed_order):
        raise ValueError("block>=2 and a square seed are required")
    order = seed_order * block
    blown = np.zeros((order, order), dtype=np.float64)
    ones = np.ones((block, block), dtype=np.float64)
    internal = ones - np.eye(block)
    for i in range(seed_order):
        for j in range(seed_order):
            rows = slice(i * block, (i + 1) * block)
            cols = slice(j * block, (j + 1) * block)
            if i == j:
                blown[rows, cols] = internal
            else:
                blown[rows, cols] = host[i, j] * ones
    return blown


def lexicographic_blowup_ceiling(seed_height: float, seed_order: int, block: int) -> float:
    """Upper bound r^2*Phi(A) + n*spherical_half_bound(r) from the blow-up estimate.

    The internal blocks contribute at most n times the spherical cap.
    """
    if seed_order < 2 or block < 2 or seed_height < 0.0:
        raise ValueError("seed_order>=2, block>=2, height>=0 required")
    return (block * block) * seed_height + seed_order * spherical_half_bound(block)


def lexicographic_blowup_floor(seed_height: float, seed_order: int, block: int) -> float:
    """Lower estimate r^2*Phi(A) - n*spherical_half_bound(r).

    It sits under the constant-block construction because m_r is at most
    the spherical cap.
    """
    if seed_order < 2 or block < 2 or seed_height < 0.0:
        raise ValueError("seed_order>=2, block>=2, height>=0 required")
    return (block * block) * seed_height - seed_order * spherical_half_bound(block)


def walsh_discrepancy_radius(order: int) -> float:
    """Hereditary discrepancy radius of the degree-two Walsh encoding.

    The encoding uses E=binom(n,2) diagonal tests. The displayed
    Akbas--Sra shape, without the universal factor kappa, is
    sqrt(E) * (1 + log(2^n / E))^2. The logarithm is natural.
    """
    if order < 2:
        raise ValueError("order>=2 required")
    edges = order * (order - 1) / 2.0
    log_ratio = order * math.log(2.0) - math.log(edges)
    return math.sqrt(edges) * (1.0 + log_ratio) ** 2


def operator_alignment(order: int, height: float, op_norm: float) -> float:
    """Prop 15.7 ratio rho(A) * ||A||_op / sqrt(n-1).

    rho is the Prop 15.6 quotient of the Boolean product by op_norm.
    Every legal operator norm, meaning op_norm^2 >= n-1 as in Prop 15.5,
    returns the same value 2*height / (n*sqrt(n-1)). A beater of a
    rho=1 conference matrix is exactly the case where this value is
    strictly less than 1.
    """
    if order < 2:
        raise ValueError("order>=2 required")
    if height < 0:
        raise ValueError("height>=0 required")
    floor = math.sqrt(order - 1)
    if op_norm < floor:
        raise ValueError("op_norm is below the Prop 15.5 floor")
    rho = boolean_operator_product(order, height) / op_norm
    return rho * op_norm / floor


def order_four_optimal_signing() -> np.ndarray:
    """Globally minimal order-four signing with one negative edge."""
    return np.array(
        [
            [0.0, 1.0, 1.0, 1.0],
            [1.0, 0.0, -1.0, 1.0],
            [1.0, -1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0, 0.0],
        ],
        dtype=np.float64,
    )


def clique_flip_norm(signing: np.ndarray) -> float:
    """Max Boolean norm over reversals of the edges inside one vertex set."""
    host = np.asarray(signing, dtype=np.float64)
    order = int(host.shape[0])
    if order < 1 or order > 8 or host.shape != (order, order):
        raise ValueError("order must be between 1 and 8")
    best = 0.0
    for mask in range(1 << order):
        flipped = host.copy()
        members = [index for index in range(order) if mask & (1 << index)]
        for left in range(len(members)):
            for right in range(left + 1, len(members)):
                i, j = members[left], members[right]
                flipped[i, j] *= -1.0
                flipped[j, i] *= -1.0
        best = max(best, phi(flipped))
    return float(best)


def ordered_skew_range_ceiling(height: float, order: int) -> float:
    """Largest ordered-skew cross term allowed by prefix energies of size at most height.

    Each of the order steps can move Q by at most twice the height, and the
    cross term is half of that total variation, so the ceiling is order*height.
    """
    if order < 2 or height < 0.0:
        raise ValueError("order>=2 and height>=0 required")
    return float(order) * height


def _antipodal_states(order: int) -> np.ndarray:
    count = 1 << (order - 1)
    states = np.ones((count, order), dtype=np.float64)
    masks = np.arange(count)
    for index in range(1, order):
        states[:, index] = np.where((masks >> (index - 1)) & 1, -1.0, 1.0)
    return states


def ordered_skew_total_variation(signing: np.ndarray) -> float:
    """Max of (1/2) sum |q_k-q_{k-1}| over index-order prefix flips."""
    host = np.asarray(signing, dtype=np.float64)
    order = int(host.shape[0])
    if order < 2 or order > 16 or host.shape != (order, order):
        raise ValueError("order must be between 2 and 16")
    current = _antipodal_states(order)
    total = np.zeros(current.shape[0], dtype=np.float64)
    for index in range(order):
        field = current @ host
        total += np.abs(-2.0 * current[:, index] * field[:, index])
        current[:, index] *= -1.0
    return float(0.5 * np.max(total))


def ordered_skew_cut_norm(signing: np.ndarray) -> float:
    """Max over Boolean states of ||R^prec y||_1."""
    host = np.asarray(signing, dtype=np.float64)
    order = int(host.shape[0])
    if order < 2 or order > 16 or host.shape != (order, order):
        raise ValueError("order must be between 2 and 16")
    upper = np.triu(host, 1)
    skew = upper - upper.T
    images = _antipodal_states(order) @ skew.T
    return float(np.max(np.sum(np.abs(images), axis=1)))


def sum_of_maxima_floor(height: float, block: int) -> float:
    """Two block norms plus the universal ±1 cross floor on that block."""
    if block < 1 or height < 0.0:
        raise ValueError("block>=1 and height>=0 required")
    return 2.0 * height + rectangle_switching_lower(block, block)


def aligned_extreme_cross_budget(height: float) -> float:
    """Largest cross term that keeps 2*height + cross inside the doubling image."""
    if height < 0.0:
        raise ValueError("height>=0 required")
    return 2.0 * (math.sqrt(2.0) - 1.0) * height


def ancestry_mass_threshold(level: float, order: int) -> float:
    """2 to the neutral increment. A compatibility tree under this mass keeps positive width."""
    if order < 1 or level < 0.0:
        raise ValueError("order>=1 and level>=0 required")
    return 2.0 ** neutral_one_step(level, order)


def even_order_stable_separation() -> int:
    """Minimum Hamming distance of distinct same-orientation stable states at even order.

    Every local field is an odd integer, so a stable margin is at least 1 and
    the separation is at least 2.
    """
    return 2


def skeleton_ball_volume(order: int, radius: int) -> int:
    """Worst-case count of other stable states inside a Hamming ball of this radius.

    The light set in the skeleton bound can be the full vertex set.
    """
    if order < 1 or radius < 0:
        raise ValueError("order>=1 and radius>=0 required")
    total = 0
    for distance in range(1, min(radius, order) + 1):
        total += math.comb(order, distance)
    return total


def row_tilt_fourth_moment_budget(coupling: float, penalty: float) -> float:
    """Per-row weighted fourth-moment ceiling 2 exp(c^2/2+4λ)(c^2+8λ)."""
    if coupling < 0.0 or penalty < 0.0:
        raise ValueError("coupling>=0 and penalty>=0 required")
    return (
        2.0
        * math.exp(0.5 * coupling * coupling + 4.0 * penalty)
        * (coupling * coupling + 8.0 * penalty)
    )


def row_tilt_window_charge(start: int, coupling: float, penalty: float) -> float:
    """Sum of the row-summed budget divided by k^{3/2} across one doubling window."""
    if start < 2:
        raise ValueError("start>=2 required")
    budget = row_tilt_fourth_moment_budget(coupling, penalty)
    total = 0.0
    for order in range(start, 2 * start):
        total += (order * budget) / (order ** 1.5)
    return total


def equipartition_product_floor(order: int) -> float:
    """Universal product floor (4/9) m^2 mu_m^2 with m = floor(order/2)."""
    if order < 2:
        raise ValueError("order>=2 required")
    half = order // 2
    moment = rademacher_absolute_moment(half)
    return (4.0 / 9.0) * (half * half) * (moment * moment)


def bootstrap_interior_root(tilt: float, weight: float) -> float:
    """Largest root of the bootstrapped interior equality at these parameters.

    The root is a lower bound for liminf alpha_n. It is not a limit.
    """
    if not 0.0 < tilt <= 1.0 or not 0.0 < weight <= 0.5:
        raise ValueError("tilt in (0,1] and weight in (0,1/2] required")
    kappa = 2.0 / math.pi
    zed = tilt * tilt / (1.0 + tilt * tilt)
    aye = 1.0 - kappa * math.asin(zed) + kappa * zed
    energy = kappa * tilt / (1.0 + tilt * tilt)
    field = math.sqrt(kappa * aye)
    beta = aye - kappa
    width = 2.0 * math.sqrt(kappa) * tilt / (1.0 + tilt * tilt)
    height = min(0.25, math.atan(math.sqrt(beta) / width) / math.pi)
    response = (1.0 - weight) ** 2 * energy + weight * (1.0 - weight) * field
    star = 2.0 / (27.0 * math.pi)

    def residual(level: float) -> float:
        left = ((1.0 + weight * weight) * level - response) * (
            level * level - star * (1.0 - height) ** 3
        )
        right = 4.0 * star * weight * weight * level * height ** 3
        return left - right

    low = 0.2
    high = 0.5
    if residual(low) >= 0.0 or residual(high) <= 0.0:
        raise ValueError("bootstrapped equality has no root in [0.2, 0.5]")
    for _ in range(80):
        mid = 0.5 * (low + high)
        if residual(mid) < 0.0:
            low = mid
        else:
            high = mid
    return 0.5 * (low + high)


def coupled_half_loss(height: float, half: int) -> float:
    """Balanced coupled loss [F - sqrt(F^2 - 4 Gamma)]/2."""
    if half < 2 or height <= 0.0:
        raise ValueError("half>=2 and height>0 required")
    product = equipartition_product_floor(half)
    discriminant = height * height - 4.0 * product
    if discriminant <= 0.0:
        raise ValueError("product floor exceeds the square of the host norm")
    return 0.5 * (height - math.sqrt(discriminant))


def two_thirds_scale(height: float) -> float:
    """H = m^{2/3}. Subadditivity of H is the zero-error form of a Fekete bound."""
    if height < 0.0:
        raise ValueError("height>=0 required")
    return height ** (2.0 / 3.0)


def step_equilibrium_level(step: float, order: int) -> float:
    """Fixed point 2K/3 of the recurrence with increment at most K sqrt(n).

    K is step/sqrt(order). The conference cap is 1/2, so a level above 1/2
    does not shrink limsup alpha_n.
    """
    if order < 1 or step < 0.0:
        raise ValueError("order>=1 and step>=0 required")
    return (2.0 / 3.0) * (step / math.sqrt(order))


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


def arcsine_saturation_level(order: int) -> float:
    """Arcsine Gaussian value L_n = n(n-1)/pi * asin(1/sqrt(n-1))."""
    if order < 2:
        raise ValueError("order>=2 required")
    return (
        order
        * (order - 1)
        / math.pi
        * math.asin(1.0 / math.sqrt(order - 1))
    )


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


def rademacher_absolute_moment(width: int) -> float:
    """Mean absolute value of a width-term Rademacher sum."""
    if width < 1:
        raise ValueError("width>=1")
    half = (width - 1) // 2
    log_moment = (
        math.log(width)
        + math.lgamma(width)
        - math.lgamma(half + 1)
        - math.lgamma(width - half)
        - (width - 1) * math.log(2.0)
    )
    return math.exp(log_moment)


def rectangle_switching_lower(rows: int, cols: int) -> float:
    """Boolean cut-norm lower bound rows * mu_cols for a sign rectangle."""
    if rows < 1 or cols < 1:
        raise ValueError("bad args")
    return float(rows) * rademacher_absolute_moment(cols)


def quadratic_variance_chernoff_threshold() -> float:
    """Ideal central-cut level from 12 alpha^2 = log 4."""
    return math.sqrt(math.log(2.0) / 6.0)


def cut_discrepancy_window(n: int) -> float:
    """Normalized slack (1 + sqrt(binom(n, 2))) / n^{3/2} from the cut-code comparison."""
    if n < 2:
        raise ValueError("n>=2")
    return (1.0 + math.sqrt(n * (n - 1) / 2.0)) / (n ** 1.5)


def best_response_mean_energy(signing: np.ndarray) -> float:
    """Mean of Q at the synchronous best response, tie value +1."""
    order = int(signing.shape[0])
    if order < 2 or order > 14 or signing.shape != (order, order):
        raise ValueError("order must be between 2 and 14")
    masks = np.arange(1 << order, dtype=np.uint32)[:, None]
    bits = (masks >> np.arange(order, dtype=np.uint32)) & 1
    states = np.where(bits == 1, 1.0, -1.0)
    fields = states @ signing
    response = np.where(fields >= 0.0, 1.0, -1.0)
    image = response @ signing
    return float(np.mean(0.5 * np.sum(response * image, axis=1)))


def pushforward_influence_ceiling(correlation: float) -> float:
    """Coarse best-response ceiling 2 + 2 sqrt(1 - |correlation|)."""
    magnitude = abs(correlation)
    if magnitude > 1.0:
        raise ValueError("correlation must lie in [-1, 1]")
    return 2.0 + 2.0 * math.sqrt(1.0 - magnitude)


def independent_skew_ideal_threshold() -> float:
    """Level sqrt(log 2 / 2) where log 4 - 4 alpha^2 changes sign."""
    return math.sqrt(math.log(2.0) / 2.0)


def independent_skew_cap_exponent() -> float:
    """First-moment rate log 4 - 1 at alpha = 1/2."""
    return math.log(4.0) - 1.0


def uniform_partition_threshold() -> float:
    """Level sqrt(log 2) / (sqrt(2) - 1) from zero near-extremal entropy cost."""
    return math.sqrt(math.log(2.0)) / (math.sqrt(2.0) - 1.0)


def directed_half_floor_constant() -> float:
    """Leading constant 1/(pi sqrt(2)) of M/sqrt(2) at the Gaussian floor."""
    return 1.0 / (math.pi * math.sqrt(2.0))


def zero_defect_spectral_bridge() -> float:
    """Leading constant 1/(2 sqrt(2)) of the zero-defect operator bound on D_to."""
    return 1.0 / (2.0 * math.sqrt(2.0))


def anticommutator_frobenius_floor(order: int) -> float:
    """Parity floor on ||AR+RA||_F^2 for sign pairs of order at least 3."""
    if order < 3:
        raise ValueError("order>=3")
    if order % 2 == 0:
        return float(2 * order * (order - 2))
    return float(2 * (order - 1) * (order - 1))


def skew_tournament(order: int) -> np.ndarray:
    """Upper-triangular +1 skew signing."""
    if order < 2:
        raise ValueError("order>=2")
    signs = np.triu(np.ones((order, order), dtype=np.float64), 1)
    return signs - signs.T


def nearly_regular_tournament(order: int) -> np.ndarray:
    """Regular tournament, or a one-vertex deletion of one when order is even."""
    if order < 2:
        raise ValueError("order>=2")
    host = order if order % 2 == 1 else order + 1
    tournament = np.zeros((host, host), dtype=np.float64)
    half = host // 2
    for source in range(host):
        for step in range(1, half + 1):
            target = (source + step) % host
            tournament[source, target] = 1.0
            tournament[target, source] = -1.0
    if order % 2 == 1:
        return tournament
    keep = np.arange(1, host)
    return tournament[np.ix_(keep, keep)]


def tournament_all_one_l1(tournament: np.ndarray) -> float:
    """L1 norm of the image of the all-ones vector."""
    return float(np.sum(np.abs(tournament.sum(axis=1))))


def directed_triangle_count(tournament: np.ndarray) -> int:
    """Number of directed 3-cycles in a tournament."""
    order = int(tournament.shape[0])
    count = 0
    for i in range(order):
        for j in range(i + 1, order):
            for k in range(j + 1, order):
                score = (
                    tournament[i, j] + tournament[j, k] + tournament[k, i]
                )
                if abs(score) == 3.0:
                    count += 1
    return count


def pentagon_signing() -> np.ndarray:
    """Order-5 signing, negative on the cycle and positive on the diagonals."""
    order = 5
    signing = np.ones((order, order), dtype=np.float64)
    np.fill_diagonal(signing, 0.0)
    for vertex in range(order):
        neighbor = (vertex + 1) % order
        signing[vertex, neighbor] = -1.0
        signing[neighbor, vertex] = -1.0
    return signing


def cosh_pressure(signing: np.ndarray, coupling: float) -> float:
    """Mean of cosh(coupling * Q / sqrt(n)) over the Boolean cube."""
    order = int(signing.shape[0])
    if order < 2 or order > 12:
        raise ValueError("order must be between 2 and 12")
    scale = coupling / math.sqrt(order)
    total = 0.0
    patterns = 1 << (order - 1)
    for mask in range(patterns):
        state = np.ones(order, dtype=np.float64)
        for index in range(order - 1):
            state[index + 1] = 1.0 if mask & (1 << index) else -1.0
        total += math.cosh(scale * form_Q(signing, state))
    return total / patterns


def product_pressure_lower_curve(order: int, coupling: float) -> float:
    """The proposed universal curve cosh(c sqrt(1-1/n))^(n/2)."""
    if order < 2:
        raise ValueError("order>=2")
    return math.cosh(coupling * math.sqrt(1.0 - 1.0 / order)) ** (order / 2.0)


def commutator_square_ceiling(order: int) -> float:
    """Elementary ceiling of Sigma(A,R)/n^4 for sign pairs."""
    if order < 3:
        raise ValueError("order>=3")
    ratio = 1.0 - 2.0 / order
    return 8.0 * (1.0 - 1.0 / order) * ratio * ratio


def arcsine_rigidity_allowance(level: float) -> float:
    """Leading allowance 64 pi (alpha - 1/pi) in the mate-rigidity inequality."""
    return 64.0 * math.pi * (level - 1.0 / math.pi)


def pressure_comparison_window(coupling: float) -> float:
    """The (log 2)/c slack between s_n(c)/c and alpha_n."""
    if coupling <= 0.0:
        raise ValueError("coupling must be positive")
    return math.log(2.0) / coupling


def critical_interpolation_defect(coupling: float) -> float:
    """Alpha-scale size c/8 of the common-temperature doubling defect."""
    if coupling <= 0.0:
        raise ValueError("coupling must be positive")
    return coupling / 8.0


def critical_block_excess(order: int, coupling: float) -> float:
    """The normalized cross term (n/2) log cosh(c/sqrt(2n))."""
    if order < 2 or coupling <= 0.0:
        raise ValueError("bad args")
    return (order / 2.0) * math.log(math.cosh(coupling / math.sqrt(2.0 * order)))


def fixed_beta_entropy_window(order: int, beta: float) -> float:
    """Normalized ceiling (log 2)/(beta sqrt(n)) on the fixed-beta pressure gap."""
    if order < 2 or beta <= 0.0:
        raise ValueError("order>=2 and beta>0 required")
    return math.log(2.0) / (beta * math.sqrt(order))


def fixed_beta_block_defect(order: int, beta: float) -> float:
    """Doubling cross term n^2 log cosh(beta), divided by n^{3/2}."""
    if order < 2 or beta <= 0.0:
        raise ValueError("order>=2 and beta>0 required")
    return math.sqrt(order) * math.log(math.cosh(beta))


def one_step_critical_pressure_sum(start: int, coupling: float) -> float:
    """Sum of the one-step cosh defects from start through one doubling."""
    if start < 2 or coupling <= 0.0:
        raise ValueError("start>=2 and coupling>0 required")
    total = 0.0
    for index in range(start, 2 * start):
        step = coupling / math.sqrt(index + 1.0)
        total += (index / (index + 1.0)) * math.log(math.cosh(step))
    return total


def fixed_beta_pressure_gap(signing: np.ndarray, beta: float) -> float:
    """beta*Phi(A) minus the log-mean of cosh(beta Q) on the half-cube.

    For every signing the minimum of that log-mean lies between
    beta*m_n - n log 2 and beta*m_n. The returned gap is the loss on this
    one signing.
    """
    order = int(signing.shape[0])
    if order < 2 or order > 16 or signing.shape != (order, order) or beta <= 0.0:
        raise ValueError("square signing of order 2..16 and beta>0 required")
    count = 1 << (order - 1)
    states = np.ones((count, order), dtype=np.float64)
    masks = np.arange(count)
    for index in range(1, order):
        states[:, index] = np.where((masks >> (index - 1)) & 1, -1.0, 1.0)
    energy = 0.5 * np.einsum("ia,ab,ib->i", states, signing, states)
    height = float(np.max(np.abs(energy)))
    mean_cosh = float(np.mean(np.cosh(beta * energy)))
    return beta * height - math.log(mean_cosh)


def tetrahedral_union_bound_level() -> float:
    """Independent-skew first-moment level sqrt(18 log 2).

    A per-state union bound on the tetrahedral profile needs this level
    before the expected number of failures drops below one.
    """
    return math.sqrt(18.0 * math.log(2.0))


def tetrahedral_budget_level() -> float:
    """Leading multiplier-three budget 3^{3/2}/2 on the normalized scale."""
    return (3.0 ** 1.5) / 2.0


def bibalanced_gaussian_headroom() -> float:
    """Limit of (T - M - N)/n^{3/2} at the Gaussian floor.

    T = (M^{2/3} + N^{2/3})^{3/2} with M and N the floors at n and 2n.
    """
    return (3.0 * math.sqrt(3.0) - 1.0 - 2.0 * math.sqrt(2.0)) / math.pi


def bibalanced_strip_cross(fraction: float) -> float:
    """Leading coefficient 4*sqrt(fraction) of the paired cross majorant."""
    if fraction < 0.0:
        raise ValueError("fraction>=0")
    return 4.0 * math.sqrt(fraction)


def bibalanced_border_coefficient(order: int) -> float:
    """Normalized size 6*sqrt(n)/n^{3/2} of the unpaired border."""
    if order < 2:
        raise ValueError("order>=2")
    return 6.0 / math.sqrt(order)


def bibalanced_composition_headroom(order: int) -> float:
    """(T - M - N)/n^{3/2} at the Gaussian floors of orders n and 2n."""
    if order < 2:
        raise ValueError("order>=2")
    small = dual_gaussian_lower_bound(order)
    large = dual_gaussian_lower_bound(2 * order)
    composed = (small ** (2.0 / 3.0) + large ** (2.0 / 3.0)) ** 1.5
    return (composed - small - large) / (order ** 1.5)


def doubling_only_oscillator(x: float, amplitude: float) -> float:
    """Scalar sharpness example for the two-ray criterion.

    G(x)=x(1+amplitude*sin(2π log2(x))) satisfies G(2x)=2G(x). For a
    sufficiently small positive amplitude the slope G(x)/x keeps oscillating,
    so exact doubling does not force the normalized sequence to converge.
    """
    if x <= 0.0 or not 0.0 < amplitude < 1.0:
        raise ValueError("x>0 and amplitude in (0,1)")
    phase = 2.0 * math.pi * math.log2(x)
    return x * (1.0 + amplitude * math.sin(phase))


def paley_gamma_envelope(order: int) -> float:
    """Gaussian upper bound on the Paley gap γ at a conference order.

    It is the difference of the spherical and Gaussian normalized floors.
    """
    if order < 2:
        raise ValueError("order>=2")
    return alpha(spherical_half_bound(order), order) - alpha(
        dual_gaussian_lower_bound(order), order
    )


def linear_undercut_gamma(prime: int) -> float:
    """Normalized gap (p-1)/n^{3/2} at the Paley order n=p^2+1."""
    if prime < 3:
        raise ValueError("prime>=3")
    order = prime * prime + 1
    return (prime - 1) / (order ** 1.5)


def min_single_skew_tetrahedral_energy(signing: np.ndarray) -> float:
    """Minimum of K_3(A; R, -R, R) over skew signings R.

    Order 2..5. The value is the exact single-skew diamond
    max |Q(x)+Q(y)-Q(z)| + |x^T R y + y^T R z + z^T R x|.
    """
    order = int(signing.shape[0])
    if order < 2 or order > 5 or signing.shape != (order, order):
        raise ValueError("order must be between 2 and 5")
    count = 1 << order
    states = np.ones((count, order), dtype=np.float64)
    masks = np.arange(count)
    for index in range(order):
        states[:, index] = np.where((masks >> index) & 1, -1.0, 1.0)
    energy = 0.5 * np.einsum("ia,ab,ib->i", states, signing, states)
    edges = [(i, j) for i in range(order) for j in range(i + 1, order)]
    best = None
    for mask in range(1 << len(edges)):
        skew = np.zeros((order, order), dtype=np.float64)
        for bit, (i, j) in enumerate(edges):
            sign = 1.0 if mask & (1 << bit) else -1.0
            skew[i, j] = sign
            skew[j, i] = -sign
        gram = states @ skew @ states.T
        joint = (
            energy[:, None, None]
            + energy[None, :, None]
            - energy[None, None, :]
        )
        cross = gram[:, :, None] + gram[None, :, :] + gram.T[:, None, :]
        value = float(np.max(np.abs(joint) + np.abs(cross)))
        if best is None or value < best:
            best = value
    return float(best)


def approximate_rank_obstruction(signing: np.ndarray, height: float) -> float:
    """How far the approximate-rank certificate misses the neutral increment.

    The certificate asks for a dimension d with d + epsilon_d < r_n.
    Zero-deficit stable states give a lower bound on epsilon_d through the
    optimal average squared distance to a d-dimensional subspace. The return
    value is the smallest such lower bound of d + epsilon_d, minus r_n.
    A positive value means the certificate fails at every dimension.
    """
    order = int(signing.shape[0])
    if order < 2 or signing.shape != (order, order):
        raise ValueError("square signing of order at least 2 required")
    count = 1 << (order - 1)
    states = np.ones((count, order), dtype=np.float64)
    masks = np.arange(count)
    for index in range(1, order):
        states[:, index] = np.where((masks >> (index - 1)) & 1, -1.0, 1.0)
    fields = states @ signing
    energy = 0.5 * np.sum(states * fields, axis=1)
    oriented = states * fields
    positive = np.all(oriented >= 0, axis=1)
    negative = np.all(oriented <= 0, axis=1)
    best = np.full(count, -1.0e18)
    best = np.where(positive, np.maximum(best, energy), best)
    best = np.where(negative, np.maximum(best, -energy), best)
    stable = best > -1.0e17
    deficit = height - best
    neutral = height * ((1.0 + 1.0 / order) ** 1.5 - 1.0)
    active = stable & (deficit + neutral <= order)
    zero = states[active & (deficit <= 0.25)]
    if zero.shape[0] == 0:
        raise ValueError("no zero-deficit stable state")
    spectrum = np.linalg.eigvalsh(zero.T @ zero)[::-1]
    total = float(np.sum(spectrum))
    closest = None
    for dim in range(order):
        tail = total - float(np.sum(spectrum[:dim]))
        average = max(tail / zero.shape[0], 0.0)
        score = dim + np.sqrt(order * average)
        if closest is None or score < closest:
            closest = score
    return float(closest - neutral)


def one_vertex_extension_minimum(signing: np.ndarray) -> float:
    """Minimum Boolean norm of a one-row extension.

    For every row a, Phi(A extended by a) equals max_x (|Q(x)| + |a·x|).
    Antipodal rows and states are represented with the first coordinate +1.
    """
    order = int(signing.shape[0])
    if order < 2 or order > 16 or signing.shape != (order, order):
        raise ValueError("order must be between 2 and 16")
    host = np.asarray(signing, dtype=np.float64)
    count = 1 << (order - 1)
    states = np.ones((count, order), dtype=np.float64)
    masks = np.arange(count)
    for index in range(1, order):
        states[:, index] = np.where((masks >> (index - 1)) & 1, -1.0, 1.0)
    energy = np.abs(0.5 * np.einsum("ia,ab,ib->i", states, host, states))
    best = None
    for start in range(0, count, 512):
        dots = np.abs(states[start : start + 512] @ states.T)
        piece = float((dots + energy).max(axis=1).min())
        if best is None or piece < best:
            best = piece
    return float(best)


def rademacher_tail(order: int, threshold: float) -> float:
    """Probability that an order-term Rademacher sum has absolute value at least threshold."""
    if order < 1 or order > 20 or threshold < 0.0:
        raise ValueError("order in 1..20 and threshold>=0 required")
    if threshold > order:
        return 0.0
    favorable = 0
    for negatives in range(order + 1):
        if abs(order - 2 * negatives) + 1e-8 >= threshold:
            favorable += math.comb(order, negatives)
    return favorable / float(2 ** order)


def partition_tail_sum(signing: np.ndarray, height: float, margin: float) -> float:
    """Sum of Rademacher tails over one representative of each antipodal pair.

    The partition-extension note says a sum strictly below 1 produces a one-row
    extension of Boolean norm strictly less than height + margin.
    """
    order = int(signing.shape[0])
    if order < 2 or order > 16 or signing.shape != (order, order):
        raise ValueError("order must be between 2 and 16")
    if height < 0.0 or margin < 0.0:
        raise ValueError("height>=0 and margin>=0 required")
    host = np.asarray(signing, dtype=np.float64)
    count = 1 << (order - 1)
    states = np.ones((count, order), dtype=np.float64)
    masks = np.arange(count)
    for index in range(1, order):
        states[:, index] = np.where((masks >> (index - 1)) & 1, -1.0, 1.0)
    energy = np.abs(0.5 * np.einsum("ia,ab,ib->i", states, host, states))
    total = 0.0
    for value in energy:
        deficit = height - float(value)
        if deficit < 0.0:
            deficit = 0.0
        total += rademacher_tail(order, deficit + margin)
    return float(total)


def exact_rank_extension_sum(signing: np.ndarray, height: float) -> float:
    """Exact-rank union bound on the active stable family.

    With d_H the rank of the active states and r the neutral increment,
    the sum is 2 * sum exp(-(D+r)^2 / (2 d_H)). A value below 1 would
    certify a one-vertex extension below r. The rank is the number of
    singular values above 1e-8.
    """
    order = int(signing.shape[0])
    if order < 2 or signing.shape != (order, order):
        raise ValueError("square signing of order at least 2 required")
    count = 1 << (order - 1)
    states = np.ones((count, order), dtype=np.float64)
    masks = np.arange(count)
    for index in range(1, order):
        states[:, index] = np.where((masks >> (index - 1)) & 1, -1.0, 1.0)
    fields = states @ signing
    energy = 0.5 * np.sum(states * fields, axis=1)
    oriented = states * fields
    positive = np.all(oriented >= 0, axis=1)
    negative = np.all(oriented <= 0, axis=1)
    best = np.full(count, -1.0e18)
    best = np.where(positive, np.maximum(best, energy), best)
    best = np.where(negative, np.maximum(best, -energy), best)
    stable = best > -1.0e17
    deficit = height - best
    neutral = height * ((1.0 + 1.0 / order) ** 1.5 - 1.0)
    active = stable & (deficit + neutral <= order)
    points = states[active]
    if points.shape[0] == 0:
        raise ValueError("no active stable state")
    rank = int(np.sum(np.linalg.svd(points, compute_uv=False) > 1.0e-8))
    if rank < 1:
        raise ValueError("active states have rank 0")
    widths = deficit[active] + neutral
    return float(np.sum(2.0 * np.exp(-(widths * widths) / (2.0 * rank))))


def spectral_extension_score(signing: np.ndarray, height: float) -> float:
    """Smallest spectral union bound on the active stable family.

    For each dimension d from 1 through n-1, project the active states onto
    the top d right singular vectors. A dimension is discarded when any
    state has margin D+r-sqrt(n)*delta at most 0. The score of a kept
    dimension is 2 * sum exp(-margin^2 / (2d)). A value below 1 would
    certify a one-vertex extension below the neutral increment. Every
    dimension discarded returns infinity.
    """
    order = int(signing.shape[0])
    if order < 2 or signing.shape != (order, order):
        raise ValueError("square signing of order at least 2 required")
    count = 1 << (order - 1)
    states = np.ones((count, order), dtype=np.float64)
    masks = np.arange(count)
    for index in range(1, order):
        states[:, index] = np.where((masks >> (index - 1)) & 1, -1.0, 1.0)
    fields = states @ signing
    energy = 0.5 * np.sum(states * fields, axis=1)
    oriented = states * fields
    positive = np.all(oriented >= 0, axis=1)
    negative = np.all(oriented <= 0, axis=1)
    best = np.full(count, -1.0e18)
    best = np.where(positive, np.maximum(best, energy), best)
    best = np.where(negative, np.maximum(best, -energy), best)
    stable = best > -1.0e17
    deficit = height - best
    neutral = height * ((1.0 + 1.0 / order) ** 1.5 - 1.0)
    active = stable & (deficit + neutral <= order)
    points = states[active]
    if points.shape[0] == 0:
        raise ValueError("no active stable state")
    widths = deficit[active] + neutral
    _, _, basis_rows = np.linalg.svd(points, full_matrices=False)
    best_score = None
    scale = math.sqrt(order)
    for dim in range(1, order):
        basis = basis_rows[:dim]
        projected = points @ basis.T @ basis
        delta = np.linalg.norm(points - projected, axis=1)
        margin = widths - scale * delta
        if np.any(margin <= 0.0):
            continue
        score = float(np.sum(2.0 * np.exp(-(margin * margin) / (2.0 * dim))))
        if best_score is None or score < best_score:
            best_score = score
    if best_score is None:
        return math.inf
    return float(best_score)


def phase_self_gluing_defect(order: int) -> float:
    """Excess of the phase-compatible rank-one lift over the factor 2*sqrt(2).

    This is the lower bound in the two-half self-gluing note: the amplification
    4 against 2*sqrt(2), charged at the Gaussian floor, minus the linear
    perturbation (2+2*sqrt(2))*order. A positive value cannot be absorbed by a
    vanishing dyadic error.
    """
    if order < 2:
        raise ValueError("order>=2")
    amplification = 4.0 - 2.0 * math.sqrt(2.0)
    perturbation = 2.0 + 2.0 * math.sqrt(2.0)
    return 2.0 * amplification * dual_gaussian_lower_bound(order) - perturbation * order


def sigma_square_mass(symmetric: np.ndarray, skew: np.ndarray) -> float:
    """Off-diagonal square mass of A^2-R^2 and AR-RA."""
    difference = symmetric @ symmetric - skew @ skew
    commutator = symmetric @ skew - skew @ symmetric
    mask = 1.0 - np.eye(symmetric.shape[0])
    return float(np.sum((difference * difference + commutator * commutator) * mask))


def anticommutator_squared_frobenius(
    symmetric: np.ndarray, skew: np.ndarray
) -> float:
    """Squared Frobenius norm of AR+RA."""
    commutator = symmetric @ skew + skew @ symmetric
    return float(np.sum(commutator * commutator))


def clifford_sphere_requirement(order: int) -> float:
    """Boolean-norm floor forced by the full complex operator bound."""
    if order < 2:
        raise ValueError("order>=2")
    return order * math.sqrt(2.0 * (order - 1)) / (2.0 * math.sqrt(2.0))


def complete_positive_signing(order: int) -> np.ndarray:
    """The Seidel matrix J-I."""
    if order < 2:
        raise ValueError("order>=2")
    return np.ones((order, order), dtype=np.float64) - np.eye(order)


def min_hermitian_op_over_skew(signing: np.ndarray) -> float:
    """Minimum operator norm of signing + iR over skew signings R."""
    order = int(signing.shape[0])
    if order < 2 or order > 4 or signing.shape != (order, order):
        raise ValueError("order must be between 2 and 4")
    edges = [(i, j) for i in range(order) for j in range(i + 1, order)]
    best = None
    base = signing.astype(np.complex128)
    for mask in range(1 << len(edges)):
        skew = np.zeros((order, order), dtype=np.complex128)
        for bit, (i, j) in enumerate(edges):
            sign = 1.0 if mask & (1 << bit) else -1.0
            skew[i, j] = 1j * sign
            skew[j, i] = -1j * sign
        spectrum = np.linalg.eigvalsh(base + skew)
        operator = float(np.max(np.abs(spectrum)))
        if best is None or operator < best:
            best = operator
    return float(best)


def expected_triangle_polynomial_radius() -> float:
    """Spectral radius of x^3 - 6x - 2."""
    roots = np.roots([1.0, 0.0, -6.0, -2.0])
    return float(np.max(np.abs(roots)))


def optimized_product_delta(order: int) -> float:
    """Largest (order-b)^2 mu_b^2/4 for the two-sided Boolean product floor."""
    if order < 2:
        raise ValueError("order>=2")
    best = 0.0
    for block in range(1, order):
        moment = rademacher_absolute_moment(block)
        value = (order - block) ** 2 * moment * moment / 4.0
        if value > best:
            best = value
    return best


def product_forced_minimum(previous_minimum: float, deleted_delta: float) -> float:
    """Positive root of x = previous + deleted_delta/x."""
    if previous_minimum < 0 or deleted_delta < 0:
        raise ValueError("negative input")
    return 0.5 * (
        previous_minimum
        + math.sqrt(previous_minimum * previous_minimum + 4.0 * deleted_delta)
    )


def phase_range_lower(order: int) -> float:
    """Elementary phase-range floor 2 max (order-b) mu_b."""
    if order < 2:
        raise ValueError("order>=2")
    best = 0.0
    for block in range(1, order):
        value = (order - block) * rademacher_absolute_moment(block)
        if value > best:
            best = value
    return 2.0 * best


def orthogonal_square_gap_minimum(previous_minimum: float, deleted: int) -> float:
    """Root of the orthogonal-completion identity with the cross term removed."""
    if deleted < 2 or previous_minimum < 0:
        raise ValueError("bad args")
    binomial = deleted * (deleted - 1) / 2.0
    return math.sqrt(previous_minimum * previous_minimum + binomial)


def max_cross_row_energy(retained: int, deleted: int) -> float:
    """Largest squared Euclidean energy of a Boolean image of a sign cross block."""
    if retained < 1 or deleted < 1:
        raise ValueError("bad args")
    return float(deleted * retained * retained)


def orthogonal_max_cross_minimum(
    previous_minimum: float, retained: int, deleted: int
) -> float:
    """Square-gap root with the cross energy replaced by its row-sum maximum."""
    base = orthogonal_square_gap_minimum(previous_minimum, deleted)
    return math.sqrt(base * base + max_cross_row_energy(retained, deleted))


def hadamard_cap_over_sqrt2_floor(n: int) -> float:
    """Conference cap divided by sqrt(2) times the Gaussian floor."""
    return spherical_half_bound(n) / (
        math.sqrt(2.0) * dual_gaussian_lower_bound(n)
    )


def hadamard_cap_exceeds_sqrt2_floor(n: int) -> bool:
    """True when the conference cap exceeds sqrt(2) times the floor plus n."""
    return spherical_half_bound(n) > (
        math.sqrt(2.0) * dual_gaussian_lower_bound(n) + n
    )


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

    Sharper than Frobenius Lipschitz (Prop 15.20: gap <= n sqrt(2k)) whenever
    2k < n sqrt(2k), i.e. k < n^2/2 (the sparse-disagreement regime).
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


def maxR_matching_level(C: np.ndarray, pairs, p: Optional[int] = None) -> int:
    """
    Exact max of R_M = Q_C - S_M on the level S_M = -p (Prop 15.32 / Γ-pairing).

    Pairing coordinates: free signs u_a = x_{i_a}, relative z_a = x_{i_a} x_{j_a}.
    Then S_M = sum_a C_{e_a} z_a and R_M = sum_{a<b} γ_ab(z) u_a u_b with
    γ ∈ {-4,-2,0,2,4}. Enumerates all z on the level and all free u (2^{n/2} each).
    Feasible for n=26 (m=13); spike criterion is maxR >= Phi - p = p(m-1).
    """
    C = np.asarray(C)
    n = int(C.shape[0])
    if n % 2 != 0:
        raise ValueError("n even required for perfect matching")
    m = n // 2
    if p is None:
        p = int(round(np.sqrt(n - 1)))
    pairs = [(int(a), int(b)) for a, b in pairs]
    if len(pairs) != m:
        raise ValueError(f"expected {m} pairs, got {len(pairs)}")
    c = np.array([int(C[i, j]) for i, j in pairs], dtype=np.int8)
    U = np.ones((1 << m, m), dtype=np.float64)
    for a in range(m):
        U[:, a] = np.where((np.arange(1 << m) >> a) & 1, -1.0, 1.0)
    best = -10**9
    C8 = np.asarray(C, dtype=np.int8)
    for mask in range(1 << m):
        z = U[mask]
        if int(c @ z) != -p:
            continue
        G = np.zeros((m, m), dtype=np.float64)
        for a in range(m):
            i, j = pairs[a]
            za = int(z[a])
            for b in range(a + 1, m):
                k, l = pairs[b]
                zb = int(z[b])
                g = (
                    int(C8[i, k])
                    + int(C8[i, l]) * zb
                    + int(C8[j, k]) * za
                    + int(C8[j, l]) * za * zb
                )
                G[a, b] = G[b, a] = g
        R = 0.5 * np.sum(U * (U @ G), axis=1)
        best = max(best, int(R.max()))
    return int(best)


def matching_score_mod4_constant(
    C: np.ndarray, Max_plus: np.ndarray, pairs
) -> int:
    """
    Return the unique S_M mod 4 on Max+ (Prop 15.32b). Raises if not constant.
    """
    mods = set()
    for y in Max_plus:
        s = 0
        for i, j in pairs:
            s += int(C[i, j]) * int(y[i]) * int(y[j])
        mods.add(s % 4)
    if len(mods) != 1:
        raise ValueError(f"S_M not constant mod 4: {mods}")
    return int(next(iter(mods)))
