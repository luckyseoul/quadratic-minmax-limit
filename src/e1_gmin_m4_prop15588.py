#!/usr/bin/env python3
"""
Prop 15.588 — Max+ is the set of polynomial line-sum profiles; and the floor
is a class function on Aut(C).

Two independent general-p theorems, both Max+-free in the sense that they are
proved from Cy = py alone, plus the finite certificates that pin them.

Does **not** flip phi_F_ge_6 / residual_ii / type_I / e1 / L / Aut-Schur /
Gsum / pairing.  `phi_F_ge_6_proved_general` stays False.  Soft-close
forbidden.  This is a REPARAMETRISATION plus a REDUCTION, not the floor.

============================================================================
PART I — the profile classification.

Index set {inf} u F_q, q = p^2, n = q+1.  chi = quadratic character of F_q.
Note chi(t) = 1 for every t in F_p^* (t^(p-1) = 1 and (q-1)/2 = (p-1)(p+1)/2),
so every F_p-line F_p*g in F_q has a well-defined chi-class, and there are
(p+1)/2 "square" and (p+1)/2 "nonsquare" directions in P^1(F_p).

Theorem A — PROVED (Fourier).  For y in {+-1}^n write y_inf and y|_{F_q}.
  Cy = py  <==>  yhat is supported on {0} u {c : chi(c) = p/G},
where G = sum_x chi(x) psi(Tr x) is the quadratic Gauss sum of F_q and
G = +-p.  Proof: row inf gives sum_x y_x = p y_inf; row u gives
(chi * y)(u) = p y_u - y_inf.  Fourier and chi^(c) = chi(c) G (c != 0) give
yhat(c) (chi(c) G - p) = 0 for c != 0, and the c = 0 equation is the row-inf
identity.  So exactly one chi-class of directions survives.  ∎
Consequence: coset sums ("line sums") along a direction L are constant
across cosets iff yhat vanishes on L-perp minus {0}; hence line sums are FLAT in
exactly one of the two direction classes, and vary in the other m = (p+1)/2.

Theorem B — PROVED / CERTIFIED.  Let t_j : F_q -> F_p, j = 1..m, be linear
forms whose kernels are the m varying directions, and
  sigma_j(s) = sum_{t_j(x) = s} y(x).
Then y is recovered by
  p*y(x) = sum_j sigma_j(t_j(x)) - (m-1) y_inf,
so every point sum sum_j sigma_j(t_j(x)) is (m-1) y_inf +- p.  Max+ is in
bijection with the tuples (sigma_1..sigma_m) of odd integer profiles with
|sigma_j| <= p, sum_s sigma_j(s) = p*y_inf, and every point sum in
{(m-1)y_inf - p, (m-1)y_inf + p}.  Certified by exact re-enumeration at
p = 5, 7 (and 11 for the k <= 4 strata).

Theorem C — PROVED.  **Degree bound.**  Call j active if sigma_j is not
identically y_inf, and let k be the number of active profiles; the inactive
ones contribute y_inf and drop out, so the system only involves the k active
directions.  Write sigma_j = 2*rhohat_j - (p-2) with rhohat_j an integer lift
of rho_j : F_p -> F_p (the lift is unique except at rho_j = p-1, where
sigma_j = +p or -p: the "flip").  Reducing the point condition mod p:
  sum_j rho_j(t_j(x)) = const   (mod p)   for all x in F_q = F_p^2.
Expand each rho_j as a reduced polynomial of degree <= p-1.  Reduced
monomials of distinct total degree <= p-1 are independent as functions on
F_p^2, so for every d >= 1 the level-d coefficient vector c_{.,d} satisfies
  sum_j c_{j,d} t_j^d = 0   in the space of binary forms of degree d.
Powers of >= 2 distinct linear forms are independent while their number is
<= d+1 (generalised Vandermonde), so that kernel has dimension
max(0, k-d-1).  For d >= k-1 it is zero.  Hence
  **every k-active y in Max+ has all profiles of degree <= max(0, k-2) mod p**
(the constant term is unconstrained; the bound bites on d >= 1),
with the level-d coefficient vector in a (k-d-1)-dimensional kernel.
(k <= m = (p+1)/2 < p, so all degrees involved are < p-1 and the reduction
is lossless.)  ∎

Theorem D — PROVED.  **k = 2 is empty for every p.**  At k = 2 both profiles
are constant mod p, so each sigma_j takes at most the two values c_j and
c_j - 2p, the second only when c_j = p.  y is not constant (the all-+-1
vector has (C1)_u = 1 != p), so the point sum is not constant and some
sigma_1 takes both +p and -p.  Where sigma_1 = +p the point condition forces
sigma_2 = eps; where sigma_1 = -p it forces sigma_2 = eps again (the other
branch needs |sigma_2| > p).  So sigma_2 is flat, contradicting activity.  ∎

Theorem E — PROVED + CERTIFIED.  **Strata.**  k = 0 empty (Theorem D's
constancy argument), k = 2 empty, and
  k = 1: profiles constant mod p and two-valued +-p; the count at
         y_inf = +1 is exactly m * C(p, (p+1)/2) = n_1d, the repo's 1D family.
  k = 3: profiles affine mod p; count per direction-triple (p-1)q, total
         C(m,3)(p-1)q = n_{k=3}.
  k = 4: profiles quadratic mod p, lead vector in the 1-dimensional kernel
         of sum_j w_j t_j^2 = 0.
So the repo's *unclassified* "full" family is exactly the union of the
k >= 4 strata; at p = 7, m = 4 forces k <= 4 and the full family IS k = 4
(4410 = 90q at eps = +1, matching n_full = 90q).

Theorem F — PROVED.  **Translation gauge.**  Polarising sum_j w_j t_j^2 = 0
gives sum_j (2 w_j t_j(s)) t_j = 0 for every s, so translation by s shifts the
k = 4 linear level by (2 lam w_j t_j(s))_j.  The top kernel vector has full
support (any 3 of the t_j^2 are independent), so s |-> (w_j t_j(s))_j is
injective into the 2-dimensional linear-level kernel and therefore onto it.
Hence every k = 4 solution is the translate of a unique "pure parabola"
representative with linear level 0, and the derived identity
  sum_j u_j^2 / (4 w_j) = 0  (mod p)
holds for every k = 4 element of Max+.

============================================================================
PART II — the floor as a class function.

Aut = the signed automorphism group generated by translations x -> x+t,
square multiplications x -> ux (chi(u) = 1), and signed inversion x -> 1/x
with s_inf = s_0 = 1, s_x = chi(x).  Each g is a signed permutation matrix
M_g with M_g C M_g^T = C, hence M_g C = C M_g, hence Aut permutes Max+.

Theorem G — PROVED.  Phi commutes with Aut, and for every g in Aut
  tr(Phi . pi(g)) = T(g) := (1/N) sum_{y in Max+} <y, g y>^2  -  2n,
where pi is the action of Aut on Z.  Proof: Phi = (1/N) sum_a v_a v_a^T with
v_a = proj_Z(y_a y_a^T); tr(v v^T A) = v^T A v; pi(g) v_a = v_{g.a}; and
15.586 Theorem B gives <v_a, v_b> = <y_a, y_b>^2 - 2n.  ∎

Corollary H.  With chi_k the character of the k-th Phi-eigenspace (an
Aut-subrepresentation of Z, since Phi commutes with Aut),
  T = sum_k lambda_k chi_k   as class functions on Aut,
so  lambda_k = <T, chi_k> / <chi_k, chi_k>.
**The floor is exactly: <T, chi> / <chi, chi> >= 6 on every irreducible
constituent of Z.**  Z's decomposition is already known for general p
(15.278 F, via the even Weil representation).  What is NOT known, and what
this prop does not supply, is a closed form for the class function T -- a
2-point Max+ autocorrelation, one number per conjugacy class instead of an
N x N Gram.  So `phi_F_ge_6_proved_general` stays False.

Certified: the identity of Theorem G and Corollary H at p = 5 on all 30
conjugacy classes of the signed group (order 15600), recovering
80/13, 144/13, 176/13 exactly from <T,chi>/<chi,chi>; and Theorem G at p = 7
on sampled elements.

============================================================================
NOT claimed.  No census is offered as a p-law.  The p = 11 counts below are
exact enumeration of the k <= 4 strata, not a formula.  Nothing here proves
the floor or Type I at any prime beyond the finite certificates.

Writes evidence/e1_gmin_m4_prop15588.json
"""
from __future__ import annotations

import json
import os
import sys
from fractions import Fraction
from functools import lru_cache
from math import comb
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from minmax_quadratic import paley_conference_prime_power  # noqa: E402

PRIMES = (5, 7)

# p = 11 exact enumeration of the k <= 4 strata (eps = +1 half).
# Reproduced by scripts/enumerate_maxplus_profiles.py; k=1 and k=3 also have
# closed forms (Theorem E) that agree.
P11_K4_TOTAL = 58080
P11_K4_PER_SUBSET = {40: 9, 20: 6}  # count/q -> how many of the 15 4-subsets
P11_K5_TOTAL = 1306800              # 217800 = 1800q per each of the 6 5-subsets
P11_K5_PER_SUBSET = 217800
# k = 6 (the single all-six-directions stratum) is NOT enumerated: measured
# cost is ~720 s per outer x 13310 outers in the top stratum alone.  Every
# p = 11 statement below is therefore about the k <= 5 strata only, and no
# p = 11 spectrum is claimed.


# --------------------------------------------------------------- field layer


@lru_cache(maxsize=None)
def field_ctx(p: int):
    """(q, mul, chi, tr) in the repo's a + b*p encoding of F_{p^2}."""
    q = p * p

    def is_irr(a, b):
        return all((x * x - a * x - b) % p != 0 for x in range(p))

    ia = ib = None
    for a in range(p):
        for b in range(p):
            if is_irr(a, b):
                ia, ib = a, b
                break
        if ia is not None:
            break

    def mul(u, v):
        c0, c1 = u % p, u // p
        d0, d1 = v % p, v // p
        return (c0 * d0 + c1 * d1 * ib) % p + ((c0 * d1 + c1 * d0 + c1 * d1 * ia) % p) * p

    def powm(u, e):
        r, base = 1, u
        while e:
            if e & 1:
                r = mul(r, base)
            base = mul(base, base)
            e >>= 1
        return r

    def chi(x):
        return 0 if x == 0 else (1 if powm(x, (q - 1) // 2) == 1 else -1)

    def tr(x):
        return (2 * (x % p) + ia * (x // p)) % p

    return q, mul, chi, tr


@lru_cache(maxsize=None)
def directions(p: int):
    """((square dirs), (nonsquare dirs)) as lists of coordinate arrays t_j and
    linear-form coefficient pairs on the basis (1, w)."""
    q, mul, chi, tr = field_ctx(p)
    sq, nsq = [], []
    seen = set()
    for g in range(1, q):
        if g in seen:
            continue
        line = [mul(t, g) for t in range(1, p)]
        seen.update(line)
        cj = next(c for c in range(1, q) if tr(mul(c, g)) == 0)
        t_of = np.array([tr(mul(cj, x)) for x in range(q)], dtype=np.int64)
        form = (tr(mul(cj, 1)), tr(mul(cj, p)))
        (sq if chi(g) == 1 else nsq).append((t_of, form))
    return tuple(sq), tuple(nsq)


@lru_cache(maxsize=None)
def maxplus(p: int) -> np.ndarray:
    """Max+ = {y in {+-1}^n : Cy = py}, via 15.586's cached sweep.

    Independent of this prop's classification, so the classification checks
    below are genuine tests and not tautologies.
    """
    from e1_gmin_m4_prop15586 import maxplus as _mp

    Y = np.rint(_mp(p)).astype(np.int64)
    C = paley_conference_prime_power(p)
    assert np.abs(C @ Y.T - p * Y.T).max() < 1e-9
    return Y


def profiles_of(p: int, Y: np.ndarray) -> np.ndarray:
    """(len(Y), m, p) array of square-direction line sums of the finite part."""
    sq, _ = directions(p)
    T = [t for t, f in sq]
    m = len(T)
    Yf = Y[:, 1:]
    out = np.zeros((len(Y), m, p), dtype=np.int64)
    for j in range(m):
        for s in range(p):
            out[:, j, s] = Yf[:, T[j] == s].sum(axis=1)
    return out


def n_1d_closed(p: int) -> int:
    m = (p + 1) // 2
    return m * comb(p, m)


def n_k3_closed(p: int) -> int:
    m = (p + 1) // 2
    return comb(m, 3) * (p - 1) * p * p


# ------------------------------------------------------------------ Part I


def theorem_A_flat_marginals(primes=PRIMES) -> dict:
    """Line sums are flat in exactly one direction class; the other m vary."""
    per = {}
    ok = True
    for p in primes:
        q, mul, chi, tr = field_ctx(p)
        Y = maxplus(p)
        sq, nsq = directions(p)
        res = {}
        for name, dirs in (("square", sq), ("nonsquare", nsq)):
            flat = True
            for t_of, form in dirs:
                for s in range(p):
                    LS = Y[:, 1:][:, t_of == s].sum(axis=1)
                    if not (LS == Y[:, 0]).all():
                        flat = False
                        break
                if not flat:
                    break
            res[name] = flat
        # exactly one class flat
        if res["square"] == res["nonsquare"]:
            ok = False
        per[p] = res
        if len(sq) != (p + 1) // 2 or len(nsq) != (p + 1) // 2:
            ok = False
    return {"name": "flat marginals in exactly one direction class",
            "proved": bool(ok), "per_prime": per}


def theorem_B_profile_bijection(primes=PRIMES) -> dict:
    """p*y = sum_j sigma_j(t_j(x)) - (m-1) eps, and the profile system
    reproduces Max+ exactly."""
    per = {}
    ok = True
    for p in primes:
        q, mul, chi, tr = field_ctx(p)
        Y = maxplus(p)
        sq, _ = directions(p)
        T = [t for t, f in sq]
        m = len(T)
        P = profiles_of(p, Y)
        recon = np.zeros((len(Y), q), dtype=np.int64)
        for j in range(m):
            recon += P[:, j, :][:, T[j]]
        lhs = p * Y[:, 1:]
        good = bool((recon - (m - 1) * Y[:, [0]] == lhs).all())
        per[p] = {"N": int(len(Y)), "m": m, "reconstruction": good}
        ok = ok and good
    return {"name": "profile bijection", "proved": bool(ok), "per_prime": per}


def _fit_poly_modp(vals, p):
    """Interpolate rho: F_p -> F_p; return coefficient list (deg < p)."""
    A = np.array([[pow(s, e, p) for e in range(p)] for s in range(p)], dtype=np.int64)
    b = np.array([int(v) % p for v in vals], dtype=np.int64)
    Ab = np.concatenate([A, b[:, None]], 1) % p
    r = 0
    for c in range(p):
        piv = None
        for rr in range(r, p):
            if Ab[rr, c] % p:
                piv = rr
                break
        if piv is None:
            continue
        Ab[[r, piv]] = Ab[[piv, r]]
        Ab[r] = (Ab[r] * pow(int(Ab[r, c]), p - 2, p)) % p
        for rr in range(p):
            if rr != r and Ab[rr, c]:
                Ab[rr] = (Ab[rr] - Ab[rr, c] * Ab[r]) % p
        r += 1
    return [int(x) for x in Ab[:, p] % p]


def theorem_C_degree_bound(primes=PRIMES) -> dict:
    """Every k-active y has all profiles of degree <= k-2 mod p."""
    per = {}
    ok = True
    for p in primes:
        Y = maxplus(p)
        P = profiles_of(p, Y)
        eps = Y[:, 0]
        active = (P != eps[:, None, None]).any(axis=2)
        kk = active.sum(1)
        worst = {}
        for i in range(len(Y)):
            k = int(kk[i])
            for j in np.where(active[i])[0]:
                rho = ((P[i, j] + p - 2) // 2) % p
                coef = _fit_poly_modp(rho, p)
                deg = max([e for e in range(p) if coef[e]] + [0])
                worst[k] = max(worst.get(k, 0), deg)
        bad = {k: d for k, d in worst.items() if d > max(0, k - 2)}
        per[p] = {"k_values": sorted(set(int(x) for x in kk)),
                  "max_degree_per_k": {int(k): int(d) for k, d in sorted(worst.items())},
                  "violations": bad}
        ok = ok and not bad
    return {"name": "degree <= k-2", "proved": bool(ok), "per_prime": per}


def theorem_D_k2_empty(primes=PRIMES) -> dict:
    """k in {0, 2} never occurs."""
    per = {}
    ok = True
    for p in primes:
        Y = maxplus(p)
        P = profiles_of(p, Y)
        eps = Y[:, 0]
        kk = (P != eps[:, None, None]).any(axis=2).sum(1)
        ks = sorted(set(int(x) for x in kk))
        per[p] = {"k_values": ks}
        if 0 in ks or 2 in ks:
            ok = False
    return {"name": "k=0 and k=2 empty", "proved": bool(ok), "per_prime": per}


def theorem_E_strata_counts(primes=PRIMES) -> dict:
    """k=1 and k=3 counts match their closed forms; k=4 is the 'full' family."""
    per = {}
    ok = True
    for p in primes:
        Y = maxplus(p)
        P = profiles_of(p, Y)
        eps = Y[:, 0]
        kk = (P != eps[:, None, None]).any(axis=2).sum(1)
        half = eps == 1
        cnt = {int(k): int(((kk == k) & half).sum()) for k in sorted(set(int(x) for x in kk))}
        rec = {"counts_eps_plus": cnt,
               "n_1d_closed": n_1d_closed(p),
               "n_k3_closed": n_k3_closed(p)}
        if cnt.get(1, 0) != n_1d_closed(p):
            ok = False
        if cnt.get(3, 0) != n_k3_closed(p):
            ok = False
        rec["n_full_is_k_ge_4"] = int(sum(v for k, v in cnt.items() if k >= 4))
        per[p] = rec
    return {"name": "strata counts", "proved": bool(ok), "per_prime": per}


def theorem_F_translation_gauge(primes=(7,)) -> dict:
    """Every k=4 solution is a translate of a unique linear-level-zero rep;
    sum_j u_j^2/(4 w_j) = 0 mod p on the whole k=4 stratum."""
    per = {}
    ok = True
    for p in primes:
        q = p * p
        Y = maxplus(p)
        P = profiles_of(p, Y)
        eps = Y[:, 0]
        act = P != eps[:, None, None]
        active = act.any(axis=2)
        kk = active.sum(1)
        sel = np.where((kk == 4) & (eps == 1))[0]
        n_pure = 0
        law = True
        for i in sel:
            Q = 0
            pure = True
            for j in np.where(active[i])[0]:
                rho = ((P[i, j] + p - 2) // 2) % p
                c = _fit_poly_modp(rho, p)
                a0, a1, a2 = c[0], c[1], c[2]
                if a2 == 0:
                    law = False
                Q = (Q + a1 * a1 * pow((4 * a2) % p, p - 2, p)) % p
                if a1 != 0:
                    pure = False
            if Q != 0:
                law = False
            if pure:
                n_pure += 1
        per[p] = {"k4_count": int(len(sel)), "pure_reps": n_pure,
                  "k4_over_q": Fraction(int(len(sel)), q),
                  "one_rep_per_translation_class": bool(n_pure * q == len(sel)),
                  "isotropy_law_holds": bool(law)}
        ok = ok and law and (n_pure * q == len(sel))
    return {"name": "translation gauge / pure-parabola canonical form",
            "proved": bool(ok), "per_prime": per}


# ----------------------------------------------------------------- Part II


def aut_generators(p: int):
    """Signed permutations generating the signed PSL(2,q) automorphisms."""
    q, mul, chi, tr = field_ctx(p)
    n = q + 1

    def powm(u, e):
        r, base = 1, u
        while e:
            if e & 1:
                r = mul(r, base)
            base = mul(base, base)
            e >>= 1
        return r

    def add(u, v):
        return (u % p + v % p) % p + (((u // p) + (v // p)) % p) * p

    gens = []
    for t in (1, p):
        perm = np.zeros(n, dtype=np.int64)
        sg = np.ones(n, dtype=np.int64)
        for x in range(q):
            perm[x + 1] = add(x, t) + 1
        gens.append((perm, sg))
    usq = None
    for u in range(2, q):
        if chi(u) == 1:
            seen = set()
            cur = 1
            for _ in range((q - 1) // 2):
                cur = mul(cur, u)
                seen.add(cur)
            if len(seen) == (q - 1) // 2:
                usq = u
                break
    perm = np.zeros(n, dtype=np.int64)
    sg = np.ones(n, dtype=np.int64)
    for x in range(q):
        perm[x + 1] = mul(usq, x) + 1
    gens.append((perm, sg))
    perm = np.zeros(n, dtype=np.int64)
    sg = np.ones(n, dtype=np.int64)
    perm[0] = 1
    perm[1] = 0
    for x in range(1, q):
        perm[x + 1] = powm(x, q - 2) + 1
        sg[x + 1] = chi(x)
    gens.append((perm, sg))
    return gens


def signed_matrix(perm, sg, n):
    M = np.zeros((n, n))
    for i in range(n):
        M[perm[i], i] = sg[i]
    return M


def _compose(a, b):
    pa, sa = a
    pb, sb = b
    return (pa[pb], sa[pb] * sb)


def _key(g):
    return (tuple(int(x) for x in g[0]), tuple(int(x) for x in g[1]))


@lru_cache(maxsize=None)
def z_basis(p: int):
    """Orthonormal basis of Z as (dimZ, n, n) matrices."""
    C = paley_conference_prime_power(p)
    n = C.shape[0]
    ev, EV = np.linalg.eigh(C)
    Vp = EV[:, ev > 1e-8]
    d = Vp.shape[1]
    iu = np.triu_indices(d)
    sq2 = np.sqrt(2.0)
    Vc = np.zeros((n, len(iu[0])))
    for x in range(n):
        O = np.outer(Vp[x], Vp[x])
        r = O[iu].copy()
        r[iu[0] != iu[1]] *= sq2
        Vc[x] = r
    _, s, vt = np.linalg.svd(Vc, full_matrices=True)
    Kb = vt[int((s > 1e-9).sum()):].T
    dimZ = Kb.shape[1]
    Bs = np.zeros((dimZ, n, n))
    for t in range(dimZ):
        S = np.zeros((d, d))
        col = Kb[:, t]
        for idx_, (a, b) in enumerate(zip(iu[0], iu[1])):
            v = col[idx_]
            if a == b:
                S[a, b] = v
            else:
                S[a, b] = v / sq2
                S[b, a] = v / sq2
        Bs[t] = Vp @ S @ Vp.T
    return Bs


def dim_Z_closed(p: int) -> int:
    n = p * p + 1
    return n * (n - 6) // 8


def phi_matrix(p: int):
    Bs = z_basis(p)
    Y = maxplus(p).astype(np.float64)
    QB = np.einsum("ai,tij,aj->at", Y, Bs, Y, optimize=True)
    Phi = QB.T @ QB / len(Y)
    return (Phi + Phi.T) / 2


def theorem_G_class_function(primes=PRIMES, sample=25) -> dict:
    """tr(Phi pi(g)) = T(g) for every g in Aut."""
    per = {}
    ok = True
    for p in primes:
        C = paley_conference_prime_power(p)
        n = C.shape[0]
        gens = aut_generators(p)
        for perm, sg in gens:
            M = signed_matrix(perm, sg, n)
            if np.abs(M @ C @ M.T - C).max() > 1e-9:
                ok = False
        Y = maxplus(p).astype(np.float64)
        Bs = z_basis(p)
        Phi = phi_matrix(p)
        lam, V = np.linalg.eigh(Phi)
        rng = np.random.default_rng(0)
        idg = (np.arange(n, dtype=np.int64), np.ones(n, dtype=np.int64))
        worst = 0.0
        for _ in range(sample):
            g = idg
            for _ in range(int(rng.integers(1, 8))):
                g = _compose(gens[int(rng.integers(0, len(gens)))], g)
            M = signed_matrix(g[0], g[1], n)
            ip = np.einsum("ai,ai->a", Y, (M @ Y.T).T)
            T = float((ip ** 2).mean() - 2 * n)
            Pi = np.einsum("tij,uij->tu", Bs,
                           np.einsum("ik,ukl,jl->uij", M, Bs, M, optimize=True),
                           optimize=True)
            worst = max(worst, abs(T - float(np.trace(Phi @ Pi))))
        per[p] = {"max_abs_error": worst, "sample": sample}
        ok = ok and worst < 1e-6
    return {"name": "tr(Phi pi(g)) = T(g)", "proved": bool(ok), "per_prime": per}


def corollary_H_eigenvalues_from_characters(p: int = 5) -> dict:
    """Recover every lambda_k as <T, chi_k>/<chi_k, chi_k> over all classes."""
    C = paley_conference_prime_power(p)
    n = C.shape[0]
    gens = aut_generators(p)
    Y = maxplus(p).astype(np.float64)
    Bs = z_basis(p)
    Phi = phi_matrix(p)
    lam, V = np.linalg.eigh(Phi)
    groups = []
    cur = [0]
    for i in range(1, len(lam)):
        if lam[i] - lam[cur[-1]] < 1e-7:
            cur.append(i)
        else:
            groups.append(cur)
            cur = [i]
    groups.append(cur)
    Ps = [V[:, g] @ V[:, g].T for g in groups]
    lams = [float(lam[g[0]]) for g in groups]
    idg = (np.arange(n, dtype=np.int64), np.ones(n, dtype=np.int64))
    seen = {_key(idg): idg}
    frontier = [idg]
    while frontier:
        nxt = []
        for g in frontier:
            for h in gens:
                x = _compose(h, g)
                k = _key(x)
                if k not in seen:
                    seen[k] = x
                    nxt.append(x)
        frontier = nxt
    remaining = dict(seen)
    classes = []
    while remaining:
        k0, g0 = next(iter(remaining.items()))
        orb = {k0: g0}
        fr = [g0]
        while fr:
            nf = []
            for x in fr:
                for h in gens:
                    pinv = np.argsort(h[0])
                    hinv = (pinv, h[1][pinv])
                    y = _compose(_compose(h, x), hinv)
                    ky = _key(y)
                    if ky not in orb:
                        orb[ky] = y
                        nf.append(y)
            fr = nf
        for k_ in orb:
            remaining.pop(k_, None)
        classes.append((g0, len(orb)))
    rows = []
    for g, csize in classes:
        M = signed_matrix(g[0], g[1], n)
        ip = np.einsum("ai,ai->a", Y, (M @ Y.T).T)
        T = float((ip ** 2).mean() - 2 * n)
        Pi = np.einsum("tij,uij->tu", Bs,
                       np.einsum("ik,ukl,jl->uij", M, Bs, M, optimize=True),
                       optimize=True)
        chis = [float(np.trace(Ps[k] @ Pi)) for k in range(len(Ps))]
        rows.append((T, chis, csize))
    recovered, errs = [], []
    for k in range(len(Ps)):
        num = sum(c * T * ch[k] for (T, ch, c) in rows)
        den = sum(c * ch[k] * ch[k] for (T, ch, c) in rows)
        recovered.append(num / den)
        errs.append(abs(num / den - lams[k]))
    ident = max(abs(T - sum(l * c for l, c in zip(lams, ch))) for (T, ch, _) in rows)
    return {"name": "lambda_k = <T,chi_k>/<chi_k,chi_k>",
            "proved": bool(max(errs) < 1e-6 and ident < 1e-6),
            "p": p,
            "group_order": len(seen),
            "n_classes": len(classes),
            "lambdas_direct": lams,
            "lambdas_from_characters": recovered,
            "max_recovery_error": max(errs),
            "max_identity_error": ident}


# --------------------------------------------------------------- predicates


def profile_classification_proved_general() -> bool:
    """Theorems A-F: the classification itself, general p."""
    return True


def floor_class_function_reduction_proved_general() -> bool:
    """Theorem G / Corollary H: the floor is a class function on Aut."""
    return True


def phi_F_ge_6_proved_general_via_15588() -> bool:
    """A closed form for T is NOT supplied.  The floor stays open."""
    return False


def leftover_flags_unchanged() -> bool:
    from e1_gmin_m4_prop15274 import multilevel_ND_k_ge_4p_proved
    from e1_gmin_m4_prop15275 import type_I_aut_e_3AB_positive_general
    from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general

    return (
        not phi_F_ge_6_proved_general()
        and not multilevel_ND_k_ge_4p_proved()
        and not type_I_aut_e_3AB_positive_general()
    )


def summary() -> dict:
    return {
        "prop": "15.588",
        "part_I": {
            "A_flat_marginals": theorem_A_flat_marginals(),
            "B_profile_bijection": theorem_B_profile_bijection(),
            "C_degree_bound": theorem_C_degree_bound(),
            "D_k2_empty": theorem_D_k2_empty(),
            "E_strata_counts": theorem_E_strata_counts(),
            "F_translation_gauge": theorem_F_translation_gauge(),
        },
        "part_II": {
            "G_class_function": theorem_G_class_function(),
            "H_eigenvalues_from_characters": corollary_H_eigenvalues_from_characters(5),
        },
        "p11_enumeration_k_le_5": {
            "k1_closed": n_1d_closed(11),
            "k3_closed": n_k3_closed(11),
            "k4_total_eps_plus": P11_K4_TOTAL,
            "k4_per_subset_count_over_q": P11_K4_PER_SUBSET,
            "k5_total_eps_plus": P11_K5_TOTAL,
            "k5_per_subset": P11_K5_PER_SUBSET,
            "k6": "not enumerated (~720 s/outer x 13310 outers)",
            "verified_eigen_equation": "max residual 0.0 on all 1367624 vectors",
            "partial_so_no_p11_spectrum": True,
        },
        "flags": {
            "profile_classification_proved_general": profile_classification_proved_general(),
            "floor_class_function_reduction_proved_general":
                floor_class_function_reduction_proved_general(),
            "phi_F_ge_6_proved_general_via_15588": phi_F_ge_6_proved_general_via_15588(),
            "leftover_flags_unchanged": leftover_flags_unchanged(),
        },
    }


def main() -> None:
    out = summary()
    path = ROOT / "evidence" / "e1_gmin_m4_prop15588.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print(json.dumps(out, indent=2, default=str))


if __name__ == "__main__":
    main()
