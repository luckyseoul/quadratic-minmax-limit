#!/usr/bin/env python3
"""
P0 structural attack: bound |m4| on |κ|=1 from Max+/boolean/C only.

Approaches (all multi-worker, F17):
  A) Combinatorial Tκ and T-action on κ-strata (Max+-free)
  B) Residual ρ = m4 - κ/p² satisfies (4p I - T)ρ = Tκ/p²
  C) Boolean + evec pointwise identities → bound |ρ| via star / 2-design
  D) Closed-form candidates for max|m4| and g_min as functions of p
  E) Exact rational m4 from column Hadamard products + frame YC=pY

Writes evidence/e1_gmin_m4_structure.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import Counter, defaultdict
from concurrent.futures import as_completed
from fractions import Fraction
from itertools import combinations
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from workers import pool, require_workers  # noqa: E402


def L_abs(p: int) -> Fraction:
    return Fraction(p - 2, 2 * p * p)


def kappa(C, S) -> int:
    a, b, c, d = S
    return int(C[a, b] * C[c, d] + C[a, c] * C[b, d] + C[a, d] * C[b, c])


def _blas1():
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[k] = "1"


# ---------------------------------------------------------------------------
# A: combinatorial Tκ on conference C (no Max+)
# ---------------------------------------------------------------------------
def _worker_Tkappa_shard(args) -> dict:
    """For each |κ|=1 quad in shard, compute (Tκ)(S) = sum_{v,r} C_vr κ(S_vr)."""
    _blas1()
    p, shard, C = args
    C = np.asarray(C, dtype=np.float64)
    n = C.shape[0]

    # Precompute κ for ALL quads needed? For neighbors only — compute on the fly.
    vals = []
    by_ext_sig = Counter()
    for S in shard:
        kap = kappa(C, S)
        if abs(kap) != 1:
            continue
        Sset = set(S)
        tkap = 0
        # also histogram external signed sums
        for r in range(n):
            if r in Sset:
                continue
            sig = tuple(int(C[v, r]) for v in S)
            by_ext_sig[sig] += 1
            # for each v, contribute C_vr * κ(S_{v→r})
            for v in S:
                others = tuple(sorted(x for x in S if x != v))
                Sp = tuple(sorted(others + (r,)))
                tkap += int(C[v, r]) * kappa(C, Sp)
        # Ext if m4=κ/p² would be Tκ/p²; residual equation uses this
        vals.append(tkap)

    if not vals:
        return {
            "n": 0,
            "max_abs_Tkappa": 0,
            "min_Tkappa": 0,
            "max_Tkappa": 0,
            "mean_Tkappa": 0,
            "unique_Tkappa": {},
            "ext_sig_sample": {},
        }
    arr = np.array(vals, dtype=int)
    uniq = Counter(int(x) for x in arr)
    return {
        "n": len(vals),
        "max_abs_Tkappa": int(np.max(np.abs(arr))),
        "min_Tkappa": int(np.min(arr)),
        "max_Tkappa": int(np.max(arr)),
        "mean_Tkappa": float(np.mean(arr)),
        "unique_Tkappa": {str(k): v for k, v in sorted(uniq.items())},
        "ext_sig_sample": {str(k): int(v) for k, v in by_ext_sig.most_common(20)},
    }


def run_Tkappa(p: int, W: int) -> dict:
    from minmax_quadratic import paley_conference_prime_power

    C = paley_conference_prime_power(p).astype(np.float64)
    n = C.shape[0]
    quads = list(combinations(range(n), 4))
    # only |κ|=1
    k1 = [q for q in quads if abs(kappa(C, q)) == 1]
    # shard
    shards = [k1[i::W] for i in range(W)]
    parts = []
    with pool(W) as ex:
        futs = [ex.submit(_worker_Tkappa_shard, (p, sh, C)) for sh in shards if sh]
        for fut in as_completed(futs):
            parts.append(fut.result())
    ntot = sum(r["n"] for r in parts)
    if ntot == 0:
        return {"p": p, "n": 0}
    max_abs = max(r["max_abs_Tkappa"] for r in parts)
    uniq: Counter = Counter()
    for r in parts:
        for k, v in r["unique_Tkappa"].items():
            uniq[k] += v
    # If m4 = κ/p² exactly, Ext = Tκ/p², |m4| = 1/p² (Wick). Residual ρ satisfies
    # (4p I - T)ρ = Tκ/p², so ρ = 0 iff Tκ=0 on the support... not true.
    # Bound: if |Tκ| ≤ B and |T| op controlled...
    thr_ext = 2.0 * (p - 4) / p
    # Wick Ext = Tκ/p²; same-sign thr on Ext
    max_wick_ext = max_abs / (p * p)
    return {
        "p": p,
        "n_kappa1": ntot,
        "max_abs_Tkappa": max_abs,
        "unique_Tkappa_counts": dict(uniq),
        "n_unique_Tkappa": len(uniq),
        "max_abs_Wick_Ext_Tkappa_over_p2": max_wick_ext,
        "thr_Ext_same_sign": thr_ext,
        "Wick_Ext_alone_le_thr": max_wick_ext <= thr_ext + 1e-12,
        "note": "Tκ is Max+-free (depends only on C). Residual ρ couples via T.",
    }


# ---------------------------------------------------------------------------
# B: pointwise boolean+evec identities (halfspace + Max+ samples)
# ---------------------------------------------------------------------------
def _worker_boolean_ids(p: int) -> dict:
    """
    For every Max+ y: Cy=py, y_i^2=1, sum y = ±(p+1).
    Derive numerical bounds on 4-wise products without claiming closed form yet.
    Also test candidate closed forms for g_min.
    """
    _blas1()
    sys.path.insert(0, str(ROOT / "src"))
    from minmax_quadratic import paley_conference_prime_power

    if p == 5:
        Y = np.load("/tmp/maxplus_p5.npy").astype(np.float64)
    elif p == 7:
        Y = np.load("/tmp/e1_p7/maxplus.npy").astype(np.float64)
    else:
        raise ValueError(p)

    C = paley_conference_prime_power(p).astype(np.float64)
    N, n = Y.shape
    # sums
    sums = Y.sum(axis=1)
    sum_ok = np.all(np.abs(np.abs(sums) - (p + 1)) < 1e-8)
    # Cy = py
    evec_ok = np.allclose(Y @ C.T, p * Y, atol=1e-6) or np.allclose(Y @ C, p * Y, atol=1e-6)
    # actually rows: C @ y = p y → (Y @ C.T) if C sym: Y @ C = p Y
    evec_ok = bool(np.allclose(Y @ C, p * Y, atol=1e-5))

    # Star: for each y,i: sum_j C_ij y_i y_j = p
    # already evec

    # Pairwise dots: E[(y·y')^2] for y' uniform = 2n for fixed Max+ y (incl antipodes?)
    # Compute distance distribution
    # Sample to save time for p=7
    rng = np.random.default_rng(p)
    take = min(N, 200)
    idx = rng.choice(N, size=take, replace=False)
    Ys = Y[idx]
    dots = Ys @ Ys.T  # take x take
    # off-diagonal unique values
    mask = ~np.eye(take, dtype=bool)
    dvals = dots[mask]
    uniq_dots = sorted(set(int(round(x)) for x in dvals))

    # Fourth moment matrix of coordinates: M_ijkl = m4
    # Candidate formulas for max|m4| on |κ|=1
    candidates = {}
    lab = L_abs(p)
    for name, val in {
        "L_abs": lab,
        "(p-2)/(p(2p+3))": Fraction(p - 2, p * (2 * p + 3)),
        "1/(2p)": Fraction(1, 2 * p),
        "1/(p+1)": Fraction(1, p + 1),
        "3/(p(2p+3))": Fraction(3, p * (2 * p + 3)),
        "(p-2)/(2p(p+1))": Fraction(p - 2, 2 * p * (p + 1)),
        "2/(p(p+1))": Fraction(2, p * (p + 1)),
        "(p+1)/(p(2p+1))": Fraction(p + 1, p * (2 * p + 1)),
    }.items():
        candidates[name] = {
            "value": float(val),
            "frac": str(val),
            "ge_L": val >= lab,  # upper bound must be ≥ actual max and ideally ≤ something...
            # for UB of |m4| we need candidate ≥ max|m4| and candidate ≤ L_abs to close
            "le_L": val <= lab,
        }

    # Actual max|m4| from prior census or recompute sample
    from e1_gmin_m4_proof import kappa_pts

    max_m4 = 0.0
    # full for p=5, sample classes for structure
    quads = list(combinations(range(n), 4))
    if p == 7:
        rng2 = np.random.default_rng(0)
        quads = [quads[i] for i in rng2.choice(len(quads), size=min(50000, len(quads)), replace=False)]
    for S in quads:
        kap = kappa(C, S)
        if abs(kap) != 1:
            continue
        a, b, c, d = S
        m = float(np.mean(Y[:, a] * Y[:, b] * Y[:, c] * Y[:, d]))
        max_m4 = max(max_m4, abs(m))

    # Which candidates are valid UBs (≥ max_m4) and ≤ L (closing)?
    closing = []
    valid_ub = []
    for name, info in candidates.items():
        if info["value"] + 1e-12 >= max_m4:
            valid_ub.append(name)
            if info["le_L"]:
                closing.append(name)

    # Try: m4 = E[prod] related to (sum y)^4 expansion
    # (sum_i y_i)^4 = sum_i y_i^4 + ... + 24 sum_{a<b<c<d} y_a y_b y_c y_d + ...
    # For each Max+ y: (sum y_i)^4 = (p+1)^4 constant!
    # Expand and average → exact linear relation among power-sum moments.
    s = p + 1
    s4 = s**4  # for each y with sum = +(p+1); antipodes have sum = -(p+1) same s^4
    # Expansion of (sum y_i)^4:
    # n * 1  (i=i=i=i)
    # + 4 * (terms type i,i,i,j) — but y_i^3 y_j = y_i y_j
    # Count partitions of 4:
    # 4: n terms → n
    # 3+1: C(4,3)=4 ways, n(n-1) pairs → 4 n(n-1) E[y_i^3 y_j]=4 n(n-1) E[y_i y_j]
    # 2+2: C(4,2)/2=3 ways? standard:
    # (sum y_i)^4 = sum_i y_i^4
    #   + 4 sum_{i≠j} y_i^3 y_j
    #   + 6 sum_{i≠j} y_i^2 y_j^2
    #   + 12 sum_{i≠j≠k} y_i^2 y_j y_k   (i distinct from j,k and j≠k)
    #   + 24 sum_{a<b<c<d} y_a y_b y_c y_d
    # With y_i^2=1: y_i^4=1, y_i^3 y_j = y_i y_j, y_i^2 y_j^2=1, y_i^2 y_j y_k = y_j y_k
    #
    # sum_i y_i^4 = n
    # 4 sum_{i≠j} y_i y_j = 4 ( (sum y)^2 - n ) = 4(s^2 - n)  pointwise!
    # 6 sum_{i≠j} 1 = 6 n(n-1)
    # 12 sum_{i≠j≠k distinct} y_j y_k ... need careful
    #
    # Pointwise identity for each Max+ y with |sum y|=p+1:
    # (p+1)^4 = n + 4((sum y)^2 - n) + 6 n(n-1) + 12 * (stuff) + 24 * sum_{quads} prod
    #
    # sum_{i≠j} y_i y_j = (sum y)^2 - n = s^2 - n  (same for ±)
    # Number of ordered (i≠j≠k all distinct) for y_i^2 y_j y_k = y_j y_k:
    # Choose i free (n-2 choices after j,k)? For fixed distinct j,k: sum_{i∉{j,k}} y_j y_k = (n-2) y_j y_k
    # Number of ordered triples i,j,k all distinct in the 12-term: the expansion coefficient for
    # sum_{i,j,k distinct} y_i^2 y_j y_k is 12? Standard multinomial...
    #
    # Use cleaner approach via elementary symmetric / power sums.
    # Let e1 = sum y_i = ±s
    # p2 = sum_{i<j} y_i y_j = (e1^2 - n)/2
    # p3 = sum_{i<j<k} y_i y_j y_k
    # p4 = sum_{a<b<c<d} y_a y_b y_c y_d
    # e1^4 expansion known:
    # e1^4 = sum y_i^4 + 4 sum_{i≠j} y_i^3 y_j + 6 sum_{i≠j} y_i^2 y_j^2
    #      + 12 sum_{i≠j≠k} y_i^2 y_j y_k + 24 p4_ordered_as_unordered
    # = n + 4(e1^2 - n) + 6 n(n-1) + 12 sum_{j≠k} y_j y_k (n-2)  [i runs outside j,k]
    #      wait: sum_{i,j,k all distinct} y_i^2 y_j y_k = sum_{j≠k} y_j y_k * (n-2)
    # = (n-2)(e1^2 - n)
    # + 24 p4
    #
    # Coefficient of the 12-term in (sum y_i)^4 is 12 for ordered distinct i,j,k with pattern i^2 j k.
    # So contribution 12 * (n-2)(e1^2 - n)
    # And 24 * p4 where p4 = sum_{a<b<c<d} prod
    #
    # Check at p=5: s=6, n=26, e1^2=36, e1^4=1296
    # n + 4(36-26) + 6*26*25 + 12*(24)*(36-26) + 24 p4
    # = 26 + 40 + 3900 + 12*24*10 + 24 p4
    # = 3966 + 2880 + 24 p4 = 6846 + 24 p4
    # 1296 = 6846 + 24 p4 ⇒ 24 p4 = 1296-6846 = -5550 ⇒ p4 = -231.25 not integer — ERROR in coefficients

    # Recompute expansion carefully using sympy-free counting:
    # (sum_i y_i)^4 = sum_{i,j,k,l} y_i y_j y_k y_l
    # Partition by equality pattern of indices.

    # Pattern 4 same: n * 1
    # Pattern 3+1: C(4;3,1)=4, number of (i≠j): n(n-1), value y_i y_j → 4 n(n-1) terms... 
    # Actually sum_{i≠j} y_i^3 y_j * 4 = 4 sum_{i≠j} y_i y_j
    # Pattern 2+2: two distinct pairs: C(4;2,2)/2! = 3, sum_{i<j} 2? 
    # sum_{i≠j} y_i^2 y_j^2 * (4!/(2!2!))/2? Standard result for ±1:
    
    # Use: e1^2 = n + 2 p2_pairs where p2_pairs = sum_{i<j} y_i y_j
    # e1^3 = sum y_i^3 + 3 sum_{i≠j} y_i^2 y_j + 6 p3
    #      = e1 + 3 sum_{i≠j} y_j + 6 p3 = e1 + 3(n-1)e1? no
    # sum_{i≠j} y_i^2 y_j = sum_{i≠j} y_j = sum_j y_j (n-1) = (n-1) e1
    # so e1^3 = e1 + 3(n-1)e1 + 6 p3 = e1(3n-2) + 6 p3
    # p3 = (e1^3 - (3n-2)e1)/6
    # For e1=±s: p3 = ±(s^3 - (3n-2)s)/6

    # e1^4 = e1 * e1^3 = e1( (3n-2)e1 + 6 p3 ) = (3n-2)e1^2 + 6 e1 p3
    # Also e1^4 = sum_{ijkl} ...
    # Relation: 24 p4 + lower = e1^4
    # Known identity:
    # p4 = (e1^4 - 6 e1^2 (n- something) ... )
    # From Wikipedia / elementary:
    # sum_{i<j<k<l} y_i y_j y_k y_l = (1/24) * (
    #   e1^4 - 6 e1^2 n + 3 n^2 + 8 e1 p3_alt ... )
    
    # Direct from power-sum / Newton for z_i=y_i with z^2=1:
    # Use generating function prod_i (1 + y_i t) = sum e_k t^k
    # e4 = p4, and for y_i=±1, the poly is determined by how many -1's.
    # If sum y = s = n - 2w where w = # of -1's, then
    # p4 = e4 of the multiset of ±1 = coefficient.
    # Number of -1's: w = (n - e1)/2
    # e_k = sum of products of k distinct y's = elementary symmetric in the ±1 list
    # = sum_{j=0}^k C(w,j) C(n-w, k-j) (-1)^j
    # For e1 = s = p+1: w = (n - s)/2 = (p^2+1 - p - 1)/2 = (p^2 - p)/2 = p(p-1)/2
    # e4 = sum_{j=0}^4 C(w,j) C(n-w, 4-j) (-1)^j

    w = p * (p - 1) // 2
    n_plus = n - w  # (n+s)/2
    assert n_plus + w == n
    e4 = 0
    from math import comb

    for j in range(5):
        if j <= w and (4 - j) <= n_plus:
            e4 += comb(w, j) * comb(n_plus, 4 - j) * ((-1) ** j)
    # Pointwise: every Max+ y with sum = +(p+1) has the SAME number of -1's, hence SAME e4 = p4!
    # Antipodal y has sum = -(p+1), w' = (n+s)/2 = n_plus, e4' = same by sign?
    # For sum = -s: w_neg = (n+s)/2 = n_plus, e4_neg = sum_j C(n_plus,j) C(w,4-j) (-1)^j
    # Not necessarily equal to e4 — but product of 4 coords: flipping all signs leaves e4 invariant (4 even).
    # Flipping all signs: each 4-product unchanged. And antipode has sum -s with w_neg # of -1's.
    # Actually antipode -y: if y has w minuses, -y has n-w minuses. e4(y)=e4(-y) since 4 even.
    # So p4 is CONSTANT on all Max+ ! 
    # Then E[p4] = p4 = e4, and also E[p4] = binom(n,4) * mean_m4_over_all_quads
    # mean over ALL 4-sets of m4 = e4 / binom(n,4)
    # But we need m4 only on |κ|=1, not global mean.

    mean_m4_all = e4 / comb(n, 4) if comb(n, 4) else 0

    # Global mean of κ: average of three pairing products
    # For conference, average C_ab C_cd over disjoint pairs known from C^2

    return {
        "p": int(p),
        "N": int(N),
        "n": int(n),
        "sum_coords_pm_p1": bool(sum_ok),
        "YC_eq_pY": bool(evec_ok),
        "uniq_pairwise_dots_sample": uniq_dots[:30],
        "n_uniq_dots_sample": len(uniq_dots),
        "w_num_minus": int(w),
        "e4_constant_on_Max": int(e4),
        "mean_m4_all_quads": float(mean_m4_all),
        "max_abs_m4_kappa1_emp": max_m4,
        "L_abs": float(lab),
        "candidates": candidates,
        "valid_upper_bounds_ge_maxm4": valid_ub,
        "closing_candidates_le_L_and_ge_maxm4": closing,
        "pointwise_e4_identity": (
            "Every Max+ vector has exactly w=p(p-1)/2 coordinates equal to -1 "
            "(when oriented with sum=+(p+1)), hence the elementary symmetric e4 "
            "is CONSTANT on Max+. This fixes the average of m4 over ALL 4-sets, "
            "not the max on |κ|=1."
        ),
    }


# ---------------------------------------------------------------------------
# C: algebraic identity from e4 + κ stratification
# ---------------------------------------------------------------------------
def _worker_kappa_weighted_m4(p: int) -> dict:
    """
    Split binom(n,4) into κ-classes; use constant e4 and known κ-distribution
    (from C only) to constrain weighted averages of m4.
    """
    _blas1()
    sys.path.insert(0, str(ROOT / "src"))
    from minmax_quadratic import paley_conference_prime_power
    from math import comb

    C = paley_conference_prime_power(p).astype(np.float64)
    n = C.shape[0]
    # κ distribution (Max+-free)
    kap_counts = Counter()
    for S in combinations(range(n), 4):
        kap_counts[kappa(C, S)] += 1

    # From e4 constant: sum_S m4(S) = e4  (since sum_y sum_S prod = sum_S N m4 = N e4,
    # and for each y, sum_S prod = e4, so sum_S m4 = e4)
    # Wait: m4(S) = (1/N) sum_y prod_S(y), sum_S m4(S) = (1/N) sum_y e4(y) = e4
    w = p * (p - 1) // 2
    n_plus = n - w
    e4 = 0
    for j in range(5):
        if j <= w and (4 - j) <= n_plus:
            e4 += comb(w, j) * comb(n_plus, 4 - j) * ((-1) ** j)

    # Load Max+ for empirical m4 by κ
    if p == 5:
        Y = np.load("/tmp/maxplus_p5.npy").astype(np.float64)
    elif p == 7:
        Y = np.load("/tmp/e1_p7/maxplus.npy").astype(np.float64)
    else:
        return {"p": p, "skip": True}

    N = len(Y)
    sum_m4 = 0.0
    m4_by_kap = defaultdict(list)
    # For p=7 full pass is fine with vectorization per κ sample
    for S in combinations(range(n), 4):
        a, b, c, d = S
        m = float(np.mean(Y[:, a] * Y[:, b] * Y[:, c] * Y[:, d]))
        k = kappa(C, S)
        sum_m4 += m
        if abs(k) == 1:
            m4_by_kap[k].append(m)
        elif abs(k) == 3:
            m4_by_kap[k].append(m)

    # Theoretical: sum_m4 should equal e4
    # Weighted: n1 * mean_m4_on_|κ|=1 * (sign structure) + n3 * mean_on_|κ|=3 = e4

    stats = {}
    for k, vals in m4_by_kap.items():
        arr = np.array(vals)
        stats[str(k)] = {
            "count": len(vals),
            "mean_m4": float(np.mean(arr)),
            "max_abs": float(np.max(np.abs(arr))),
            "min": float(np.min(arr)),
            "max": float(np.max(arr)),
        }

    n1 = kap_counts[1] + kap_counts[-1]
    n3 = kap_counts[3] + kap_counts[-3]
    # If m4 had constant absolute value M on |κ|=1 with sign(m4)=sign(κ), 
    # and M3 on |κ|=3 with sign(m4)=sign(κ), then
    # sum_S m4 = M*(n_{κ=1} - n_{κ=-1}) + M3*(n3pos - n3neg) = e4
    # And n_{κ=1} may equal n_{κ=-1} by symmetry → contribution 0 from |κ|=1 if odd function of κ!

    return {
        "p": int(p),
        "kappa_counts": {str(k): int(v) for k, v in sorted(kap_counts.items())},
        "e4": int(e4),
        "sum_m4_empirical": sum_m4,
        "sum_m4_vs_e4_err": abs(sum_m4 - e4),
        "m4_stats_by_kappa": stats,
        "n_kappa_abs1": int(n1),
        "n_kappa_abs3": int(n3),
        "symmetry_note": (
            "If m4 is odd in κ (m4(-S_signs) = -m4) and κ-counts are symmetric, "
            "the e4 identity constrains |κ|=3 more than |κ|=1 max."
        ),
    }


# ---------------------------------------------------------------------------
# D: prove bound via residual z on two coordinates (conditional)
# ---------------------------------------------------------------------------
def _worker_residual_z(p: int) -> dict:
    """
    Fix edge {i,j}, condition on (y_i,y_j). Residual z_k = y_k - μ_k.
    Prop 15.50: μ = Σ_{*S} Σ_SS^{-1} y_S. Bound E[z_k z_l] and 4th moments.
    Target: show |E[y_a y_b y_c y_d]| ≤ L_abs when |κ|=1.
    """
    _blas1()
    sys.path.insert(0, str(ROOT / "src"))
    from minmax_quadratic import paley_conference_prime_power

    if p == 5:
        Y = np.load("/tmp/maxplus_p5.npy").astype(np.float64)
    elif p == 7:
        Y = np.load("/tmp/e1_p7/maxplus.npy").astype(np.float64)
    else:
        raise ValueError(p)

    C = paley_conference_prime_power(p).astype(np.float64)
    N, n = Y.shape
    # Σ = I + C/p
    # For S={0,1}, Σ_SS = [[1, c/p],[c/p, 1]] with c=C_01
    # For disjoint k,l with various κ involving 0,1,k,l

    # Pick a fixed pair i=0,j=1
    i, j = 0, 1
    cij = float(C[i, j])
    # Conditional mean for k: μ_k = (C_ki y_i + C_kj y_j)/(p + something)?
    # Σ_{kS} Σ_SS^{-1} y_S
    # Σ_SS^{-1} = p^2/(p^2-1) * [[1, -c/p],[-c/p, 1]] wait det = 1 - 1/p^2 = (p^2-1)/p^2
    # Actually det(Σ_SS) = 1 - (cij/p)^2 = (p^2-1)/p^2 since cij=±1
    # Σ_SS^{-1} = p^2/(p^2-1) * [[1, -cij/p], [-cij/p, 1]]
    
    # Empirical: for each Max+ y, compute z and check z_k^2 + 2 μ_k z_k = 1 - μ_k^2
    # (boolean identity y_k^2=1)

    # Direct bound attempt: E[y_i y_j y_k y_l] = E[ y_i y_j (μ_k + z_k)(μ_l + z_l) ]
    # = E[y_i y_j μ_k μ_l] + E[y_i y_j μ_k z_l] + E[y_i y_j z_k μ_l] + E[y_i y_j z_k z_l]
    # μ linear in y_i,y_j so first term is Wick-like; cross terms may vanish; last is residual.

    # Compute for all k,l disjoint from {i,j} the empirical m4 vs formula
    errs = []
    max_m4_k1 = 0.0
    n_k1 = 0
    # μ coefficients: μ_k = α_k y_i + β_k y_j
    inv_scale = (p * p) / (p * p - 1.0)
    for k in range(n):
        if k in (i, j):
            continue
        for l in range(k + 1, n):
            if l in (i, j):
                continue
            kap = kappa(C, (i, j, k, l))
            if abs(kap) != 1:
                continue
            n_k1 += 1
            # Σ_{k,{i,j}} = [C_ki/p, C_kj/p]
            # μ_k = [C_ki/p, C_kj/p] Σ_SS^{-1} [y_i, y_j]^T
            rowk = np.array([C[k, i] / p, C[k, j] / p])
            rowl = np.array([C[l, i] / p, C[l, j] / p])
            SSinv = inv_scale * np.array([[1.0, -cij / p], [-cij / p, 1.0]])
            # For each y:
            yi, yj = Y[:, i], Y[:, j]
            yk, yl = Y[:, k], Y[:, l]
            # μ_k(y) = rowk @ SSinv @ [yi,yj]
            coef_k = rowk @ SSinv  # shape (2,)
            coef_l = rowl @ SSinv
            mu_k = coef_k[0] * yi + coef_k[1] * yj
            mu_l = coef_l[0] * yi + coef_l[1] * yj
            z_k = yk - mu_k
            z_l = yl - mu_l
            m4 = float(np.mean(yi * yj * yk * yl))
            # Expansion
            t_mumu = float(np.mean(yi * yj * mu_k * mu_l))
            t_muz = float(np.mean(yi * yj * mu_k * z_l))
            t_zmu = float(np.mean(yi * yj * z_k * mu_l))
            t_zz = float(np.mean(yi * yj * z_k * z_l))
            pred = t_mumu + t_muz + t_zmu + t_zz
            errs.append(abs(pred - m4))
            max_m4_k1 = max(max_m4_k1, abs(m4))
            # boolean residual identity: y^2=1 ⇒ z^2 + 2μz = 1-μ^2
            # check mean abs violation
    # Bound |t_zz|: this is the hard residual
    # At p=5 Fréchet was too weak; measure how large t_zz is vs budget

    # Recompute with storing t_zz for |κ|=1
    tzz_same = []
    tzz_all = []
    mumu_all = []
    for k in range(n):
        if k in (i, j):
            continue
        for l in range(k + 1, n):
            if l in (i, j):
                continue
            kap = int(C[i, j] * C[k, l] + C[i, k] * C[j, l] + C[i, l] * C[j, k])
            if abs(kap) != 1:
                continue
            rowk = np.array([C[k, i] / p, C[k, j] / p])
            rowl = np.array([C[l, i] / p, C[l, j] / p])
            SSinv = inv_scale * np.array([[1.0, -cij / p], [-cij / p, 1.0]])
            coef_k = rowk @ SSinv
            coef_l = rowl @ SSinv
            yi, yj = Y[:, i], Y[:, j]
            mu_k = coef_k[0] * yi + coef_k[1] * yj
            mu_l = coef_l[0] * yi + coef_l[1] * yj
            z_k = Y[:, k] - mu_k
            z_l = Y[:, l] - mu_l
            t_mumu = float(np.mean(yi * yj * mu_k * mu_l))
            t_zz = float(np.mean(yi * yj * z_k * z_l))
            m4 = float(np.mean(yi * yj * Y[:, k] * Y[:, l]))
            tzz_all.append(t_zz)
            mumu_all.append(t_mumu)
            # same-sign of m4 with κ
            if np.sign(m4) * np.sign(kap) > 0:
                tzz_same.append(t_zz)

    lab = float(L_abs(p))
    return {
        "p": int(p),
        "n_kappa1_from_edge_01": n_k1,
        "expansion_max_err": float(max(errs)) if errs else None,
        "max_abs_m4": max_m4_k1,
        "L_abs": lab,
        "max_abs_tzz": float(np.max(np.abs(tzz_all))) if tzz_all else None,
        "max_abs_tzz_same_sign_m4": float(np.max(np.abs(tzz_same))) if tzz_same else None,
        "max_abs_mumu": float(np.max(np.abs(mumu_all))) if mumu_all else None,
        "mean_abs_tzz": float(np.mean(np.abs(tzz_all))) if tzz_all else None,
        "note": (
            "m4 = E[y_i y_j μ_k μ_l] + cross + E[y_i y_j z_k z_l]. "
            "μ is Max+-free (depends on C,y_i,y_j only). "
            "Closing L requires |t_mumu + t_zz + cross| ≤ L_abs."
        ),
    }


def main():
    W = require_workers(min_workers=4)
    print(f"m4_structure: W={W}", flush=True)
    out = {"workers": W, "parts": {}}

    # A: Tκ at p=5,7 (combinatorial, multi-worker)
    for p in (5, 7):
        print(f"Tκ census p={p}...", flush=True)
        out["parts"][f"Tkappa_p{p}"] = run_Tkappa(p, W)
        print(f"  {out['parts'][f'Tkappa_p{p}']}", flush=True)

    # B,C,D in process pool
    with pool(W) as ex:
        futs = {
            ex.submit(_worker_boolean_ids, 5): "bool5",
            ex.submit(_worker_boolean_ids, 7): "bool7",
            ex.submit(_worker_kappa_weighted_m4, 5): "kw5",
            ex.submit(_worker_kappa_weighted_m4, 7): "kw7",
            ex.submit(_worker_residual_z, 5): "z5",
            ex.submit(_worker_residual_z, 7): "z7",
        }
        for fut in as_completed(futs):
            tag = futs[fut]
            r = fut.result()
            out["parts"][tag] = r
            print(f"  {tag}: keys={list(r.keys())[:8]}...", flush=True)

    # Algebra synthesis
    lab5 = L_abs(5)
    # Candidate proof sketch from e4 + T
    out["algebra_synthesis"] = {
        "master": "m4=κ/p²+Ext/(4p), Ext=Tm4, σ_sum=4κ",
        "e4_constant": (
            "Every Max+ y has exactly w=p(p-1)/2 minus signs (oriented sum=p+1), "
            "so e4=∑_{a<b<c<d} y_a y_b y_c y_d is CONSTANT on Max+. "
            "Hence ∑_S m4(S) = e4 (exact, Max+ structure)."
        ),
        "Tkappa_role": (
            "If ρ=m4-κ/p², then (4pI-T)ρ = Tκ/p². "
            "Tκ is pure conference combinatorics. "
            "Bound follows if ‖(4pI-T)^{-1}‖ is controlled on the |κ|=1 subspace "
            "and Tκ is known."
        ),
        "L_abs_p5": str(lab5),
        "closing_status": "OPEN",
    }

    # Check if any closing candidate appeared
    closing = []
    for tag in ("bool5", "bool7"):
        r = out["parts"].get(tag, {})
        closing.extend(r.get("closing_candidates_le_L_and_ge_maxm4") or [])
    out["closing_candidates_found"] = closing

    # Structural theorem candidates that ARE proved in this run
    proved = []
    for p in (5, 7):
        b = out["parts"].get(f"bool{p}" if False else f"bool{p}", {})
    # fix tags
    for p, tag in ((5, "bool5"), (7, "bool7")):
        b = out["parts"].get(tag, {})
        if b.get("sum_coords_pm_p1") and b.get("YC_eq_pY"):
            proved.append(f"p={p}: sum=±(p+1) and YC=pY on all Max+")
        if b.get("e4_constant_on_Max") is not None:
            proved.append(f"p={p}: e4={b['e4_constant_on_Max']} constant (formula w=p(p-1)/2)")
    for p, tag in ((5, "kw5"), (7, "kw7")):
        k = out["parts"].get(tag, {})
        if k.get("sum_m4_vs_e4_err") is not None and k["sum_m4_vs_e4_err"] < 1e-6:
            proved.append(f"p={p}: ∑_S m4(S)=e4 certified err={k['sum_m4_vs_e4_err']:.2e}")

    out["proved_this_run"] = proved
    out["status"] = (
        "Structural lemmas: e4 constant on Max+ (boolean weight w=p(p-1)/2); "
        "∑ m4 = e4; Tκ combinatorial spectrum computed; residual z expansion identity. "
        "OPEN: convert these into |m4|≤L_abs on |κ|=1 for all primes p≥5 "
        "(no per-p Max+ census as proof)."
    )

    path = ROOT / "evidence" / "e1_gmin_m4_structure.json"
    path.write_text(json.dumps(out, indent=2))
    print(out["status"], flush=True)
    print(f"wrote {path}", flush=True)


if __name__ == "__main__":
    main()
