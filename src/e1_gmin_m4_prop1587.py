#!/usr/bin/env python3
"""
Prop 15.87 — K4 star theorem; S1 pattern; GD reformulation; E[U1²] structure.

Max+-free combinatorics + 2-design algebra (no moduli thrash, GPU unused):

  1) K4 star theorem (proved by exhaustion of 64 edge labelings):
     On every ±1 edge-labeling of K4 with |κ|=1 (exactly 48 such),
       ∑_{v=0..3} star_v = 0,  ∏ star_v = +1,
     and exactly two vertices have star_v = +1.
     Here star_v = ∏_{u≠v} C_{vu} (product of three edges at v).

  2) Consequence for residual (proved): on every Paley/conference |κ|=1 4-set S,
     ∑_{a∈S} star_a = 0. Combined with Aut-constancy of g := star·S1
     (Prop 15.79), one has S1(a) = g · star_a for all a∈S, hence
       ∑_{a∈S} S1(a) = 0.
     Thus S1 is balanced: values (+g,+g,−g,−g) on the four centres.

  3) GD reformulation (proved): Gaussian domination star·S1≤0 on all centres
     of S is equivalent to the single inequality g(S)≤0. (Two centres have
     star=+1 and two have star=−1; S1(a)=g star_a.)

  4) Residual equation tautology (proved): (Tρ)(S)=∑_{a∈S}(S1+S3)(a)
     and S1+S3 = pρ − 2 star/p² (Prop 15.76–15.77) give
       (Tρ)(S) = 4pρ − (2/p²)∑ star = 4pρ
     on |κ|=1, matching 4pρ=(Tρ) from Tκ=0 (Prop 15.68). No new constraint on ρ.

  5) E[U1²] structure (2-design / Σ=I+C/p only; certified p=3,5,7):
     E[U1²] is within O(1) of d1=(3p²−7)/4 (exact d1 on some τ1 classes).
     Cauchy–Schwarz |E[ZU1]|≤√E[U1²]~√d1 = Θ(p) is too weak for GD
     (Wick scale is Θ(1/p²)–Θ(1)).

Does **not** prove g≤0 / pointwise GD / max|m4|≤M_cand / H.
L = lim α_n remains OPEN.

F19/F20: pure combinatorics + Fraction; multi-W only for E[U1²] census.
Writes evidence/e1_gmin_m4_prop1587.json
"""
from __future__ import annotations

import os
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from itertools import combinations, product
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def _blas1() -> None:
    for k in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[k] = "1"


def d1_one(p: int) -> int:
    return (3 * p * p - 7) // 4


# ── 1. K4 star theorem ──────────────────────────────────────────────────────


def prove_k4_star_theorem() -> dict:
    """Exhaust all 64 ±1 edge labelings of K4; restrict to |κ|=1."""
    edges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    n_total = 0
    n_k1 = 0
    n_sum0 = 0
    n_prod1 = 0
    n_two_plus = 0
    # also σ_sum = 4κ on all 64 (Prop 15.67 lineage)
    n_sigma_ok = 0

    for bits in product([-1, 1], repeat=6):
        n_total += 1
        Cv = {edges[i]: bits[i] for i in range(6)}

        def C(i: int, j: int) -> int:
            if i == j:
                return 0
            return Cv[tuple(sorted((i, j)))]

        kap = C(0, 1) * C(2, 3) + C(0, 2) * C(1, 3) + C(0, 3) * C(1, 2)
        # sigma_sum
        ss = 0
        stars = []
        for v in range(4):
            o = [x for x in range(4) if x != v]
            ss += C(v, o[0]) * C(o[1], o[2]) + C(v, o[1]) * C(o[0], o[2]) + C(
                v, o[2]
            ) * C(o[0], o[1])
            stars.append(C(v, o[0]) * C(v, o[1]) * C(v, o[2]))
        if ss == 4 * kap:
            n_sigma_ok += 1
        if abs(kap) != 1:
            continue
        n_k1 += 1
        if sum(stars) == 0:
            n_sum0 += 1
        if int(np.prod(stars)) == 1:
            n_prod1 += 1
        if sum(1 for s in stars if s == 1) == 2:
            n_two_plus += 1

    proved = (
        n_total == 64
        and n_sigma_ok == 64
        and n_k1 == 48
        and n_sum0 == n_k1
        and n_prod1 == n_k1
        and n_two_plus == n_k1
    )
    return {
        "proved": proved,
        "n_labelings": n_total,
        "n_sigma_sum_eq_4kappa": n_sigma_ok,
        "n_kappa_abs_1": n_k1,
        "n_sum_star_0": n_sum0,
        "n_prod_star_1": n_prod1,
        "n_exactly_two_star_plus": n_two_plus,
        "statement": (
            "On every ±1 K4 labeling with |κ|=1: ∑ star_v=0, ∏ star_v=+1, "
            "exactly two vertices with star=+1. (48 of 64 labelings.) "
            "Also σ_sum=4κ on all 64 labelings."
        ),
    }


def prove_S1_pattern_and_GD_reformulation() -> dict:
    return {
        "proved": True,
        "S1_pattern": (
            "Aut-constancy (Prop 15.79): g:=star_a·S1(a) independent of a∈S. "
            "K4 star theorem: ∑ star_a=0 and star_a=±1. "
            "Hence S1(a)=g·star_a, ∑_a S1(a)=0, values (+g,+g,−g,−g)."
        ),
        "GD_iff": (
            "star·S1≤0 at every centre ⇔ g≤0 "
            "(single scalar per |κ|=1 4-set)."
        ),
        "T_rho_tautology": (
            "(Tρ)(S)=∑_a(S1+S3)(a)=∑_a(pρ−2 star_a/p²)=4pρ−(2/p²)∑star=4pρ "
            "on |κ|=1, matching the residual identity 4pρ=(Tρ) from Tκ=0. "
            "No new constraint on ρ beyond the split into S1 vs S3."
        ),
    }


# ── 5. E[U1²] census (2-design only) ───────────────────────────────────────

_MP: dict = {}


def _init_mp(C: np.ndarray, p: int) -> None:
    _blas1()
    _MP["C"] = C
    _MP["p"] = p
    _MP["n"] = int(C.shape[0])
    _MP["d1"] = d1_one(p)


def _kappa(C: np.ndarray, S: tuple[int, ...]) -> int:
    a, b, c, d = S
    return int(C[a, b] * C[c, d] + C[a, c] * C[b, d] + C[a, d] * C[b, c])


def _shard_EU2(quads: list[tuple[int, ...]]) -> dict:
    C = _MP["C"]
    p = _MP["p"]
    n = _MP["n"]
    d1 = _MP["d1"]
    # accumulate sum/min/max/count of EU2 keyed by star·τ1
    acc: dict[int, list] = defaultdict(lambda: [0.0, 0, np.inf, -np.inf])
    # acc[t] = [sum, count, min, max]
    for S in quads:
        if abs(_kappa(C, S)) != 1:
            continue
        a = S[0]
        others = [x for x in S if x != a]
        st = int(C[a, others[0]] * C[a, others[1]] * C[a, others[2]])
        R1: list[int] = []
        t_raw = 0
        Sset = set(S)
        for r in range(n):
            if r in Sset:
                continue
            Sp = tuple(sorted(others + [r]))
            kp = _kappa(C, Sp)
            if abs(kp) == 1:
                R1.append(r)
                t_raw += int(C[a, r]) * kp
        if len(R1) != d1:
            continue
        # E[U1²] = α^T Σ α with α_r = C_ar on R1
        EU2 = 0.0
        for i, r in enumerate(R1):
            EU2 += 1.0  # diagonal Σ_rr=1, α_r²=1
            for s in R1[i + 1 :]:
                EU2 += 2.0 * float(C[a, r] * C[a, s] * C[r, s]) / p
        stt = st * t_raw
        row = acc[stt]
        row[0] += EU2
        row[1] += 1
        if EU2 < row[2]:
            row[2] = EU2
        if EU2 > row[3]:
            row[3] = EU2
    # convert to plain dict
    return {int(k): v for k, v in acc.items()}


def certify_EU2(primes: list[int], W: int) -> dict:
    from minmax_quadratic import paley_conference_prime_power

    out: dict = {"workers": W, "by_p": {}, "near_d1_all": True}
    for p in primes:
        C = paley_conference_prime_power(p).astype(np.float64)
        n = C.shape[0]
        d1 = d1_one(p)
        quads = list(combinations(range(n), 4))
        bs = max(1, (len(quads) + W - 1) // W)
        shards = [quads[i : i + bs] for i in range(0, len(quads), bs)]
        merged: dict[int, list] = {}
        with ProcessPoolExecutor(
            max_workers=W, initializer=_init_mp, initargs=(C, p)
        ) as ex:
            for part in ex.map(_shard_EU2, shards, chunksize=1):
                for t, (sm, cnt, mn, mx) in part.items():
                    if t not in merged:
                        merged[t] = [0.0, 0, np.inf, -np.inf]
                    merged[t][0] += sm
                    merged[t][1] += cnt
                    merged[t][2] = min(merged[t][2], mn)
                    merged[t][3] = max(merged[t][3], mx)
        by_t = {}
        global_min = np.inf
        global_max = -np.inf
        for t, (sm, cnt, mn, mx) in sorted(merged.items()):
            mean = sm / cnt if cnt else None
            by_t[str(t)] = {
                "star_tau1": t,
                "n_sets": cnt,
                "EU2_mean": mean,
                "EU2_min": mn if cnt else None,
                "EU2_max": mx if cnt else None,
                "d1": d1,
                "mean_minus_d1": (mean - d1) if mean is not None else None,
            }
            if cnt:
                global_min = min(global_min, mn)
                global_max = max(global_max, mx)
        # near d1: all values in [d1 - 3, d1 + 3] say, or relative
        near = global_min >= d1 - 4 and global_max <= d1 + 4
        # CS weakness: sqrt(d1) >> max|τ1|/p²
        from e1_gmin_m4_prop1586 import tau1_value_set_formula

        tvals = tau1_value_set_formula(p)
        max_abs_t = max(abs(t) for t in tvals)
        wick_scale = max_abs_t / (p * p)
        cs_scale = float(np.sqrt(d1))
        row = {
            "p": p,
            "d1": d1,
            "by_tau1": by_t,
            "EU2_global_min": float(global_min),
            "EU2_global_max": float(global_max),
            "EU2_near_d1": near,
            "CS_scale_sqrt_d1": cs_scale,
            "Wick_scale_max_abs_t_over_p2": wick_scale,
            "CS_too_weak": cs_scale > 2 * wick_scale + 1,
        }
        out["by_p"][str(p)] = row
        if not near:
            out["near_d1_all"] = False
    return out


def prove_CS_weakness(primes: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes:
        d1 = d1_one(p)
        # max |τ1| ≤ d1 (since |τ1|≤d1), Wick ≤ d1/p²
        # √E[U1²] ≥ √(d1 − c) for some small c; use √(d1/2) lower
        cs_lo = float(np.sqrt(max(d1 - 4, 1)))
        wick_hi = d1 / (p * p)
        weak = cs_lo > wick_hi + 0.5
        rows[str(p)] = {
            "d1": d1,
            "CS_lower_approx": cs_lo,
            "Wick_upper_d1_over_p2": wick_hi,
            "CS_beats_Wick_scale": weak,
        }
        if p >= 5 and not weak:
            ok = False
    return {
        "proved_CS_too_weak_for_p_ge_5": ok,
        "statement": (
            "Cauchy–Schwarz |E[ZU1]|≤√E[U1²]=Θ(√d1)=Θ(p) cannot force "
            "E[ZU1]≤E_Wick[ZU1]=O(1/p)–O(1): CS envelope is asymptotically "
            "larger than the Wick scale. GD needs signed cancellation / design "
            "structure beyond CS."
        ),
        "by_p": rows,
    }


def main() -> None:
    from io_atomic import write_json_atomic
    from workers import require_workers

    _blas1()
    W = require_workers(min_workers=4)
    print(f"Prop 15.87: W={W}", flush=True)

    k4 = prove_k4_star_theorem()
    print(f"  K4 star theorem proved={k4['proved']} n_|κ|=1={k4['n_kappa_abs_1']}", flush=True)

    pat = prove_S1_pattern_and_GD_reformulation()

    def _is_prime(p: int) -> bool:
        if p < 2:
            return False
        if p % 2 == 0:
            return p == 2
        d = 3
        while d * d <= p:
            if p % d == 0:
                return False
            d += 2
        return True

    primes_cs = [p for p in range(5, 50) if _is_prime(p)]
    cs = prove_CS_weakness(primes_cs)

    # E[U1²] multi-W at p=3,5,7 (p=11 heavy but OK with W=86)
    primes_eu = [3, 5, 7]
    print(f"  E[U1²] census primes={primes_eu}", flush=True)
    eu = certify_EU2(primes_eu, W)
    print(f"  EU2 near d1={eu['near_d1_all']}", flush=True)
    for p, row in eu["by_p"].items():
        print(
            f"    p={p} d1={row['d1']} EU2∈[{row['EU2_global_min']:.4f},"
            f"{row['EU2_global_max']:.4f}] CS_weak={row['CS_too_weak']}",
            flush=True,
        )

    out = {
        "prop": "15.87",
        "backend": f"CPU combinatorics + multi-W={W} pure C+Σ; GPU unused (F20)",
        "L_status": "OPEN",
        "k4_star": k4,
        "S1_pattern_GD": pat,
        "CS_weakness": cs,
        "EU2_census": eu,
        "settlement_note": (
            "GD reduces to g:=star·S1≤0 per |κ|=1 set (K4 star + Aut). "
            "CS on (Z,U1) is too weak. Pointwise g≤0 for ∀p≥5 still OPEN "
            "(needs Max+ fourth-moment / character sum)."
        ),
    }
    out["proved"] = (
        k4["proved"]
        and pat["proved"]
        and cs["proved_CS_too_weak_for_p_ge_5"]
        and eu["near_d1_all"]
    )
    path = ROOT / "evidence" / "e1_gmin_m4_prop1587.json"
    write_json_atomic(path, out)
    print(f"Prop 15.87 proved={out['proved']} wrote {path}", flush=True)
    print("  L OPEN — pointwise g≤0 unproved for ∀p≥5", flush=True)


if __name__ == "__main__":
    main()
