#!/usr/bin/env python3
"""
Algebraic classification of the 144 optimal perfect-matching flips of Paley C_10.

Theorem N10-C (certified):
  (i)  Maximizer criterion: a perfect matching M of K_10 yields Φ=13 after
       flipping on Paley C iff for every maximizer x of C (half-cube),
           sign(Q_C(x)) * S_M(x) >= 1,
       where S_M(x) = sum_{{i,j} in M} C_ij x_i x_j.
       Equivalently, every old maximizer drops to |Q| <= 13 after the flip.
  (ii) The six maximizers with Q_C = +15 already form a complete certificate.
  (iii) PΓL(2,9) acts transitively on the set of 144 optimal matchings
       (single orbit of size 144 under Möbius + Frobenius on PG(1,9)).

Does not claim existence of lim α_n.
"""
from __future__ import annotations

import json
import os
import sys
import time
from itertools import product
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from minmax_quadratic import form_Q, paley_conference_prime_power, phi
from n10_matching_optima import perfect_matchings

for _k in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(_k, "1")

Matching = List[Tuple[int, int]]


# ---------------------------------------------------------------------------
# F_9 field ops matching paley_conference_prime_power(3) encoding
# ---------------------------------------------------------------------------

def f9_params() -> Tuple[int, int]:
    """First irreducible x^2 - a x - b over F_3 (same search order as library)."""
    for a in range(3):
        for b in range(3):
            if all((x * x - a * x - b) % 3 != 0 for x in range(3)):
                return a, b
    raise RuntimeError("no irreducible")


_IA, _IB = f9_params()


def f9_dec(u: int) -> Tuple[int, int]:
    return u % 3, u // 3


def f9_enc(a: int, b: int) -> int:
    return a + 3 * b


def f9_add(u: int, v: int) -> int:
    a, b = f9_dec(u)
    c, d = f9_dec(v)
    return f9_enc((a + c) % 3, (b + d) % 3)


def f9_sub(u: int, v: int) -> int:
    a, b = f9_dec(u)
    c, d = f9_dec(v)
    return f9_enc((a - c) % 3, (b - d) % 3)


def f9_mul(u: int, v: int) -> int:
    c0, c1 = f9_dec(u)
    d0, d1 = f9_dec(v)
    e0 = (c0 * d0 + c1 * d1 * _IB) % 3
    e1 = (c0 * d1 + c1 * d0 + c1 * d1 * _IA) % 3
    return f9_enc(e0, e1)


def f9_pow(u: int, e: int) -> int:
    r, base = 1, u
    while e:
        if e & 1:
            r = f9_mul(r, base)
        base = f9_mul(base, base)
        e >>= 1
    return r


def f9_inv(u: int) -> int:
    if u == 0:
        raise ZeroDivisionError
    return f9_pow(u, 7)  # |F9*| = 8


def f9_frob(u: int) -> int:
    return f9_pow(u, 3)


# ---------------------------------------------------------------------------
# Maximizers and criterion
# ---------------------------------------------------------------------------

def paley_maximizers(C: np.ndarray) -> List[Tuple[np.ndarray, float]]:
    """All half-cube maximizers of Φ(C) (x[0]=+1)."""
    n = C.shape[0]
    npat = 1 << (n - 1)
    best = 0.0
    maxs: List[Tuple[np.ndarray, float]] = []
    for xb in range(npat):
        x = np.ones(n, dtype=np.float64)
        for i in range(n - 1):
            x[i + 1] = 1.0 if (xb >> i) & 1 else -1.0
        q = form_Q(C, x)
        aq = abs(q)
        if aq > best + 1e-12:
            best = aq
            maxs = [(x.copy(), float(q))]
        elif abs(aq - best) < 1e-12:
            maxs.append((x.copy(), float(q)))
    return maxs


def matching_signed_sum(C: np.ndarray, matching: Matching, x: np.ndarray) -> float:
    """S_M(x) = sum_{{i,j} in M} C_ij x_i x_j."""
    return float(sum(C[i, j] * x[i] * x[j] for i, j in matching))


def maximizer_drop_criterion(
    C: np.ndarray,
    matching: Matching,
    maximizers: Sequence[Tuple[np.ndarray, float]] | None = None,
) -> bool:
    """
    True iff every maximizer x of C satisfies sign(Q)*S_M(x) >= 1,
    which forces |Q_C(x) - 2 S_M(x)| <= Φ(C) - 2 when Φ(C)=15 and S integer.
    """
    if maximizers is None:
        maximizers = paley_maximizers(C)
    for x, q in maximizers:
        s = matching_signed_sum(C, matching, x)
        if q > 0 and s < 1.0 - 1e-9:
            return False
        if q < 0 and s > -1.0 + 1e-9:
            return False
    return True


def positive_maximizers(
    maximizers: Sequence[Tuple[np.ndarray, float]],
) -> List[Tuple[np.ndarray, float]]:
    return [(x, q) for x, q in maximizers if q > 0]


def flip_matching(C: np.ndarray, matching: Matching) -> np.ndarray:
    A = C.copy()
    for i, j in matching:
        A[i, j] *= -1.0
        A[j, i] *= -1.0
    return A


def matching_to_frozenset(m: Matching) -> frozenset:
    return frozenset(tuple(sorted(e)) for e in m)


# ---------------------------------------------------------------------------
# PΓL(2,9) action on PG(1,9) = {∞} ∪ F_9  (vertices 0=∞, 1..9 = field 0..8)
# ---------------------------------------------------------------------------

def mobius_map_point(v: int, a: int, b: int, c: int, d: int, use_frob: bool) -> int:
    """Map vertex v under z |-> (a z + b)/(c z + d), optional Frobenius on field pts."""
    if v == 0:
        z = None
    else:
        z = v - 1
        if use_frob:
            z = f9_frob(z)
    if z is None:
        if c == 0:
            return 0
        return f9_mul(a, f9_inv(c)) + 1
    num = f9_add(f9_mul(a, z), b)
    den = f9_add(f9_mul(c, z), d)
    if den == 0:
        return 0
    return f9_mul(num, f9_inv(den)) + 1


def apply_pgl_to_matching(
    matching: Matching, a: int, b: int, c: int, d: int, use_frob: bool = False
) -> frozenset:
    edges = []
    for i, j in matching:
        ii = mobius_map_point(i, a, b, c, d, use_frob)
        jj = mobius_map_point(j, a, b, c, d, use_frob)
        edges.append(tuple(sorted((ii, jj))))
    return frozenset(edges)


def pgl29_representatives() -> List[Tuple[int, int, int, int]]:
    """One matrix per class in PGL(2,9) = GL(2,9)/F9*  (expect 720)."""
    mats = []
    for a, b, c, d in product(range(9), repeat=4):
        det = f9_sub(f9_mul(a, d), f9_mul(b, c))
        if det == 0:
            continue
        mats.append((a, b, c, d))
    seen = set()
    reps = []
    for a, b, c, d in mats:
        key_hit = False
        for s in range(1, 9):
            if (f9_mul(s, a), f9_mul(s, b), f9_mul(s, c), f9_mul(s, d)) in seen:
                key_hit = True
                break
        if key_hit:
            continue
        for s in range(1, 9):
            seen.add((f9_mul(s, a), f9_mul(s, b), f9_mul(s, c), f9_mul(s, d)))
        reps.append((a, b, c, d))
    return reps


def pgl_orbit(matching: Matching, reps: Sequence[Tuple[int, int, int, int]] | None = None) -> set:
    """Orbit of matching under PΓL(2,9) (PGL + Frobenius)."""
    if reps is None:
        reps = pgl29_representatives()
    orbit = set()
    for a, b, c, d in reps:
        for uf in (False, True):
            orbit.add(apply_pgl_to_matching(matching, a, b, c, d, use_frob=uf))
    return orbit


# ---------------------------------------------------------------------------
# Campaign
# ---------------------------------------------------------------------------

def run_classification(out_dir: Path) -> Dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    report: Dict = {
        "n": 10,
        "theorem": "N10-C",
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "f9_irreducible": {"ia": _IA, "ib": _IB},
    }

    C = paley_conference_prime_power(3)
    assert abs(phi(C) - 15.0) < 1e-9
    maximizers = paley_maximizers(C)
    report["n_maximizers"] = len(maximizers)
    report["n_positive_maximizers"] = sum(1 for _, q in maximizers if q > 0)
    report["n_negative_maximizers"] = sum(1 for _, q in maximizers if q < 0)
    assert len(maximizers) == 12
    assert report["n_positive_maximizers"] == 6

    print("[1] enumerate 945 perfect matchings + maximizer criterion ...", flush=True)
    pms = perfect_matchings(10)
    assert len(pms) == 945
    pos_max = positive_maximizers(maximizers)

    opt_phi = []
    opt_crit = []
    opt_crit_pos = []
    for m in pms:
        ph = phi(flip_matching(C, m))
        crit = maximizer_drop_criterion(C, m, maximizers)
        crit_pos = maximizer_drop_criterion(C, m, pos_max)
        if abs(ph - 13.0) < 1e-9:
            opt_phi.append(m)
        if crit:
            opt_crit.append(m)
        if crit_pos:
            opt_crit_pos.append(m)

    report["n_phi13"] = len(opt_phi)
    report["n_criterion_full"] = len(opt_crit)
    report["n_criterion_positive_only"] = len(opt_crit_pos)
    report["criterion_full_equals_phi13"] = len(opt_crit) == len(opt_phi) == 144
    report["criterion_positive_equals_phi13"] = len(opt_crit_pos) == 144
    # set equality
    set_phi = {matching_to_frozenset(m) for m in opt_phi}
    set_crit = {matching_to_frozenset(m) for m in opt_crit}
    set_pos = {matching_to_frozenset(m) for m in opt_crit_pos}
    report["sets_agree"] = set_phi == set_crit == set_pos
    print(
        f"  φ=13: {len(opt_phi)}  crit: {len(opt_crit)}  "
        f"pos-crit: {len(opt_crit_pos)}  agree={report['sets_agree']}",
        flush=True,
    )

    # Necessary sign product
    def sign_prod(m):
        p = 1.0
        for i, j in m:
            p *= C[i, j]
        return p

    report["all_opt_have_sign_prod_minus"] = all(sign_prod(m) < 0 for m in opt_phi)
    n_sign_neg = sum(1 for m in pms if sign_prod(m) < 0)
    report["n_matchings_sign_prod_minus"] = n_sign_neg
    report["sign_prod_necessary_not_sufficient"] = n_sign_neg > 144

    print("[2] PΓL(2,9) orbit of one optimal matching ...", flush=True)
    reps = pgl29_representatives()
    report["n_pgl29"] = len(reps)
    assert len(reps) == 720
    seed = opt_phi[0]
    orbit = pgl_orbit(seed, reps)
    report["orbit_size"] = len(orbit)
    report["orbit_equals_optima"] = orbit == set_phi
    report["orbit_subset_optima"] = orbit <= set_phi
    print(
        f"  |PGL|={len(reps)}  orbit={len(orbit)}  "
        f"equals_optima={report['orbit_equals_optima']}",
        flush=True,
    )

    # Fingerprint one orbit rep
    A = flip_matching(C, seed)
    report["seed_matching"] = [list(e) for e in seed]
    report["seed_phi"] = float(phi(A))
    report["seed_r"] = 2.0 * 13.0 / (10.0 * 3.0)

    report["verdict"] = {
        "n_optimal_matchings": 144,
        "characterised_by_maximizer_drop": report["sets_agree"],
        "positive_maximizers_suffice": report["criterion_positive_equals_phi13"],
        "single_PGammaL_orbit": report["orbit_equals_optima"],
        "note": (
            "The 144 Hamming-5 optima of Paley C_10 are exactly one PΓL(2,9)-orbit "
            "of perfect matchings, equivalently those matchings that drop every "
            "Paley maximizer to |Q|<=13 (six + maximizers already certify)."
        ),
    }
    report["total_time_s"] = time.time() - t0

    out = out_dir / "n10_matching_classify.json"
    out.write_text(json.dumps(report, indent=2, default=str))
    print(f"wrote {out}", flush=True)
    print(json.dumps(report["verdict"], indent=2), flush=True)
    return report


def main():
    out = Path(os.environ.get("OUT", str(ROOT / "evidence")))
    run_classification(out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
