#!/usr/bin/env python3
"""
Prop 15.192 — Gsum row algebra (Max+-free); Aut_e averaging for dual-eq
ker-box; residual (i) still OPEN for general p≥5.

Does **not** flip gsum_disj_lb_proved_general. Soft-close forbidden.

SETUP
  Gsum = E₊[f fᵀ] + E₋[f fᵀ] on edges of K_n, f_e = C_e y_i y_j.
  Dual-eq normal form (15.182): x = α1 − 2e_* + κ, α=6/(pn),
  κ ∈ ker(Gsum), 1ᵀκ=0, κ_e=2−α, −α≤κ_g≤1−α (g≠e).

PROVED Max+-free
  A. Gsum diagonal. Gsum_ee = E₊[f_e²] + E₋[f_e²] = 1+1 = 2.

  B. Adjacent Gsum vanishes.  Already 15.170/15.189: for edges ab,ac,
     Gsum_ab,ac = C_ab C_ac (E₊[y_b y_c] + E₋[y_b y_c]) = 0 from π=±C/p.

  C. Row sum. For every edge e,
        ∑_g Gsum_eg = n.
     Proof. ∑_g f_g = ½ yᵀ C y = ±(p n)/2 on Max± (yᵀCy = ±p n).
     E₊[f_e] = 1/p, E₋[f_e] = −1/p (from π).
     E₊[f_e ∑_g f_g] = (p n/2)(1/p) = n/2;
     E₋[f_e ∑_g f_g] = (−p n/2)(−1/p) = n/2;
     row sum = n.  ∎
     (Uses only π and the star-sum identity ∑_g f_g = ±pn/2; both Max+-free
     for normalised conference / Paley π from 15.189.)

  D. Average disj Gsum.  Off-diag row sum = n−2; adjacent contribute 0;
     #disj partners of e is binom(n−2,2)=(n−2)(n−3)/2; hence
        avg_{g disj e} Gsum_eg = 2/(n−3).  ∎

  E. Aut_e averaging (reduction).  Let Aut_e be the stabiliser in Aut(C) of
     the unordered freeness edge e.  Aut_e preserves Max± setwise and hence
     Gsum.  The dual-eq constraints (Gsum κ=0, 1ᵀκ=0, box, κ_e fixed) are
     Aut_e-invariant and convex.  Therefore if any κ is dual-eq feasible,
     its Aut_e-average is Aut_e-invariant and still dual-eq feasible, with
     the same κ_e.  So dual-eq is possible iff it is possible in the
     Aut_e-invariant subspace.  ∎
     (Group Aut_e is nontrivial for Paley: multiplications by squares fixing
     {∞,0}, Frobenius, and inversion swapping the two endpoints.)

  F. Scheme-ker still blocked (15.190): pure scheme-ker max κ_e =
     α(n−2)/2 = 3(n−2)/(pn) < 2−α for all primes p≥3.

CERTIFIED (not general)
  G. Aut_e-invariant dual-eq ker-box (HiGHS on orbit variables + full Max±
     Gsum):
        p=3: max κ_e = 14/5 = 2.8 ≥ need 9/5  (FEASIBLE; matches full ker)
        p=5: max κ_e = 369/455 ≈ 0.811 < need 127/65  (INFEASIBLE)
        p=7: max κ_e = 11736/19775 ≈ 0.593 < need 347/175  (INFEASIBLE)
     Full ker-box max matches Aut_e-inv max at these p (averaging sharp).
  H. Ratio (full max κ_e) / scheme_max ≈ 3.5 (p=3), 41/28 (p=5), 163/113 (p=7).
     For p≥5 the ratio is < 3/2, and
        (3/2)·3(n−2)/(pn) < 2−α   for all primes p≥5
     (⇔ 9(n−2) < 4 p n − 12 with n=p²+1; cubic 4p³−9p²+4p−3>0 on [5,∞)).
     So a Max+-free proof of max κ_e ≤ (3/2) scheme_max for p≥5 would close
     residual (i) via dual-eq ker-box.

OPEN (blocks residual i / predicate flip)
  I. Prove Max+-free that for all primes p≥5
        max{κ_e : Gsum κ=0, 1ᵀκ=0, −α≤κ_g≤1−α (g≠e)} < 2−α
     e.g. by max κ_e ≤ (3/2)·3(n−2)/(pn), or any stricter Aut_e-invariant
     bound.  Or prove |μ|≤2/n (or ≤1/(2p)) on |κ|=1 (15.188–191).

Until I: gsum_disj_lb_proved_general False; residual_i False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15192.json
"""
from __future__ import annotations

import json
import sys
from collections import deque
from fractions import Fraction
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    is_prime,
    n_of,
)
from e1_gmin_m4_prop15190 import (  # noqa: E402
    alpha_particular,
    scheme_ker_max_kappa_e,
)


def theorem_gsum_diag_two() -> dict:
    return {
        "proved": True,
        "theorem": "Gsum_ee = 2 for every edge (E±[f_e²]=1).",
    }


def theorem_gsum_row_sum_n() -> dict:
    return {
        "proved": True,
        "theorem": (
            "∑_g Gsum_eg = n for every edge e. "
            "From ∑_g f_g = ±(p n)/2 on Max± and E±[f_e]=±1/p."
        ),
        "depends_on": ["pi_equals_C_over_p", "star_sum_features"],
    }


def theorem_avg_disj_gsum() -> dict:
    primes = [p for p in range(3, 60) if is_prime(p)]
    rows = {}
    ok = True
    for p in primes:
        n = n_of(p)
        avg = Fraction(2, n - 3)
        # check matches (n-2) / binom(n-2,2)
        n_disj = (n - 2) * (n - 3) // 2
        avg2 = Fraction(n - 2, n_disj)
        rows[str(p)] = {"avg": str(avg), "match": avg == avg2}
        if avg != avg2:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "avg_{g disj e} Gsum_eg = 2/(n−3). "
            "Off-diag row sum n−2, adjacent 0, #disj=binom(n−2,2)."
        ),
        "depends_on": ["gsum_row_sum_n", "adjacent_zero"],
        "by_p_sample": {k: rows[k] for k in ("3", "5", "7", "11") if k in rows},
    }


def theorem_aut_e_averaging() -> dict:
    return {
        "proved": True,
        "theorem": (
            "If Aut_e ≤ Aut(C) stabilises the freeness edge e, then dual-eq "
            "feasibility is equivalent to Aut_e-invariant dual-eq feasibility "
            "(convex Aut_e-average preserves Gsum κ=0, mass, box, and κ_e)."
        ),
    }


def theorem_three_halves_scheme_beats_need() -> dict:
    """(3/2)*scheme_max < 2-α for all primes p≥5."""
    primes = [p for p in range(5, 80) if is_prime(p)]
    rows = {}
    ok = True
    for p in primes:
        n = n_of(p)
        sch = Fraction(3 * (n - 2), p * n)
        bound = sch * Fraction(3, 2)
        need = 2 - alpha_particular(p)
        rows[str(p)] = {
            "three_halves_scheme": str(bound),
            "need": str(need),
            "strict_lt": bound < need,
        }
        if not (bound < need):
            ok = False
    # Algebra: 9(n-2) < 4pn - 12 with n=p²+1
    # 9p²-9 < 4p³+4p-12 ⇔ 0 < 4p³-9p²+4p-3
    def cubic(p: int) -> int:
        return 4 * p**3 - 9 * p**2 + 4 * p - 3

    cubic_ok = all(cubic(p) > 0 for p in primes)
    return {
        "proved": ok and cubic_ok,
        "theorem": (
            "(3/2)·3(n−2)/(pn) < 2−α for all primes p≥5 "
            "(⇔ 4p³−9p²+4p−3>0)."
        ),
        "by_p_sample": {k: rows[k] for k in ("5", "7", "11", "13") if k in rows},
        "note": (
            "If max dual-eq κ_e ≤ (3/2) scheme_max for all p≥5, residual (i) "
            "closes. Ratio certified <3/2 at p=5,7; open in general."
        ),
    }


def certify_aut_e_kerbox(p: int = 5) -> dict:
    """Census Aut_e-invariant dual-eq ker-box max κ_e."""
    import numpy as np
    from scipy.optimize import linprog
    from minmax_quadratic import paley_conference_prime_power
    from e1_bitight_infeas import boolean_evecs

    C = paley_conference_prime_power(p).astype(float)
    n = C.shape[0]
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

    def mul(u: int, v: int) -> int:
        c0, c1 = u % p, u // p
        d0, d1 = v % p, v // p
        return ((c0 * d0 + c1 * d1 * ib) % p) + (
            (c0 * d1 + c1 * d0 + c1 * d1 * ia) % p
        ) * p

    def powm(u: int, e: int) -> int:
        r, base = 1, u
        while e:
            if e & 1:
                r = mul(r, base)
            base = mul(base, base)
            e >>= 1
        return r

    squares = {powm(x, 2) for x in range(1, q)}
    gens = []

    def perm_from(f_fin, inf_image: int = 0):
        g = np.zeros(n, dtype=int)
        g[0] = inf_image
        for x in range(q):
            g[x + 1] = f_fin(x) + 1
        return g

    for s in squares:
        gens.append(perm_from(lambda x, s=s: mul(s, x)))
    gens.append(perm_from(lambda x: powm(x, p)))
    ginv = np.zeros(n, dtype=int)
    ginv[0] = 1
    ginv[1] = 0
    for x in range(1, q):
        ginv[x + 1] = powm(x, q - 2) + 1
    gens.append(ginv)

    seen = {tuple(range(n))}
    dq = deque([np.arange(n)])
    while dq:
        g = dq.popleft()
        for h in gens:
            nh = g[h]
            t = tuple(nh.tolist())
            if t not in seen:
                seen.add(t)
                dq.append(nh)
    group = [np.array(t) for t in seen]

    edges = list(combinations(range(n), 2))
    e2i = {e: i for i, e in enumerate(edges)}
    e_idx = e2i[(0, 1)]
    visited: set[int] = set()
    orbits: list[list[int]] = []
    for ei, ed in enumerate(edges):
        if ei in visited:
            continue
        orb: list[int] = []
        for g in group:
            a, b = sorted((int(g[ed[0]]), int(g[ed[1]])))
            ej = e2i[(a, b)]
            if ej not in visited:
                visited.add(ej)
                orb.append(ej)
        orbits.append(orb)

    cache = Path("/tmp/grok-goal-03c774661dc8/implementer/e1_p7")
    if p == 7 and (cache / "maxplus.npy").is_file():
        Mp = np.load(cache / "maxplus.npy").astype(float)
        Mm = np.load(cache / "maxminus.npy").astype(float)
    elif p > 5:
        return {"p": p, "skipped": True, "reason": "no Max± cache for p>5"}
    else:
        Mp = boolean_evecs(C, float(p), 0).astype(float)
        Mm = boolean_evecs(C, -float(p), 0).astype(float)

    Fp = np.array([[C[i, j] * y[i] * y[j] for i, j in edges] for y in Mp])
    Fm = np.array([[C[i, j] * y[i] * y[j] for i, j in edges] for y in Mm])
    Gsum = (Fp.T @ Fp) / len(Mp) + (Fm.T @ Fm) / len(Mm)

    alpha = float(alpha_particular(p))
    need = 2 - alpha
    n_orb = len(orbits)
    O = np.zeros((len(edges), n_orb))
    for t, orb in enumerate(orbits):
        for ei in orb:
            O[ei, t] = 1.0
    e_orb = next(t for t, orb in enumerate(orbits) if e_idx in orb)
    A = np.vstack([Gsum @ O, O.sum(0, keepdims=True)])
    _, S, Vt = np.linalg.svd(A, full_matrices=True)
    rnk = int(np.sum(S > 1e-7))
    null = Vt[rnk:].T
    d = null.shape[1]
    if d == 0:
        return {
            "p": p,
            "max_kappa_e": 0.0,
            "need": need,
            "blocked": True,
            "aut_e_order": len(group),
            "n_orbits": n_orb,
            "ker_dim": 0,
        }
    Glist = []
    h = []
    for t in range(n_orb):
        if t == e_orb:
            continue
        Glist.append(null[t])
        h.append(1 - alpha)
        Glist.append(-null[t])
        h.append(alpha)
    G = np.asarray(Glist)
    hh = np.asarray(h)
    res = linprog(
        -null[e_orb],
        A_ub=G,
        b_ub=hh,
        bounds=[(None, None)] * d,
        method="highs",
    )
    max_ke = float(-res.fun) if res.success else None
    sch = float(scheme_ker_max_kappa_e(p))
    return {
        "p": p,
        "max_kappa_e": max_ke,
        "max_kappa_e_frac": str(Fraction(max_ke).limit_denominator(100000))
        if max_ke is not None
        else None,
        "scheme_max": sch,
        "need": need,
        "ratio_to_scheme": (max_ke / sch) if max_ke else None,
        "blocked": max_ke is not None and max_ke < need - 1e-9,
        "aut_e_order": len(group),
        "n_orbits": n_orb,
        "ker_dim": d,
        "proved_general": False,
    }


def hinge_status_192() -> dict:
    return {
        "residual_ii_closed": True,
        "residual_i_closed": False,
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "gsum_row_algebra_proved": True,
        "aut_e_averaging_proved": True,
        "three_halves_scheme_lt_need_proved": True,
        "max_ke_le_three_halves_scheme_general": False,
        "open": (
            "Prove max dual-eq κ_e ≤ (3/2)·scheme_max (or any bound <2−α) "
            "for all p≥5 Max+-free; or |μ|≤2/n on |κ|=1."
        ),
    }


def main() -> dict:
    diag = theorem_gsum_diag_two()
    row = theorem_gsum_row_sum_n()
    avg = theorem_avg_disj_gsum()
    aut = theorem_aut_e_averaging()
    th = theorem_three_halves_scheme_beats_need()
    cert = {str(p): certify_aut_e_kerbox(p) for p in (3, 5)}
    c7 = certify_aut_e_kerbox(7)
    if not c7.get("skipped"):
        cert["7"] = c7
    out = {
        "title": (
            "Prop 15.192 Gsum row algebra; Aut_e dual-eq reduction; "
            "residual (i) OPEN"
        ),
        "proved": {
            "gsum_diag_two": diag["proved"],
            "gsum_row_sum_n": row["proved"],
            "avg_disj_gsum": avg["proved"],
            "aut_e_averaging": aut["proved"],
            "three_halves_scheme_lt_need": th["proved"],
            "max_ke_bound_general": False,
            "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
            "residual_i_closed": False,
            "E1_closed": e1_closed_general(),
            "L_closed": False,
        },
        "algebra": {
            "diag": diag,
            "row_sum": row,
            "avg_disj": avg,
            "aut_e": aut,
            "three_halves": th,
        },
        "certified": {"aut_e_kerbox": cert},
        "hinge": hinge_status_192(),
        "F3": (
            "Aut_e reduction + (3/2) scheme bound would close residual (i). "
            "Census blocks p=5,7; general max_ke bound open. Soft-close forbidden."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15192.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.192 Gsum row algebra; Aut_e reduction; residual (i) OPEN")
    print(f"  row sum n proved={row['proved']}; avg disj=2/(n-3) proved={avg['proved']}")
    print(f"  Aut_e averaging proved={aut['proved']}")
    print(f"  (3/2)scheme < need for p≥5 proved={th['proved']}")
    for k, v in cert.items():
        if v.get("skipped"):
            print(f"  Aut_e kerbox p={k}: skipped")
            continue
        print(
            f"  Aut_e kerbox p={k}: max_ke={v.get('max_kappa_e_frac')} "
            f"blocked={v.get('blocked')} ratio={v.get('ratio_to_scheme')}"
        )
    print(f"  gsum={gsum_disj_lb_proved_general()}; e1={e1_closed_general()}")
    print(f"  open: {out['hinge']['open']}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
