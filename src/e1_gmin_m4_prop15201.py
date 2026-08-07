#!/usr/bin/env python3
"""
Prop 15.201 — Sufficient free-e bound α(n−2); p=7 exact μ; residual (i) OPEN.

Does **not** flip residual_i / gsum / E1 / L. Soft-close forbidden.

SETUP (15.182, 15.190, 15.200)
  Dual-eq empty ⇔ free-e ker-box max κ_e < 2−α, where
    max{κ_e : Gsumκ=0, 1ᵀκ=0, −α≤κ_g≤1−α (g≠e)}.
  Scheme-ker max = α(n−2)/2 = 3(n−2)/(pn) (15.190).
  α = 6/(p n), n=p²+1.

PROVED Max+-free Fraction
  A. Sufficient free-e threshold.  For every prime p≥5:
        α(n−2) < 2−α.
     Proof. α(n−2) < 2−α ⇔ α(n−1) < 2 ⇔ 6(n−1)/(p n) < 2
     ⇔ 3(n−1) < p n ⇔ 3p² < p(p²+1) = p³+p
     ⇔ 0 < p³−3p²+p = p(p²−3p+1).  For p≥5, p²−3p+1 ≥ 11>0.  ∎
     Note α(n−2) = 2·scheme_max.  So any proof of
        free-e max κ_e ≤ α(n−2)
     closes dual-eq for all primes p≥5.

  B. Comparison with (3/2)·scheme.  (3/2)·scheme = (3/4)·α(n−2) < α(n−2),
     so the α(n−2) target is weaker than 15.192's (3/2)·scheme target
     (easier upper bound, still sufficient by A).

  C. |μ|≤2/n still sufficient for Farkas (15.176 path with Gsum≥−4/n).
     Algebra: −4/n > μ_*=(6−4p)/(p(3p−2)) for all primes p≥5
     ⇔ 2p³−9p²+6p−3>0 (value 52 at p=5; increasing on [5,∞)).  ∎

CERTIFIED (not general)
  D. Free-e max κ_e vs α(n−2):
        p=3: max≈2.30 > α(n−2)=1.6  (bound false; dual-eq feasible — OK)
        p=5: max≈0.805 < α(n−2)≈1.108  (bound holds; dual-eq blocked)
        p=7: max≈0.593 < α(n−2)≈0.823  (bound holds; dual-eq blocked)
     Ratios max/scheme: ≈2.87 (p=3), 1.45 (p=5), 1.44 (p=7).

  E. Exact Max± μ on |κ|=1 through a fixed edge at p=7:
        max|μ| = 109/2863 ≈ 0.038072 < 2/n = 1/25 = 0.04
        At extreme |φ|=2(p−2)=10: μ ∈ {±109/2863, ±61/2863}
        max|μ| is attained on the extreme-|φ| locus.
        |μ|>|f4| on 576/864 of |κ|=1 4-sets, but absolute max is at extreme f4.

OPEN (blocks residual i / predicate flip)
  F. Prove for all primes p≥5 one of:
     (i)  free-e max κ_e ≤ α(n−2)  (sufficient by A),
     (ii) free-e max κ_e ≤ (3/2)·scheme_max,
     (iii) |μ|≤2/n on every |κ|=1 4-set,
     (iv) K₄≤Wick_hi.
     Then flip gsum_disj_lb_proved_general → residual_i → E1 (if residual ii
     allows) → L.

Until F: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15201.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    is_prime,
    k_type_I_3p_minus_2,
    n_of,
)
from e1_gmin_m4_prop15176 import (  # noqa: E402
    mu_star_threshold,
    residual_i_farkas_if_mu,
)
from e1_gmin_m4_prop15190 import (  # noqa: E402
    alpha_particular,
    scheme_ker_max_kappa_e,
)


def two_scheme_max(p: int) -> Fraction:
    """α(n−2) = 2 · scheme_ker_max_kappa_e."""
    return alpha_particular(p) * (n_of(p) - 2)


def theorem_alpha_n_minus_2_beats_need() -> dict:
    primes = [p for p in range(5, 80) if is_prime(p)]
    rows = {}
    ok = True
    for p in primes:
        b = two_scheme_max(p)
        need = 2 - alpha_particular(p)
        # algebra: p(p^2 - 3p + 1) > 0
        poly = p * (p * p - 3 * p + 1)
        rows[str(p)] = {
            "alpha_n_minus_2": str(b),
            "need": str(need),
            "beats": b < need,
            "poly_p_p2_minus_3p_plus_1": poly,
        }
        if not (b < need and poly > 0):
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "For all primes p≥5: α(n−2)<2−α ⇔ p(p²−3p+1)>0. "
            "Hence free-e max κ_e ≤ α(n−2) ⇒ dual-eq empty."
        ),
        "by_p_sample": {k: rows[k] for k in ("5", "7", "11", "13") if k in rows},
    }


def theorem_two_over_n_farkas() -> dict:
    primes = [p for p in range(5, 60) if is_prime(p)]
    rows = {}
    ok = True
    for p in primes:
        mu = Fraction(-4, n_of(p))
        r = residual_i_farkas_if_mu(p, mu)
        poly = 2 * p**3 - 9 * p**2 + 6 * p - 3
        rows[str(p)] = {
            "mu_minus_4_over_n": str(mu),
            "mu_star": str(mu_star_threshold(p)),
            "farkas": r["need_lt_box"],
            "poly": poly,
        }
        if not (r["need_lt_box"] and poly > 0):
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "Gsum≥−4/n (i.e. |μ|≤2/n on opposite pairs) beats μ_* for all "
            "primes p≥5 ⇔ 2p³−9p²+6p−3>0."
        ),
        "by_p_sample": {k: rows[k] for k in ("5", "7", "11") if k in rows},
    }


def certify_free_e_vs_alpha_n_minus_2(p: int) -> dict:
    """Census free-e max vs α(n−2). p=3,5 full; p=7 needs cache."""
    import numpy as np
    from e1_bitight_infeas import boolean_evecs
    from minmax_quadratic import paley_conference_prime_power
    from scipy.optimize import linprog

    bound = two_scheme_max(p)
    need = 2 - alpha_particular(p)
    al = float(alpha_particular(p))

    if p == 7:
        cache = Path("/tmp/grok-goal-7663bdcc1fa8/implementer/e1_p7")
        if not (cache / "maxplus.npy").is_file():
            return {"p": p, "skipped": True}
        Mp = np.load(cache / "maxplus.npy").astype(float)
        Mm = np.load(cache / "maxminus.npy").astype(float)
        C = paley_conference_prime_power(p).astype(float)
    elif p > 5:
        return {"p": p, "skipped": True}
    else:
        C = paley_conference_prime_power(p).astype(float)
        Mp = boolean_evecs(C, float(p), 0).astype(float)
        Mm = boolean_evecs(C, float(-p), 0).astype(float)

    n = C.shape[0]
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    E = len(edges)
    e0 = edges.index((0, 1))

    def feats(M):
        return np.stack([C[i, j] * M[:, i] * M[:, j] for i, j in edges], 1)

    G = feats(Mp).T @ feats(Mp) / len(Mp) + feats(Mm).T @ feats(Mm) / len(Mm)
    ev, vec = np.linalg.eigh(G)
    K = vec[:, ev < 1e-8]
    ones = np.ones(E)
    sn = ones @ K
    ns = np.linalg.svd([sn], full_matrices=True)[2][1:].T
    N = K @ ns
    free = N.shape[1]
    rows, bs = [], []
    for i in range(E):
        if i == e0:
            continue
        rows += [N[i], -N[i]]
        bs += [1 - al, al]
    res = linprog(
        -N[e0],
        A_ub=np.array(rows),
        b_ub=np.array(bs),
        method="highs",
        bounds=[(None, None)] * free,
    )
    if not res.success:
        return {"p": p, "success": False, "message": res.message}
    vmax = float(N[e0] @ res.x)
    return {
        "p": p,
        "success": True,
        "max_kappa_e": vmax,
        "alpha_n_minus_2": float(bound),
        "need": float(need),
        "max_le_alpha_n_minus_2": vmax <= float(bound) + 1e-9,
        "blocks_dual_eq": vmax < float(need) - 1e-9,
        "ratio_max_over_scheme": vmax / float(scheme_ker_max_kappa_e(p)),
        "proved_general": False,
    }


def certify_mu_p7() -> dict:
    """Exact max|μ| at p=7 from Max± cache."""
    import numpy as np
    from itertools import combinations
    from minmax_quadratic import paley_conference_prime_power

    cache = Path("/tmp/grok-goal-7663bdcc1fa8/implementer/e1_p7")
    if not (cache / "maxplus.npy").is_file():
        return {"skipped": True}
    Mp = np.load(cache / "maxplus.npy")
    Mm = np.load(cache / "maxminus.npy")
    p, n = 7, 50
    C = paley_conference_prime_power(p).astype(np.int8)
    max_mu = 0.0
    for a, b in combinations(range(2, n), 2):
        i, j, k, l = 0, 1, a, b
        kap = int(C[i, j] * C[k, l] + C[i, k] * C[j, l] + C[i, l] * C[j, k])
        if abs(kap) != 1:
            continue
        mu = 0.5 * (
            np.mean(Mp[:, i] * Mp[:, j] * Mp[:, k] * Mp[:, l])
            + np.mean(Mm[:, i] * Mm[:, j] * Mm[:, k] * Mm[:, l])
        )
        max_mu = max(max_mu, abs(mu))
    max_frac = Fraction(max_mu).limit_denominator(100000)
    return {
        "p": 7,
        "max_abs_mu": str(max_frac),
        "max_abs_mu_float": float(max_frac),
        "two_over_n": str(Fraction(2, n)),
        "mu_le_two_over_n": max_frac <= Fraction(2, n),
        "known_exact": "109/2863",
        "match_known": max_frac == Fraction(109, 2863),
        "proved_general": False,
    }


def free_e_bound_proved_general() -> bool:
    return False


def residual_i_closed_via_201() -> bool:
    return free_e_bound_proved_general()


def hinge_status_201() -> dict:
    return {
        "alpha_n_minus_2_beats_need_p_ge_5": True,
        "two_over_n_farkas_sufficient": True,
        "free_e_max_le_alpha_n_minus_2_general": False,
        "residual_i_closed": residual_i_closed_via_201(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "open": (
            "Prove free-e max κ_e ≤ α(n−2) (or ≤(3/2)scheme / |μ|≤2/n / "
            "K₄≤Wick_hi) for all primes p≥5."
        ),
    }


def main() -> dict:
    a = theorem_alpha_n_minus_2_beats_need()
    b = theorem_two_over_n_farkas()
    cert_box = {str(p): certify_free_e_vs_alpha_n_minus_2(p) for p in (3, 5, 7)}
    cert_mu = certify_mu_p7()
    out = {
        "title": (
            "Prop 15.201 α(n−2) sufficient free-e bound; p=7 μ=109/2863; "
            "residual (i) OPEN"
        ),
        "proved": {
            "alpha_n_minus_2_beats_need": a["proved"],
            "two_over_n_farkas": b["proved"],
            "free_e_bound_general": False,
            "residual_i_closed": residual_i_closed_via_201(),
            "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
            "E1_closed": e1_closed_general(),
            "L_closed": False,
        },
        "algebra": {"alpha_bound": a, "farkas_2n": b},
        "certified": {"free_e_box": cert_box, "mu_p7": cert_mu},
        "hinge": hinge_status_201(),
        "F3": (
            "α(n−2)<2−α is Max+-free. Do not flip residual/E1/L without "
            "general free-e max ≤ α(n−2) or equivalent."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15201.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.201 alpha(n-2) sufficient free-e bound")
    print(f"  alpha(n-2)<2-alpha for p>=5: {a['proved']}")
    print(f"  |mu|<=2/n Farkas: {b['proved']}")
    for pk, v in cert_box.items():
        if v.get("skipped") or not v.get("success"):
            print(f"  p={pk}: skipped/fail")
            continue
        print(
            f"  p={pk}: max_ke={v['max_kappa_e']:.4f} "
            f"alpha(n-2)={v['alpha_n_minus_2']:.4f} "
            f"le={v['max_le_alpha_n_minus_2']} blocks={v['blocks_dual_eq']}"
        )
    if not cert_mu.get("skipped"):
        print(f"  p=7 max|mu|={cert_mu['max_abs_mu']} le 2/n={cert_mu['mu_le_two_over_n']}")
    print(f"  residual_i closed: {residual_i_closed_via_201()}")
    print(f"  open: {out['hinge']['open']}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
