#!/usr/bin/env python3
"""
Prop 15.191 — Cy-expansion of m₄ on |κ|=1; f4 majorant path; residual (i) OPEN.

Does **not** flip gsum_disj_lb_proved_general. Soft-close forbidden.

SETUP
  C symmetric conference, Cy=±p y on Max±.  Four-set S, |κ(S)|=1,
  φ(S)=∑_{r∉S} ∏_{i∈S} C_{ri},  μ₄=½(E₊+E₋)[∏_{i∈S} y_i].
  Target: |μ₄|≤1/(2p) (⇔ Gsum≥−1/p) or the sufficient |μ₄|≤2/n.

PROVED Max+-free
  A. Derangement permanent of C[S,S] on |κ|=1.
     For every zero-diagonal symmetric ±1 K4 labelling with |κ|=1,
        ∑_{σ derangement of S} ∏_{i∈S} C_{i,σ(i)} = 1.
     (Exhaustion of 64 edge-labelings; |κ|=3 gives 9.)  ∎

  B. Star-sum vanishes on |κ|=1.
     ∑_{s∈S} ∏_{i∈S\\{s}} C_{is} = 0 whenever |κ|=1.
     (Same 64-exhaustion; |κ|=3 gives ±4.)  ∎

  C. Size-(2,2) contribution to E[∏_{i∈S}(Cy)_i] equals −3φ.
     Partition domain S into two pairs (A,B); map A→u, B→v (u≠v);
     moment E[y_u² y_v²]=1.  For each ordered pair-partition,
        ∑_{u≠v}(∏_{i∈A} C_{iu})(∏_{j∈B} C_{jv})
          = (∑_u ∏_A C_{iu})(∑_v ∏_B C_{jv}) − ∑_u ∏_{i∈S} C_{iu}
          = 0 − φ(S)
     by C²=(n−1)I (off-diagonal of C² vanishes).  Six ordered partitions
     double-count each (2,2)-map once each way ⇒ total weight −3φ.  ∎
     Holds for any symmetric conference (no Max+).

  D. Size-(3,1) contribution vanishes on |κ|=1 (Paley + π).
     Singleton s→v, triple A→u (u≠v); Max+ moment π_{uv}=C_{uv}/p (15.189).
     C² reduction yields a multiple of the star-sum (B), hence 0 on |κ|=1.
     Max− with π=−C/p: same star-sum factor ⇒ 0.  ∎

  E. Size-1 maps contribute φ (moment 1).  Combined with (C)+(D):
        size1 + size2 = φ − 3φ = −2φ
     on every |κ|=1 four-set.  ∎

  F. f4 majorant (Fraction, Max+-free under |φ|≤2(p−2) from 15.186).
     Define f4(κ,φ):=(4κ−φ)/(p n).  On |κ|=1 with |φ|≤2(p−2):
        |f4| ≤ (4 + 2(p−2))/(p n) = 2p/(p n) = 2/n.
     And 2/n = 2/(p²+1) ≤ 1/(2p) for all primes p≥5
     (⇔ 4p ≤ p²+1 ⇔ p²−4p+1≥0, true on [5,∞)).  ∎
     Hence |μ|≤|f4| on |κ|=1 ⇒ |μ|≤2/n≤1/(2p) ⇒ Gsum≥−1/p
     ⇒ residual-(i) Farkas (15.176).

CERTIFIED (not general)
  G. Size-3 contribution to the Cy-expansion is 0 at p=3,5 (full enum of
     image-size-3 maps; unconstrained C² sum already 0; collision corrections
     cancel on |κ|=1).  Sample-consistent at p=7.
  H. p^4 m₄⁺ = −2φ + Q₄(m₄⁺) at p=5 on |κ|=1 sample (Q₄ = injection /
     Per contribution), and m₄⁺=m₄⁻ there, so p^4 μ = −2φ + Q₄(μ).
  I. Pointwise |μ| vs |f4| / 2/n on |κ|=1:
        p=3: |μ|≰|f4| and |μ|≰2/n (max|μ|=1/3; out of residual scope);
        p=5: μ=f4 pointwise on each (κ,φ) class; max|μ|=3/65 < 2/n=1/13;
        p=7: |μ|≰|f4| on many classes (f4 is **not** a pointwise majorant);
             but max|μ|=109/2863 < 2/n=1/25 < 1/14.
     Viable sufficient target remains |μ|≤2/n (not |μ|≤|f4|).
  J. Full dual-eq ker-box empty at p=5,7 (15.190 census); feasible at p=3.

OPEN (blocks residual i / predicate flip)
  K. Prove |μ₄| ≤ 2/n (or ≤1/(2p)) on |κ|=1 for all primes p≥5 (Max+-free),
     **or** full dual-eq ker-box empty for all p≥5.
     Then flip gsum_disj_lb_proved_general → residual-(i) → E1 → L.
     Note: |μ|≤|f4| fails at p=7; do not use f4 as a general majorant.

Until K: gsum_disj_lb_proved_general False; residual_i False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15191.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from itertools import combinations, permutations, product
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    is_prime,
    n_of,
)
from e1_gmin_m4_prop15177 import mu4_threshold  # noqa: E402
from e1_gmin_m4_prop15186 import theorem_phi_bound_general  # noqa: E402
from e1_gmin_m4_prop15188 import two_over_n  # noqa: E402


def kappa_of_labeling(C6: tuple[int, ...]) -> int:
    """κ from 6 upper-triangle ±1 edge labels of K4 (order 01,02,03,12,13,23)."""
    c01, c02, c03, c12, c13, c23 = C6
    return c01 * c23 + c02 * c13 + c03 * c12


def matrix_from_labeling(C6: tuple[int, ...]):
    C = [[0] * 4 for _ in range(4)]
    e = 0
    for i in range(4):
        for j in range(i + 1, 4):
            C[i][j] = C[j][i] = C6[e]
            e += 1
    return C


def derangement_permanent(C) -> int:
    total = 0
    for sig in permutations(range(4)):
        if any(sig[i] == i for i in range(4)):
            continue
        w = 1
        for i in range(4):
            w *= C[i][sig[i]]
        total += w
    return total


def star_sum(C) -> int:
    total = 0
    for s in range(4):
        pr = 1
        for i in range(4):
            if i != s:
                pr *= C[i][s]
        total += pr
    return total


def theorem_derangement_perm_on_kappa1() -> dict:
    """Max+-free 64-exhaust: derang_perm ≡ 1 on |κ|=1."""
    by = {}
    ok = True
    for bits in product((-1, 1), repeat=6):
        C = matrix_from_labeling(bits)
        k = kappa_of_labeling(bits)
        d = derangement_permanent(C)
        by[(k, d)] = by.get((k, d), 0) + 1
        if abs(k) == 1 and d != 1:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "On every zero-diag symmetric ±1 K4 with |κ|=1, "
            "the derangement permanent of C[S,S] equals 1. "
            "(|κ|=3 ⇒ permanent 9.) Exhaustion of 64 labelings."
        ),
        "multiplicities": {f"kappa={k},derang={d}": c for (k, d), c in sorted(by.items())},
    }


def theorem_star_sum_vanishes_on_kappa1() -> dict:
    """Max+-free 64-exhaust: star-sum ≡ 0 on |κ|=1."""
    by = {}
    ok = True
    for bits in product((-1, 1), repeat=6):
        C = matrix_from_labeling(bits)
        k = kappa_of_labeling(bits)
        st = star_sum(C)
        by[(k, st)] = by.get((k, st), 0) + 1
        if abs(k) == 1 and st != 0:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "On every zero-diag symmetric ±1 K4 with |κ|=1, "
            "∑_{s∈S} ∏_{i≠s} C_{is} = 0. Exhaustion of 64 labelings."
        ),
        "multiplicities": {f"kappa={k},star={st}": c for (k, st), c in sorted(by.items())},
    }


def theorem_size22_contribution() -> dict:
    return {
        "proved": True,
        "theorem": (
            "Any symmetric conference C²=(n−1)I: the (2,2)-image contribution "
            "to E[∏_{i∈S}(Cy)_i] equals −3 φ(S) (six ordered pair-partitions "
            "each give −φ, each map counted twice; moment 1). No Max+ used."
        ),
    }


def theorem_size31_vanishes_on_kappa1() -> dict:
    return {
        "proved": True,
        "theorem": (
            "On |κ|=1, (3,1)-image contribution with π_±=±C/p reduces by C² "
            "to a multiple of the star-sum, hence vanishes by B. Uses Paley π "
            "(15.189) and star-sum (B)."
        ),
        "depends_on": ["star_sum_kappa1", "pairwise_pi"],
    }


def theorem_size1_plus_size2() -> dict:
    return {
        "proved": True,
        "theorem": (
            "On |κ|=1 (Paley): size1+size2 contribution to E[∏(Cy)_i] equals "
            "−2 φ(S).  (size1=φ; size2=−3φ from (C)+(D).)"
        ),
        "depends_on": ["size22", "size31", "size1"],
    }


def f4_formula(p: int, kappa: int, phi: int) -> Fraction:
    """(4κ − φ)/(p n)."""
    return Fraction(4 * kappa - phi, p * n_of(p))


def theorem_f4_majorant() -> dict:
    """|f4| ≤ 2/n ≤ 1/(2p) for p≥5 under |φ|≤2(p−2), |κ|=1."""
    primes = [p for p in range(5, 80) if is_prime(p)]
    rows = {}
    ok = True
    for p in primes:
        n = n_of(p)
        # worst |4κ−φ| = 4 + 2(p−2) = 2p
        worst = Fraction(2 * p, p * n)  # = 2/n
        two_n = Fraction(2, n)
        thr = mu4_threshold(p)
        le_two_n = worst == two_n
        le_thr = two_n <= thr
        rows[str(p)] = {
            "f4_maj": str(worst),
            "two_over_n": str(two_n),
            "one_over_2p": str(thr),
            "f4_maj_eq_two_n": le_two_n,
            "two_n_le_thr": le_thr,
        }
        if not (le_two_n and le_thr):
            ok = False
    # algebraic check 2/n ≤ 1/(2p) ⇔ 4p ≤ p²+1
    alg_ok = all(4 * p <= p * p + 1 for p in primes)
    return {
        "proved": ok and alg_ok and theorem_phi_bound_general()["proved"],
        "theorem": (
            "On |κ|=1 with |φ|≤2(p−2): |4κ−φ|/(pn) ≤ 2/n ≤ 1/(2p) for all "
            "primes p≥5.  Hence |μ|≤|f4| ⇒ residual-(i) Farkas via 15.176."
        ),
        "by_p_sample": {k: rows[k] for k in ("5", "7", "11", "13") if k in rows},
        "two_over_n_le_half_over_p": alg_ok,
        "depends_on": ["phi_bound_15_186"],
    }


def implication_two_n_bound_closes_residual_i() -> dict:
    return {
        "proved_implication": True,
        "bound_proved_general": False,
        "theorem": (
            "IF |μ₄|≤2/n on every |κ|=1 four-set for all primes p≥5, "
            "THEN |μ₄|≤1/(2p) ⇒ Gsum≥−1/p ⇒ residual-(i) dual-eq Farkas "
            "⇒ gsum_disj_lb_proved_general may flip. "
            "(|μ|≤|f4| is sufficient when it holds but fails at p=7.)"
        ),
    }


def certify_mu_le_f4(p: int = 5) -> dict:
    """Census: |μ| ≤ |f4| on |κ|=1."""
    import numpy as np
    from e1_bitight_infeas import boolean_evecs
    from minmax_quadratic import paley_conference_prime_power

    C = paley_conference_prime_power(p).astype(float)
    n = C.shape[0]
    if p == 7:
        cache = Path("/tmp/grok-goal-03c774661dc8/implementer/e1_p7")
        mp_path, mm_path = cache / "maxplus.npy", cache / "maxminus.npy"
        if not (mp_path.is_file() and mm_path.is_file()):
            return {"p": p, "skipped": True, "reason": "no Max± cache"}
        Mp = np.load(mp_path).astype(float)
        Mm = np.load(mm_path).astype(float)
    else:
        Mp = boolean_evecs(C, float(p), 0).astype(float)
        Mm = boolean_evecs(C, -float(p), 0).astype(float)

    viol = 0
    n_k1 = 0
    max_mu = 0.0
    max_ratio = 0.0
    for S in combinations(range(n), 4):
        a, b, c, d = S
        kap = C[a, b] * C[c, d] + C[a, c] * C[b, d] + C[a, d] * C[b, c]
        if abs(abs(kap) - 1) > 1e-9:
            continue
        n_k1 += 1
        ph = sum(C[r, a] * C[r, b] * C[r, c] * C[r, d] for r in range(n))
        mu = 0.5 * (
            float(np.mean(Mp[:, a] * Mp[:, b] * Mp[:, c] * Mp[:, d]))
            + float(np.mean(Mm[:, a] * Mm[:, b] * Mm[:, c] * Mm[:, d]))
        )
        f4 = float((4 * kap - ph) / (p * n))
        max_mu = max(max_mu, abs(mu))
        if abs(f4) > 1e-15:
            max_ratio = max(max_ratio, abs(mu) / abs(f4))
        if abs(mu) > abs(f4) + 1e-9:
            viol += 1
    two_n = 2.0 / n
    thr = 1.0 / (2 * p)
    return {
        "p": p,
        "n_kappa1": n_k1,
        "max_abs_mu": max_mu,
        "two_over_n": two_n,
        "one_over_2p": thr,
        "max_abs_mu_over_f4": max_ratio,
        "violations_mu_gt_f4": viol,
        "mu_le_f4": viol == 0,
        "mu_le_two_n": max_mu <= two_n + 1e-12,
        "mu_le_thr": max_mu <= thr + 1e-12,
        "proved_general": False,
    }


def hinge_status_191() -> dict:
    return {
        "residual_ii_closed": False,  # full open (15.193); affine branch closed by 15.179
        "residual_i_closed": False,
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "f4_envelope_le_two_n_proved": True,
        "mu_le_two_n_proved_general": False,
        "mu_le_f4_fails_at_p7": True,
        "open": (
            "Prove |μ|≤2/n (or ≤1/(2p)) on |κ|=1 for all p≥5 Max+-free; "
            "or full dual-eq ker-box empty for all p≥5. "
            "Do not use |μ|≤|f4| (false at p=7)."
        ),
    }


def main() -> dict:
    der = theorem_derangement_perm_on_kappa1()
    star = theorem_star_sum_vanishes_on_kappa1()
    s22 = theorem_size22_contribution()
    s31 = theorem_size31_vanishes_on_kappa1()
    s12 = theorem_size1_plus_size2()
    f4m = theorem_f4_majorant()
    impl = implication_two_n_bound_closes_residual_i()
    cert = {str(p): certify_mu_le_f4(p) for p in (3, 5)}
    c7 = certify_mu_le_f4(7)
    if not c7.get("skipped"):
        cert["7"] = c7
    out = {
        "title": (
            "Prop 15.191 Cy-expansion on |κ|=1; 2/n target; residual (i) OPEN"
        ),
        "proved": {
            "derangement_perm_kappa1": der["proved"],
            "star_sum_kappa1": star["proved"],
            "size22_minus_3phi": s22["proved"],
            "size31_vanishes_kappa1": s31["proved"],
            "size1_plus_size2_minus_2phi": s12["proved"],
            "f4_envelope_le_two_n_le_half_over_p": f4m["proved"],
            "mu_le_two_n_general": False,
            "mu_le_f4_general": False,
            "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
            "residual_i_closed": False,
            "E1_closed": e1_closed_general(),
            "L_closed": False,
        },
        "algebra": {
            "derangement_perm": der,
            "star_sum": star,
            "size22": s22,
            "size31": s31,
            "size1_plus_size2": s12,
            "f4_envelope": f4m,
            "implication": impl,
        },
        "certified": {"mu_vs_f4_and_two_n": cert},
        "hinge": hinge_status_191(),
        "F3": (
            "Viable target |μ|≤2/n (≤1/(2p) for p≥5). |μ|≤|f4| fails at p=7. "
            "Soft-close forbidden; census is not a general proof."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15191.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.191 Cy-expansion; 2/n target; residual (i) OPEN")
    print(f"  derang_perm=1 on |κ|=1 proved={der['proved']}")
    print(f"  star_sum=0 on |κ|=1 proved={star['proved']}")
    print(f"  size1+size2=-2φ proved={s12['proved']}")
    print(f"  f4 envelope ≤2/n≤1/(2p) proved={f4m['proved']}")
    for k, v in cert.items():
        if v.get("skipped"):
            print(f"  census p={k}: skipped")
            continue
        print(
            f"  census p={k}: max|μ|={v['max_abs_mu']:.6f} "
            f"μ≤|f4|={v['mu_le_f4']} μ≤2/n={v['mu_le_two_n']} "
            f"f4_viol={v['violations_mu_gt_f4']}"
        )
    print(f"  gsum={gsum_disj_lb_proved_general()}; e1={e1_closed_general()}")
    print(f"  open: {out['hinge']['open']}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
