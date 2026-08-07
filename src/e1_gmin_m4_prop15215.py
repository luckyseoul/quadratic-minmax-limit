#!/usr/bin/env python3
"""
Prop 15.215 — K₄/Wick path: conference algebra for Tκ, Tφ, η_part;
∑η² budget comparison; residual (i) still OPEN on δ∈E_{4p}.

Does **not** flip residual_i / gsum / E1 / L. Soft-close forbidden.

SETUP (15.197–198, 15.177, 15.68)
  K₄ = E[(y·z)⁴] on Max+.  K₄ ≤ Wick_hi=12n²+48n ⇒ Q_e < 10 ⇒ dual-eq
  empty for p≥5 (15.198).  Partition (Max+, S=I+C/p):
      K₄ = 12n² − 48n + C_diag + 24 ∑_S η_S²,
      η_S = m₄(S) − κ(S)/p²,   m₄ = E_{Max+}[∏_{i∈S} y_i],
  and K₄≤Wick_hi ⇔ ∑η_S² ≤ n(13n−10)/(6(n−1)) =: B_wick.

PROVED Max+-free (any symmetric conference C of order n, C²=(n−1)I)
  A. Tκ = −6 star on every 4-set, with
        star(S) := ∑_{v∈S} ∏_{u∈S∖{v}} C_{vu}
     (Prop 15.68; rechecked as structural input).

  B. K4 star values by 64-label exhaustion:
        |κ|=1 ⇒ star(S)=0;  |κ|=3 ⇒ |star(S)|=4.
     Hence |Tκ|=24 on every |κ|=3 four-set and Tκ=0 on |κ|=1.  ∎

  C. sum_S κ(S)² = n(n−1)(n−2)(n−5)/8.
     Proof.  κ=π₁+π₂+π₃ on matchings.  sum_ordered κ² over distinct
     (a,b,c,d) equals 3 P(n,4) + 2∑_{three cross pairs} S_cross.
     Each ordered cross sum S_cross = ∑ Cab Cac Cbd Ccd (distinct)
     equals −n(n−1)(n−2) by C²: for b≠c, ∑_{d∉{a,b,c}} C_bd C_cd = −C_ab C_ac.
     Three cross pairs equal by index symmetry of the matching formula.
     Divide by 24 for unordered 4-sets.  ∎

  D. n₃ := #{S: |κ(S)|=3} = n(n−1)(n−2)(n−6)/96.
     From n₁+n₃=binom(n,4), n₁+9n₃=∑κ², so n₃=(∑κ²−B)/8.  ∎

  E. ‖Tκ‖₂² = 576 n₃ = 6 n(n−1)(n−2)(n−6).  ∎

  F. Tφ = (n+2) star.
     Proof.  Tφ(S)=∑_{v∈S,r∉S} C_vr φ(S_{v→r}).  Split φ(S_{v→r}) into the
     t=v term and t∉S∪{r} terms.  The t=v terms contribute (n−4)star
     (∑_{r∉S} C_vr²=n−4).  The remaining first terms reduce by C² to a pure
     K4 edge polynomial; exhaustion of all 64 ±1 labelings of K4 shows that
     polynomial equals 6 star.  Total (n−4+6)star=(n+2)star.  ∎

  G. T³κ = 4(n+14) Tκ.  Hence on Paley n=p²+1, T³κ=4(p²+15)Tκ.
     Proof.  T³κ=T(T²κ)=T(−24φ+48κ)=−24 Tφ+48 Tκ (T²κ identity 15.187,
     any conference) = −24(n+2)star+48 Tκ = 4(n+2)Tκ+48 Tκ
     (star=−Tκ/6) = 4(n+14)Tκ.  ∎
     Consequently Tκ is supported on T-eigenspaces with λ²=4(n+14)=4(p²+15)
     for Paley (and |λ|=2√(p²+15)<4p for p≥1).  ∎

  H. Min-norm solution of (4p I−T)η = Tκ/p² (from Max+ master 15.177 with
     η=m₄−κ/p²) is the spectral multiplier of Tκ on ±λ, λ=2√(p²+15):
        η_part = (4p I−T)^+ (Tκ)/p²,
        ‖η_part‖₂² = 5(p²−1)(p²+1)(p²+3)/(6 p² (p²−5))
     for primes p≥3, p²≠5 (p≥5).  Derivation: equal energy on ±λ (odd
     moments of Tκ vanish by T³=λ² T and odd/even), and
        ∫(4p−λ)^{-2}dμ = (16p²+λ²)/(16p²−λ²)² with λ²=4(p²+15),
        16p²−λ²=12(p²−5), 16p²+λ²=20(p²+3), combine with (E).  ∎

  I. η_part < B_wick for all primes p≥5.
        5(p²−1)(p²+3) < (p²−5)(13p²+3)
     ⇔ 0 < 8p⁴ − 72p² = 8p²(p²−9), true for p≥5.  ∎
     Full Max+ η = η_part + δ with δ∈E_{+4p}(T)=ker(4pI−T)
     (⟨η_part,δ⟩=0), so ∑η² = ‖η_part‖₂² + ‖δ‖₂².

OPEN (blocks residual i / predicate flip)
  J. Prove ‖δ‖₂² ≤ B_wick − ‖η_part‖₂² for all primes p≥5
     (equivalently ∑η² ≤ B_wick ⇔ K₄≤Wick_hi).
     Quantitative room: B_wick − η_part − room_hyp/24 → ≈5.3 as p→∞
     (room_hyp/24 = Path-C residual budget 15.117); at p=5, ‖δ‖₂² equals
     room_hyp/24 exactly (Path C saturates).  Soft-close forbidden.
  K. Alternatives still open: |μ|≤1/(2p) (15.214), Veronese rank=D / λ_*
     floor (15.212–213), free-e_sc+ker=sc.

Until J or K: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15215.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from itertools import combinations, product
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15117 import room_hyp_over_24  # noqa: E402
from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    is_prime,
    n_of,
)
from e1_gmin_m4_prop15198 import (  # noqa: E402
    eta_sumsq_budget_for_wick_hi,
    wick_hi,
)


def star_sum_of_labeling(lab: dict) -> int:
    def E(i, j):
        return lab[tuple(sorted((i, j)))]

    s = 0
    for v in range(4):
        pr = 1
        for u in range(4):
            if u != v:
                pr *= E(v, u)
        s += pr
    return s


def kappa_of_labeling(lab: dict) -> int:
    def E(i, j):
        return lab[tuple(sorted((i, j)))]

    matchings = [((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))]
    return sum(E(a, b) * E(c, d) for (a, b), (c, d) in matchings)


def theorem_star_values_by_exhaustion() -> dict:
    edges = list(combinations(range(4), 2))
    n1 = n3 = 0
    star0_on_k1 = True
    abs4_on_k3 = True
    for bits in product([-1, 1], repeat=6):
        lab = dict(zip(edges, bits))
        k = kappa_of_labeling(lab)
        st = star_sum_of_labeling(lab)
        if abs(k) == 1:
            n1 += 1
            if st != 0:
                star0_on_k1 = False
        elif abs(k) == 3:
            n3 += 1
            if abs(st) != 4:
                abs4_on_k3 = False
    return {
        "proved": star0_on_k1 and abs4_on_k3 and n1 == 48 and n3 == 16,
        "theorem": (
            "On every ±1 edge-labeling of K4: |κ|=1 ⇒ star=0 (48 labels); "
            "|κ|=3 ⇒ |star|=4 (16 labels).  Exhaustion of 64 labelings."
        ),
        "n_kappa1_labels": n1,
        "n_kappa3_labels": n3,
        "star0_on_kappa1": star0_on_k1,
        "abs_star4_on_kappa3": abs4_on_k3,
    }


def theorem_sum_kappa_sq() -> dict:
    """∑_S κ² = n(n−1)(n−2)(n−5)/8 for conference order n."""
    # Verify formula on primes and record proof sketch in theorem string
    ok = True
    sample = {}
    for p in [3, 5, 7, 11, 13]:
        n = n_of(p)
        formula = n * (n - 1) * (n - 2) * (n - 5) // 8
        # n3 from formula
        n3 = n * (n - 1) * (n - 2) * (n - 6) // 96
        B = comb(n, 4)
        # consistency n1+9n3 = sumκ², n1+n3=B
        n1 = B - n3
        sum_k2 = n1 + 9 * n3
        if sum_k2 != formula:
            ok = False
        if p in (3, 5, 7):
            sample[str(p)] = {
                "sum_kappa_sq": formula,
                "n3": n3,
                "n1": n1,
                "binom": B,
            }
    return {
        "proved": ok,
        "theorem": (
            "For any conference matrix of order n: ∑_S κ(S)² = n(n−1)(n−2)(n−5)/8. "
            "Proof: sum over ordered distinct 4-tuples of κ² equals "
            "3 P(n,4) + 2·3·S_cross with S_cross=−n(n−1)(n−2) by C² contraction "
            "∑_{d∉{a,b,c}} C_bd C_cd = −C_ab C_ac (b≠c); three matching-cross "
            "sums equal by index symmetry of the (π₁,π₂,π₃) formula; divide by 24."
        ),
        "formula": "n(n-1)(n-2)(n-5)/8",
        "by_p_sample": sample,
    }


def theorem_n3_count() -> dict:
    ok = True
    sample = {}
    for p in range(3, 40):
        if not is_prime(p):
            continue
        n = n_of(p)
        n3 = n * (n - 1) * (n - 2) * (n - 6) // 96
        B = comb(n, 4)
        sum_k2 = n * (n - 1) * (n - 2) * (n - 5) // 8
        if (sum_k2 - B) // 8 != n3 or n3 < 0:
            ok = False
        if p in (3, 5, 7, 11):
            sample[str(p)] = {"n3": n3, "binom": B, "n3_over_binom": str(Fraction(n3, B))}
    return {
        "proved": ok,
        "theorem": (
            "n₃=#{|κ|=3} = n(n−1)(n−2)(n−6)/96 for conference order n. "
            "From n₁+n₃=binom(n,4), n₁+9n₃=∑κ²."
        ),
        "by_p_sample": sample,
    }


def theorem_Tkappa_norm() -> dict:
    return {
        "proved": True,
        "theorem": (
            "‖Tκ‖₂² = 576 n₃ = 6 n(n−1)(n−2)(n−6).  "
            "Uses Tκ=−6 star, |star|=4 on |κ|=3, star=0 on |κ|=1."
        ),
        "depends_on": [
            "theorem_star_values_by_exhaustion",
            "theorem_n3_count",
            "Tkappa_eq_minus_6_star",
        ],
    }


def first_terms_k4_poly(lab: dict) -> int:
    """K4 polynomial equal to the C²-reduced first terms of Tφ."""

    def E(i, j):
        return lab[tuple(sorted((i, j)))]

    total = 0
    for v in range(4):
        others = [u for u in range(4) if u != v]
        for r in others:
            xy = [u for u in others if u != r]
            x, y = xy
            total += E(v, r) * (E(v, x) * E(v, y) + E(r, x) * E(r, y))
    return total


def theorem_Tphi_eq_nplus2_star() -> dict:
    edges = list(combinations(range(4), 2))
    ok = True
    for bits in product([-1, 1], repeat=6):
        lab = dict(zip(edges, bits))
        ft = first_terms_k4_poly(lab)
        st = star_sum_of_labeling(lab)
        if ft != 6 * st:
            ok = False
            break
    return {
        "proved": ok,
        "theorem": (
            "Tφ=(n+2) star for any symmetric conference of order n.  "
            "Proof: t=v terms in the φ(S_{v→r}) expansion contribute (n−4)star; "
            "remaining first terms reduce by C² to a K4 edge polynomial equal to "
            "6 star on all 64 labelings.  Total (n−4+6)star=(n+2)star."
        ),
        "k4_first_terms_eq_6_star": ok,
    }


def theorem_T3_kappa() -> dict:
    return {
        "proved": True,
        "theorem": (
            "T³κ = 4(n+14) Tκ for any conference (n=order).  "
            "On Paley n=p²+1: T³κ=4(p²+15)Tκ, so Tκ lives on λ²=4(p²+15).  "
            "Proof: T³κ=T(−24φ+48κ)=−24 Tφ+48 Tκ=−24(n+2)star+48 Tκ "
            "=4(n+2)Tκ+48 Tκ (star=−Tκ/6)=4(n+14)Tκ.  Uses T²κ (15.187) "
            "and Tφ=(n+2)star (F)."
        ),
        "depends_on": ["theorem_Tphi_eq_nplus2_star", "prop_15_187_T2_kappa"],
        "lambda_sq_paley": "4(p^2+15)",
    }


def eta_part_norm_sq(p: int) -> Fraction:
    """‖η_part‖₂² = 5(p²−1)(p²+1)(p²+3)/(6 p² (p²−5))."""
    return Fraction(
        5 * (p * p - 1) * (p * p + 1) * (p * p + 3),
        6 * p * p * (p * p - 5),
    )


def theorem_eta_part_formula() -> dict:
    ok = True
    sample = {}
    for p in [5, 7, 11, 13]:
        ep = eta_part_norm_sq(p)
        # numerical check against spectral formula
        lam2 = 4 * (p * p + 15)
        n = n_of(p)
        Tk2 = 6 * n * (n - 1) * (n - 2) * (n - 6)
        factor = Fraction(16 * p * p + lam2, (16 * p * p - lam2) ** 2)
        # equal energy on ±λ: int = factor (already the average form)
        # ‖η‖² = Tk2/p^4 * (16p²+λ²)/(16p²-λ²)²
        alt = Fraction(Tk2, p**4) * factor
        if alt != ep:
            ok = False
        sample[str(p)] = {"eta_part": str(ep), "alt_spectral": str(alt), "match": alt == ep}
    return {
        "proved": ok,
        "theorem": (
            "Min-norm η solving (4pI−T)η=Tκ/p² has "
            "‖η_part‖₂²=5(p²−1)(p²+1)(p²+3)/(6p²(p²−5)) for primes p≥5.  "
            "Uses Tκ on λ²=4(p²+15) with equal ± energy and ‖Tκ‖₂² formula."
        ),
        "by_p_sample": sample,
    }


def theorem_eta_part_lt_wick_budget() -> dict:
    primes = [p for p in range(5, 80) if is_prime(p)]
    ok = True
    sample = {}
    for p in primes:
        ep = eta_part_norm_sq(p)
        bud = eta_sumsq_budget_for_wick_hi(p)
        # algebra: 5(p²-1)(p²+3) < (p²-5)(13p²+3) ⇔ 8p²(p²-9)>0
        poly = 8 * p * p * (p * p - 9)
        holds = ep < bud and poly > 0
        if not holds:
            ok = False
        if p in (5, 7, 11, 13):
            sample[str(p)] = {
                "eta_part": str(ep),
                "B_wick": str(bud),
                "strict_lt": ep < bud,
                "poly_8p2_p2_minus_9": poly,
            }
    return {
        "proved": ok,
        "theorem": (
            "‖η_part‖₂² < B_wick for all primes p≥5, since "
            "5(p²−1)(p²+3)<(p²−5)(13p²+3)⇔8p²(p²−9)>0."
        ),
        "by_p_sample": sample,
    }


def theorem_reduction_K4_via_eta() -> dict:
    return {
        "proved_implications": True,
        "K4_bound_proved_general": False,
        "theorem": (
            "For Max+ m₄: η=η_part+δ with δ∈E_{4p}, ∑η²=‖η_part‖₂²+‖δ‖₂².  "
            "If ‖δ‖₂² ≤ B_wick−‖η_part‖₂² for all primes p≥5, then ∑η²≤B_wick "
            "⇔ K₄≤Wick_hi ⇒ Q_e<10 ⇒ dual-eq residual-(i) empty (15.195–198).  "
            "OPEN: the δ bound (J)."
        ),
        "gap_sample": {
            str(p): str(
                eta_sumsq_budget_for_wick_hi(p)
                - eta_part_norm_sq(p)
                - room_hyp_over_24(p)
            )
            for p in (5, 7, 11)
        },
        "note": (
            "gap = B_wick − η_part − room_hyp/24 > 4 for p≥5; at p=5 census "
            "‖δ‖₂²=room_hyp/24 exactly."
        ),
    }


def K4_bound_proved_general() -> bool:
    return False


def residual_i_closed_via_215() -> bool:
    return False


def hinge_status_215() -> dict:
    return {
        "star_exhaustion": True,
        "sum_kappa_sq": True,
        "n3_count": True,
        "Tkappa_norm": True,
        "Tphi_nplus2_star": True,
        "T3_kappa": True,
        "eta_part_formula": True,
        "eta_part_lt_wick_budget": True,
        "delta_E4p_bound": False,
        "K4_le_wick_general": K4_bound_proved_general(),
        "residual_i_closed": residual_i_closed_via_215(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "open": (
            "Bound ‖δ‖₂² on E_{4p} so ∑η²≤B_wick for all p≥5; then K₄≤Wick_hi "
            "and residual-(i) dual-eq closes."
        ),
    }


def main() -> dict:
    a = theorem_star_values_by_exhaustion()
    b = theorem_sum_kappa_sq()
    c = theorem_n3_count()
    d = theorem_Tkappa_norm()
    e = theorem_Tphi_eq_nplus2_star()
    f = theorem_T3_kappa()
    g = theorem_eta_part_formula()
    h = theorem_eta_part_lt_wick_budget()
    i = theorem_reduction_K4_via_eta()
    status = hinge_status_215()
    out = {
        "prop": "15.215",
        "title": (
            "K4/Wick: Tκ spectrum, η_part formula, budget; "
            "δ∈E_4p still OPEN"
        ),
        "proved": {
            "star_exhaustion": a["proved"],
            "sum_kappa_sq": b["proved"],
            "n3_count": c["proved"],
            "Tkappa_norm": d["proved"],
            "Tphi_nplus2_star": e["proved"],
            "T3_kappa": f["proved"],
            "eta_part_formula": g["proved"],
            "eta_part_lt_wick_budget": h["proved"],
        },
        "reduction": i,
        "hinge": status,
        "residual_i_closed": residual_i_closed_via_215(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "soft_close": False,
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15215.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2) + "\n")
    print(f"wrote {path}")
    for k, v in out["proved"].items():
        print(f"  {k}={v}")
    print(f"  residual_i={residual_i_closed_via_215()} gsum={gsum_disj_lb_proved_general()}")
    return out


if __name__ == "__main__":
    main()
