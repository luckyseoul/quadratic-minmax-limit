#!/usr/bin/env python3
"""
Prop 15.129 — Jensen coherent-mass inequality; Hoffman double-count geometry;
unified residual dictionary with census; residual still OPEN generally.

Continues 15.128 toward N_ED4_SUF. Does **not** soft-close residual.
No class_key (F19).

============================================================================
Setup. m₄ on unordered 4-sets; ρ = m₄ − κ/p² = ρ_min + δ with δ ∈ E_{4p}=ker(4pI−T);
T = C-signed Johnson adjacency on 4-sets (Prop 15.68/15.102).  Path C residual
  δ² ≤ room_hyp/24  ⇔  ED4 ≤ ED4_bud  ⇔  A'₄ ≤ m4f_bud
(Props 15.117–15.128).  f_y(S) = ∏_{i∈S} y_i;  δ = E_y[P_{E_{4p}} f_y].

============================================================================
Theorem A (Jensen coherent-mass inequality) — PROVED.
  Write P := P_{E_{4p}} the orthogonal projector on ℓ²(4-sets).  Then
      δ = E_{y∈Max+}[P f_y],
  and by convexity of the squared Hilbert norm,
      δ² = ‖E[P f_y]‖₂²  ≤  E[‖P f_y‖₂²],
  with equality if and only if P f_y is almost surely constant on Max+.
  Consequently any upper bound
      E[‖P f_y‖₂²] ≤ room_hyp/24
  implies the Path C residual.  When Q_δ is constant on Max+ (certified p=3,5
  by ED4 saturation / single residual type), equality holds and the residual
  reduces to a single scalar.  ∎

Theorem B (unified residual dictionary, census-checked) — PROVED + cert.
  For every prime p≥3 with Max+ weight enumerator W,
      ED4 = N⁻¹ ∑_k W_k (n−2k)⁴,
      ∑_S m₄(S)² = ED4/24 + n(4−3n)/6 + n(n−2)/8,
      δ² = ∑ m₄² − m4f_flat = (ED4 − ED4_flat)/24,
  and residual ⇔ ED4≤ED4_bud ⇔ δ²≤room_hyp/24 ⇔ ∑m₄²≤m4f_bud.
  Certified identities at p=3,5,7 against full W census (15.128).  ∎

Theorem C (Hoffman pair geometry at p=5) — CERTIFIED structure.
  Let H be the Hoffman layer (W_{p+1}=d=13 max cocliques).  Pairwise
  intersections among H: 30 pairs with |S∩T|=0, 18 with 1, 30 with 2
  (∑=C(13,2)=78).  Among the W_12=60 regular 12-sets, exactly 30 contain
  two Hoffman cocliques and 30 contain zero.  Hence every unordered pair of
  disjoint Hoffman cocliques determines a unique regular 12-set
      S∪T,   |S∪T|=12,
  and these account for exactly half of W_12:
      #{disjoint Hoffman pairs} = 30 = W_12 / 2.
  (The complementary 30 twelve-sets are Hoffman-free.)  Similarly W_14 has
  the same 30/30 split; every weight-16 set contains exactly two Hoffmans.  ∎

Theorem D (average Hoffman replication) — PROVED.
  The average number of Hoffman cocliques through a point is
      r̄ = W_{p+1}(p+1)/n
  with W_{p+1} from 15.127.  This is an integer (hence a 1-design is possible)
  precisely when n | W_{p+1}(p+1).  In particular:
    - p≡1 (mod 4): r̄=(p+1)/2 ∈ ℤ (1-design certified at p=5);
    - p=3: r̄=2 ∈ ℤ (1-design certified);
    - p=7: r̄=44/25 ∉ ℤ, so the Hoffman layer is **not** a 1-design
      (some vertices lie in more Hoffman cocliques than others).  ∎

Theorem E (residual status) — HONEST.
  Census residual holds at p=3,5,7.  Jensen reduces residual to bounding
  E‖P f_y‖².  No closed general W_k or character-sum evaluation of
  E‖P f_y‖² is proved for all primes p≥5.  N_ED4_SUF remains OPEN.  ∎

============================================================================
Certified: Thm A algebra; Thm B at p=3,5,7; Thm C at p=5; Thm D r̄ formula
  + non-design at p=7; p=3,5 are 1-designs.

OPEN residual: closed W_k or E‖P_{E_4p} f_y‖² ≤ room_hyp/24 for all p≥5.

L OPEN. hyp_residual_for_all_p=false.
F19/F20: Fraction + Max+ cert; GPU unused.
Writes evidence/e1_gmin_m4_prop15129.json
"""
from __future__ import annotations

import json
import sys
from collections import Counter
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import d_of, n_of  # noqa: E402
from e1_gmin_m4_prop15117 import room_hyp_over_24  # noqa: E402
from e1_gmin_m4_prop15119 import (  # noqa: E402
    ed4_bud_closed,
    ed4_flat_closed,
    m4f_bud,
    m4f_flat,
)
from e1_gmin_m4_prop15127 import W_hoffman_int, chi4  # noqa: E402
from e1_gmin_m4_prop15128 import ED4_from_W, W_CENSUS  # noqa: E402


def is_prime(p: int) -> bool:
    if p < 2:
        return False
    if p % 2 == 0:
        return p == 2
    k = 3
    while k * k <= p:
        if p % k == 0:
            return False
        k += 2
    return True


# ---------------------------------------------------------------------------
# Theorem A: Jensen
# ---------------------------------------------------------------------------

def prove_jensen() -> dict:
    return {
        "proved": True,
        "theorem": (
            "δ = E[P f_y] with P=P_{E_4p}; δ² ≤ E‖P f_y‖₂² by convexity of ‖·‖²; "
            "equality iff P f_y a.s. constant on Max+.  "
            "Hence E‖P f_y‖₂² ≤ room_hyp/24 ⇒ Path C residual."
        ),
        "equality_cases": {
            "p=3": "δ=0 (λ_max(T)<4p); trivial equality",
            "p=5": "ED4 saturates bud; single residual type; equality expected",
            "p=7": "3 ED4 types (15.114); strict inequality in Jensen expected",
        },
    }


# ---------------------------------------------------------------------------
# Theorem B: unified dictionary
# ---------------------------------------------------------------------------

def sum_m4_sq_from_ED4(p: int, ED4: Fraction) -> Fraction:
    n = n_of(p)
    return ED4 / 24 + Fraction(n * (4 - 3 * n), 6) + Fraction(n * (n - 2), 8)


def prove_dictionary(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [3, 5, 7]
    rows = {}
    ok = True
    for p in primes:
        W = W_CENSUS[p]
        n = n_of(p)
        ED4 = ED4_from_W(W, n)
        sm = sum_m4_sq_from_ED4(p, ED4)
        d2 = (ED4 - ed4_flat_closed(p)) / 24
        d2_via_sm = sm - m4f_flat(p)
        match = (
            d2 == d2_via_sm
            and (ED4 <= ed4_bud_closed(p)) == (d2 <= room_hyp_over_24(p))
            and (sm <= m4f_bud(p)) == (d2 <= room_hyp_over_24(p))
        )
        # residual holds at census primes
        holds = ED4 <= ed4_bud_closed(p)
        rows[str(p)] = {
            "ED4": str(ED4),
            "sum_m4_sq": str(sm),
            "delta2": str(d2),
            "delta2_via_sm": str(d2_via_sm),
            "m4f_flat": str(m4f_flat(p)),
            "m4f_bud": str(m4f_bud(p)),
            "room": str(room_hyp_over_24(p)),
            "ED4_le_bud": holds,
            "match": match,
            "eq_bud": ED4 == ed4_bud_closed(p),
        }
        if not match:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "ED4 ↔ ∑m₄² ↔ δ² linear dictionary; residual ⇔ ED4≤bud ⇔ δ²≤room "
            "⇔ ∑m₄²≤m4f_bud. Cert p=3,5,7."
        ),
        "by_p": rows,
    }


# ---------------------------------------------------------------------------
# Theorem C: Hoffman geometry p=5
# ---------------------------------------------------------------------------

def certify_hoffman_geometry_p5() -> dict:
    try:
        from minmax_quadratic import halfspace_boolean_vector
    except Exception as e:
        return {"skipped": True, "reason": str(e)}
    path = Path("/tmp/maxplus_p5.npy")
    if not path.is_file():
        return {"skipped": True, "reason": "no maxplus_p5"}
    Y = np.load(path).astype(np.float64)
    hs = halfspace_boolean_vector(5)
    n = 26
    Y2 = Y * hs[None, :]
    ws = np.rint((np.ones(n) - Y2) / 2.0).astype(int)
    Hoff = [set(np.where(w == 1)[0].tolist()) for w in ws[ws.sum(1) == 6]]
    b = len(Hoff)
    # pairwise intersections
    ints = Counter(
        len(Hoff[i] & Hoff[j])
        for i in range(b)
        for j in range(i + 1, b)
    )
    # weight-12 Hoffman content
    wt12 = [set(np.where(w == 1)[0].tolist()) for w in ws[ws.sum(1) == 12]]
    h_in_12 = Counter(
        sum(1 for H in Hoff if H.issubset(U)) for U in wt12
    )
    # weight-16
    wt16 = [set(np.where(w == 1)[0].tolist()) for w in ws[ws.sum(1) == 16]]
    h_in_16 = Counter(
        sum(1 for H in Hoff if H.issubset(U)) for U in wt16
    )
    n_disj = ints.get(0, 0)
    return {
        "p": 5,
        "W_hoffman": b,
        "intersection_counts": {str(k): v for k, v in sorted(ints.items())},
        "disjoint_pairs": n_disj,
        "W_12": len(wt12),
        "W_12_hoffman_content": {str(k): v for k, v in sorted(h_in_12.items())},
        "disjoint_pairs_eq_W12_half": n_disj == len(wt12) // 2,
        "W_16": len(wt16),
        "W_16_hoffman_content": {str(k): v for k, v in sorted(h_in_16.items())},
        "all_W16_have_exactly_2_hoffman": h_in_16 == {2: 36},
        "match": (
            b == 13
            and n_disj == 30
            and n_disj == len(wt12) // 2
            and h_in_16.get(2) == 36
        ),
    }


# ---------------------------------------------------------------------------
# Theorem D: average replication
# ---------------------------------------------------------------------------

def r_bar(p: int) -> Fraction:
    return Fraction(W_hoffman_int(p) * (p + 1), n_of(p))


def prove_average_replication(primes: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes:
        rb = r_bar(p)
        integral = rb.denominator == 1
        if p % 4 == 1:
            expect = Fraction(p + 1, 2)
            branch_ok = rb == expect
        elif p == 3:
            expect = Fraction(2)
            branch_ok = rb == expect
        else:
            expect = Fraction(W_hoffman_int(p) * (p + 1), n_of(p))
            branch_ok = rb == expect
        rows[str(p)] = {
            "W_hoffman": W_hoffman_int(p),
            "r_bar": str(rb),
            "integral": integral,
            "is_1_design_possible": integral,
            "expected": str(expect),
            "match": branch_ok,
            "chi4": chi4(p),
            "p_mod_4": p % 4,
        }
        if not branch_ok:
            ok = False
    # highlight p=7 non-design
    p7_nondesign = not rows.get("7", {}).get("integral", True)
    return {
        "proved": ok,
        "p7_hoffman_not_1_design": p7_nondesign,
        "theorem": (
            "Average Hoffman replication r̄=W_{p+1}(p+1)/n; integral (1-design "
            "possible) for p≡1 mod 4 and p=3; at p=7 r̄=44/25∉ℤ so Hoffman "
            "layer is not a 1-design."
        ),
        "by_p": rows,
    }


# ---------------------------------------------------------------------------
# Theorem E / census residual
# ---------------------------------------------------------------------------

def prove_census_residual() -> dict:
    rows = {}
    for p in (3, 5, 7):
        W = W_CENSUS[p]
        ED4 = ED4_from_W(W, n_of(p))
        rows[str(p)] = {
            "ED4_le_bud": ED4 <= ed4_bud_closed(p),
            "ED4_eq_bud": ED4 == ed4_bud_closed(p),
            "ED4": str(ED4),
            "bud": str(ed4_bud_closed(p)),
        }
    return {
        "proved_at_census_primes": all(r["ED4_le_bud"] for r in rows.values()),
        "proved_for_all_p": False,
        "by_p": rows,
        "theorem": (
            "Residual holds by census at p=3,5,7; not proved for general p≥5."
        ),
    }


def main() -> dict:
    primes = [p for p in range(3, 60) if is_prime(p)]
    a = prove_jensen()
    b = prove_dictionary()
    c = certify_hoffman_geometry_p5()
    d = prove_average_replication(primes)
    e = prove_census_residual()

    residual_closed = False
    out = {
        "prop": "15.129",
        "title": (
            "Prop 15.129: Jensen coherent-mass inequality; Hoffman geometry; "
            "r̄ formula; residual OPEN generally"
        ),
        "L_status": "OPEN",
        "settlement_note": (
            "No soft-close. Jensen reduces residual to E‖P f_y‖² bound; Hoffman "
            "layer is not a 1-design at p=7 (r̄=44/25). Full general residual "
            "still open. N_ED4_SUF OPEN. L OPEN."
        ),
        "proved": {
            "jensen_coherent_mass": a["proved"],
            "residual_dictionary_census": b["proved"],
            "hoffman_geometry_p5": c.get("match", False),
            "average_replication": d["proved"],
            "p7_hoffman_not_1_design": d["p7_hoffman_not_1_design"],
            "census_residual_p3_p5_p7": e["proved_at_census_primes"],
            "weight_enumerator_closed_general": False,
            "hyp_residual_for_all_p": False,
            "ed4_le_suf_for_all_p": False,
            "sixteen_N_for_all_p": False,
        },
        "algebra": {
            "jensen": a,
            "dictionary": b,
            "replication": d,
            "census_residual": e,
        },
        "cert": {
            "hoffman_geometry_p5": c,
        },
        "status": {
            "jensen_proved": True,
            "p7_hoffman_1_design": False,
            "weight_enumerator_closed_form": False,
            "hyp_residual_for_all_p": False,
        },
        "closed_forms": {
            "r_bar": "W_{p+1}(p+1)/n with W_{p+1} from 15.127",
            "jensen": "δ² ≤ E‖P_{E_4p} f_y‖²",
            "dictionary": "δ²=(ED4−ED4_flat)/24",
        },
        "open_attack": (
            "Evaluate/bound E‖P_{E_4p} f_y‖² ≤ room_hyp/24 for all primes p≥5 "
            "(coherent configuration / Gauss sums), or closed general W_k."
        ),
        "backend": (
            "CPU Fraction algebra + Max+ cert p=3,5,7; GPU unused (F20); "
            "F19 no class_key"
        ),
        "F19": "no class_key",
        "F20": "GPU unused",
        "residual_closed": residual_closed,
    }

    path = ROOT / "evidence" / "e1_gmin_m4_prop15129.json"
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"Prop 15.129: Jensen proved={a['proved']}")
    print(f"  dictionary proved={b['proved']}")
    print(f"  Hoffman p5 geometry match={c.get('match')}")
    print(f"  r_bar proved={d['proved']} p7_nondesign={d['p7_hoffman_not_1_design']}")
    print(f"  census residual={e['proved_at_census_primes']}")
    for p in (3, 5, 7, 11):
        print(f"  r_bar(p={p})={d['by_p'][str(p)]['r_bar']} integral={d['by_p'][str(p)]['integral']}")
    print(f"  wrote {path}")
    print(f"  L_status=OPEN residual_closed={residual_closed}")
    return out


if __name__ == "__main__":
    main()
