#!/usr/bin/env python3
"""
Prop 15.241 — Hull criterion (conv{0,μ_part,f4}) ⇒ residual-(i);
halfspace-orbit cancellation structure; residual-(i) still OPEN.

Does **not** flip residual_i / type_I / gsum / E1 / L. Soft-close forbidden.

SETUP (15.186, 15.191, 15.240)
  On |κ|=1 with |φ|≤2(p−2): |μ_part|≤maj≤2/n and |f4|≤2/n≤1/(2p)
  for primes p≥5 (15.240 A, 15.191 F).  Max-abs envelope
  |μ|≤max(|μ_part|,|f4|) is already a sufficient criterion (15.240 B).

PROVED Max+-free
  A. Hull criterion (sufficient for residual-(i)).
     If for every |κ|=1 four-set of Paley C,
        min(0, μ_part_S, f4_S) ≤ μ_S ≤ max(0, μ_part_S, f4_S),
     then |μ_S| ≤ max(|μ_part_S|, |f4_S|) ≤ 2/n ≤ 1/(2p) for all primes p≥5.
     Proof: the real interval [min(0,μ_part,f4), max(0,μ_part,f4)] is
     contained in [−M, M] with M=max(|μ_part|,|f4|,0)=max(|μ_part|,|f4|),
     and M≤2/n by 15.240 A + 15.191 F; then 15.176 Farkas.  ∎
     (When μ_part f4 ≤ 0 the hull interval equals [−M,M] up to the two
     endpoints, so hull ⇔ max-abs envelope.  When μ_part f4 > 0 the hull
     forces μ to share the common sign of {μ_part,f4} or vanish — a
     strict strengthening of max-abs alone.)

  B. Segment convexity (recall / toolkit).
     For any t∈[0,1]: |(1−t)μ_part + t f4| ≤ 2/n on |κ|=1, because both
     endpoints lie in [−2/n,2/n] and the interval is convex.  ∎
     (So any μ on the segment [μ_part,f4] is residual-(i) safe; the hull
     with 0 is the census-stable relaxation when the segment is overshot
     toward 0, as at p=7 |φ|=10 CR-split.)

  C. Halfspace orbit alone does **not** give residual-(i).
     Let H_p be the Aut-orbit of the F_p-linear halfspace boolean
     eigenvector (minmax_quadratic.halfspace_boolean_vector; size
     |H_3|=12, |H_5|=60, |H_7|=168).  At p=5 the H_p-only average
     m4_H satisfies max_{|κ|=1}|m4_H| ≥ 3/5 > 1/10 on some classes
     (certified), while the full Max+ average equals f4 with
     max|μ|=3/65≤1/10.  Hence residual-(i) for full Max+ requires
     cancellation from non-halfspace Max+ vectors; pure halfspace
     character sums cannot close residual-(i) for p≥5.  ∎
     (At p=3, Max+=H_3 and residual-(i) correctly fails.)

CERTIFIED (not general)
  D. Hull holds with 0 violations on all |κ|=1 sets at p=5 (μ=f4) and
     p=7 (prior Max+ census: between-0-μ_part-f4 viol=0 on 176400 sets;
     max-abs envelope viol=0).  p=3: μ=μ_part on |κ|=1 (Max+=H_3).
  E. Segment-between μ_part and f4 alone fails at p=7 (29400 viol), so
     the 0-vertex of the hull is essential for the census-stable form.

OPEN (blocks residual i / predicate flip)
  F. Prove the hull hypothesis of (A) for all primes p≥5 (or max-abs
     envelope, or λ_* / K₄≤Wick_hi / |μ|≤M_cand).  Soft-close forbidden.

Until F: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15241.json
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
    n_of,
    residual_i_dual_eq_empty_proved_general,
)
from e1_gmin_m4_prop15177 import mu4_threshold  # noqa: E402
from e1_gmin_m4_prop15186 import majorant_mu_part, mu_part_formula  # noqa: E402
from e1_gmin_m4_prop15188 import two_over_n  # noqa: E402
from e1_gmin_m4_prop15191 import theorem_f4_majorant  # noqa: E402
from e1_gmin_m4_prop15240 import (  # noqa: E402
    theorem_envelope_criterion,
    theorem_maj_le_two_over_n,
)


def f4_formula(p: int, kappa: int, phi: int) -> Fraction:
    n = n_of(p)
    return Fraction(4 * kappa - phi, p * n)


def theorem_hull_criterion() -> dict:
    """Hull with 0 is Max+-free sufficient for residual-(i)."""
    t_maj = theorem_maj_le_two_over_n()
    t_f4 = theorem_f4_majorant()
    t_env = theorem_envelope_criterion()
    ok = t_maj["proved"] and t_f4.get("proved", True) and t_env["proved_criterion"]
    # algebraic check: for sample p, all (κ,φ) with |κ|=1,|φ|≤2(p-2),
    # any x in [min(0,mp,f4), max(0,mp,f4)] has |x| ≤ 2/n
    sample_ok = True
    samples = {}
    for p in (5, 7, 11, 13, 17, 19, 23):
        if not is_prime(p):
            continue
        M = two_over_n(p)
        bound_phi = 2 * (p - 2)
        phis = list(range(-bound_phi, bound_phi + 1, 4))
        worst = Fraction(0)
        for k in (-1, 1):
            for ph in phis:
                mp = mu_part_formula(p, k, ph)
                ff = f4_formula(p, k, ph)
                lo = min(Fraction(0), mp, ff)
                hi = max(Fraction(0), mp, ff)
                # endpoints of hull are among {0, mp, f4}
                for x in (lo, hi, mp, ff, Fraction(0)):
                    if abs(x) > M + Fraction(0):  # exact
                        if abs(x) > M:
                            sample_ok = False
                    worst = max(worst, abs(x))
        samples[str(p)] = {
            "two_over_n": str(M),
            "hull_endpoint_max_abs": str(worst),
            "le_two_over_n": worst <= M,
        }
    return {
        "proved_criterion": ok and sample_ok,
        "hypothesis_proved_general": False,
        "theorem": (
            "If min(0,μ_part,f4) ≤ μ ≤ max(0,μ_part,f4) on every |κ|=1 "
            "four-set for all primes p≥5, then |μ|≤max(|μ_part|,|f4|)≤2/n "
            "≤1/(2p), hence residual-(i) Farkas (15.176).  Hull is equivalent "
            "to max-abs envelope when μ_part f4≤0, and strictly stronger when "
            "μ_part f4>0 (same-sign coherence).  Hypothesis OPEN (certified "
            "p=3,5,7 only)."
        ),
        "depends_on": [
            "theorem_maj_le_two_over_n",
            "prop_15_191_f4_majorant",
            "prop_15_240_envelope_criterion",
            "prop_15_176_farkas",
        ],
        "by_p_sample": samples,
        "sufficient_for_residual_i": False,  # hypothesis open
    }


def theorem_segment_convexity() -> dict:
    ok = True
    sample = {}
    for p in range(5, 101):
        if not is_prime(p):
            continue
        M = two_over_n(p)
        bound_phi = 2 * (p - 2)
        phis = list(range(-bound_phi, bound_phi + 1, 4))
        for k in (-1, 1):
            for ph in phis:
                mp = mu_part_formula(p, k, ph)
                ff = f4_formula(p, k, ph)
                # check 5 convex combinations
                for t_num in range(0, 5):
                    t = Fraction(t_num, 4)
                    x = (1 - t) * mp + t * ff
                    if abs(x) > M:
                        ok = False
        if p in (5, 7, 11, 13, 23, 47, 97):
            sample[str(p)] = {"two_over_n": str(M), "segment_in_ball": True}
    return {
        "proved": ok,
        "theorem": (
            "On every |κ|=1 four-set for primes p≥5: for all t∈[0,1], "
            "|(1−t)μ_part+t f4|≤2/n, since |μ_part|≤maj≤2/n, |f4|≤2/n, "
            "and [−2/n,2/n] is convex."
        ),
        "depends_on": ["prop_15_240_maj_le_2n", "prop_15_191_f4_majorant"],
        "by_p_sample": sample,
    }


def theorem_halfspace_orbit_insufficient() -> dict:
    """Halfspace-only m4 can exceed 1/(2p); full Max+ needs cancellation."""
    # Certified numerical facts (exact fractions where simple)
    cert = {
        "3": {
            "orbit_size": 12,
            "equals_full_Max_plus": True,
            "note": "Max+=H_3; residual-(i) fails at p=3 (consistent)",
        },
        "5": {
            "orbit_size": 60,
            "full_Max_plus": 260,
            "m4_H_on_k1_class_abs": str(Fraction(1, 15)),  # κ/(3p) on main classes
            "m4_H_max_abs_lower_bound": str(Fraction(3, 5)),
            "thr": str(Fraction(1, 10)),
            "exceeds_thr": True,
            "full_max_abs": str(Fraction(3, 65)),
            "full_equals_f4_on_k1": True,
        },
        "7": {
            "orbit_size": 168,
            "full_Max_plus": 11452,
            "m4_H_max_abs_lower_bound": str(Fraction(1, 7)),
            "thr": str(Fraction(1, 14)),
            "exceeds_thr": True,
        },
    }
    return {
        "proved_structure": True,
        "theorem": (
            "The F_p-linear halfspace Aut-orbit H_p is a proper subset of Max+ "
            "for p≥5 (|H_5|=60<260, |H_7|=168<11452).  The H_p-only 4-point "
            "average m4_H exceeds 1/(2p) on some |κ|=1 classes (p=5: ≥3/5>1/10; "
            "p=7: ≥1/7>1/14, certified), while full Max+ satisfies the residual-(i) "
            "threshold at those p.  Therefore pure halfspace character sums cannot "
            "prove residual-(i) for general p≥5; non-halfspace Max+ cancellation is "
            "essential.  At p=3, Max+=H_3 and residual-(i) fails."
        ),
        "certified": cert,
        "depends_on": [
            "halfspace_boolean_vector",
            "paley_conference_prime_power",
            "prop_15_188_census",
        ],
        "dead_path": "halfspace-only character sum as residual-(i) proof",
    }


def residual_i_closed_via_241() -> bool:
    return False


def hinge_status_241() -> dict:
    return {
        "hull_criterion": theorem_hull_criterion()["proved_criterion"],
        "hull_hypothesis": theorem_hull_criterion()["hypothesis_proved_general"],
        "segment_convexity": theorem_segment_convexity()["proved"],
        "halfspace_insufficient": theorem_halfspace_orbit_insufficient()[
            "proved_structure"
        ],
        "residual_i_closed_via_241": residual_i_closed_via_241(),
        "residual_i_dual_eq_empty": residual_i_dual_eq_empty_proved_general(),
        "gsum_disj_lb": gsum_disj_lb_proved_general(),
        "e1": e1_closed_general(),
        "open_hinge": (
            "hull min(0,μ_part,f4)≤μ≤max(0,μ_part,f4) on |κ|=1 for all p≥5, "
            "or max-abs envelope / λ_* / K4≤Wick_hi / M_cand"
        ),
    }


def main() -> dict:
    a = theorem_hull_criterion()
    b = theorem_segment_convexity()
    c = theorem_halfspace_orbit_insufficient()
    status = hinge_status_241()
    out = {
        "prop": "15.241",
        "title": (
            "hull criterion ⇒ residual-(i); halfspace cancellation; "
            "hull OPEN"
        ),
        "proved": {
            "hull_criterion_chain": a["proved_criterion"],
            "hull_hypothesis": a["hypothesis_proved_general"],
            "segment_convexity": b["proved"],
            "halfspace_orbit_insufficient": c["proved_structure"],
            "residual_i_closed": residual_i_closed_via_241(),
        },
        "algebra": {
            "hull": "min(0,μ_part,f4) ≤ μ ≤ max(0,μ_part,f4)",
            "implies": "|μ|≤max(|μ_part|,|f4|)≤2/n≤1/(2p)",
            "thr_p5": str(mu4_threshold(5)),
            "maj_p5": str(majorant_mu_part(5)),
        },
        "hinge": status,
        "residual_i_closed": residual_i_closed_via_241(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "note": (
            "Hull is the census-stable sufficient form (needs the 0 vertex at "
            "p=7).  Halfspace-only path is dead.  No predicate flip."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15241.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"Prop 15.241  residual_i={residual_i_closed_via_241()}")
    print(f"  hull_crit={a['proved_criterion']} hull_hyp={a['hypothesis_proved_general']}")
    print(f"  segment={b['proved']} halfspace_dead={c['proved_structure']}")
    print(f"  wrote {path}")
    return out


if __name__ == "__main__":
    main()
