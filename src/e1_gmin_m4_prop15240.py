#!/usr/bin/env python3
"""
Prop 15.240 — maj ≤ 2/n; envelope criterion ⇒ |μ|≤2/n≤1/(2p);
residual-(i) still OPEN on the envelope hypothesis.

Does **not** flip residual_i / type_I / gsum / E1 / L. Soft-close forbidden.

SETUP (15.186, 15.191, 15.205, 15.239)
  On |κ|=1 with |φ|≤2(p−2): |μ_part|≤maj=(p²+4p−9)/(p²(p²−5)),
  |f4|≤2/n with f4=(4κ−φ)/(p n), and 2/n≤1/(2p) for primes p≥5 (15.191).
  Census: μ lies in the interval with endpoints among {0, μ_part, f4}
  at p=3,5,7 (“envelope”).  Free-e path still needs ker=sc (15.207–213).

PROVED Max+-free
  A. maj ≤ 2/n for every prime p≥5.
        (p²+4p−9)/(p²(p²−5)) ≤ 2/(p²+1)
     ⇔ (p²+4p−9)(p²+1) ≤ 2 p² (p²−5)
     ⇔ 0 ≤ p⁴ − 4p³ − 2p² − 4p + 9
     ⇔ 0 ≤ (p−1)(p³ − 3p² − 5p − 9).
     p−1>0.  Cubic c(p)=p³−3p²−5p−9: c(5)=16>0 and
     c'(p)=3p²−6p−5, c'(5)=40>0, c''>0 on [5,∞) so c'>0 and c>0.  ∎

  B. Envelope criterion (sufficient for residual-(i)).
     If for every |κ|=1 four-set of Paley C,
        |μ_S| ≤ max(|μ_part_S|, |f4_S|),
     then |μ_S| ≤ 2/n ≤ 1/(2p) for all primes p≥5.
     Proof: |μ_part|≤maj≤2/n by (A) and 15.186; |f4|≤2/n by 15.191 F;
     2/n≤1/(2p) by 15.191 F.  Then Gsum≥−1/p (15.176) ⇒ residual-(i)
     dual-eq empty once wired.  ∎
     (Equivalent form: μ ∈ conv{μ_part, f4, −μ_part, −f4} is stronger than
     needed; the max-abs comparison is the direct L∞ statement.)

  C. Free-e reduction recall (no new algebra).
     free-e_sc ≤ r(p)·scheme_max with r≤3/2 (15.205) and
     (3/2)·scheme < 2−α (15.192) ⇒ free-e_sc < 2−α for all p≥5
     **when** the free-e max equals free-e_sc i.e. ker=sc.
     ker=sc ⇔ G₊≻0 on 𝒲₊₊⁰ ⇔ E[q²]≥λ_*=8(n−6)/n ‖B‖² (15.207–213).
     cost_D≤8/p + ker=sc is an alternate free-e close (15.210).  ∎

CERTIFIED
  D. Envelope |μ|≤max(|μ_part|,|f4|) on all |κ|=1 at p=5 (μ=f4) and by
     prior Max± census at p=3,7 (15.191 I / session handoff).
  E. maj≤2/n as Fraction for all primes 5≤p≤500.

OPEN (blocks residual i / predicate flip)
  F. Prove the envelope |μ|≤max(|μ_part|,|f4|) on |κ|=1 for all primes
     p≥5, **or** E[q²]≥λ_* on 𝒲₊₊⁰ (ker=sc), **or** |μ|≤M_cand / K₄ thr
     / Path C directly.  Soft-close forbidden.

Until F: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15240.json
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
from e1_gmin_m4_prop15186 import majorant_mu_part  # noqa: E402
from e1_gmin_m4_prop15188 import two_over_n  # noqa: E402
from e1_gmin_m4_prop15191 import theorem_f4_majorant  # noqa: E402


def maj_le_two_over_n_diff(p: int) -> Fraction:
    """(p−1)(p³−3p²−5p−9) / positive den — sign of 2/n − maj."""
    # 2/n - maj = N/D with N = (p-1)(p^3-3p^2-5p-9) up to positive factors
    return Fraction(
        (p - 1) * (p**3 - 3 * p * p - 5 * p - 9),
        p * p * (p * p - 5) * (p * p + 1),
    )


def theorem_maj_le_two_over_n() -> dict:
    ok = True
    sample = {}
    for p in range(5, 501):
        if not is_prime(p):
            continue
        maj = majorant_mu_part(p)
        twn = two_over_n(p)
        d = maj_le_two_over_n_diff(p)
        if maj > twn or d != twn - maj or d <= 0:
            ok = False
            break
        if p in (5, 7, 11, 13, 17, 19, 23, 31, 47, 101):
            sample[str(p)] = {
                "maj": str(maj),
                "two_over_n": str(twn),
                "slack": str(d),
            }
    # cubic positivity proof data
    cubic_ok = all(
        (p**3 - 3 * p * p - 5 * p - 9) > 0
        for p in range(5, 501)
        if is_prime(p)
    )
    return {
        "proved": ok and cubic_ok,
        "theorem": (
            "For every prime p≥5: maj=(p²+4p−9)/(p²(p²−5)) ≤ 2/n=2/(p²+1).  "
            "Equivalent to (p−1)(p³−3p²−5p−9)≥0; cubic is 16 at p=5 with "
            "positive derivative 3p²−6p−5 on [5,∞)."
        ),
        "by_p_sample": sample,
        "depends_on": ["prop_15_186_majorant", "prop_15_188_two_over_n"],
    }


def theorem_envelope_criterion() -> dict:
    t_maj = theorem_maj_le_two_over_n()
    t_f4 = theorem_f4_majorant()
    # 2/n <= 1/(2p) is in theorem_f4_majorant / 15.191
    ok_chain = t_maj["proved"] and t_f4.get("proved", True)
    return {
        "proved_criterion": ok_chain,
        "hypothesis_proved_general": False,
        "theorem": (
            "If |μ|≤max(|μ_part|,|f4|) on every |κ|=1 four-set for all primes "
            "p≥5, then |μ|≤2/n≤1/(2p) (by maj≤2/n, |f4|≤2/n, 2/n≤1/(2p)), "
            "hence disj Gsum≥−1/p and residual-(i) Farkas applies.  "
            "The pointwise envelope hypothesis remains OPEN in general "
            "(certified p=3,5,7 only)."
        ),
        "depends_on": [
            "theorem_maj_le_two_over_n",
            "prop_15_191_f4_majorant",
            "prop_15_176_farkas",
        ],
        "sufficient_for_residual_i": False,  # hypothesis open
    }


def theorem_free_e_sc_budget_recall() -> dict:
    return {
        "proved_structure": True,
        "ker_sc_proved_general": False,
        "theorem": (
            "free-e_sc ≤ r(p)·scheme with r≤3/2 and (3/2)·scheme<2−α for all "
            "p≥5 ⇒ free-e_sc<2−α.  Dual-eq empty for all p≥5 once ker=sc "
            "(G₊≻0 / E[q²]≥λ_* on 𝒲₊₊⁰).  ker=sc still OPEN."
        ),
        "depends_on": ["prop_15_192", "prop_15_205", "prop_15_207", "prop_15_211"],
    }


def residual_i_closed_via_240() -> bool:
    return False


def hinge_status_240() -> dict:
    return {
        "maj_le_two_over_n": theorem_maj_le_two_over_n()["proved"],
        "envelope_criterion": theorem_envelope_criterion()["proved_criterion"],
        "envelope_hypothesis": theorem_envelope_criterion()[
            "hypothesis_proved_general"
        ],
        "free_e_sc_structure": theorem_free_e_sc_budget_recall()[
            "proved_structure"
        ],
        "residual_i_closed_via_240": residual_i_closed_via_240(),
        "residual_i_dual_eq_empty": residual_i_dual_eq_empty_proved_general(),
        "gsum_disj_lb": gsum_disj_lb_proved_general(),
        "e1": e1_closed_general(),
        "open_hinge": (
            "envelope |μ|≤max(|μ_part|,|f4|) or ker=sc (λ_*) or |μ|≤M_cand/K4/Path C"
        ),
    }


def main() -> dict:
    a = theorem_maj_le_two_over_n()
    b = theorem_envelope_criterion()
    c = theorem_free_e_sc_budget_recall()
    status = hinge_status_240()
    out = {
        "prop": "15.240",
        "title": (
            "maj≤2/n; envelope criterion ⇒ residual-(i); envelope OPEN"
        ),
        "proved": {
            "maj_le_two_over_n": a["proved"],
            "envelope_criterion_chain": b["proved_criterion"],
            "envelope_hypothesis": b["hypothesis_proved_general"],
            "residual_i_closed": residual_i_closed_via_240(),
        },
        "algebra": {
            "maj_le_2n_slack": "(p-1)(p^3-3p^2-5p-9)/(p^2(p^2-5)(p^2+1))",
            "maj_p5": str(majorant_mu_part(5)),
            "two_over_n_p5": str(two_over_n(5)),
            "thr_p5": str(mu4_threshold(5)),
        },
        "hinge": status,
        "residual_i_closed": residual_i_closed_via_240(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "note": (
            "Envelope is a clean sufficient Max+-free target reduced to "
            "2/n by maj≤2/n.  Hypothesis not proved for general p.  "
            "No predicate flip."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15240.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"Prop 15.240  residual_i={residual_i_closed_via_240()}")
    print(f"  maj_le_2n={a['proved']} envelope_hyp={b['hypothesis_proved_general']}")
    print(f"  wrote {path}")
    return out


if __name__ == "__main__":
    main()
