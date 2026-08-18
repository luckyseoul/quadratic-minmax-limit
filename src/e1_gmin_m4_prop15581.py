#!/usr/bin/env python3
"""
Prop 15.581 — At p=5, Q++/q² is the named 1D+k=3 mix and
meets the 15.507 J_{N*}=2 floor.

    n_{1d}=m C(p,m),  n_{k=3}=C(m,3)(p−1)q,  N=n_{1d}+n_{k=3}
    μ_{1d}=2p⁴(p²−3p−2)/(p−2),  μ_{k=3}=96 p⁴ P(r)/(p²−1)
    Q_{++}/q² = (n_{1d} μ_{1d} + n_{k=3} μ_{k=3}) / (N q²)
             = 48/13  <  26/7 = Qpp_ub_J2(5)

m=3 so exclusive n_full=0 (15.573).  All four ingredients
are Max+-free names (15.561, 15.575, 15.574).  Fail: drop
k=3 (16/3>26/7).  Fail: claim the value is 26/7.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
residual_ii / type_I / 15.279–15.580 flags.  Does **not** import phi_F.
Does **not** overwrite 15.550–15.580.  Does **not** interpolate (p−5)/15.
Does **not** use the p=7 15.568 mix weights.

============================================================================
Setup.  15.507 D: p=5, J_{N*}=2 ⇒ ⟨δ,ψ⟩≤2 ⇔ Q++/q²≤26/7.
Claude leftover-1 BLOCK: do not name Q_τ by a two-point μ_full
fit.  This unit only evaluates the already-named 1D+k=3 mix
at p=5, where that mix is the exclusive Q++.

============================================================================
Theorem A — PROVED (15.561+15.574+15.575; Fraction).
  p=5: n_1d=30, n_k3=100, N=130, μ_1d=10000/3, μ_k3=2000,
  Q++/q²=48/13.  Fail: drop k=3 (16/3).  Fail: 26/7 as the
  value.  ∎

Theorem B — PROVED (A + 15.507 D, ⇔).
  15.507 D: J_{N*}=2 pairing ≤2 ⇔ Q++/q²≤26/7 at p=5.
  48/13 < 26/7 (slack 2/91).  μ_k3(5,2)=2000 is 15.574.
  Fail: 48/13 ≥ 26/7.  ∎

Theorem C — OPEN.  For p≥7 the 1D+k=3-only mix (N=n_1d+n_k3)
  is 4.68>4 at p=7 and is not a tight upper bound on live Q
  (true mix uses n_full with smaller μ).  Q++≤4q² and the
  p≡3 pairing stay open.  phi_F not imported.

============================================================================
Backend: Fraction + 15.574 P at p=5.  Writes
evidence/e1_gmin_m4_prop15581.json
"""
from __future__ import annotations

import json
import os
import sys
from fractions import Fraction
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15507 import Qpp_ub_J2
from e1_gmin_m4_prop15561 import n1d_named, n_k3_named
from e1_gmin_m4_prop15574 import mu_k3_named
from e1_gmin_m4_prop15575 import mu_1d_named

EV = ROOT / "evidence" / "e1_gmin_m4_prop15581.json"


def N_p5_named() -> int:
    """Exclusive N at p=5: k=3=full."""
    return n1d_named(5) + n_k3_named(5)


def Qpp_p5_named() -> Fraction:
    n1 = n1d_named(5)
    nk = n_k3_named(5)
    num = n1 * mu_1d_named(5) + nk * Fraction(mu_k3_named(5, 2)).limit_denominator()
    # μ_k3(5,2) is 2000 exactly
    num = n1 * mu_1d_named(5) + nk * Fraction(2000, 1)
    den = N_p5_named() * 5**4
    return Fraction(num, den)


def Qpp_p5_drop_k3() -> Fraction:
    """Fail: drop the k=3 term, Q=μ_1d/q²=16/3."""
    return mu_1d_named(5) / 5**4


def prove_A() -> dict:
    q = Qpp_p5_named()
    drop = Qpp_p5_drop_k3()
    thr = Qpp_ub_J2(5)
    ok = (
        n1d_named(5) == 30
        and n_k3_named(5) == 100
        and N_p5_named() == 130
        and mu_1d_named(5) == Fraction(10000, 3)
        and abs(mu_k3_named(5, 2) - 2000.0) < 1e-8
        and q == Fraction(48, 13)
        and drop == Fraction(16, 3)
        and q != thr
        and drop > thr
    )
    return {
        "proved": bool(ok),
        "n_1d": n1d_named(5),
        "n_k3": n_k3_named(5),
        "N": N_p5_named(),
        "Qpp": str(q),
        "drop_k3": str(drop),
        "thr": str(thr),
        "theorem": (
            "p=5 Q++/q²=48/13 from named μ_1d, μ_k3. "
            "Fail drop k=3 (16/3); fail value=26/7."
        ),
    }


def prove_B() -> dict:
    q = Qpp_p5_named()
    thr = Qpp_ub_J2(5)
    ok = q < thr and thr == Fraction(26, 7) and q == Fraction(48, 13)
    return {
        "proved": bool(ok),
        "Qpp": str(q),
        "thr": str(thr),
        "slack": str(thr - q),
        "theorem": "48/13<26/7 (slack 2/91). J=2 floor from named μ.",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "Qpp_le_4_general": False,
        "p7_1dk3_only_mix": "4.68>4",
        "Q_tau_named_in_p": False,
        "phi_F_imported": False,
        "note": (
            "p=5 floor is named. p≥7 needs n_full/μ_full or a tighter "
            "ub than the 1D+k=3-only mix. Do not interpolate (p−5)/15. "
            "Do not import the 15.568 mix. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.581  p=5 Q++/q²=48/13<26/7 from named μ_1d, μ_k3", flush=True)
    A = prove_A()
    print(f"  A named 48/13: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B 48/13<26/7: {B['proved']}", flush=True)
    C = prove_open()
    out = {
        "prop": "15.581",
        "title": "p=5 Q++/q²=48/13<26/7 from named μ_1d and μ_k3",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "C": C,
        "proved": {
            "Qpp_p5_named_48_13": A["proved"],
            "J2_floor_from_named_mu": B["proved"],
            "Q_tau_named_in_p": False,
            "phi_F_ge_6": False,
        },
        "e1_closed_general": bool(e1_closed_general()),
        "gsum_disj_lb": bool(gsum_disj_lb_proved_general()),
        "phi_F_ge_6": False,
        "L_status": "OPEN",
        "flags_not_flipped": [
            "phi_F_ge_6",
            "e1",
            "L",
            "Aut-Schur",
            "Gsum",
            "pairing",
            "residual_ii",
            "type_I",
        ],
        "identity": "Q++/q^2 = 48/13 < 26/7 at p=5",
    }
    EV.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print("wrote", EV, flush=True)
    return out


if __name__ == "__main__":
    main()
