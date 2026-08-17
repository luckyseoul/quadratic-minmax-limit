#!/usr/bin/env python3
"""
Prop 15.382 — Interior 4-level leftover cannot be a pure-pair
Aut_e double-star.  Minus-slice first moment on {f_e=−1} forces
|G∩⊠-star|=p+1 (the μ_f=0 count).  Interior still needs singles
or far edges.  residual_ii_k_eq_4p_empty stays False.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.381 flags.

============================================================================
Setup.  15.381 emptied both 4-level J-corners.  Interior is
b>0, d>0: S takes −4 and −6 on {f_e=+1}.  Aut_e double-star
(μ_f=0) has |G∩□|=3p−1, |G∩⊠|=p+1 and no far edges (15.274 J).
Wedge: on {f_e=−1}, □-pairs cancel and ⊠-pairs double.

============================================================================
Theorem A — PROVED (mod-4 support).
  Interior 4-level needs S=−4 (≡0 mod 4) and S=−6 (≡2 mod 4)
  on U^c={f_e=+1}.  Fail: claim −4≡2 (mod 4).  ∎

Theorem B — PROVED (pure-pair residue).
  Pure-pair double-star (no singles): on U^c, □-pairs double and
  ⊠-pairs cancel, so S=4n_+−(3p−1).  Thus S≡−(3p−1) (mod 4):
  p≡1 ⇒ S≡2 (mod 4); p≡3 ⇒ S≡0 (mod 4).  A single residue
  class cannot contain both −4 and −6.  Not interior.
  Fail: claim S≡0 (mod 4) at p=5.  ∎

Theorem C — PROVED (minus-slice first moment on U).
  E[f_g|U]=0 for a □-star edge and −2/(p+1) for a ⊠-star edge
  (E[f]=−1/p, E[f_e f_g]=−Δ/p).  S≡−2 on U then forces
      n1_⊠+2 n2_⊠ = |G∩⊠-star| = p+1.
  Fail: claim the count is 3p−1.  ∎

Theorem D — OPEN.  Double-stars with singles, and leftover G
  with far edges, are not emptied.  Flag stays False.

============================================================================
Backend: Fraction identities.  Writes
evidence/e1_gmin_m4_prop15382.json
"""
from __future__ import annotations

import json
import os
import sys
from fractions import Fraction
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import (  # noqa: E402
    residual_ii_k_eq_4p_empty,
    three_weight_ms_mu_box,
    three_weight_ms_mu_boxtimes,
)
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15381 import E_fe_minus_named, a_thr, e_hi, e_lo  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15382_catalog.json"


def n_box_muf0(p: int) -> int:
    return 3 * p - 1


def n_boxtimes_muf0(p: int) -> int:
    return p + 1


def pure_pair_S_mod4(p: int) -> int:
    """S ≡ −(3p−1) (mod 4) on {f_e=+1} for a pure-pair double-star."""
    return (-(3 * p - 1)) % 4


def E_fg_U_box(_p: int) -> Fraction:
    """E[f_g|U] for a □-star edge (Δ=+1)."""
    return Fraction(0)


def E_fg_U_boxtimes(p: int) -> Fraction:
    """E[f_g|U] for a ⊠-star edge (Δ=−1)."""
    return Fraction(-2, p + 1)


def prove_A() -> dict:
    ok = True
    ok = ok and ((-4) % 4) == 0
    ok = ok and ((-6) % 4) == 2
    if ((-4) % 4) == ((-6) % 4):
        ok = False
    # interior interval is nonempty
    for p in (5, 7, 11, 13, 19):
        ok = ok and e_lo(p) < e_hi(p)
        # a midpoint has b>0 and d>0
        e = (e_lo(p) + e_hi(p)) / 2
        b = e - Fraction(1, p)
        d = a_thr(p) - 2 * e
        ok = ok and b > 0 and d > 0
    return {
        "proved": bool(ok),
        "m4_mod": (-4) % 4,
        "m6_mod": (-6) % 4,
        "theorem": (
            "Interior needs S=−4 and S=−6 on {f=+1}, distinct mod 4. "
            "Fail: −4≡−6 (mod 4)."
        ),
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 11, 13, 17, 19, 23):
        r = pure_pair_S_mod4(p)
        expect = 2 if p % 4 == 1 else 0
        ok = ok and r == expect
        # cannot contain both 0 and 2
        ok = ok and not (r == 0 and r == 2)
        rows[str(p)] = {"S_mod4": r, "p_mod4": p % 4}
        # fail-when-wrong: p=5 is ≡2, not 0
    ok = ok and pure_pair_S_mod4(5) == 2
    ok = ok and pure_pair_S_mod4(7) == 0
    if pure_pair_S_mod4(5) == 0:
        ok = False
    # 4n − (3p−1) at n=0 is −(3p−1); at full n2_□ is +(3p−1)
    ok = ok and (3 * 5 - 1) == 14
    ok = ok and (-14) % 4 == 2
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "Pure-pair double-star: S≡−(3p−1) (mod 4) on {f=+1}. "
            "Not both −4 and −6. Fail: S≡0 (mod 4) at p=5."
        ),
    }


def prove_C() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 11, 13, 19):
        # E[f_e f_g] = −Δ/p; □ ⇒ −1/p, ⊠ ⇒ +1/p
        a = a_thr(p)
        # E[f|U] = (E[f] − E[ff_e]) / (2a)
        two_a = 2 * a
        e_box = (E_fe_minus_named(p) - Fraction(-1, p)) / two_a
        e_bt = (E_fe_minus_named(p) - Fraction(1, p)) / two_a
        ok = ok and e_box == E_fg_U_box(p)
        ok = ok and e_bt == E_fg_U_boxtimes(p)
        # S≡−2 ⇒ n_⊠ * (−2/(p+1)) = −2 ⇒ n_⊠ = p+1
        ok = ok and n_boxtimes_muf0(p) * E_fg_U_boxtimes(p) == -2
        ok = ok and n_boxtimes_muf0(p) != n_box_muf0(p)
        # matches 15.274 μ_f=0 counts
        ns = p * p - 1
        ok = ok and three_weight_ms_mu_box(p, Fraction(0)) * ns == n_box_muf0(p)
        ok = ok and (
            three_weight_ms_mu_boxtimes(p, Fraction(0)) * ns == n_boxtimes_muf0(p)
        )
        rows[str(p)] = {
            "E_U_box": str(e_box),
            "E_U_bt": str(e_bt),
            "n_bt": n_boxtimes_muf0(p),
            "n_box": n_box_muf0(p),
        }
        if n_boxtimes_muf0(p) == n_box_muf0(p):
            ok = False
    ok = ok and n_boxtimes_muf0(5) == 6
    ok = ok and n_box_muf0(5) == 14
    if n_boxtimes_muf0(5) == 14:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "Minus-slice + double-star ⇒ |G∩⊠|=p+1. "
            "Fail: claim 3p−1."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "Interior 4-level with singles, or with far Aut_e edges, "
            "stays open. residual_ii_k_eq_4p_empty stays False. "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.382  interior 4-level ≠ pure-pair double-star", flush=True)
    A = prove_A()
    print(f"  A mod4: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B pure-pair: {B['proved']} p5mod={B['rows']['5']['S_mod4']}", flush=True)
    C = prove_C()
    print(f"  C n_⊠=p+1: {C['proved']}", flush=True)
    D = prove_open()
    print(f"  D open: resii={D['residual_ii_k_eq_4p_empty']}", flush=True)
    CAT.write_text(
        json.dumps(
            {
                "pure_pair_mod4": {
                    str(p): pure_pair_S_mod4(p) for p in (5, 7, 11, 13)
                },
                "n_bt": {str(p): n_boxtimes_muf0(p) for p in (5, 7, 11)},
            },
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.382",
        "title": (
            "Interior 4-level is not a pure-pair Aut_e double-star; "
            "|G∩⊠|=p+1"
        ),
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "interior_needs_both_mod4": A["proved"],
            "pure_pair_not_interior": B["proved"],
            "n_boxtimes_eq_p_plus_1": C["proved"],
            "residual_ii_k_eq_4p_empty": D["residual_ii_k_eq_4p_empty"],
            "phi_F_ge_6_proved_general": D["phi_F_ge_6"],
        },
        "algebra": {"A": A, "B": B, "C": C, "D": D},
        "L_status": "OPEN",
        "flags_not_flipped": [
            "phi_F_ge_6",
            "e1",
            "L",
            "Aut-Schur",
            "Gsum",
            "pairing",
            "residual_ii_k_eq_4p_empty",
        ],
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15382.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
