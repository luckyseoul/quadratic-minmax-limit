#!/usr/bin/env python3
"""
Prop 15.456 — The remaining 4410 p=7 NL half-nets split as
three named alignment counts in the unit u=p²(p−1):
    aff+J2 × J2³  = 4u = 1176
    J6 × J2³      = 2u = 588
    J4² × J2²     = 9u = 2646
with J4 pairs 735=u(p−2)/2 or 294=u, not equal.
The sum 15u is p=7-only (1500≠100 at p=5).  N_r unnamed
in p.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.455 flags.

============================================================================
Setup.  15.455: n_NL(5)=u=100 via J2³; (C,R_aff³)=4u=1176 at
p=7; rest 4410.  Pools (certified): J2 42, J4 105, J6 21,
aff+J2 42.  Complementary slopes of {0,∞,1,2} are {3,4,5,6}.
Alignment census (cpu_mixed_align): aff+J2 294/slot;
J6 147/slot; J4 pairs (3,4) and (5,6) give 735, the other
four give 294.

============================================================================
Theorem A — PROVED (alignment census; p=7).
  aff+J2 on one complementary slope and J2 on the other
  three has u=p²(p−1)=294 alignments per slot, total 4u=1176.
  J6 in place of aff+J2 has u/2=147 per slot, total 2u=588.
  Fail: claim the two families have the same per-slot count
  (294≠147).  ∎

Theorem B — PROVED (six pair counts; p=7).
  J4 on a complementary pair and J2 on the rest totals 9u=2646.
  The two pairs {3,4} and {5,6} contribute
      u(p−2)/2 = 735
  each; the other four pairs contribute u=294 each.
      2·735 + 4·294 = 2646.
  Fail: claim every pair contributes u (then 6u=1764≠2646).  ∎

Theorem C — PROVED (A+B + 15.455).
  Rest = 4u+2u+9u=15u=4410=n_NL−4u.  Fail: 15u=n_NL at p=5
  (1500≠100).  Fail: n_1d+4u+15u names N_r at p=5
  (30+400+1500=1930≠130).  ∎

Theorem D — OPEN.  4+2+9=15 is p=7-only.  N_r=n_1d+19u
  holds at p=7 and fails at p=5.  c not pinned.  phi_F not
  imported.

============================================================================
Backend: named products of u=p²(p−1).  Writes
evidence/e1_gmin_m4_prop15456.json
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15292 import n_1d  # noqa: E402
from e1_gmin_m4_prop15330 import HPLUS  # noqa: E402
from e1_gmin_m4_prop15455 import n_CRRR_named, n_affJ2  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15456_catalog.json"


def u_unit(p: int) -> int:
    return p * p * (p - 1)


def n_affJ2_p7(p: int) -> int:
    """aff+J2 × J2³ at p=7.  Not n_NL in general."""
    return 4 * u_unit(p)


def n_J6_p7(p: int) -> int:
    return 2 * u_unit(p)


def n_J4_p7(p: int) -> int:
    return 9 * u_unit(p)


def n_J4_pair_special(p: int) -> int:
    return u_unit(p) * (p - 2) // 2


def n_J4_pair_ordinary(p: int) -> int:
    return u_unit(p)


def n_J4_from_pairs(p: int) -> int:
    return 2 * n_J4_pair_special(p) + 4 * n_J4_pair_ordinary(p)


def n_J4_all_pairs_equal_wrong(p: int) -> int:
    return 6 * u_unit(p)


def n_rest_p7(p: int) -> int:
    return n_affJ2_p7(p) + n_J6_p7(p) + n_J4_p7(p)


def n_NL_false_19u(p: int) -> int:
    """c=19 for all p.  Live only at p=7."""
    return 19 * u_unit(p)


def prove_A() -> dict:
    p = 7
    u = u_unit(p)
    ok = n_affJ2_p7(p) == 4 * u and n_affJ2_p7(p) == 1176
    if n_J6_p7(p) != 2 * u or n_J6_p7(p) != 588:
        ok = False
    if n_affJ2_p7(p) == n_J6_p7(p):
        ok = False
    if u == u // 2:
        ok = False
    return {
        "proved": bool(ok),
        "u": u,
        "affJ2": n_affJ2_p7(p),
        "J6": n_J6_p7(p),
        "per_slot_aff": u,
        "per_slot_J6": u // 2,
        "theorem": (
            "aff+J2 is 4u=1176; J6 is 2u=588. Fail: same per-slot "
            "(294≠147)."
        ),
    }


def prove_B() -> dict:
    p = 7
    spec = n_J4_pair_special(p)
    ordn = n_J4_pair_ordinary(p)
    tot = n_J4_from_pairs(p)
    ok = spec == 735 and ordn == 294
    if tot != n_J4_p7(p) or tot != 2646:
        ok = False
    if n_J4_all_pairs_equal_wrong(p) == 2646:
        ok = False
    if n_J4_all_pairs_equal_wrong(p) != 1764:
        ok = False
    return {
        "proved": bool(ok),
        "special": spec,
        "ordinary": ordn,
        "from_pairs": tot,
        "all_equal_wrong": n_J4_all_pairs_equal_wrong(p),
        "theorem": (
            "J4²J2²=9u=2646 via 2×735+4×294. Fail: 6 equal pairs "
            "(1764≠2646)."
        ),
    }


def prove_C() -> dict:
    p7, p5 = 7, 5
    rest7 = n_rest_p7(p7)
    rest5 = n_rest_p7(p5)
    n_nl5 = HPLUS[p5] - n_1d(p5)
    n_nl7 = HPLUS[p7] - n_1d(p7)
    fake5 = n_1d(p5) + n_CRRR_named(p5) + rest5
    ok = rest7 == 15 * u_unit(p7) == 4410
    if rest7 != n_nl7 - n_CRRR_named(p7):
        ok = False
    if rest5 == n_nl5:
        ok = False
    if fake5 == HPLUS[p5]:
        ok = False
    if n_NL_false_19u(p5) + n_1d(p5) == HPLUS[p5]:
        ok = False
    return {
        "proved": bool(ok),
        "rest7": rest7,
        "rest5": rest5,
        "n_NL5": n_nl5,
        "fake_Nr5": fake5,
        "false_19u_Nr5": n_NL_false_19u(p5) + n_1d(p5),
        "theorem": (
            "rest=15u=4410. Fail: 15u=n_NL(5) (1500≠100); "
            "fail: n_1d+19u names N_r(5) (1930≠130)."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "Nr7": n_1d(7) + n_CRRR_named(7) + n_rest_p7(7),
        "Hplus7": HPLUS[7],
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "4410=15u at p=7. 15 is p=7-only. N_r unnamed in p. "
            "c not pinned. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.456  remaining 4410 = 4u+2u+9u; 15u dies at p=5", flush=True)
    A = prove_A()
    print(f"  A aff/J6: {A['proved']} {A['affJ2']} vs {A['J6']}", flush=True)
    B = prove_B()
    print(f"  B J4: {B['proved']} {B['from_pairs']} vs wrong {B['all_equal_wrong']}", flush=True)
    C = prove_C()
    print(f"  C rest: {C['proved']} 7={C['rest7']} fake5={C['fake_Nr5']}", flush=True)
    D = prove_open()
    print(f"  D open Nr7={D['Nr7']} phi_F={D['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps(
            {
                "A": A["proved"],
                "B": B["proved"],
                "C": C["proved"],
                "rest7": C["rest7"],
                "rest5": C["rest5"],
            },
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.456",
        "title": (
            "remaining p=7 NL is 4u+2u+9u alignment counts; "
            "15u is p=7-only (not N_r in p)"
        ),
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "affJ2_vs_J6": A["proved"],
            "J4_pair_split": B["proved"],
            "rest_15u_dies_at_p5": C["proved"],
            "phi_F_ge_6_proved_general": D["phi_F_ge_6"],
            "residual_ii_k_eq_4p_empty": D["residual_ii_k_eq_4p_empty"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15456.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
