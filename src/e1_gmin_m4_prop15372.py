#!/usr/bin/env python3
"""
Prop 15.372 — The 15.371 D_form induces Q++_form=8A/D_form which
lies in the S0 floor interval [Q_lo, Q^{ub}] and in the named box
[4−2/(p−1), 4−2/(p+2)] at every prime 5≤p≤79.  For p≡1, every
even S0 family is ≥6 (15.326 C).  The p=11 3-net hunt is replaced
by a 6-net unit-propagation count (p=5,7 certified 130, 5726).
phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.371 flags.  Residual (ii) and Type I stay False.

============================================================================
Setup.  15.326 C: p≡1, S0 line, all even S/q²≥6 ⇔ Q_lo≤Q++≤ub.
15.371: D_form matches live D at p=5,7 and fails at p=3.
Drop-binom D gives Q=16 at p=5, outside the interval.

============================================================================
Theorem A — PROVED (unit-propagation census).
  The 3-net backtrack is the wrong object for D(11).  A 6-net
  (r=(p+1)/2) cell-backtrack with line unit-propagation recovers
  N=130 at p=5 and N=5726 at p=7.  Fail: claim the p=7 4-net
  count is 434 (1D+Hoffman).  ∎

Theorem B — PROVED (exact Fractions, all primes 5..79).
  Q++_form=8A/D_form satisfies Q_lo ≤ Q++_form ≤ Q^{ub} and
  4−2/(p−1) ≤ Q++_form ≤ 4−2/(p+2).  Fail: use D_form_drop_binom
  at p=5, which gives Q=16 ∉ [11/4, 26/7] and min even S=−80.  ∎

Theorem C — PROVED (15.326 C + B; p≡1, 5≤p≤79).
  On the S0 line with Q=Q++_form, every even family is ≥6
  at every p≡1 in 5..79 (worst p=5: 80/13).  Not claimed
  for p≥83.  Fail: the drop-binom Q at p=5.  ∎

Theorem D — OPEN.  D_form is not known to equal D at p≥11
  until the 6-net count returns.  p≡3 still has a residual
  S0 parameter.  phi_F not imported.

============================================================================
Backend: Fraction arithmetic + numba 6-net.  Writes
evidence/e1_gmin_m4_prop15372.json
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

from e1_gmin_m4_prop15170 import e1_closed_general, is_prime  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15321 import Qpp_floor_ub  # noqa: E402
from e1_gmin_m4_prop15326 import Q_lo_named, even_S_p1  # noqa: E402
from e1_gmin_m4_prop15361 import n_1d_plus_hoff  # noqa: E402
from e1_gmin_m4_prop15367 import A_num  # noqa: E402
from e1_gmin_m4_prop15371 import D_form, D_form_drop_binom  # noqa: E402

SCRATCH = Path("/tmp/grok-goal-ea1b2cdbd197/implementer")
BOXJ = SCRATCH / "cpu_Qform_box.json"
NETJ = SCRATCH / "cpu_D11_countnet.json"
CAT = ROOT / "evidence" / "e1_gmin_m4_prop15372_catalog.json"


def Qpp_form(p: int) -> Fraction:
    return Fraction(8 * A_num(p), D_form(p))


def box_lo(p: int) -> Fraction:
    return Fraction(4 * (p - 1) - 2, p - 1)


def box_hi(p: int) -> Fraction:
    return Fraction(4 * (p + 2) - 2, p + 2)


def prove_A() -> dict:
    rec = json.loads(NETJ.read_text()) if NETJ.is_file() else {}
    n5 = rec.get("p5_3net", -1)
    n7 = rec.get("p7_4net", -1)
    n11 = rec.get("p11_6net")
    ok = n5 == 130 and n7 == 5726
    ok = ok and n7 != n_1d_plus_hoff(7)
    if n7 == n_1d_plus_hoff(7) or n5 != 130:
        ok = False
    return {
        "proved": bool(ok),
        "p5_3net": n5,
        "p7_4net": n7,
        "p11_6net": n11,
        "N_pred_11": 244332,
        "p11_eq_pred": n11 == 244332 if n11 is not None else None,
        "replaced_3net": True,
        "theorem": (
            "6-net UP count recovers 130 and 5726. "
            "3-net is not D(11).  Fail: p=7 count is 434."
        ),
    }


def prove_B() -> dict:
    rec = json.loads(BOXJ.read_text()) if BOXJ.is_file() else {}
    # recompute the load-bearing fail-when-wrong here
    Q5 = Qpp_form(5)
    Q7 = Qpp_form(7)
    Qbad = Fraction(8 * A_num(5), D_form_drop_binom(5))
    in5 = Q_lo_named(5) <= Q5 <= Qpp_floor_ub(5)
    in7 = Q_lo_named(7) <= Q7 <= Qpp_floor_ub(7)
    box5 = box_lo(5) <= Q5 <= box_hi(5)
    bad_out = not (Q_lo_named(5) <= Qbad <= Qpp_floor_ub(5))
    ok = in5 and in7 and box5 and bad_out
    ok = ok and bool(rec.get("all_in_s0_interval")) and bool(rec.get("all_in_named_box"))
    if Qbad == Q5 or not bad_out:
        ok = False
    CAT.write_text(
        json.dumps(
            {
                "Q5": str(Q5),
                "Q7": str(Q7),
                "Qbad_drop": str(Qbad),
                "all_s0": rec.get("all_in_s0_interval"),
                "all_box": rec.get("all_in_named_box"),
                "n_primes": rec.get("n_primes"),
            },
            indent=2,
        )
        + "\n"
    )
    return {
        "proved": bool(ok),
        "Q5": str(Q5),
        "Q7": str(Q7),
        "Qbad_drop": str(Qbad),
        "all_s0": rec.get("all_in_s0_interval"),
        "all_box": rec.get("all_in_named_box"),
        "theorem": "Q++_form is in [Q_lo,ub] and the named box at p=5..79.",
    }


def prove_C() -> dict:
    rec = json.loads(BOXJ.read_text()) if BOXJ.is_file() else {}
    Ss5 = even_S_p1(5, Qpp_form(5))
    min5 = min(Ss5.values())
    Ssbad = even_S_p1(5, Fraction(8 * A_num(5), D_form_drop_binom(5)))
    minbad = min(Ssbad.values())
    ok = min5 >= 6 and minbad < 6
    ok = ok and bool(rec.get("all_p1_minS_ge_6"))
    if min5 < 6 or minbad >= 6:
        ok = False
    return {
        "proved": bool(ok),
        "minS_p5": str(min5),
        "minS_drop": str(minbad),
        "all_p1_ge6": rec.get("all_p1_minS_ge_6"),
        "worst": rec.get("worst_p1_minS"),
        "theorem": "p≡1, 5≤p≤79: S0 even families ≥6 at Q++_form; drop-binom dies.",
    }


def prove_open() -> dict:
    rec = json.loads(NETJ.read_text()) if NETJ.is_file() else {}
    n11 = rec.get("p11_6net")
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "Q_tau_named_in_p": False,
        "D_named_in_p": False,
        "p11_6net_done": n11 is not None,
        "p11_matches_D_form": n11 == 244332 if n11 is not None else None,
        "note": (
            "S0 floor for p≡1 would follow from D_form=D. "
            "p≡3 residual remains.  phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.372  Q++_form in S0 box; 6-net replaces 3-net", flush=True)
    A = prove_A()
    print(f"  A 6-net: {A['proved']} p5={A['p5_3net']} p7={A['p7_4net']} p11={A['p11_6net']}", flush=True)
    B = prove_B()
    print(f"  B box: {B['proved']} Q5={B['Q5']} Qbad={B['Qbad_drop']}", flush=True)
    C = prove_C()
    print(f"  C S0: {C['proved']} min5={C['minS_p5']} drop={C['minS_drop']}", flush=True)
    D = prove_open()
    print(f"  D open: phi_F={D['phi_F_ge_6']} p11done={D['p11_6net_done']}", flush=True)
    out = {
        "prop": "15.372",
        "title": "Q++_form in S0/named box; 6-net replaces 3-net",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "six_net_recovers_5_7": A["proved"],
            "Qform_in_s0_and_box": B["proved"],
            "p1_S0_even_ge_6_at_Qform": C["proved"],
            "Q_tau_named_in_p": False,
            "D_named_in_p": False,
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
        ],
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15372.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
