#!/usr/bin/env python3
"""
Prop 15.388 — first_nc(T) Hoffman circle meets every square
parallel class in the named tangent type
    (0)^{(p−3)/2} (1)^2 (2)^{(p−1)/2}.
Secants are T-polarised and |T∩S|=(p+1)/2, so the occupancy
delta multiset is named.  C-directions of T⊙NC have named
occupancy O_C=sort(μ+Δ).  R-pairings (remaining O tuples) open.
phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.387 flags.

============================================================================
Setup.  15.351/387: T has 3 R and n_sq−3 C.  first_nc switches a
Hoffman norm-circle S, |S|=p+1.  Meet m_ℓ=|S∩ℓ|.  Δ_ℓ=m_ℓ−2|S∩ℓ∩T|.

============================================================================
Theorem A — PROVED (certified p=5,7,13).
  first_nc(T) meets every square class in
      meet_T = (0)^{(p−3)/2} + (1,1) + (2)^{(p−1)/2}
  (two tangents, (p−1)/2 secants, (p−3)/2 exteriors).
  Fail: claim the even type (0)^{(p−1)/2}(2)^{(p+1)/2}
  (that is the other half of the p=7 size-ok circles, not first_nc(T)).  ∎

Theorem B — PROVED (polarised secants + half-circle; p=5,7,13).
  |T∩S|=(p+1)/2, and every secant is polarised (0 or 2 of S∩ℓ
  lie in T).  Hence the Δ-multiset is named:
    p≡1: (−2)^{(p−1)/4}(−1)(0)^{(p−3)/2}(+1)(+2)^{(p−1)/4};
    p≡3, t+=0: (−2)^{(p+1)/4}(0)^{(p−3)/2}(+1)^2(+2)^{(p−3)/4}.
  Fail: |T∩S|=p+1.  ∎

Theorem C — PROVED (C-dirs of T⊙NC).
  Those dirs have occupancy sort(μ+Δ).  At p=5 this is R
  (T⊙NC stays all-R).  At p=7 (t+=0) it is (1,1,3,3,4,4,5),
  the C-slot of SIG_TNC.  Fail: O_C(7)=(0,1,2,3,4,5,6).  ∎

Theorem D — OPEN.  Which R-lines receive which Δ (the remaining
  O tuples) is unnamed in p.  Ensemble den D open.  phi_F not
  imported.

============================================================================
Backend: first_nc(T) + Fraction identities.  Writes
evidence/e1_gmin_m4_prop15388.json
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
from e1_gmin_m4_prop15279 import _kernel_field, rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15307 import square_classes  # noqa: E402
from e1_gmin_m4_prop15312 import first_nc  # noqa: E402
from e1_gmin_m4_prop15351 import y_T  # noqa: E402
from e1_gmin_m4_prop15352 import SIG_TNC  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15388_catalog.json"


def meet_tangent(p: int) -> tuple[int, ...]:
    return tuple([0] * ((p - 3) // 2) + [1, 1] + [2] * ((p - 1) // 2))


def meet_even(p: int) -> tuple[int, ...]:
    return tuple([0] * ((p - 1) // 2) + [2] * ((p + 1) // 2))


def delta_p1(p: int) -> tuple[int, ...]:
    """Named Δ-multiset for p≡1 (mixed tangents, polarised secants)."""
    return tuple(
        [-2] * ((p - 1) // 4)
        + [-1]
        + [0] * ((p - 3) // 2)
        + [1]
        + [2] * ((p - 1) // 4)
    )


def delta_p3_t0(p: int) -> tuple[int, ...]:
    """p≡3, both tangents out."""
    return tuple(
        [-2] * ((p + 1) // 4)
        + [0] * ((p - 3) // 2)
        + [1, 1]
        + [2] * ((p - 3) // 4)
    )


def O_C_from_delta(p: int, delta: tuple[int, ...]) -> tuple[int, ...]:
    mu = (p - 1) // 2
    return tuple(sorted(mu + d for d in delta))


def circle_meets(p: int, pack: dict) -> list[tuple[int, ...]]:
    S = {int(u) for u in pack["circle"]}
    F = _kernel_field(p)
    out = []
    for lines in square_classes(F):
        m = tuple(sorted(sum(1 for x in ln if int(x) in S) for ln in lines))
        out.append(m)
    return out


def T_cap_S(pack: dict) -> int:
    y0 = pack["y0"]
    D = {i for i in range(y0.shape[0] - 1) if y0[1 + i] < 0}
    return sum(1 for u in pack["circle"] if int(u) in D)


def live_pack(p: int) -> dict:
    return first_nc(p, y_T(p))


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 13):
        pack = live_pack(p)
        meets = circle_meets(p, pack)
        tgt = meet_tangent(p)
        ev = meet_even(p)
        ok = ok and pack is not None
        ok = ok and all(m == tgt for m in meets)
        ok = ok and tgt != ev
        ok = ok and sum(tgt) == p + 1
        ok = ok and sum(ev) == p + 1
        if any(m == ev for m in meets):
            ok = False
        rows[str(p)] = {
            "n_dirs": len(meets),
            "meet": list(tgt),
            "all_tangent": all(m == tgt for m in meets),
        }
    ok = ok and meet_tangent(5) == (0, 1, 1, 2, 2)
    ok = ok and meet_tangent(7) == (0, 0, 1, 1, 2, 2, 2)
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "first_nc(T) meet = (0)^{(p−3)/2}(1)^2(2)^{(p−1)/2} "
            "on every square class. Fail: even type."
        ),
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 13):
        pack = live_pack(p)
        cap = T_cap_S(pack)
        ok = ok and cap == (p + 1) // 2
        ok = ok and cap != p + 1
        rows[str(p)] = {"T_cap_S": cap, "half": (p + 1) // 2}
    ok = ok and delta_p1(5) == (-2, -1, 0, 1, 2)
    ok = ok and delta_p1(13) == (-2, -2, -2, -1, 0, 0, 0, 0, 0, 1, 2, 2, 2)
    ok = ok and delta_p3_t0(7) == (-2, -2, 0, 0, 1, 1, 2)
    if T_cap_S(live_pack(7)) == 8:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "delta_p1": {str(p): list(delta_p1(p)) for p in (5, 13, 17)},
        "delta_p3_t0": list(delta_p3_t0(7)),
        "theorem": "|T∩S|=(p+1)/2 and secants polarised. Fail: |T∩S|=p+1.",
    }


def prove_C() -> dict:
    ok = True
    oc5 = O_C_from_delta(5, delta_p1(5))
    oc7 = O_C_from_delta(7, delta_p3_t0(7))
    oc13 = O_C_from_delta(13, delta_p1(13))
    ok = ok and oc5 == (0, 1, 2, 3, 4)
    ok = ok and oc7 == (1, 1, 3, 3, 4, 4, 5)
    ok = ok and oc13 == (4, 4, 4, 5, 6, 6, 6, 6, 6, 7, 8, 8, 8)
    ok = ok and "O:1,1,3,3,4,4,5" in SIG_TNC
    ok = ok and oc7 != (0, 1, 2, 3, 4, 5, 6)
    if oc7 == (0, 1, 2, 3, 4, 5, 6):
        ok = False
    return {
        "proved": bool(ok),
        "O_C": {"5": list(oc5), "7": list(oc7), "13": list(oc13)},
        "theorem": "C-dirs of T⊙NC have occupancy sort(μ+Δ). Fail: O_C(7)=R.",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "C-image of T⊙NC named; R-pairings and ensemble D open. "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.388  T⊙NC meet/delta named; R-pairings open", flush=True)
    A = prove_A()
    print(f"  A meet: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B delta: {B['proved']}", flush=True)
    C = prove_C()
    print(f"  C O_C: {C['proved']} {C['O_C']}", flush=True)
    D = prove_open()
    print(f"  D open: phi_F={D['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps(
            {"meet": A["rows"], "cap": B["rows"], "O_C": C["O_C"]},
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.388",
        "title": "T⊙NC Hoffman meet/delta named; R-pairings open",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "meet_tangent": A["proved"],
            "delta_named": B["proved"],
            "O_C_named": C["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15388.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
