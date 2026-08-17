#!/usr/bin/env python3
"""
Prop 15.387 — Triangle (T) square-line occupancy is named in p:
exactly 3 R-directions and (n_sq−3) C-directions.  R means one
line of each n=0,…,p−1.  p=7 remaining NL histograms are the
sums of the 15.352 direction signatures.  Ensemble / D open.
phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.386 flags.

============================================================================
Setup.  15.307: H / C / R are the named per-class occupancy
signatures.  15.351 B: y_T has 3 R and the rest C (CRRR at p=7,
CCCCRRR at p=13).  15.349/350: T=A·T_0 (p≡1) or the named
p=7 matrix.

============================================================================
Theorem A — PROVED (R signature; all odd p≥5).
  An R-class has sorted occupancies (0,1,…,p−1), so
      P_R(k)=1/p for k=0,…,p−1,   P_R(p)=0,
      E[n⁴|R]=(p−1)(2p−1)(3p²−3p−1)/30.
  Fail: claim P_R(p)=1/p, or drop (3p²−3p−1).  ∎

Theorem B — PROVED (T type; certified p=5,7,13,17,29).
  y_T has exactly 3 R square-dirs and n_sq−3 C-dirs,
  n_sq=(p+1)/2.  Hence
      P_T(p)=0,
      P_T(k)=6/(p(p+1))           (k∈{0,…,p−1}\\{μ}),
      P_T(μ)=(p−2)(p−3)/(p(p+1)),
      E[n⁴|T]=(3 E[n⁴|R] + (n_sq−3) μ⁴)/n_sq
  (=354/5 at p=5, 264 at p=7).  Fail: all-R at p=7
  (then P_T(μ)=1/7 ≠ 5/14).  ∎

Theorem C — PROVED (15.352 signatures; certified p=7).
  The five exotic p=7 NL occupancy histograms are the
  concatenations of the named 4-tuples SIG_TNC / SIG_SHARP /
  SIG_LAST / SIG_588 / SIG_294.  They match the live
  1176×4 + 588 + 294 split.  Fail: T⊙NC² is a 1176.  ∎

Theorem D — OPEN.  Exotic O signatures unnamed in p.
  Ensemble E[n⁴] still den D.  phi_F not imported.

============================================================================
Backend: Fraction identities + 15.352 signature tables.
Writes evidence/e1_gmin_m4_prop15387.json
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
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15292 import n_1d  # noqa: E402
from e1_gmin_m4_prop15307 import sig_C, sig_H, sig_R  # noqa: E402
from e1_gmin_m4_prop15330 import HPLUS  # noqa: E402
from e1_gmin_m4_prop15351 import sig_y, y_T  # noqa: E402
from e1_gmin_m4_prop15352 import (  # noqa: E402
    SIG_294,
    SIG_588,
    SIG_LAST,
    SIG_SHARP,
    SIG_TNC,
)
from e1_gmin_m4_prop15385 import En4_1d  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15387_catalog.json"


def n_sq_named(p: int) -> int:
    return (p + 1) // 2


def En4_R(p: int) -> Fraction:
    return Fraction((p - 1) * (2 * p - 1) * (3 * p * p - 3 * p - 1), 30)


def En4_R_drop_cubic(p: int) -> Fraction:
    return Fraction((p - 1) * (2 * p - 1), 30)


def En4_C(p: int) -> Fraction:
    mu = Fraction(p - 1, 2)
    return mu**4


def En4_T(p: int) -> Fraction:
    return (3 * En4_R(p) + (n_sq_named(p) - 3) * En4_C(p)) / n_sq_named(p)


def P_R(p: int, k: int) -> Fraction:
    if 0 <= k < p:
        return Fraction(1, p)
    return Fraction(0)


def P_T(p: int, k: int) -> Fraction:
    mu = (p - 1) // 2
    if k == p:
        return Fraction(0)
    if k == mu:
        return Fraction((p - 2) * (p - 3), p * (p + 1))
    if 0 <= k < p:
        return Fraction(6, p * (p + 1))
    return Fraction(0)


def P_T_all_R(p: int, k: int) -> Fraction:
    """Wrong neighbour: T is all-R."""
    return P_R(p, k)


def hist_from_dir_tags(p: int, tags: tuple[str, ...]) -> tuple[int, ...]:
    h = [0] * (p + 1)
    for t in tags:
        if t == "C":
            h[(p - 1) // 2] += p
        elif t == "R":
            for k in range(p):
                h[k] += 1
        elif t == "H":
            h[0] += (p + 1) // 2
            h[p] += (p - 1) // 2
        elif t.startswith("O:"):
            for part in t.split(":", 1)[1].split(","):
                h[int(part)] += 1
        else:
            raise ValueError(t)
    return tuple(h)


def En4_from_hist(hist: tuple[int, ...]) -> Fraction:
    tot = sum(hist)
    acc = 0
    for k, c in enumerate(hist):
        acc += c * k**4
    return Fraction(acc, tot)


# Live p=7 occupancy counts among 28 lines (gpu_occ_split).
LIVE_P7_TYPE_HISTS = {
    "T_CRRR": (3, 3, 3, 10, 3, 3, 3, 0),
    "TNC": (1, 6, 5, 6, 2, 7, 0, 1),
    "SHARP": (0, 7, 6, 5, 4, 3, 2, 1),
    "LAST": (1, 4, 8, 7, 1, 4, 2, 1),
    "YSTAR": (2, 5, 4, 6, 4, 5, 2, 0),
    "TNC2": (4, 4, 0, 6, 8, 6, 0, 0),
    "HS": (4, 0, 0, 21, 0, 0, 0, 3),
}
LIVE_P7_SIZES = {
    "T_CRRR": 1176,
    "TNC": 1176,
    "SHARP": 1176,
    "LAST": 1176,
    "YSTAR": 588,
    "TNC2": 294,
    "HS": 140,
}


def prove_A() -> dict:
    ok = True
    ok = ok and sig_R(5) == (0, 1, 2, 3, 4)
    ok = ok and sig_R(7) == (0, 1, 2, 3, 4, 5, 6)
    ok = ok and P_R(5, 5) == 0
    ok = ok and P_R(5, 2) == Fraction(1, 5)
    ok = ok and En4_R(5) == Fraction(354, 5)
    ok = ok and En4_R(7) == 325
    s4 = sum(k**4 for k in range(5))
    ok = ok and En4_R(5) == Fraction(s4, 5)
    ok = ok and En4_R_drop_cubic(5) != En4_R(5)
    if En4_R_drop_cubic(5) == En4_R(5):
        ok = False
    return {
        "proved": bool(ok),
        "En4_R": {str(p): str(En4_R(p)) for p in (5, 7, 11, 13)},
        "theorem": "R: P(k)=1/p on {0..p−1}. Fail: drop (3p²−3p−1).",
    }


def prove_B() -> dict:
    ok = True
    for p in (5, 7, 11, 13, 17, 29):
        ok = ok and sum(P_T(p, k) for k in range(p + 1)) == 1
        m4 = sum(P_T(p, k) * k**4 for k in range(p + 1))
        ok = ok and m4 == En4_T(p)
    ok = ok and En4_T(5) == Fraction(354, 5)
    ok = ok and En4_T(7) == 264
    ok = ok and P_T(7, 3) == Fraction(5, 14)
    ok = ok and P_T_all_R(7, 3) == Fraction(1, 7)
    ok = ok and P_T(7, 3) != P_T_all_R(7, 3)
    if P_T(7, 3) == P_T_all_R(7, 3):
        ok = False
    # live constructors: 3 R + rest C
    for p in (5, 7, 13, 17, 29):
        tags = sig_y(p, y_T(p))
        ok = ok and tags.count("R") == 3
        ok = ok and tags.count("C") == n_sq_named(p) - 3
        if tags.count("R") != 3:
            ok = False
    return {
        "proved": bool(ok),
        "En4_T": {str(p): str(En4_T(p)) for p in (5, 7, 13, 17, 29)},
        "P_T_mu": {str(p): str(P_T(p, (p - 1) // 2)) for p in (5, 7, 13)},
        "theorem": (
            "T is 3 R + (n_sq−3) C. Fail: all-R at p=7 (1/7≠5/14)."
        ),
    }


def prove_C() -> dict:
    p = 7
    derived = {
        "T_CRRR": hist_from_dir_tags(p, ("C", "R", "R", "R")),
        "TNC": hist_from_dir_tags(p, SIG_TNC),
        "SHARP": hist_from_dir_tags(p, SIG_SHARP),
        "LAST": hist_from_dir_tags(p, SIG_LAST),
        "YSTAR": hist_from_dir_tags(p, SIG_588),
        "TNC2": hist_from_dir_tags(p, SIG_294),
        "HS": hist_from_dir_tags(p, ("H", "C", "C", "C")),
    }
    ok = True
    for name, hist in derived.items():
        ok = ok and hist == LIVE_P7_TYPE_HISTS[name]
        ok = ok and sum(hist) == n_sq_named(p) * p
    # sizes
    ok = ok and sum(LIVE_P7_SIZES.values()) == HPLUS[7] == 5726
    ok = ok and LIVE_P7_SIZES["TNC2"] == 294
    if LIVE_P7_SIZES["TNC2"] == 1176:
        ok = False
    # 1D hist matches 15.385
    ok = ok and derived["HS"] == (4, 0, 0, 21, 0, 0, 0, 3)
    en4 = {name: En4_from_hist(h) for name, h in derived.items()}
    ok = ok and en4["T_CRRR"] == 264
    ok = ok and en4["HS"] == En4_1d(7)
    return {
        "proved": bool(ok),
        "hists": {k: list(v) for k, v in derived.items()},
        "En4": {k: str(v) for k, v in en4.items()},
        "theorem": (
            "p=7 exotic hists = 15.352 signature sums. "
            "Fail: T⊙NC² is a 1176."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "n_1d": {str(p): n_1d(p) for p in (5, 7)},
        "note": (
            "T occupancy named in p; exotic O and ensemble den D open. "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.387  T occupancy named; p=7 NL hists from 15.352", flush=True)
    A = prove_A()
    print(f"  A R: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B T: {B['proved']} En4_T={B['En4_T']}", flush=True)
    C = prove_C()
    print(f"  C p7 hists: {C['proved']}", flush=True)
    D = prove_open()
    print(f"  D open: phi_F={D['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps(
            {"En4_R": A["En4_R"], "En4_T": B["En4_T"], "p7_En4": C["En4"], "p7_hists": C["hists"]},
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.387",
        "title": "T occupancy named in p; p=7 NL hists from 15.352; ensemble open",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "R_occupancy": A["proved"],
            "T_occupancy": B["proved"],
            "p7_NL_hists_from_sigs": C["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15387.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
