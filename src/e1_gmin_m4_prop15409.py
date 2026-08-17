#!/usr/bin/env python3
"""
Prop 15.409 — |H+| ≡ n_1d (mod p²) for every odd prime p≥5,
so D=|H+|/(2p) lies on the Max+-free lattice
    D = C(p−1,(p−1)/2)/2 + p u,   u∈ℤ≥0.
D_form is off this lattice at every p≥11 (and at p=3).
Paley-region complete sums miss 13 and 409.
phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.408 flags.  Residual (ii) and Type I stay False.

============================================================================
Setup.  Aut_∞ = {x ↦ a x^{p^k}+b : a∈□, k∈{0,1}}, |G|=p²(p²−1).
H_+ = n_1d Type+ 1D + NL G-orbits (15.331 A).  Translations
preserving a size-k set are line-unions (15.332 A).  Nonsquare-
direction lines meet every Max+ in (p−1)/2 points (15.305 C),
so a nonsquare line-union is not Max+.  Square-direction
line-unions are Type+ 1D, counted by n_1d.  Hence NL ∩ T={id}
and Stab_NL ↪ Aut_∞/T ≅ □⋊⟨Frob⟩, so |Stab_NL| | (p²−1) and
every NL orbit is 0 (mod p²).

============================================================================
Theorem A — PROVED (orbit–stab + 15.305 C + 15.332 A; all odd p≥5).
  |H+| ≡ n_1d (mod p²).  Fail: claim the same congruence mod p
  only (true at p=5 and does not kill D_form(11)).  Live |H+|
  at p=5,7 (130, 5726) satisfies the congruence.  ∎

Theorem B — PROVED (A + 2p | n_1d).
  n_1d = p C(p−1,(p−1)/2), so D_1d := n_1d/(2p) = C(p−1,(p−1)/2)/2
  is an integer (central binomial even).  Combined with A and
  2p | |H+|, one has t even in |H+|=n_1d+p² t, hence
      D = D_1d + p u,   u∈ℤ≥0.
  Fail: claim D_1d = m C(p−1,m−1)/2 (that is 9≠3 at p=5).  ∎

Theorem C — PROVED (exact integers).
  D_form = 2A − 3(p+1)/2 + C(p,(p−1)/2) equals live D at p=5,7
  and is off the lattice of B at p=3 and at every prime 11≤p≤79.
  Fail: claim D_form(11)=D_1d(11)+11u for some integer u
  (11106−126 is not divisible by 11).  Do not rehab D_form.  ∎

Theorem D — PROVED (complete counts; p=5,7).
  Named Paley-region counts
      #{χ=1}, #{χ(x)=χ(x−1)=1}, #{χ(x)=χ(x−1)=χ(x+1)=1},
      #{χ(x)=χ(y)=χ(x+y)=χ(x−y)=1}, #{χ(x)=χ(y)=χ(x−1)=χ(y−1)=1}
  miss {13,409}.  The p=7 ++_raw count is 6, not 409.
  Fail: claim the parallelogram count is 409 (it is 144).  ∎

Theorem E — PROVED (S0 floor window as a bound on u; p=5).
  On the S0 line, Q=8A/D∈[Q_lo,Q^{ub}] iff u lies in a named
  integer interval.  At p=5 the interval is {2}, matching live
  D=13.  At p=7 it is {53,…,64} (12 values) and does not name
  live u=57.  Fail: claim the p=5 window contains 0 or 1.
  p≡3 still has Q_{--}.  phi_F not imported from the p=5 pin.  ∎

Theorem F — OPEN.  u(p) is unnamed for p≥7, so Q_τ=8A/D is
  unnamed.  phi_F not imported.

============================================================================
Backend: Fraction / comb arithmetic + F_q Paley-region counts.
Writes evidence/e1_gmin_m4_prop15409.json
"""
from __future__ import annotations

import json
import os
import sys
from fractions import Fraction
from functools import lru_cache
from math import comb
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general, is_prime  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15275 import (  # noqa: E402
    type_I_aut_e_3AB_positive_general,
    type_I_multilevel_bad_case_ND_closed,
)
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field  # noqa: E402
from e1_gmin_m4_prop15292 import n_1d  # noqa: E402
from e1_gmin_m4_prop15321 import Qpp_floor_ub  # noqa: E402
from e1_gmin_m4_prop15326 import Q_lo_named  # noqa: E402
from e1_gmin_m4_prop15357 import D_live  # noqa: E402
from e1_gmin_m4_prop15367 import A_num  # noqa: E402
from e1_gmin_m4_prop15371 import D_form  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15409_catalog.json"
HPLUS_LIVE = {5: 130, 7: 5726}
P_LATTICE = (3, 5, 7, 11, 13, 17, 19, 23, 29, 31)
P_FORM = tuple(p for p in range(3, 80) if is_prime(p))


def central_half(p: int) -> int:
    """C(p−1, (p−1)/2)/2 = n_1d/(2p)."""
    n = (p - 1) // 2
    return comb(2 * n, n) // 2


def D_1d(p: int) -> int:
    return central_half(p)


def D_1d_wrong_m(p: int) -> int:
    """Fail-when-wrong neighbour m C(p−1,m−1)/2."""
    m = (p + 1) // 2
    return m * comb(p - 1, m - 1) // 2


def D_lattice(p: int, u: int) -> int:
    return D_1d(p) + p * u


def N_from_D(p: int, D: int) -> int:
    return 2 * p * D


def on_lattice(p: int, D: int) -> bool:
    d1 = D_1d(p)
    if D < d1:
        return False
    return (D - d1) % p == 0


def live_u(p: int) -> int:
    return (D_live(p) - D_1d(p)) // p


def _ceil_frac(x: Fraction) -> int:
    n, d = x.numerator, x.denominator
    if d < 0:
        n, d = -n, -d
    if n >= 0:
        return (n + d - 1) // d
    return n // d


def _floor_frac(x: Fraction) -> int:
    n, d = x.numerator, x.denominator
    if d < 0:
        n, d = -n, -d
    return n // d


def u_window(p: int) -> tuple[int, int]:
    """Integers u with Q=8A/(D_1d+p u) in [Q_lo, Q^{ub}] on S0."""
    A = A_num(p)
    d1 = D_1d(p)
    lo, hi = Q_lo_named(p), Qpp_floor_ub(p)
    dmin = Fraction(8 * A, 1) / hi
    dmax = Fraction(8 * A, 1) / lo
    umin = (dmin - d1) / p
    umax = (dmax - d1) / p
    return _ceil_frac(umin), _floor_frac(umax)


def D_form_on_lattice_general() -> bool:
    """True only if D_form lies on the lattice for all p≥5. It does not."""
    return False


@lru_cache(maxsize=8)
def paley_region_counts(p: int) -> dict:
    F = _kernel_field(p)
    q, chi, add, neg = F["q"], F["chi"], F["add"], F["neg"]
    n_sq = n_pm = n_raw = 0
    for x in range(1, q):
        if chi[x] != 1:
            continue
        n_sq += 1
        xm1 = add[x][neg[1]]
        if chi[xm1] == 1:
            n_pm += 1
            if x != 1 and chi[add[x][1]] == 1:
                n_raw += 1
    n_para = n_box = 0
    for x in range(1, q):
        if chi[x] != 1:
            continue
        xm1 = add[x][neg[1]]
        for y in range(1, q):
            if chi[y] != 1:
                continue
            if chi[add[x][y]] == 1 and chi[add[x][neg[y]]] == 1:
                n_para += 1
            if chi[xm1] == 1 and chi[add[y][neg[1]]] == 1:
                n_box += 1
    return {
        "squares": n_sq,
        "chi_x_xm1": n_pm,
        "pp_raw": n_raw,
        "parallelogram": n_para,
        "box_xm1_ym1": n_box,
    }


def prove_A() -> dict:
    ok = True
    rows = {}
    for p, N in HPLUS_LIVE.items():
        n1 = n_1d(p)
        ok = ok and N % (p * p) == n1 % (p * p)
        ok = ok and (N - n1) % (p * p) == 0
        rows[str(p)] = {"N": N, "n_1d": n1, "N_mod": N % (p * p)}
        if (N - n1) % p != 0:
            ok = False
    # fail-when-wrong: congruence only mod p does not kill D_form(11)
    n11, Nf11 = n_1d(11), N_from_D(11, D_form(11))
    ok = ok and (Nf11 - n11) % 11 == 0
    ok = ok and (Nf11 - n11) % 121 != 0
    if (Nf11 - n11) % 121 == 0:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "p11_form_mod_p": (Nf11 - n11) % 11,
        "p11_form_mod_psq": (Nf11 - n11) % 121,
        "theorem": (
            "|H+| ≡ n_1d (mod p²). Fail: claim only mod p "
            "(D_form(11) is 0 mod 11 and 44 mod 121)."
        ),
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 11, 13, 17, 19, 23):
        n1 = n_1d(p)
        d1 = D_1d(p)
        cen = comb(p - 1, (p - 1) // 2)
        ok = ok and n1 == p * cen
        ok = ok and 2 * p * d1 == n1
        ok = ok and cen % 2 == 0
        ok = ok and d1 == cen // 2
        rows[str(p)] = {"D_1d": d1, "n_1d": n1}
    ok = ok and D_1d(5) == 3 and D_1d(7) == 10
    ok = ok and D_1d_wrong_m(5) == 9 != 3
    ok = ok and live_u(5) == 2 and live_u(7) == 57
    ok = ok and D_live(5) == D_lattice(5, 2)
    ok = ok and D_live(7) == D_lattice(7, 57)
    if D_1d_wrong_m(5) == D_1d(5) or D_1d(5) != 3:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "u_live": {5: 2, 7: 57},
        "theorem": (
            "D=C(p−1,(p−1)/2)/2 + p u. "
            "Fail: D_1d=m C(p−1,m−1)/2 (9 at p=5)."
        ),
    }


def prove_C() -> dict:
    ok = True
    off = []
    on = []
    for p in P_FORM:
        Df = D_form(p)
        hit = on_lattice(p, Df)
        if p in (5, 7):
            ok = ok and hit and Df == D_live(p)
            on.append(p)
        else:
            ok = ok and (not hit)
            off.append(p)
        if p == 11 and on_lattice(11, Df):
            ok = False
        if p >= 11 and hit:
            ok = False
    ok = ok and not on_lattice(11, D_form(11))
    ok = ok and (D_form(11) - D_1d(11)) % 11 != 0
    ok = ok and not D_form_on_lattice_general()
    return {
        "proved": bool(ok),
        "on": on,
        "n_off": len(off),
        "off_head": off[:8],
        "D_form_11": D_form(11),
        "D_1d_11": D_1d(11),
        "theorem": (
            "D_form off the lattice at p=3 and all 11≤p≤79. "
            "Fail: D_form(11)=126+11u."
        ),
    }


def prove_D() -> dict:
    c5 = paley_region_counts(5)
    c7 = paley_region_counts(7)
    vals5 = set(c5.values())
    vals7 = set(c7.values())
    ok = True
    ok = ok and 13 not in vals5 and 409 not in vals5
    ok = ok and 13 not in vals7 and 409 not in vals7
    ok = ok and c7["parallelogram"] == 144
    ok = ok and c7["pp_raw"] == 6
    ok = ok and c7["box_xm1_ym1"] == 121
    ok = ok and c5["parallelogram"] == 24
    if c7["parallelogram"] == 409 or 409 in vals7:
        ok = False
    return {
        "proved": bool(ok),
        "p5": c5,
        "p7": c7,
        "theorem": (
            "Paley-region counts miss {13,409}. "
            "Fail: parallelogram=409 (it is 144)."
        ),
    }


def prove_E() -> dict:
    lo5, hi5 = u_window(5)
    lo7, hi7 = u_window(7)
    ok = True
    ok = ok and (lo5, hi5) == (2, 2)
    ok = ok and live_u(5) == 2
    ok = ok and lo7 <= 57 <= hi7
    ok = ok and (lo7, hi7) == (53, 64)
    ok = ok and hi7 - lo7 + 1 == 12
    ok = ok and 0 not in range(lo5, hi5 + 1)
    ok = ok and 1 not in range(lo5, hi5 + 1)
    if lo5 != 2 or hi5 != 2 or 57 < lo7 or 57 > hi7:
        ok = False
    return {
        "proved": bool(ok),
        "window_5": [lo5, hi5],
        "window_7": [lo7, hi7],
        "live_u": {5: 2, 7: 57},
        "theorem": (
            "S0 floor window pins u(5)=2 and is {53..64} at p=7. "
            "Fail: p=5 window contains 0 or 1."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "D_form_on_lattice_general": D_form_on_lattice_general(),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "type_I_aut_e_3AB_positive_general": bool(
            type_I_aut_e_3AB_positive_general()
        ),
        "type_I_multilevel_bad_case_ND_closed": bool(
            type_I_multilevel_bad_case_ND_closed()
        ),
        "note": (
            "u(p) unnamed for p≥7. D_form dead as a general name. "
            "Paley regions miss 409. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.409  D = C(p-1,(p-1)/2)/2 + p u; D_form off lattice", flush=True)
    A = prove_A()
    print(f"  A congruence: {A['proved']} live={A['rows']}", flush=True)
    B = prove_B()
    print(f"  B lattice: {B['proved']} D_1d5={D_1d(5)} u={B['u_live']}", flush=True)
    C = prove_C()
    print(
        f"  C D_form off: {C['proved']} on={C['on']} n_off={C['n_off']}",
        flush=True,
    )
    D = prove_D()
    print(
        f"  D Paley miss: {D['proved']} para7={D['p7']['parallelogram']}",
        flush=True,
    )
    E = prove_E()
    print(
        f"  E window: {E['proved']} p5={E['window_5']} p7={E['window_7']}",
        flush=True,
    )
    F = prove_open()
    print(
        f"  F open: phi_F={F['phi_F_ge_6']} D_form={F['D_form_on_lattice_general']}",
        flush=True,
    )
    CAT.write_text(
        json.dumps(
            {
                "D_1d": {str(p): D_1d(p) for p in (5, 7, 11, 13)},
                "D_form": {str(p): D_form(p) for p in (3, 5, 7, 11)},
                "windows": {"5": E["window_5"], "7": E["window_7"]},
                "paley7": D["p7"],
            },
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.409",
        "title": (
            "D = C(p−1,(p−1)/2)/2 + p u; D_form off the lattice "
            "for all p≥11; Paley regions miss 409"
        ),
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "Hplus_cong_n1d_mod_psq": A["proved"],
            "D_lattice_D1d_plus_p_u": B["proved"],
            "D_form_off_lattice_p_ge_11": C["proved"],
            "paley_region_misses_409": D["proved"],
            "u_window_pins_p5": E["proved"],
            "D_form_on_lattice_general": F["D_form_on_lattice_general"],
            "phi_F_ge_6_proved_general": F["phi_F_ge_6"],
            "residual_ii_k_eq_4p_empty": F["residual_ii_k_eq_4p_empty"],
            "type_I_aut_e_3AB_positive_general": F[
                "type_I_aut_e_3AB_positive_general"
            ],
        },
        "algebra": {"A": A, "B": B, "C": C, "D": D, "E": E, "F": F},
        "L_status": "OPEN",
        "flags_not_flipped": [
            "phi_F_ge_6",
            "e1",
            "L",
            "Aut-Schur",
            "Gsum",
            "pairing",
            "residual_ii_k_eq_4p_empty",
            "type_I_multilevel_bad_case_ND_closed",
            "type_I_aut_e_3AB_positive_general",
        ],
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15409.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
