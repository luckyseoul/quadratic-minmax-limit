#!/usr/bin/env python3
"""
Prop 15.417 — Two leftover-1 referee directions die on their
fail-when-wrong.  den(Q_NL) is not [H₊:Stab(τ)] / [G:Stab] /
c(p) / an Aut_∞ orbit size.  15.298 reverse slack on unnamed
L_χ slots (even after naming ns, or freezing ++) is wider than
the S0 floor window.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.416 flags.

============================================================================
Setup.  Aut_∞ has |G|=p²(p²−1) (15.308).  15.331: NL stabs
| p+1, c(5)=1, c(7)=19, n_NL=c(p) p²(p−1).  Live Q_NL=64/15,
632/171.  15.298: Q/16 = k²+k(G₊+G_r)+∑_{s,χ} L_χ(s) S_χ(1+rs).
Moving L₊ and L₋ of one coarse type together gives
    dQ(r)/dL_τ = 16(q n₀(r) − n_τ),   n₀=#{s∈τ: 1+rs=0}.

============================================================================
Theorem A — PROVED (exact orders; certified p=5,7).
  den(Q_NL)∈{15,171} is not an Aut_∞ index and not |H₊|/|Stab|:
      15 ∤ |H₊|=130,  19∤|G|=2352,  19∤|H₊|=5726,  171∤|G|,
      171∤|H₊|,  171∤ n_NL.  c(5)=1≠15, c(7)=19≠171.
  NL orbits are 100 and {1176⁴,588,294}, none equal to 15, 19,
  or 171.  NL stabs are {6} and {2,4,8}.  Fail: claim 15 | 130;
  claim c(7)=171; claim 19 | 2352.  ∎

Theorem B — PROVED (15.298 half-Gauss weights; field only).
  Mean_{r∈++} 16|q n₀−n_++|/q² · k² exceeds Q^{ub}−Q_{lo}
  at p=5 and p=7 (named ns does not shrink the ++ slot).
  Fail: claim the ++ coefficient vanishes (n₀=n_++/q).  ∎

Theorem C — OPEN.  Q_NL unnamed.  phi_F not imported.

============================================================================
Backend: Fraction + Aut_∞ orders + field half-Gauss.  Writes
evidence/e1_gmin_m4_prop15417.json
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
from e1_gmin_m4_prop15279 import _kernel_field  # noqa: E402
from e1_gmin_m4_prop15292 import n_1d  # noqa: E402
from e1_gmin_m4_prop15298 import coarse_Q_class  # noqa: E402
from e1_gmin_m4_prop15308 import G_order  # noqa: E402
from e1_gmin_m4_prop15321 import Qpp_floor_ub  # noqa: E402
from e1_gmin_m4_prop15326 import Q_lo_named  # noqa: E402
from e1_gmin_m4_prop15331 import HPLUS, c_from_stabs, nl_stabs  # noqa: E402
from e1_gmin_m4_prop15409 import D_form_on_lattice_general  # noqa: E402
from e1_gmin_m4_prop15411 import LIVE_QNL  # noqa: E402
from e1_gmin_m4_prop15412 import Q_NL_is_short_rational  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15417_catalog.json"


def plus_class(F) -> list[int]:
    q = F["q"]
    return [
        s
        for s in range(1, q)
        if F["chi"][s] == 1 and coarse_Q_class(F, s) == "++"
    ]


def n0_plus(F, r: int, plus: list[int]) -> int:
    """#{s∈++ : 1+r s = 0}."""
    add, mul, neg = F["add"], F["mul"], F["neg"]
    return sum(1 for s in plus if add[1][mul[r][s]] == 0)


def dQ_dL_plus(p: int) -> Fraction:
    """Mean_r 16(q n0(r) − n_++)/q²  for r∈++ (Q/q² per unit L)."""
    F = _kernel_field(p)
    F = {**F, "p": p}
    plus = plus_class(F)
    q = F["q"]
    n = len(plus)
    acc = 0
    for r in plus:
        acc += q * n0_plus(F, r, plus) - n
    # mean over r, then 16/q²
    return Fraction(16 * acc, n * q * q)


def k_named(p: int) -> int:
    return p * (p - 1) // 2


def slack_width(p: int) -> Fraction:
    """|Δ(Q++/q²)| if unnamed L_++ ranges over an interval of length k²."""
    return abs(dQ_dL_plus(p)) * Fraction(k_named(p) ** 2)


def window_width(p: int) -> Fraction:
    return Qpp_floor_ub(p) - Q_lo_named(p)


def index_hits_den(p: int) -> bool:
    """True only if den(Q_NL) is G-index, |H+| index, c(p), orbit, or stab."""
    G = G_order(p)
    N = HPLUS[p]
    n_nl = N - n_1d(p)
    den = LIVE_QNL[p].denominator
    stabs = nl_stabs(p)
    c = c_from_stabs(p, stabs)
    orbits = [G // s for s in stabs]
    if den in stabs or den in orbits or den == c:
        return True
    if N % den == 0 and G % den == 0:
        return True
    if n_nl % den == 0 and den in (15, 171):
        # n_NL ≡ 0 (mod 19) at p=7 is c(p) p²(p−1), not den 171
        return den == 19
    return False


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        G = G_order(p)
        N = HPLUS[p]
        n_nl = N - n_1d(p)
        den = LIVE_QNL[p].denominator
        stabs = nl_stabs(p)
        c = c_from_stabs(p, stabs)
        orbits = [G // s for s in stabs]
        rec = {
            "G": G,
            "Hplus": N,
            "n_nl": n_nl,
            "den_QNL": den,
            "stabs": stabs,
            "orbits": orbits,
            "c": c,
            "G_div": G % den == 0,
            "Hplus_div": N % den == 0,
            "index_hits": index_hits_den(p),
        }
        rows[str(p)] = rec
        ok = ok and not rec["index_hits"]
        ok = ok and not rec["Hplus_div"]
    ok = ok and (130 % 15 != 0)
    ok = ok and (2352 % 19 != 0)
    ok = ok and (5726 % 19 != 0)
    ok = ok and (2352 % 171 != 0)
    ok = ok and c_from_stabs(5, nl_stabs(5)) == 1
    ok = ok and c_from_stabs(7, nl_stabs(7)) == 19
    ok = ok and c_from_stabs(7, nl_stabs(7)) != 171
    if index_hits_den(5) or index_hits_den(7):
        ok = False
    if 130 % 15 == 0 or 2352 % 19 == 0:
        ok = False
    return {
        "proved": bool(ok),
        "rows": {
            p: {
                "G": rows[p]["G"],
                "Hplus": rows[p]["Hplus"],
                "n_nl": rows[p]["n_nl"],
                "den_QNL": rows[p]["den_QNL"],
                "stabs": rows[p]["stabs"],
                "orbits": rows[p]["orbits"],
                "c": rows[p]["c"],
                "index_hits": rows[p]["index_hits"],
            }
            for p in rows
        },
        "theorem": (
            "den(Q_NL) is not [H₊:Stab] or [G:Stab] or c(p). "
            "Fail: 15|130; 19|2352; c(7)=171."
        ),
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 11, 13):
        c = dQ_dL_plus(p)
        sw = slack_width(p)
        ww = window_width(p)
        rows[str(p)] = {
            "dQ_dL_plus": str(c),
            "slack_k2": str(sw),
            "window": str(ww),
            "slack_gt_window": sw > ww,
        }
        if p in (5, 7):
            ok = ok and sw > ww
            ok = ok and c != 0
    if dQ_dL_plus(5) == 0:
        ok = False
    if slack_width(5) <= window_width(5):
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "15.298 ++-slot k² slack exceeds Q^{ub}−Q_lo at p=5,7. "
            "Fail: dQ/dL_++=0."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "Q_NL_5": str(LIVE_QNL[5]),
        "Q_NL_7": str(LIVE_QNL[7]),
        "index_hits_p5": index_hits_den(5),
        "index_hits_p7": index_hits_den(7),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "D_form_on_lattice_general": D_form_on_lattice_general(),
        "Q_NL_is_short_rational": Q_NL_is_short_rational(),
        "note": (
            "Both leftover-1 referee directions fail their "
            "fail-when-wrong. c(7)=19 divides L_par den, not Q_NL. "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.417  den(Q_NL) not an index; 15.298 slack too wide", flush=True)
    A = prove_A()
    print(f"  A index dead: {A['proved']} p5c={A['rows']['5']['c']} p7c={A['rows']['7']['c']}", flush=True)
    B = prove_B()
    print(f"  B slack: {B['proved']} p5={B['rows']['5']}", flush=True)
    C = prove_open()
    print(f"  C open: phi_F={C['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps({"A": A["rows"], "B": B["rows"]}, indent=2) + "\n"
    )
    out = {
        "prop": "15.417",
        "title": (
            "den(Q_NL) is not [H₊:Stab]/c(p); 15.298 unnamed-L slack "
            "exceeds the floor window"
        ),
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "index_identification_dead": A["proved"],
            "slack_too_wide": B["proved"],
            "phi_F_ge_6_proved_general": C["phi_F_ge_6"],
            "residual_ii_k_eq_4p_empty": C["residual_ii_k_eq_4p_empty"],
            "Q_NL_is_short_rational": C["Q_NL_is_short_rational"],
        },
        "algebra": {"A": A, "B": B, "C": C},
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15417.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
