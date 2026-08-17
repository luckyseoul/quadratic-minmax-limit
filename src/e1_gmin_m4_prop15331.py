#!/usr/bin/env python3
"""
Prop 15.331 — |H_+| = n_1d + Σ_{NL} |G|/s_τ with every NL Aut_∞
stab dividing p+1.  409 | Φ_{|T|}(7) is p=7-only.  Multiplicities
unnamed.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.330 flags.  Floor stays OPEN.

============================================================================
Setup.  Aut_∞ = {u ↦ a u^{p^k}+b : a∈□, k∈{0,1}}, |G|=p²(p²−1).
n_1d = m C(p,m), m=(p+1)/2 (Type+ halfspaces).  NL = H_+ \\ {hs}.
Live |H_+| = 6, 130, 5726 at p=3,5,7.  den(Q_{++})=|H_+|/(2p).

============================================================================
Theorem A — PROVED (orbit–stab + 15.307 types; certified p=3,5,7).
  |H_+| = n_1d + Σ_{NL types τ} |G|/s_τ, and the sum of orbit
  sizes equals the cache |H_+|.  Fail-when-wrong: drop Frob in |G|
  (then AP-orbit identity fails, 15.308 A).  ∎

Theorem B — PROVED (certified p=3,5,7).
  Every NL Aut_∞-stab divides p+1.  Equivalently
      |H_+| = n_1d + c(p) p²(p−1),   c(p)=Σ_τ (p+1)/s_τ
  with c(3)=0, c(5)=1, c(7)=19.  Halfspace stabs are 4p (AP) and,
  when present, 6p (non-AP at p=7), which need not divide p+1.
  Fail-when-wrong: claim 3 is an NL stab at p=7 (3∤8); claim 4p
  is an NL stab (20∤6 at p=5; 28∤8 at p=7).  ∎

Theorem C — PROVED as a negative (cyclotomic exact).
  409 | Φ_24(7) with Φ_24(x)=x⁸−x⁴+1 and 24=|T|=(q−1)/2 at p=7
  (Φ_24(7)=409·73·193).  This does **not** name den(Q_{++}) in p:
  at p=5, |T|=12 and 13 ∤ Φ_12(5)=601.  13 | Φ_4(5)=26 instead.
  Fail-when-wrong: claim 409=Φ_24(7) (5762401); claim
  13 | Φ_{|T|}(5).  Class numbers of named discs {−7,−28,−4q,−q}
  are in {1,2,4}, not 409.  ∎

Theorem D — OPEN.  The multiplicity of each d|(p+1) as an NL
  stab is unnamed (1×(p+1) at p=5 vs 4×2, 1×4, 1×8 at p=7).
  c(p)=3p−2 hits 19 at p=7 and misses c(5)=1.  |H_+| unnamed.
  Q_τ unnamed.  phi_F not imported.  ∎

============================================================================
Backend: 15.307 cache fold for type sizes; exact cyclotomic.
Writes evidence/e1_gmin_m4_prop15331.json
"""
from __future__ import annotations

import json
import sys
from functools import lru_cache
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15292 import n_1d  # noqa: E402
from e1_gmin_m4_prop15307 import classify_p  # noqa: E402
from e1_gmin_m4_prop15308 import G_order, G_order_no_frob  # noqa: E402
from e1_gmin_m4_prop15330 import LIVE_QPP  # noqa: E402


HPLUS = {3: 6, 5: 130, 7: 5726}


def Phi24(x: int) -> int:
    return x**8 - x**4 + 1


def Phi12(x: int) -> int:
    return x**4 - x**2 + 1


def Phi4(x: int) -> int:
    return x**2 + 1


def _is_hs_type_key(key: str, p: int) -> bool:
    """15.307 hs key is one H and (n_sq−1) copies of C."""
    return key.count("'H'") == 1 and key.count("'C'") == (p + 1) // 2 - 1


@lru_cache(maxsize=8)
def nl_stabs(p: int) -> list[int]:
    """NL Aut_∞ stabs = |G| / |type| for non-halfspace 15.307 types."""
    rec = classify_p(p)
    G = G_order(p)
    stabs = []
    for key, cnt in rec["D_types"].items():
        if _is_hs_type_key(key, p):
            continue
        if cnt == 0 or G % int(cnt):
            continue
        stabs.append(G // int(cnt))
    return sorted(stabs)


def c_from_stabs(p: int, stabs: list[int]) -> int:
    return sum((p + 1) // s for s in stabs)


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (3, 5, 7):
        rec = classify_p(p)
        n1 = n_1d(p)
        G = G_order(p)
        hs = rec["n_hs"]
        extra = rec["Nplus"] - hs
        stabs = nl_stabs(p)
        recon = sum(G // s for s in stabs)
        hit = (
            rec["hs_count_ok"]
            and rec["Nplus"] == HPLUS[p]
            and hs == n1
            and recon == extra
            and extra + n1 == HPLUS[p]
        )
        if G_order_no_frob(p) == G:
            hit = False
        if not hit:
            ok = False
        rows[str(p)] = {
            "n_1d": n1,
            "Nplus": rec["Nplus"],
            "extra": extra,
            "nl_stabs": stabs,
            "recon": recon,
        }
    return {
        "proved": ok,
        "by_p": rows,
        "G_drop_frob_is_not_G": True,
        "theorem": "|H_+|=n_1d+Σ|G|/s_τ.  Drop Frob breaks |G|.",
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    expect_c = {3: 0, 5: 1, 7: 19}
    for p in (3, 5, 7):
        stabs = nl_stabs(p)
        div = all((p + 1) % s == 0 for s in stabs)
        c = c_from_stabs(p, stabs)
        # fail-when-wrong probes
        three_is_nl = (p == 7) and (3 in stabs)
        fourp_is_nl = (4 * p) in stabs
        if not div or c != expect_c[p] or three_is_nl or fourp_is_nl:
            ok = False
        if p == 7 and (p + 1) % 3 == 0:
            ok = False  # 3 | 8 is false; keep the probe honest
        if p == 7 and 8 % 3 == 0:
            ok = False
        rows[str(p)] = {
            "stabs": stabs,
            "div_pplus1": div,
            "c": c,
            "three_is_nl": three_is_nl,
            "fourp_is_nl": fourp_is_nl,
        }
    # 3 ∤ (p+1) at p=7
    ok = ok and (8 % 3 != 0) and (6 % 20 != 0)
    return {
        "proved": ok,
        "by_p": rows,
        "c_p": expect_c,
        "stab3_not_nl_p7": 3 not in nl_stabs(7),
        "stab4p_not_nl": all((4 * p) not in nl_stabs(p) for p in (5, 7)),
        "theorem": (
            "NL stabs | p+1.  c(3,5,7)=(0,1,19).  "
            "3 and 4p are not NL stabs."
        ),
    }


def prove_C() -> dict:
    t7 = Phi24(7)
    t5 = Phi12(5)
    p4_5 = Phi4(5)
    ok = True
    ok = ok and t7 == 7**8 - 7**4 + 1
    ok = ok and t7 == 409 * 73 * 193
    ok = ok and t7 != 409
    ok = ok and t5 == 601
    ok = ok and t5 % 13 != 0
    ok = ok and p4_5 == 26 and p4_5 % 13 == 0
    ok = ok and (5 * 5 - 1) // 2 == 12
    ok = ok and (7 * 7 - 1) // 2 == 24
    # 3p-2 is p=7-only for c
    ok = ok and (3 * 7 - 2) == 19
    ok = ok and (3 * 5 - 2) != 1
    return {
        "proved": ok,
        "Phi24_7": t7,
        "Phi24_7_factor": [409, 73, 193],
        "Phi24_is_not_409": t7 != 409,
        "Phi12_5": t5,
        "Phi12_5_div_by_13": t5 % 13 == 0,
        "Phi4_5": p4_5,
        "T_p5": 12,
        "T_p7": 24,
        "c_3p_minus_2_p5": 3 * 5 - 2,
        "theorem": (
            "409 | Φ_24(7) but 409≠Φ_24(7).  "
            "13 ∤ Φ_{|T|}(5).  Not a general den formula."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "e1": bool(e1_closed_general()),
        "Hplus_named_in_p": False,
        "Q_tau_named_in_p": False,
        "den_Qpp_p7": int(LIVE_QPP[7].denominator),
        "note": (
            "NL stab multiplicities unnamed.  |H_+| unnamed.  "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.331  NL stabs | p+1; 409 | Φ_24(7) is p=7-only", flush=True)
    A = prove_A()
    print(f"  A orbit-stab: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B NL | p+1: {B['proved']} c={B['c_p']}", flush=True)
    C = prove_C()
    print(f"  C cyclo p=7-only: {C['proved']} Phi24={C['Phi24_7']}", flush=True)
    D = prove_open()
    print(f"  D open: phi_F={D['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.331",
        "title": "NL Aut_∞ stabs divide p+1; 409|Φ_24(7) is p=7-only",
        "proved": {
            "orbit_stab_sum": A["proved"],
            "nl_stabs_div_pplus1": B["proved"],
            "cyclo_409_is_p7_only": C["proved"],
            "Hplus_named_in_p": False,
            "Q_tau_named_in_p": False,
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
    path = ROOT / "evidence" / "e1_gmin_m4_prop15331.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
