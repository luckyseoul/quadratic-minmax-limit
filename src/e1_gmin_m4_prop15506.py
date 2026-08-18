#!/usr/bin/env python3
"""
Prop 15.506 — S0 is Max+-free: T acts regularly on Ω, so
    ∑_{r∈T} Q(r)/q² = 2(q−1),     ∑_{off} (4 − Q/q²) = 8.
The p=7 three-type A-shifts that recover D=409 are not a p≡3 law
(D=77097/7 ∉ ℤ at p=11).  Ensemble Q_τ unnamed.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.505 flags.

============================================================================
Setup.  Ω={ξ : χ̂(ξ)=p}, χ̂(ξ)=∑_x χ(x) e(Tr(ξ x)).  T=□ ⊂ F_q^×.
Q(r)=|Ω|⁻¹ ∑_{ξ∈Ω} u(ξ) u(rξ), u=|ẑ|².  15.296 B: ∑_Ω u = q(q−1)
on every Max+.  15.491 certified the two sums on the p=5,7 caches.

============================================================================
Theorem A — PROVED (quadratic Gauss sum; certified p=5,7,11,13,17).
  χ̂(ξ)=χ(ξ)⁻¹ χ̂(1) and χ̂(1)=±p, so |Ω|=(q−1)/2.  For r∈T,
  χ̂(rξ)=χ̂(ξ), hence T·Ω=Ω.  The action is free, therefore regular
  (and transitive).  Fail: |Ω|=p+1; claim Ω=T at p=5 (χ̂(1)=−5).  ∎

Theorem B — PROVED (A + 15.296 B).
  ∑_{r∈T} u(rξ) = ∑_Ω u = q(q−1) for every ξ∈Ω, so
      ∑_{r∈T} Q(r) = 2 q² (q−1).
  Fail: 2 q² (q+1).  Recovers 15.491 at p=5,7 from the named
  two/three-type live Q values (no new census).  ∎

Theorem C — PROVED (B + Wick Q(±1)=8q²).
  q=p²≡1 (mod 8) ⇒ −1∈T, so |off|=(q−1)/2−2.
  ∑_{r∈T\\{±1}} (4 − Q(r)/q²) = 8.  Fail: 4 or 16.  ∎

Theorem D — PROVED (S0 linear algebra; not a D-name).
  The p=7 identities Q++=8A/D, Q_N*=8(A−n_pp)/D, Q_{--N1}=8(A−k)/D
  plus (C) recover D=409.  The same shifts at p=11 give D=77097/7,
  hence |H+|=2p D=1696134/7 ∉ ℤ.  Fail: claim they are a p≡3 law
  or that 77097/7 names |H+|/(2p).  ∎

Theorem E — OPEN.  S0 does not name the off-values.  Two-type 8A/D
  + (C) still forces D∉ℤ at p=13 (15.502).  phi_F not imported.  ∎

============================================================================
Backend: Gauss on F_{p²} + Fraction S0.  Writes
evidence/e1_gmin_m4_prop15506.json
"""
from __future__ import annotations

import json
import os
from fractions import Fraction
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import _e_tr, _kernel_field
from e1_gmin_m4_prop15290 import omega_set
from e1_gmin_m4_prop15330 import LIVE_QPP
from e1_gmin_m4_prop15355 import n_pp_named
from e1_gmin_m4_prop15367 import A_num, k_named
from e1_gmin_m4_prop15490 import LIVE_Q_MM_NS

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "evidence" / "e1_gmin_m4_prop15506.json"

CERT_OMEGA = (5, 7, 11, 13, 17)
LIVE_QMM = {7: Fraction(1376, 409)}


def n_Omega_named(p: int) -> int:
    return (p * p - 1) // 2


def n_Omega_wrong(p: int) -> int:
    """Fail-when-wrong: |Ω|=p+1."""
    return p + 1


def n_off(p: int) -> int:
    return n_Omega_named(p) - 2


def n_mm_named(p: int) -> int:
    return 0 if p % 4 == 1 else (p + 1) // 2


def n_ns_named(p: int) -> int:
    return n_off(p) - n_pp_named(p) - n_mm_named(p)


def sum_T_Q_over_q2_named(p: int) -> int:
    return 2 * (p * p - 1)


def sum_T_Q_over_q2_wrong(p: int) -> int:
    """Fail-when-wrong: 2(q+1)."""
    return 2 * (p * p + 1)


def live_Qn(p: int) -> Fraction:
    if p == 5:
        return LIVE_Q_MM_NS[5]
    return LIVE_Q_MM_NS[7]


def sum_T_Q_from_live_types(p: int) -> Fraction:
    """∑_T Q/q² from named class sizes and live coarse Q."""
    s = Fraction(16)  # Q(±1)/q² = 8+8
    s += n_pp_named(p) * LIVE_QPP[p]
    s += n_ns_named(p) * live_Qn(p)
    if p == 7:
        s += n_mm_named(p) * LIVE_QMM[7]
    return s


def sum_off_delta_from_live_types(p: int) -> Fraction:
    s = n_pp_named(p) * (4 - LIVE_QPP[p])
    s += n_ns_named(p) * (4 - live_Qn(p))
    if p == 7:
        s += n_mm_named(p) * (4 - LIVE_QMM[7])
    return s


def D_threetype_npp_k(p: int) -> Fraction:
    """S0 + Q++=8A/D, Qn=8(A−n_pp)/D, Qmm=8(A−k)/D."""
    A = A_num(p)
    NUM = (
        n_pp_named(p) * A
        + n_ns_named(p) * (A - n_pp_named(p))
        + n_mm_named(p) * (A - k_named(p))
    )
    return Fraction(2 * NUM, n_off(p) - 2)


def omega_row(p: int) -> dict:
    F = _kernel_field(p)
    Om = omega_set(F)
    T = [r for r in range(1, F["q"]) if F["chi"][r] == 1]
    g1 = sum(F["chi"][x] * _e_tr(F, x) for x in range(1, F["q"]))
    xi0 = Om[0]
    orb = {F["mul"][r][xi0] for r in T}
    return {
        "n_Omega": len(Om),
        "named": n_Omega_named(p),
        "G1": float(g1.real),
        "Omega_eq_T": set(Om) == set(T),
        "transitive": orb == set(Om),
    }


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in CERT_OMEGA:
        r = omega_row(p)
        rows[str(p)] = r
        if r["n_Omega"] != n_Omega_named(p):
            ok = False
        if r["n_Omega"] == n_Omega_wrong(p):
            ok = False
        if not r["transitive"]:
            ok = False
    if rows["5"]["Omega_eq_T"]:
        ok = False
    if not rows["7"]["Omega_eq_T"]:
        ok = False
    if abs(rows["5"]["G1"] + 5) > 1e-6:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": "|Ω|=(q−1)/2; T acts regularly. Fail |Ω|=p+1; Ω=T at p=5.",
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        got = sum_T_Q_from_live_types(p)
        named = Fraction(sum_T_Q_over_q2_named(p))
        wrong = Fraction(sum_T_Q_over_q2_wrong(p))
        rows[str(p)] = {"got": str(got), "named": str(named)}
        if got != named:
            ok = False
        if got == wrong:
            ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": "∑_T Q/q² = 2(q−1). Fail 2(q+1).",
    }


def prove_C() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        got = sum_off_delta_from_live_types(p)
        rows[str(p)] = str(got)
        if got != 8:
            ok = False
        if got in (4, 16):
            ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": "∑_off δ = 8. Fail 4 or 16.",
    }


def prove_D() -> dict:
    d7 = D_threetype_npp_k(7)
    d11 = D_threetype_npp_k(11)
    ok = d7 == 409 and d11.denominator != 1
    if d11 == 409 or d11.denominator == 1:
        ok = False
    return {
        "proved": bool(ok),
        "D7": str(d7),
        "D11": str(d11),
        "Hplus_p11": str(2 * 11 * d11),
        "theorem": (
            "p=7 three-type A-shifts recover D=409; at p=11 D=77097/7 "
            "makes |H+|=2pD=1696134/7 not an integer. Not a p≡3 law."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "Q_tau_named_in_p": False,
        "phi_F_imported": False,
        "note": (
            "S0 is Max+-free but does not name the off Q values. "
            "p=7 three-type shifts die at p=11. phi_F not imported."
        ),
    }


def main() -> dict:
    print("15.506 S0 Max+-free; p=7 three-type shifts die at p=11", flush=True)
    A, B, C, D, E = prove_A(), prove_B(), prove_C(), prove_D(), prove_open()
    out = {
        "prop": "15.506",
        "title": "S0 is Max+-free; p=7 A-shifts are not a p≡3 law",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "C": C,
        "D": D,
        "E": E,
        "e1_closed_general": bool(e1_closed_general()),
        "gsum_disj_lb": bool(gsum_disj_lb_proved_general()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "flags_not_flipped": [
            "phi_F_ge_6",
            "e1",
            "L",
            "Aut-Schur",
            "Gsum",
            "pairing",
        ],
        "L_status": "OPEN",
    }
    EV.parent.mkdir(parents=True, exist_ok=True)
    EV.write_text(json.dumps(out, indent=2) + "\n")
    print(
        f"A={A['proved']} B={B['proved']} C={C['proved']} D={D['proved']} "
        f"phi_F={out['phi_F_ge_6']}",
        flush=True,
    )
    return out


if __name__ == "__main__":
    main()
