#!/usr/bin/env python3
"""
Prop 15.429 — The 15.298 nonsquare Gauss coefficient is
    C_ns,χ(r) = −(q−1)/4
for both χ and every r∈T\\{±1}.  Field-only (Max+-free).
With L_ns=μ₊μ₋ the ns+base block of Q/q² is
    (q−1)(4p+5−p²)/(2p²)
and vanishes at p=5.  This does not name Q_NL or u in p.
phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.428 flags.

============================================================================
Setup.  15.298: Q(r)/16 = k²+k(G₊+G_r)+∑_{s,χ} L_χ(s) S_χ(1+rs).
NL first moment is μ± (15.428 fold), so the base is
A0=k(k+p)=p²(q−1)/4 and 16 A0/q²=4(q−1)/q.
L_ns=μ₊μ₋=p²(p−1)(p−3)/16 both χ (NL p=5,7).
S_χ(β)=(q−1)/2 if β=0, else (−1+χ(d)χ(β)G)/2, G=+p on Ω.

============================================================================
Theorem A — PROVED (half-Gauss + Jacobi; Max+-free; all odd p).
  χ(−1)=1 on F_{p²}, so 1+rs=0 forces s=−r^{−1} square when
  r is square: no β=0 term on ns.  1_ns=(1−χ)/2.
  ∑_{s≠0} χ(1+rs)=−1 and ∑_{s≠0} χ(s(1+rs))=−χ(r), hence
  ∑_{s ns} χ(1+rs)=(χ(r)−1)/2=0 on T.  Therefore
      C_ns,χ(r)=−(#ns)/2=−(q−1)/4
  both χ, every r∈T\\{±1}.  Fail: C_ns=0 (0≠−6 at p=5).  ∎

Theorem A2 — PROVED (Gauss of μ; Max+-free; all odd p).
  μ= a+bχ with a=p(p−2)/4, b=p/4.  G_μ(scale)=−a+b G(χ,α).
  On Ω, G(χ,ξ)=p; r square ⇒ G(χ,ξr)=p, so G_μ(r)=p/2.
  Hence G_μ(1)+G_μ(r)=p and A0=k(k+p)=p²(q−1)/4.
  Fail: claim G_μ(1)+G_μ(r)=0.  ∎

Theorem B — PROVED (A + named L_ns + named A0).
  The ns+base contribution to Q/q² is
      4(q−1)/q + (16/q²)·2·μ₊μ₋·(−(q−1)/4)
    = (q−1)(4p+5−p²)/(2p²).
  Equals 0 at p=5 and −384/49 at p=7.
  Fail: claim this equals Q_NL (0≠64/15; −384/49≠632/171).  ∎

Theorem C — OPEN.  Square-s slots still carry Q_NL.
  u and c(p) unnamed.  phi_F not imported.

============================================================================
Backend: field tables only.  Writes
evidence/e1_gmin_m4_prop15429.json
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
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _e_tr, _kernel_field  # noqa: E402
from e1_gmin_m4_prop15290 import omega_set  # noqa: E402
from e1_gmin_m4_prop15298 import half_gauss, mu_minus, mu_plus  # noqa: E402
from e1_gmin_m4_prop15356 import Q_1d_pp_named  # noqa: E402
from e1_gmin_m4_prop15409 import D_form_on_lattice_general  # noqa: E402
from e1_gmin_m4_prop15411 import LIVE_QNL  # noqa: E402
from e1_gmin_m4_prop15412 import Q_NL_is_short_rational  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15429_catalog.json"

CERT = (3, 5, 7, 11, 13, 17, 19, 23, 29, 31)


def C_ns_named(p: int) -> Fraction:
    q = p * p
    return Fraction(-(q - 1), 4)


def A0_named(p: int) -> Fraction:
    """k(k+p)=p²(q−1)/4 from G_μ(1)+G_μ(r)=p."""
    q = p * p
    return Fraction(p * p * (q - 1), 4)


def C_ns_zero_wrong(_p: int) -> Fraction:
    return Fraction(0)


def Gchi_plus(F) -> float:
    q = F["q"]
    xi = omega_set(F)[0]
    return float(
        sum(F["chi"][x] * _e_tr(F, F["mul"][xi][x]) for x in range(1, q)).real
    )


def C_ns_live(F, r: int, cd: int) -> Fraction:
    q = F["q"]
    G = Gchi_plus(F)
    acc = 0.0
    for s in range(1, q):
        if F["chi"][s] != -1:
            continue
        beta = F["add"][1][F["mul"][r][s]]
        acc += half_gauss(F, beta, cd, G, q)
    return Fraction(acc).limit_denominator(20)


def base_ns_Q(p: int) -> Fraction:
    """ns+base block of Q/q²."""
    q = p * p
    return Fraction((q - 1) * (4 * p + 5 - p * p), 2 * p * p)


def L_ns_named(p: int) -> Fraction:
    return mu_plus(p) * mu_minus(p)


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in CERT:
        F = _kernel_field(p)
        tgt = C_ns_named(p)
        rs = [r for r in F["squares"] if r not in (1, F["neg1"])]
        sample = rs[:4]
        hits = 0
        n = 0
        for r in sample:
            for cd in (1, -1):
                n += 1
                if C_ns_live(F, r, cd) == tgt:
                    hits += 1
                else:
                    ok = False
        if C_ns_zero_wrong(p) == tgt:
            ok = False
        if C_ns_zero_wrong(5) == C_ns_named(5):
            ok = False
        rows[str(p)] = {
            "C_ns": str(tgt),
            "n_r": len(rs),
            "sample_hits": f"{hits}/{n}",
        }
    if C_ns_named(5) != Fraction(-6):
        ok = False
    if C_ns_named(7) != Fraction(-12):
        ok = False
    return {
        "proved": bool(ok),
        "by_p": rows,
        "zero_wrong": "0",
        "theorem": (
            "C_ns,χ=−(q−1)/4 on T\\{±1}, both χ, p=3..31. "
            "Fail: C_ns=0."
        ),
    }


def prove_B() -> dict:
    b5 = base_ns_Q(5)
    b7 = base_ns_Q(7)
    ok = b5 == 0
    ok = ok and b7 == Fraction(-384, 49)
    ok = ok and L_ns_named(5) == Fraction(25, 2)
    ok = ok and L_ns_named(7) == Fraction(147, 2)
    ok = ok and b5 != LIVE_QNL[5]
    ok = ok and b7 != LIVE_QNL[7]
    if b5 == LIVE_QNL[5] or b7 == LIVE_QNL[7]:
        ok = False
    return {
        "proved": bool(ok),
        "p5": str(b5),
        "p7": str(b7),
        "QNL5": str(LIVE_QNL[5]),
        "QNL7": str(LIVE_QNL[7]),
        "theorem": (
            "ns+base=(q−1)(4p+5−p²)/(2p²). 0 at p=5, −384/49 at p=7. "
            "Fail: claim this is Q_NL."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "Q_NL_5": str(LIVE_QNL[5]),
        "Q_1d_5": str(Q_1d_pp_named(5)),
        "Q_NL_gt_Q1d_p5": LIVE_QNL[5] > Q_1d_pp_named(5),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "D_form_on_lattice_general": D_form_on_lattice_general(),
        "Q_NL_is_short_rational": Q_NL_is_short_rational(),
        "note": (
            "C_ns named in p. Square-s slots still carry Q_NL. "
            "u unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.429  C_ns=−(q−1)/4 named", flush=True)
    A = prove_A()
    print(f"  A C_ns: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B ns+base: {B['proved']} p5={B['p5']} p7={B['p7']}", flush=True)
    C = prove_open()
    print(f"  C open: QNL>Q1d={C['Q_NL_gt_Q1d_p5']} phi_F={C['phi_F_ge_6']}", flush=True)
    CAT.write_text(json.dumps({"A": "-(q-1)/4", "B": {"p5": B["p5"], "p7": B["p7"]}}, indent=2) + "\n")
    out = {
        "prop": "15.429",
        "title": "C_ns=−(q−1)/4; ns+base is not Q_NL; u still open",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "C_ns_named": A["proved"],
            "base_ns_not_QNL": B["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15429.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
