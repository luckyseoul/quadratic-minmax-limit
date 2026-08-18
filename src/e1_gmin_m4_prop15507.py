#!/usr/bin/env python3
"""
Prop 15.507 — For p≡1 the even pairing is a named function of one
mass M=n_pp(4−Q++/q²).  The worst branch J_{N*}=2 satisfies
    ⟨δ,ψ⟩≤2  ⇔  Q++/q² ≤ Qpp_floor_ub(p)
(26/7 at p=5).  Fail: 16/|N*|≤2 at p=5 (that is 4).  Does not
prove the inequality Max+-free.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.506 flags.

============================================================================
Setup.  p≡1≥5: --N1 empty (15.355).  T\\{±1}=++ ⊔ N*,
n_pp=2(p−2), |N*|=(p−1)(p−3)/2.  S0 (15.506): n_pp δ++ + |N*| δ_N*=8.
Even ψ∉{1,χ}: ∑_T ψ=0 ⇒ J_++ + J_{N*}=−2.  J_{N*} named (15.497).
M=n_pp δ++ , m=8−M.

============================================================================
Theorem A — PROVED (trace on T; certified p=5,13,17).
  J_++(ψ_k)=−2−J_{N*}(k) for even k∉{0,(q−1)/2}.  Fail: flip to
  −2+J_{N*} (misses live J_++=−4 on the generic k at p=5).  ∎

Theorem B — PROVED (A + S0).
  ⟨δ,ψ⟩ = J_++ M/n_pp + J_{N*} m/|N*|.  At p=5, J_{N*}=2, live
  M=24/13 recovers 24/13.  Fail: the flipped J_++ of A.  ∎

Theorem C — PROVED (B, J_{N*}=2 branch; p≡1).
  ⟨δ,ψ⟩=16/|N*| − M(4/n_pp+2/|N*|).  This is ≤2 for all M≥0
  iff 16/|N*|≤2 iff p≥13.  At p=5 one has 16/4=4>2, so the
  automatic bound fails.  Fail: claim 16/|N*|≤2 at p=5.  ∎

Theorem D — PROVED (C + live Q at p=5).
  At p=5, J_{N*}=2 pairing ≤2 ⇔ M≥12/7 ⇔ Q++/q²≤26/7
  =Qpp_floor_ub(5).  Live 48/13<26/7.  Fail: claim the threshold
  is 4−2/(p−2)=10/3, or claim live Q++ equals 26/7.  ∎

Theorem E — OPEN.  Q++≤Qpp_floor_ub is not proved Max+-free.
  Two-type Q-constancy for p≥13 is 15.298 B (census p≤7).
  δ≥0 on off is not proved.  phi_F not imported.  ∎

============================================================================
Backend: Fraction algebra + field character sums (ProcessPool).
Writes evidence/e1_gmin_m4_prop15507.json
"""
from __future__ import annotations

import json
import os
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import _kernel_field, _psi_vec
from e1_gmin_m4_prop15298 import coarse_Q_class
from e1_gmin_m4_prop15321 import Qpp_floor_ub
from e1_gmin_m4_prop15330 import LIVE_QPP
from e1_gmin_m4_prop15355 import n_pp_named
from e1_gmin_m4_prop15497 import J_Nstar_named, even_ks

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "evidence" / "e1_gmin_m4_prop15507.json"

CERT = (5, 13, 17)
W = 86


def n_nstar_p1(p: int) -> int:
    return (p - 1) * (p - 3) // 2


def J_pp_named(p: int, k: int) -> int:
    return -2 - J_Nstar_named(p, k)


def J_pp_flipped(p: int, k: int) -> int:
    """Fail-when-wrong: flip the trace sign."""
    return -2 + J_Nstar_named(p, k)


def M_from_Qpp(p: int, qpp: Fraction) -> Fraction:
    return n_pp_named(p) * (4 - qpp)


def live_M(p: int) -> Fraction:
    return M_from_Qpp(p, LIVE_QPP[p])


def pair_from_M(p: int, jn: int, M: Fraction) -> Fraction:
    npp = n_pp_named(p)
    nns = n_nstar_p1(p)
    jpp = -2 - jn
    return jpp * M / npp + jn * (8 - M) / nns


def pair_J2(p: int, M: Fraction) -> Fraction:
    """J_{N*}=2 branch."""
    return pair_from_M(p, 2, M)


def auto_J2_ub(p: int) -> Fraction:
    """16/|N*|. ≤2 iff p≥13."""
    return Fraction(16, n_nstar_p1(p))


def M_need_J2(p: int) -> Fraction:
    """pair_J2≤2 ⇔ M≥ this (the coeff of M is negative)."""
    npp = n_pp_named(p)
    nns = n_nstar_p1(p)
    # 16/nns - M(4/npp+2/nns) ≤ 2
    # 16/nns - 2 ≤ M(4/npp+2/nns)
    lhs = auto_J2_ub(p) - 2
    coeff = Fraction(4, npp) + Fraction(2, nns)
    if coeff == 0:
        return Fraction(0)
    return lhs / coeff


def Qpp_ub_J2(p: int) -> Fraction:
    """Q++/q² ≤ this ⇔ pair_J2≤2 (when coeff>0)."""
    return 4 - M_need_J2(p) / n_pp_named(p)


def Qpp_lo_wrong(p: int) -> Fraction:
    """Fail-when-wrong: 4−2/(p−2), not the J=2 threshold."""
    return 4 - Fraction(2, p - 2)


def pp_set(F) -> list[int]:
    """Merged ++ at p≡1: T\\{±1}\\N* (sub ∪ all N1)."""
    out = []
    for r in F["squares"]:
        if r in (1, F["neg1"]):
            continue
        if coarse_Q_class(F, int(r)) != "N*":
            out.append(int(r))
    return out


def J_pp_live(F, k: int, pp: list[int]) -> float:
    psi = _psi_vec(F, k)
    return float(sum(psi[r].real for r in pp))


def _jpp_job(p: int) -> dict:
    F = _kernel_field(p)
    pp = pp_set(F)
    ok = len(pp) == n_pp_named(p)
    flip_miss = False
    for k in even_ks(p):
        got = J_pp_live(F, k, pp)
        named = J_pp_named(p, k)
        if abs(got - named) > 1e-8:
            ok = False
        if abs(got - J_pp_flipped(p, k)) > 1e-6:
            flip_miss = True
    if not flip_miss:
        ok = False
    return {
        "p": p,
        "n_pp": len(pp),
        "named_n_pp": n_pp_named(p),
        "ok": ok,
        "generic_Jpp": J_pp_named(p, 2),
    }


def prove_A() -> dict:
    with ProcessPoolExecutor(max_workers=min(W, len(CERT))) as ex:
        rows = list(ex.map(_jpp_job, CERT))
    ok = all(r["ok"] for r in rows)
    if J_pp_named(5, 2) != -4:
        ok = False
    if J_pp_flipped(5, 2) == -4:
        ok = False
    return {
        "proved": bool(ok),
        "rows": {str(r["p"]): r for r in rows},
        "theorem": "J_++=−2−J_{N*}. Fail flip sign.",
    }


def prove_B() -> dict:
    M = live_M(5)
    got = pair_from_M(5, 2, M)
    bad = J_pp_flipped(5, 2) * M / n_pp_named(5) + 2 * (8 - M) / n_nstar_p1(5)
    ok = got == Fraction(24, 13) and M == Fraction(24, 13)
    if bad == Fraction(24, 13):
        ok = False
    return {
        "proved": bool(ok),
        "M5": str(M),
        "pair": str(got),
        "flipped": str(bad),
        "theorem": "p=5 J=2 live M recovers 24/13. Fail flipped J_++.",
    }


def prove_C() -> dict:
    ok = auto_J2_ub(5) == 4 and auto_J2_ub(5) > 2
    ok = ok and auto_J2_ub(13) <= 2
    if auto_J2_ub(5) <= 2:
        ok = False
    return {
        "proved": bool(ok),
        "ub5": str(auto_J2_ub(5)),
        "ub13": str(auto_J2_ub(13)),
        "theorem": "16/|N*|≤2 iff p≥13. Fail at p=5 (4>2).",
    }


def prove_D() -> dict:
    qub = Qpp_ub_J2(5)
    live = LIVE_QPP[5]
    ok = qub == Fraction(26, 7) and qub == Qpp_floor_ub(5)
    ok = ok and live < qub
    ok = ok and live != qub
    if Qpp_lo_wrong(5) == qub:
        ok = False
    if live == Qpp_lo_wrong(5):
        ok = False
    # pairing at live M
    if pair_J2(5, live_M(5)) != Fraction(24, 13):
        ok = False
    if pair_J2(5, Fraction(12, 7)) != 2:
        ok = False
    return {
        "proved": bool(ok),
        "Qpp_ub_J2": str(qub),
        "floor_ub": str(Qpp_floor_ub(5)),
        "wrong_lo": str(Qpp_lo_wrong(5)),
        "live": str(live),
        "theorem": "p=5 J=2 floor ⇔ Q++≤26/7=Qpp_floor_ub. Fail 10/3.",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "Q_tau_named_in_p": False,
        "phi_F_imported": False,
        "note": (
            "Pairing at p≡1 is named in one mass. The J=2 bound is "
            "Q++≤Qpp_floor_ub, not proved Max+-free. phi_F not imported."
        ),
    }


def main() -> dict:
    print("15.507 p≡1 pairing is one mass; J=2 iff Q++≤Qpp_floor_ub", flush=True)
    A, B, C, D, E = prove_A(), prove_B(), prove_C(), prove_D(), prove_open()
    out = {
        "prop": "15.507",
        "title": "p≡1 pairing is one mass; J=2 ⇔ Q++≤Qpp_floor_ub",
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
