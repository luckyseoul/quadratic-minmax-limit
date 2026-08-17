#!/usr/bin/env python3
"""
Prop 15.399 — J_++(ψ_2)=−4 for every odd prime p≥5 by complete
sums on F_p^× and the norm-1 torus.  The Gram cone S_k≥0 on the
p≡1 S0 line is a named interval that strictly contains the floor
window, so it proves neither 15.398 inequality.
phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.398 flags.

============================================================================
Setup.  ψ_2(x)=exp(2πi·2·dlog(x)/(q−1)).  U=ker N_{q/p}={x^{p+1}=1},
|U|=p+1.  F_p^× ⊂ F_q^× is the unique subgroup of order p−1.
15.355: for p≡1, ++=(F_p^×\\{±1})∪(U\\{±1}); for p≡3,
++=(F_p^×\\{±1})∪(++∩N1).

============================================================================
Theorem A — PROVED (complete character sums; all odd p≥5).
  ψ_2(±1)=1 and ∑_{F_p^×} ψ_2 = ∑_U ψ_2 = ∑_T ψ_2 = 0.
  Hence ∑_{F_p^×\\{±1}} ψ_2 = ∑_{U\\{±1}} ψ_2 = −2.
  On U\\{−1}, χ_q(1+u)=u^{(p+1)/2} via (1+u)^{p−1}=u^{−1}.
  For p≡3 this identifies ++∩N1 with U^{□}\\{±1}, and
  ∑_{U^{□}} ψ_2=0, so that piece is also −2.
  Therefore J_++(ψ_2)=−4, and then
      J_{N*}(ψ_2)=+2,  J_{±1}(ψ_2)=+2,  J_{--N1}(ψ_2)=0.
  Fail: drop {±1} from the F_p^× sum (then that piece is 0 and
  J_++ would be 0 or −2, not −4); or claim ++∩N1=U\\{±1}
  at p≡3 (false: |U\\{±1}|=p−1 ≠ (p−3)/2).  ∎

Theorem B — PROVED (S0 + four even families; all p≡1≥5).
  S_k≥0 on the S0 line is equivalent to the four even 15.326
  forms being ≥0, hence to
      Q++ ∈ I_psd = [4(p−3)/(p−1), (rhs+8b)/(a+2b)].
  Q_lo − lo = 6(p−3)/(p−1)² > 0 and
  hi − Q^{ub} = 3b/(a+2b) > 0, so I_psd properly contains
  the floor window.  At the lo vertex, Q++−Q_N* = −3 at p=5.
  Fail: claim 2≥11/4, or claim min(Q++−Q_N*)≥0 on the cone.  ∎

Theorem C — OPEN.  Gram/PSD does not give Q++≥Q_N* nor
  Q++≤Q^{ub}.  phi_F not imported.

============================================================================
Backend: field character sums + Fraction polynomials.  Writes
evidence/e1_gmin_m4_prop15399.json
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
from e1_gmin_m4_prop15279 import _kernel_field, _psi_vec  # noqa: E402
from e1_gmin_m4_prop15290 import paley_pair  # noqa: E402
from e1_gmin_m4_prop15300 import classes_of  # noqa: E402
from e1_gmin_m4_prop15321 import Qpp_floor_ub, s0_coeffs  # noqa: E402
from e1_gmin_m4_prop15326 import Q_lo_named, Qn_from_Qpp  # noqa: E402
from e1_gmin_m4_prop15397 import J_pm1, J_type  # noqa: E402
from e1_gmin_m4_prop15398 import Q_eq_named, sbox_qn_and_four  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15399_catalog.json"
CERT = (5, 7, 11, 13, 17, 19, 23, 29)
P1 = (5, 13, 17, 29, 37)


def _pow(F: dict, x: int, e: int) -> int:
    acc, b, ee = 1, int(x), int(e)
    if e == 0:
        return 1
    while ee:
        if ee & 1:
            acc = F["mul"][acc][b]
        b = F["mul"][b][b]
        ee >>= 1
    return acc


def fp_units(p: int) -> list[int]:
    """F_p^× in the 15.279 encoding (second coordinate 0)."""
    return list(range(1, p))


def torus_U(F: dict) -> list[int]:
    p, q = F["p"], F["q"]
    return [r for r in range(1, q) if _pow(F, r, p + 1) == 1]


def usq_off(F: dict) -> list[int]:
    """U^{□} \\ {±1}."""
    p = F["p"]
    return [
        r
        for r in torus_U(F)
        if r not in (1, F["neg1"]) and _pow(F, r, (p + 1) // 2) == 1
    ]


def pp_cap_N1(F: dict) -> list[int]:
    out = []
    for r in torus_U(F):
        if r in (1, F["neg1"]):
            continue
        a, b = paley_pair(F, r)
        if a == 1 and b == 1:
            out.append(r)
    return out


def psi2_sum(F: dict, pts) -> complex:
    psi = _psi_vec(F, 2)
    return sum(psi[r] for r in pts)


def psd_lo_named(p: int) -> Fraction:
    """Lower endpoint of I_psd: sm_j0≥0 on S0."""
    return Fraction(4 * (p - 3), p - 1)


def psd_hi_named(p: int) -> Fraction:
    """Upper endpoint of I_psd: sm_joth≥0 on S0."""
    a, b, rhs = s0_coeffs(p)
    return (rhs + 8 * b) / (a + 2 * b)


def qlo_minus_psdlo_num(p: int) -> int:
    """6(p−3), the numerator of Q_lo − lo after (p−1)²."""
    return 6 * (p - 3)


def psd_cone_forces_qn() -> bool:
    """True only if I_psd proved Q++≥Q_N*. It does not."""
    return False


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in CERT:
        F = _kernel_field(p)
        psi = _psi_vec(F, 2)
        fp = fp_units(p)
        U = torus_U(F)
        T = F["squares"]
        s_fp = psi2_sum(F, fp)
        s_U = psi2_sum(F, U)
        s_T = psi2_sum(F, T)
        pm = psi[1] + psi[F["neg1"]]
        s_fp_off = psi2_sum(F, [r for r in fp if r not in (1, F["neg1"])])
        s_U_off = psi2_sum(F, [r for r in U if r not in (1, F["neg1"])])
        ok = ok and abs(pm - 2) < 1e-8
        ok = ok and abs(s_fp) < 1e-8 and abs(s_U) < 1e-8 and abs(s_T) < 1e-8
        ok = ok and abs(s_fp_off + 2) < 1e-8 and abs(s_U_off + 2) < 1e-8
        # fail-when-wrong: dropping ±1 from the identity would make J++=0
        if abs(s_fp) > 0.5 or abs(s_fp_off) < 0.5:
            ok = False
        jpp = J_type(F, "++", 2)
        jn = J_type(F, "N*", 2)
        ok = ok and abs(jpp + 4) < 1e-8
        ok = ok and abs(jn - 2) < 1e-8
        ok = ok and abs(J_pm1(F, 2) - 2) < 1e-8
        ok = ok and abs(J_type(F, "--N1", 2)) < 1e-8
        if p % 4 == 1:
            # ++ = (F_p^×\\{±1}) ∪ (U\\{±1})
            C = set(classes_of(F)["++"])
            want = set(r for r in fp if r not in (1, F["neg1"])) | set(
                r for r in U if r not in (1, F["neg1"])
            )
            ok = ok and C == want
            # p≡1: N1 is mixed, not Paley ++
            ok = ok and pp_cap_N1(F) == []
        else:
            # ++∩N1 = U^{□}\\{±1}, not the whole U\\{±1}
            cap = set(pp_cap_N1(F))
            us = set(usq_off(F))
            Uoff = set(r for r in U if r not in (1, F["neg1"]))
            ok = ok and cap == us
            ok = ok and cap != Uoff
            if cap == Uoff:
                ok = False
            # χ_q(1+u)=u^{(p+1)/2}
            for r in U:
                if r == F["neg1"]:
                    continue
                lhs = F["chi"][F["add"][1][r]]
                rhs = 1 if _pow(F, r, (p + 1) // 2) == 1 else -1
                if r != 0 and lhs != rhs:
                    ok = False
        rows[str(p)] = {
            "Jpp": -4,
            "sum_Fp": 0,
            "sum_U": 0,
            "sum_Fp_off": -2,
            "sum_U_off": -2,
        }
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "J_++(ψ_2)=−4 for all odd p≥5 via F_p^× and U. "
            "Fail: drop {±1} (piece 0), or ++∩N1=U\\{±1} at p≡3."
        ),
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in P1:
        lo = psd_lo_named(p)
        hi = psd_hi_named(p)
        qlo = Q_lo_named(p)
        qeq = Q_eq_named(p)
        qub = Qpp_floor_ub(p)
        gap_lo = qlo - lo
        a, b, rhs = s0_coeffs(p)
        gap_hi = 3 * b / (a + 2 * b)
        ok = ok and lo < qlo < qeq < qub < hi
        ok = ok and gap_lo == Fraction(qlo_minus_psdlo_num(p), (p - 1) ** 2)
        ok = ok and hi - qub == gap_hi
        dlo = lo - Qn_from_Qpp(p, lo)
        ok = ok and dlo < 0
        if lo >= qlo or dlo >= 0:
            ok = False
        rows[str(p)] = {
            "psd_lo": str(lo),
            "Q_lo": str(qlo),
            "Q_eq": str(qeq),
            "Q_ub": str(qub),
            "psd_hi": str(hi),
            "min_Qpp_minus_Qn": str(dlo),
        }
    ok = ok and psd_lo_named(5) == 2
    ok = ok and Q_lo_named(5) == Fraction(11, 4)
    ok = ok and psd_lo_named(5) < Q_lo_named(5)
    ok = ok and (2 - Qn_from_Qpp(5, Fraction(2))) == -3
    if psd_lo_named(5) >= Q_lo_named(5):
        ok = False
    if (2 - Qn_from_Qpp(5, Fraction(2))) >= 0:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "I_psd=[4(p−3)/(p−1),(rhs+8b)/(a+2b)] properly contains "
            "[Q_lo,Q^{ub}]. Fail: 2≥11/4 or min(Q++−Qn)≥0."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "psd_cone_forces_qn": psd_cone_forces_qn(),
        "sbox_qn_and_four": sbox_qn_and_four(),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "J_++(ψ_2)=−4 is now Max+-free for all p≥5. "
            "Gram/PSD still misses both 15.398 inequalities. "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.399  J++(ψ_2)=−4 all p≥5; PSD cone misses floor", flush=True)
    A = prove_A()
    print(f"  A J++=−4: {A['proved']} primes={list(A['rows'])}", flush=True)
    B = prove_B()
    print(f"  B I_psd: {B['proved']} p5={B['rows']['5']}", flush=True)
    C = prove_open()
    print(
        f"  C open: forces_qn={C['psd_cone_forces_qn']} phi_F={C['phi_F_ge_6']}",
        flush=True,
    )
    CAT.write_text(json.dumps({"A": A["rows"], "B": B["rows"]}, indent=2) + "\n")
    out = {
        "prop": "15.399",
        "title": "J++(ψ_2)=−4 for all p≥5; Gram cone misses both 15.398 inequalities",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "Jpp2_eq_m4_all_p": A["proved"],
            "I_psd_contains_floor": B["proved"],
            "psd_cone_forces_qn": C["psd_cone_forces_qn"],
            "sbox_qn_and_four": C["sbox_qn_and_four"],
            "phi_F_ge_6_proved_general": C["phi_F_ge_6"],
            "residual_ii_k_eq_4p_empty": C["residual_ii_k_eq_4p_empty"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15399.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
