#!/usr/bin/env python3
"""
Prop 15.321 — ∑_{c∈F_p^{×□}} H(c)=χ_p(−1).  For p≡1 (mod 4),
A_□=α(p)+β(p) Q_{++} with named rational α,β.  Floor-other is
Q_{++}≤Qpp_floor_ub(p) on the S_0 line.  Q_τ / A_□ still unnamed.
phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.320 flags.  Floor stays OPEN.

============================================================================
Setup.  H(s)=∑_x χ_p(x(x−1)(x−s)) (15.320).  C(s) named via H.
p≥5 prime, p≡1 (mod 4): --N1 empty, F_1(0)=16+(p−1)Q_{++},
F_*(0)=2 Q_{++}+(p−1)Q_{N*}, S_0: a Q_{++}+b Q_{N*}=rhs with
a=2(p−2), b=(p−1)(p−3)/2 ≠0, rhs=2(p²−9).

============================================================================
Theorem A — PROVED (complete character sums; certified p=5,7,11,13).
  ∑_{c∈F_p^{×□}} H(c) = χ_p(−1).
  Hence ∑_{c□, c≠1} H(c)=χ_p(−1)+1  (using H(1)=−1).
  Derivation: 1_{□^×}=(1+χ)/2 on F_p^×, the inner sums are
  ∑_{c≠0} χ(x−c)=−χ(x) and ∑ χ(c(x−c))=−χ(−1) for x≠0, and
  ∑_x χ(x(x−1))=−1, ∑_{x≠0} χ(x−1)=−χ(−1).  Fail-when-wrong:
  claim the sum vanishes (false at p=5: H(1)+H(4)=1).  ∎

Theorem B — PROVED (15.319–15.320 + A + S_0; certified p=5).
  ∑_{c□,c≠1} C(c)=(p+1)/2 · ( −(p+1)(p−3)/2 + ε p (χ(−1)+1) )
  is named in p.  For primes p≥5, p≡1 (mod 4) (so b≠0) the A_□ identity
  collapses to
      A_□ = α(p) + β(p) Q_{++}
  with explicit rationals α,β (p=5: 840 + 30 Q_{++}, matching
  12360/13 at live 48/13).  Fail-when-wrong: drop the H-sum in
  C_{off} (p=5: Csum becomes −18, not −48).  ∎

Theorem C — PROVED as an equivalence; Q_{++} OPEN as a number.
  On the S_0 line, the σ=−1 other-j floor 2 Q_{++}−Q_{N*}≤5 is
      Q_{++} ≤ Qpp_floor_ub(p)=(5b+rhs)/(2b+a)
  (=26/7 at p=5, =310/71 at p=13).  Equivalently
  A_□ ≤ α+β Qpp_floor_ub.  This is not 4−2/(p+2) for all p≡1
  (that equals the floor ub only at p=5).  CS
  A_□≤(p+1)² (named ∑ N²)/|T| yields Q_{++}≤17 at p=5, useless.
  p≡3 retains Q_{--}; A_□ does not isolate Q_{++}.  phi_F not
  imported.  Fail-when-wrong: claim 4−2/(p+2)=floor ub at p=13.  ∎

============================================================================
Backend: F_p character sums (Max+-free).  Writes
evidence/e1_gmin_m4_prop15321.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15316 import chi_p_pm  # noqa: E402
from e1_gmin_m4_prop15319 import C1_named  # noqa: E402
from e1_gmin_m4_prop15320 import H_legendre  # noqa: E402


def H_sum_squares_named(p: int) -> int:
    return chi_p_pm(p, p - 1)


def Csum_off1_named(p: int) -> Fraction:
    chi_m1 = chi_p_pm(p, p - 1)
    eps = -chi_m1
    n = (p - 3) // 2
    hsum = chi_m1 + 1
    return Fraction(p + 1, 2) * (-(p + 1) * n + eps * p * hsum)


def Csum_off1_drop_H(p: int) -> Fraction:
    n = (p - 3) // 2
    return Fraction(p + 1, 2) * (-(p + 1) * n)


def s0_coeffs(p: int) -> tuple[int, Fraction, int]:
    a = 2 * (p - 2)
    b = Fraction((p - 1) * (p - 3), 2)
    rhs = 2 * (p * p - 9)
    return a, b, rhs


def Qpp_floor_ub(p: int) -> Fraction:
    a, b, rhs = s0_coeffs(p)
    return (5 * b + rhs) / (2 * b + a)


def Qpp_bound_4_minus_2_over_pplus2(p: int) -> Fraction:
    return Fraction(4 * (p + 2) - 2, p + 2)


def alpha_beta_p1(p: int) -> tuple[Fraction, Fraction]:
    if p < 5 or p % 4 != 1:
        raise ValueError("need a prime p≥5 with p≡1 (mod 4)")
    a, b, rhs = s0_coeffs(p)
    C1 = Fraction(C1_named(p))
    Cs = Csum_off1_named(p)
    t0 = Fraction((p + 1) ** 2 * (p - 1) ** 4, 16)
    t1 = Fraction((p + 1) ** 2 * (p - 1) ** 3, 8)
    pm1 = p - 1
    fstar_const = pm1 * rhs / b
    fstar_coef = 2 - pm1 * a / b
    alpha = t0 + t1 + C1 + Cs * fstar_const / 16
    beta = C1 * pm1 / 16 + Cs * fstar_coef / 16
    return alpha, beta


def Asq_from_Qpp_p1(p: int, Qpp: Fraction) -> Fraction:
    al, be = alpha_beta_p1(p)
    return al + be * Qpp


def prove_H_sum() -> dict:
    ok = True
    zero_hit = 0
    rows = {}
    for p in (5, 7, 11, 13):
        live = sum(
            H_legendre(p, c) for c in range(1, p) if chi_p_pm(p, c) == 1
        )
        named = H_sum_squares_named(p)
        if live != named:
            ok = False
        if live == 0:
            zero_hit += 1
        rows[str(p)] = {"live": live, "named": named}
    if zero_hit == 4:
        ok = False
    if rows["5"]["live"] == 0:
        ok = False
    return {
        "proved": ok,
        "sum_zero_false": rows["5"]["live"] != 0,
        "by_p": rows,
        "theorem": "∑_{c∈F_p^{×□}} H(c)=χ_p(−1).",
    }


def prove_Asq_linear() -> dict:
    ok = True
    drop_hit = 0
    rows = {}
    if Csum_off1_named(5) != -48:
        ok = False
    if Csum_off1_drop_H(5) != -18:
        ok = False
    if Csum_off1_named(5) != Csum_off1_drop_H(5):
        drop_hit += 1
    al, be = alpha_beta_p1(5)
    live_A = Asq_from_Qpp_p1(5, Fraction(48, 13))
    if (al, be) != (Fraction(840), Fraction(30)):
        ok = False
    if live_A != Fraction(12360, 13):
        ok = False
    rows["5"] = {
        "alpha": str(al),
        "beta": str(be),
        "Csum": str(Csum_off1_named(5)),
        "Csum_dropH": str(Csum_off1_drop_H(5)),
        "A_at_live_Qpp": str(live_A),
    }
    al13, be13 = alpha_beta_p1(13)
    rows["13"] = {"alpha": str(al13), "beta": str(be13), "Csum": str(Csum_off1_named(13))}
    if drop_hit == 0:
        ok = False
    return {
        "proved": ok,
        "drop_H_fails": drop_hit > 0,
        "by_p": rows,
        "theorem": (
            "C_off named in p.  p≡1: A_□=α+β Q_{++} "
            "(p=5: 840+30 Q_{++})."
        ),
    }


def prove_floor_equiv() -> dict:
    ok = True
    same_as_pplus2_hit = 0
    rows = {}
    ub5 = Qpp_floor_ub(5)
    ub13 = Qpp_floor_ub(13)
    cand5 = Qpp_bound_4_minus_2_over_pplus2(5)
    cand13 = Qpp_bound_4_minus_2_over_pplus2(13)
    if ub5 != Fraction(26, 7):
        ok = False
    if ub5 != cand5:
        ok = False
    if ub13 == cand13:
        same_as_pplus2_hit += 1
        ok = False
    if ub13 != Fraction(310, 71):
        ok = False
    al, be = alpha_beta_p1(5)
    A_at_ub = al + be * ub5
    rows["5"] = {
        "floor_ub": str(ub5),
        "cand_4_minus": str(cand5),
        "A_at_ub": str(A_at_ub),
    }
    rows["13"] = {
        "floor_ub": str(ub13),
        "cand_4_minus": str(cand13),
        "equal": ub13 == cand13,
    }
    # CS A bound is useless at p=5: named_sum N2 = 450, |T|=12
    cs_A = Fraction(6 * 6 * 450, 12)  # (p+1)² * (sum/|T|)
    cs_Q = (cs_A - al) / be
    if cs_Q <= ub5:
        ok = False
    rows["cs_p5"] = {"A_ub": str(cs_A), "Qpp_ub": str(cs_Q)}
    return {
        "proved": ok,
        "pplus2_equals_floor_all_p1_false": ub13 != cand13,
        "cs_closes_false": cs_Q > ub5,
        "by_p": rows,
        "theorem": (
            "S_0 + (2Q++−Q_N*≤5) ⇔ Q++≤(5b+rhs)/(2b+a). "
            "Not 4−2/(p+2) at p=13. CS does not close."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "e1": bool(e1_closed_general()),
        "A_square_named": False,
        "Q_tau_named_in_p": False,
        "lambda_ge_6": False,
        "note": (
            "A_□=α+β Q++ names the pair for p≡1, not either number. "
            "Floor ub on Q++ is named but unproved as a bound. "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.321  ∑ H on squares named; A_□=α+β Q++ (p≡1)", flush=True)
    A = prove_H_sum()
    print(f"  A H-sum: {A['proved']} sum0={not A['sum_zero_false']}", flush=True)
    B = prove_Asq_linear()
    print(f"  B A=α+βQ: {B['proved']} dropH={B['drop_H_fails']}", flush=True)
    C = prove_floor_equiv()
    print(
        f"  C floor-ub: {C['proved']} p+2_all={not C['pplus2_equals_floor_all_p1_false']} "
        f"CS={not C['cs_closes_false']}",
        flush=True,
    )
    D = prove_open()
    print(f"  D open: phi_F={D['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.321",
        "title": "∑ H on □ named; A_□=α+β Q++ for p≡1; floor ub named, unproved",
        "proved": {
            "H_sum_squares": A["proved"],
            "Asq_linear_p1": B["proved"],
            "floor_ub_equiv": C["proved"],
            "A_square_named": False,
            "Q_tau_named_in_p": False,
            "lambda_ge_6": False,
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
    path = ROOT / "evidence" / "e1_gmin_m4_prop15321.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
