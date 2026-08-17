#!/usr/bin/env python3
"""
Prop 15.396 — Dual form of the floor: λ(ψ)=λ_mean+⟨ΔQ,J_ψ⟩
with named λ_mean=8(q−1)/(q−5)>6.  The Q≥0 and Q∈[0,4]
rearrangement duals miss the threshold −2(q+11)/(q−5).
phi_F is imported only if a dual closed; it does not.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.395 flags.  Not a 15.395 occupancy-CS bound.

============================================================================
Setup.  15.290: λ=S_□/q², Q(±1)=8q², even ψ(±1)=1.
S0=2(q−1)=∑_T (Q/q²) gives the off-average
    Q_avg=4(q−9)/(q−5)   (15.300 S0 + Wick ±1).
∑_{T\\{±1}} ψ=−2.

============================================================================
Theorem A — PROVED (S0 + Wick ±1; all odd p≥5).
  Q_avg=4(q−9)/(q−5).  Fail: claim the Wick neighbour
  4(q−1)/(q+1) (equals 48/13 at p=5 but is not S0).  ∎

Theorem B — PROVED (character split; all odd p≥5).
  λ(ψ)=8(q−1)/(q−5) + ⟨ΔQ,J_ψ⟩  for even ψ with ψ(±1)=1,
  ΔQ=Q/q²−Q_avg on T\\{±1}.  The mean mode
  λ_mean=8(q−1)/(q−5)>6.  Fail: claim λ_mean=6 at p=5
  (48/5≠6), or drop λ_mean.  ∎

Theorem C — PROVED (rearrangement dual; all odd p≥5).
  Floor ⇔ ⟨ΔQ,J⟩≥−2(q+11)/(q−5).  Q≥0 bang-bang gives
  only ≥−8(q−9)/(q−5), and 8(q−9)>2(q+11) (p=5: 128>72).
  Q∈[0,4] bang-bang with two ψ=1 off-points misses by 6.
  Fail: claim 8(q−9)≤2(q+11) at p=5.  ∎

Theorem D — OPEN.  No dual certificate yet forces ⟨ΔQ,J⟩
  above the threshold.  dual_ge_6() is False; phi_F not imported.

============================================================================
Backend: Fraction identities.  Writes
evidence/e1_gmin_m4_prop15396.json
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
from e1_gmin_m4_prop15299 import Q_pp_tempting  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15396_catalog.json"
P_ODD = (5, 7, 11, 13, 17, 19, 29, 37)


def q_of(p: int) -> int:
    return p * p


def Q_avg_S0(p: int) -> Fraction:
    """Off-average of Q/q² from S0=2(q−1) and Wick Q(±1)/q²=8."""
    q = q_of(p)
    return Fraction(4 * (q - 9), q - 5)


def Q_avg_wick_neigh(p: int) -> Fraction:
    """15.365 neighbour 4(q−1)/(q+1), not S0."""
    return Q_pp_tempting(p)


def lambda_mean(p: int) -> Fraction:
    """16 − 2 Q_avg = 8(q−1)/(q−5)."""
    q = q_of(p)
    return Fraction(8 * (q - 1), q - 5)


def thresh_cov(p: int) -> Fraction:
    """Floor ⇔ ⟨ΔQ,J⟩ ≥ this (negative)."""
    q = q_of(p)
    return Fraction(-2 * (q + 11), q - 5)


def dual_Qge0(p: int) -> Fraction:
    """Q≥0 rearrangement lower bound on ⟨ΔQ,J⟩."""
    q = q_of(p)
    return Fraction(-8 * (q - 9), q - 5)


def dual_box04_psi1(p: int) -> Fraction:
    """Q∈[0,4] bang-bang with two off-points at ψ=1."""
    # −8 − 4(1+1) + 2 Q_avg = −16 + 8(q−9)/(q−5)
    q = q_of(p)
    return Fraction(-16) + Fraction(8 * (q - 9), q - 5)


def dual_ge_6() -> bool:
    """True only if a named rearrangement dual meets thresh_cov."""
    for p in P_ODD:
        if dual_Qge0(p) < thresh_cov(p):
            return False
        if dual_box04_psi1(p) < thresh_cov(p):
            return False
    return True


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in P_ODD:
        qa = Q_avg_S0(p)
        qw = Q_avg_wick_neigh(p)
        q = q_of(p)
        # S0 check: 16 + Q_avg * ((q-5)/2) = 2(q-1)
        n_off = Fraction(q - 5, 2)
        s0 = 16 + qa * n_off
        ok = ok and s0 == 2 * (q - 1)
        rows[str(p)] = {"Q_avg": str(qa), "Q_wick": str(qw), "S0": str(s0)}
        if p == 5:
            ok = ok and qa == Fraction(16, 5)
            ok = ok and qw == Fraction(48, 13)
            ok = ok and qa != qw
    if Q_avg_S0(5) == Q_avg_wick_neigh(5):
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": "Q_avg=4(q−9)/(q−5) from S0. Fail: 4(q−1)/(q+1).",
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in P_ODD:
        lm = lambda_mean(p)
        qa = Q_avg_S0(p)
        ok = ok and lm == 16 - 2 * qa
        ok = ok and lm > 6
        rows[str(p)] = {"lambda_mean": str(lm), "thresh": str(thresh_cov(p))}
    ok = ok and lambda_mean(5) == Fraction(48, 5)
    ok = ok and lambda_mean(5) != 6
    if lambda_mean(5) == 6:
        ok = False
    # live p=5 min λ=80/13 < mean, so cov is negative
    ok = ok and Fraction(80, 13) < lambda_mean(5)
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "λ=8(q−1)/(q−5)+⟨ΔQ,J⟩ and the mean mode is >6. "
            "Fail: λ_mean=6 at p=5."
        ),
    }


def prove_C() -> dict:
    ok = True
    rows = {}
    for p in P_ODD:
        th = thresh_cov(p)
        d0 = dual_Qge0(p)
        d4 = dual_box04_psi1(p)
        ok = ok and th < 0
        ok = ok and d0 < th
        ok = ok and d4 < th
        rows[str(p)] = {
            "thresh": str(th),
            "dual_Qge0": str(d0),
            "dual_box04": str(d4),
        }
    # p=5: 8(q-9)=128, 2(q+11)=72
    ok = ok and 8 * (25 - 9) > 2 * (25 + 11)
    if 8 * (25 - 9) <= 2 * (25 + 11):
        ok = False
    # Q∈[0,4] miss is exactly 6: dual - thresh = -6
    ok = ok and dual_box04_psi1(5) - thresh_cov(5) == -6
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "Q≥0 and Q∈[0,4] rearrangement duals miss "
            "−2(q+11)/(q−5). Fail: 8(q−9)≤2(q+11) at p=5."
        ),
    }


def prove_open() -> dict:
    ge6 = dual_ge_6()
    return {
        "proved": False,
        "dual_ge_6": ge6,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "Mean mode is above 6; rearrangement duals miss the cov "
            "threshold. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.396  dual λ=λ_mean+⟨ΔQ,J⟩; rearrangement misses", flush=True)
    A = prove_A()
    print(f"  A Q_avg: {A['proved']} p5={A['rows']['5']['Q_avg']}", flush=True)
    B = prove_B()
    print(f"  B mean: {B['proved']} p5 λ={B['rows']['5']['lambda_mean']}", flush=True)
    C = prove_C()
    print(f"  C dual miss: {C['proved']} p5 th={C['rows']['5']['thresh']}", flush=True)
    D = prove_open()
    print(f"  D open: dual_ge6={D['dual_ge_6']} phi_F={D['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps({"A": A["rows"], "B": B["rows"], "C": C["rows"]}, indent=2) + "\n"
    )
    out = {
        "prop": "15.396",
        "title": "Dual floor λ=λ_mean+⟨ΔQ,J⟩; Q≥0 dual misses; phi_F unimported",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "Q_avg_from_S0": A["proved"],
            "lambda_mean_gt_6": B["proved"],
            "rearrangement_dual_misses": C["proved"],
            "dual_ge_6": D["dual_ge_6"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15396.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
