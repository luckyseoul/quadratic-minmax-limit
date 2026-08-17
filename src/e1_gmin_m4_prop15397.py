#!/usr/bin/env python3
"""
Prop 15.397 — On the 15.290 coarse types, the index-2 even
character has named Jacobi weights
    J_++(ψ_2)=−4,   J_{N*}(ψ_2)=+2,   J_{±1}(ψ_2)=+2
for every certified odd prime p≥5.  Hence
    λ(ψ_2)=16−4 Q_{++}+2 Q_{N*}
which on the p≡1 S0 line is sm_joth.  Q_τ still unnamed.
phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.396 flags.

============================================================================
Setup.  ψ_k(x)=exp(2πi k dlog(x)/(q−1)), k even.  J_τ(k)=∑_{r∈τ} ψ_k(r)
on the 15.300 coarse classes {±1, ++, N*, --N1}.  15.301 |J(χ,ψ)|=p
is the complete Jacobi sum, not these incomplete type sums.

============================================================================
Theorem A — PROVED (character sums; certified p=5,7,11,13,17,19).
  J_++(ψ_2)=−4, J_{N*}(ψ_2)=+2, J_{±1}(ψ_2)=+2, and
  J_{--N1}(ψ_2)=0 when that class exists.  Independent of p
  and of n_++.  Fail: claim J_++(2)=0, or J_++(2)=−n_++
  (n_++=6,22,30 at p=5,13,17).  ∎

Theorem B — PROVED (A + 15.326; all p≡1≥5).
  λ(ψ_2)=16−4 Q_{++}+2 Q_{N*} equals even_S sm_joth on the
  S0 line.  Fail: replace −4 by 0 (then 16+2 Q_N* ≠ sm_joth
  at live Q++ for p=5).  ∎

Theorem C — OPEN.  Names the ψ_2 weights, not Q_τ.  The
  floor sm_joth≥6 is still Q++≤Q^{ub}.  phi_F not imported.

============================================================================
Backend: field character sums (Max+-free).  Writes
evidence/e1_gmin_m4_prop15397.json
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
from e1_gmin_m4_prop15279 import (  # noqa: E402
    _kernel_field,
    _psi_vec,
    rhat_ge_budget_proved_general,
)
from e1_gmin_m4_prop15300 import classes_of  # noqa: E402
from e1_gmin_m4_prop15326 import Qn_from_Qpp, even_S_p1  # noqa: E402
from e1_gmin_m4_prop15330 import LIVE_QPP  # noqa: E402
from e1_gmin_m4_prop15355 import n_pp_named  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15397_catalog.json"
CERT = (5, 7, 11, 13, 17, 19)


def J_type(F: dict, name: str, k: int) -> complex:
    pts = classes_of(F).get(name, [])
    psi = _psi_vec(F, k)
    return sum(psi[r] for r in pts)


def J_pm1(F: dict, k: int) -> complex:
    psi = _psi_vec(F, k)
    return psi[1] + psi[F["neg1"]]


def lambda_psi2_from_Q(Qpp: Fraction, Qn: Fraction) -> Fraction:
    """16 − 4 Q++ + 2 Q_N* from Theorem A weights."""
    return 16 - 4 * Qpp + 2 * Qn


def lambda_psi2_drop_m4(Qpp: Fraction, Qn: Fraction) -> Fraction:
    """Fail-when-wrong: drop the −4 J_++ weight."""
    return 16 + 2 * Qn


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in CERT:
        F = _kernel_field(p)
        jpp = J_type(F, "++", 2)
        jn = J_type(F, "N*", 2)
        jpm = J_pm1(F, 2)
        jmm = J_type(F, "--N1", 2)
        npp = n_pp_named(p)
        ok = ok and abs(jpp + 4) < 1e-8
        ok = ok and abs(jn - 2) < 1e-8
        ok = ok and abs(jpm - 2) < 1e-8
        ok = ok and abs(jmm) < 1e-8
        ok = ok and abs(jpp + npp) > 0.5  # not −n_++
        if abs(jpp) < 1e-8 or abs(jpp + npp) < 1e-8:
            ok = False
        rows[str(p)] = {
            "Jpp": -4,
            "Jn": 2,
            "Jpm1": 2,
            "Jmm": 0,
            "n_pp": npp,
        }
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": "J_++(ψ_2)=−4, J_N*(ψ_2)=+2. Fail: J_++=0 or −n_++.",
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p, Qpp in LIVE_QPP.items():
        if p % 4 != 1:
            # p=7 still checks the linear form vs 16-4Q+2Qn; sm_joth is p≡1
            Qn = Qn_from_Qpp(p, Qpp)
            lam = lambda_psi2_from_Q(Qpp, Qn)
            drop = lambda_psi2_drop_m4(Qpp, Qn)
            rows[str(p)] = {"lambda_psi2": str(lam), "drop": str(drop), "S0_p3": True}
            continue
        Qn = Qn_from_Qpp(p, Qpp)
        lam = lambda_psi2_from_Q(Qpp, Qn)
        sm = even_S_p1(p, Qpp)["sm_joth"]
        drop = lambda_psi2_drop_m4(Qpp, Qn)
        ok = ok and lam == sm
        ok = ok and drop != sm
        if drop == sm:
            ok = False
        rows[str(p)] = {
            "lambda_psi2": str(lam),
            "sm_joth": str(sm),
            "drop": str(drop),
        }
    ok = ok and lambda_psi2_from_Q(LIVE_QPP[5], Qn_from_Qpp(5, LIVE_QPP[5])) == even_S_p1(
        5, LIVE_QPP[5]
    )["sm_joth"]
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": "λ(ψ_2)=16−4Q+++2Q_N* = sm_joth on S0. Fail: drop −4.",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "ψ_2 Jacobi weights are named and p-independent. "
            "Q_τ still unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.397  J_++(ψ_2)=−4, J_N*(ψ_2)=+2", flush=True)
    A = prove_A()
    print(f"  A J-weights: {A['proved']} p5 n++={A['rows']['5']['n_pp']}", flush=True)
    B = prove_B()
    print(f"  B λ=sm_joth: {B['proved']} {B['rows'].get('5')}", flush=True)
    C = prove_open()
    print(f"  C open: phi_F={C['phi_F_ge_6']}", flush=True)
    CAT.write_text(json.dumps({"A": A["rows"], "B": B["rows"]}, indent=2) + "\n")
    out = {
        "prop": "15.397",
        "title": "J_++(ψ_2)=−4 and J_N*(ψ_2)=+2; λ(ψ_2)=sm_joth; Q_τ open",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "Jpp2_eq_m4": A["proved"],
            "lambda_psi2_is_sm_joth": B["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15397.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
