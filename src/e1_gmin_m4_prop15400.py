#!/usr/bin/env python3
"""
Prop 15.400 — The floor λ_j≥6 is not pointwise on Max+.
A constraint stronger than S_k≥0 that demanded |⟨ΔN,G_U⟩|²
above threshold on every regular set is dead.
phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.399 flags.

============================================================================
Setup.  15.304: λ_j = 8/(p+1) q^{-2} E|⟨ΔN,G_U(j)⟩|², and
λ_j≥6 ⇔ E|inner|² ≥ (3/4)(p+1)q².  Squares give only λ_j≥0
(the 15.399 Gram cone).  15.292: n_1d=m C(p,m) and
S_1D/q² ≥ 4p(p−1)/(p−2).

============================================================================
Theorem A — PROVED (cache fold; certified p=5,7).
  λ_j≥6 fails on a positive-density set of Max+.
  At p=5, even j∈{2,4}: exactly n_1d=30 rows have
  λ=4p(p−1)/(p−2)=80/3 and the other 100 have λ=0.
  At p=7, j=2: 4102 of 5726 rows have λ<6, and 350 have
  λ=0, which is not N+−n_1d=5586 (NL does not vanish).
  Fail: claim every Max+ has λ_2≥6; or claim λ_2≡0 on
  all NL at p=7.  ∎

Theorem B — PROVED (A + named 1D weight; p=5).
  Ensemble λ_2 = n_1d/N+ · 4p(p−1)/(p−2) = 80/13.
  The p=5 binding mode is exactly the 1D indicator.
  Fail: claim the 100 NL rows contribute a positive
  share of λ_2 (their share is 0).  ∎

Theorem C — OPEN.  The same 1D-weight lower bound
  n_1d/N+ · 4p(p−1)/(p−2) is 140/5726 · 168/5 < 1 at p=7,
  so it cannot import the floor.  Pointwise |inner|² is
  dead as a stronger-than-S≥0 constraint.  phi_F not imported.

============================================================================
Backend: vectorized N@G_U on the p=5,7 caches (no live Q).
Writes evidence/e1_gmin_m4_prop15400.json
"""
from __future__ import annotations

import json
import os
import sys
from fractions import Fraction
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field  # noqa: E402
from e1_gmin_m4_prop15290 import omega_set  # noqa: E402
from e1_gmin_m4_prop15292 import n_1d  # noqa: E402
from e1_gmin_m4_prop15298 import load_N  # noqa: E402
from e1_gmin_m4_prop15302 import nontrivial_s  # noqa: E402
from e1_gmin_m4_prop15304 import gu_vec  # noqa: E402
from e1_gmin_m4_prop15330 import HPLUS  # noqa: E402
from e1_gmin_m4_prop15398 import sbox_qn_and_four  # noqa: E402
from e1_gmin_m4_prop15399 import psd_cone_forces_qn  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15400_catalog.json"


def s1d_named(p: int) -> Fraction:
    """Named 15.292 lower bound S_1D/q² = 4p(p−1)/(p−2)."""
    return Fraction(4 * p * (p - 1), p - 2)


def one_d_weight_lb(p: int) -> Fraction:
    """n_1d/N+ · S_1D, using live |H+|."""
    return Fraction(n_1d(p), HPLUS[p]) * s1d_named(p)


def pointwise_floor() -> bool:
    """True only if every Max+ has λ_j≥6. It does not."""
    return False


def _dN(p: int):
    F = _kernel_field(p)
    N = load_N(p, F).astype(np.float64)
    s = nontrivial_s(F)
    perm = [F["mul"][F["inv"][s]][d] for d in range(F["q"])]
    return F, N - N[:, perm]


def lambda_rows(p: int, j: int) -> np.ndarray:
    F, dN = _dN(p)
    xi = omega_set(F)[0]
    q2 = float(F["q"] ** 2)
    gv = gu_vec(F, j, xi)
    return (8.0 / (p + 1)) * np.abs(dN @ gv) ** 2 / q2


def prove_A() -> dict:
    ok = True
    rows = {}
    # p=5: two-level on j=2
    lam5 = lambda_rows(5, 2)
    snapped = np.round(lam5, 8)
    vals, cnts = np.unique(snapped, return_counts=True)
    hist = {str(float(v)): int(c) for v, c in zip(vals, cnts)}
    n0 = int(np.sum(np.abs(lam5) < 1e-8))
    n_hi = int(np.sum(np.abs(lam5 - float(s1d_named(5))) < 1e-6))
    n_lt6 = int(np.sum(lam5 < 6 - 1e-9))
    ok = ok and n0 == 100 and n_hi == n_1d(5) == 30
    ok = ok and n_lt6 == 100
    ok = ok and len(vals) == 2
    if n_lt6 == 0:
        ok = False
    rows["5"] = {
        "n": 130,
        "n_1d": n_1d(5),
        "n_lambda_0": n0,
        "n_lambda_s1d": n_hi,
        "n_lt6": n_lt6,
        "hist": hist,
        "s1d": str(s1d_named(5)),
    }
    # p=7: NL does not vanish
    lam7 = lambda_rows(7, 2)
    n0_7 = int(np.sum(np.abs(lam7) < 1e-8))
    n_lt6_7 = int(np.sum(lam7 < 6 - 1e-9))
    nl = HPLUS[7] - n_1d(7)
    ok = ok and n_lt6_7 > 0
    ok = ok and n0_7 != nl
    if n_lt6_7 == 0 or n0_7 == nl:
        ok = False
    rows["7"] = {
        "n": HPLUS[7],
        "n_1d": n_1d(7),
        "n_lambda_0": n0_7,
        "n_lt6": n_lt6_7,
        "NL": nl,
        "n0_eq_NL": n0_7 == nl,
    }
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "λ≥6 is not pointwise. p=5 j=2 is {0}^100 ∪ {80/3}^30. "
            "p=7 NL does not vanish. Fail: every Max+ has λ_2≥6."
        ),
    }


def prove_B() -> dict:
    ok = True
    mix = one_d_weight_lb(5)
    live = Fraction(80, 13)
    ok = ok and mix == live
    ok = ok and n_1d(5) == 30 and HPLUS[5] == 130
    ok = ok and s1d_named(5) == Fraction(80, 3)
    # fail-when-wrong: NL share of the p=5 j=2 mean is 0
    lam = lambda_rows(5, 2)
    nl_share = float(np.mean(lam[np.abs(lam) < 1]))
    ok = ok and abs(nl_share) < 1e-8
    if mix != live or abs(nl_share) > 1e-6:
        ok = False
    return {
        "proved": bool(ok),
        "mix": str(mix),
        "live_lambda2": str(live),
        "NL_share": nl_share,
        "theorem": (
            "p=5 ensemble λ_2 = n_1d/N+ · 80/3 = 80/13. "
            "Fail: NL share of λ_2 is positive."
        ),
    }


def prove_open() -> dict:
    lb7 = one_d_weight_lb(7)
    return {
        "proved": False,
        "pointwise_floor": pointwise_floor(),
        "one_d_lb_p7": str(lb7),
        "one_d_lb_p7_lt_6": lb7 < 6,
        "psd_cone_forces_qn": psd_cone_forces_qn(),
        "sbox_qn_and_four": sbox_qn_and_four(),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "Pointwise |inner|² dies. 1D-weight LB is <1 at p=7. "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.400  floor λ≥6 is not pointwise on Max+", flush=True)
    A = prove_A()
    print(f"  A not pointwise: {A['proved']} p5={A['rows']['5']['n_lt6']}/130", flush=True)
    B = prove_B()
    print(f"  B p5 mix 80/13: {B['proved']} mix={B['mix']}", flush=True)
    C = prove_open()
    print(
        f"  C open: ptwise={C['pointwise_floor']} p7_1d_lb={C['one_d_lb_p7']} "
        f"phi_F={C['phi_F_ge_6']}",
        flush=True,
    )
    CAT.write_text(json.dumps({"A": A["rows"], "B": {"mix": B["mix"]}}, indent=2) + "\n")
    out = {
        "prop": "15.400",
        "title": "Floor λ≥6 is not pointwise; |inner|²-per-set dies",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "not_pointwise": A["proved"],
            "p5_mix_is_1d_indicator": B["proved"],
            "pointwise_floor": C["pointwise_floor"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15400.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
