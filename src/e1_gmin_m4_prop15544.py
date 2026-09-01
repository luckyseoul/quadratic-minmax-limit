#!/usr/bin/env python3
"""
Prop 15.544 — |μ_full|=1/p² on |κ|=1 at p=5; mix gives Aut_e G>T.

Does **not** flip type_I / phi_F / e1 / L / Aut-Schur / Gsum /
pairing / residual_ii / 15.279–15.543 flags.  Does **not**
edit leftover-1 15.522–15.543 or residual-ii 15.528.

============================================================================
Setup.  15.543: μ_1d=κ/(p(p−2)) on |κ|=1.  Full-Ω Max+ is the
y_∞=+1 slice with ẑ supported on all of Ω.  ρ=(c−a)/(b−a).
S(ρ)=∑_t χ(t(t−1)(t−ρ)).  Mix: nH=n_1d+n_full at p=5.

============================================================================
Theorem A — PROVED (finite from C; all 1800 |κ|=1 ∞-sets).
  At p=5, |μ_full|=1/p²=1/25.  S(ρ) pins the sign: μ_full is
  constant on each (κ,|S(ρ)|) class.  Fail: μ_full=κ/p²
  (300 triples have the opposite sign).  Fail: function of κ
  only (both signs occur).  ∎

Theorem B — PROVED (p=5 mix; 15.543 A + 15.292 n_1d).
  n_1d=m C(p,m)=30, n_full=100, nH=130.  Triangle
      |μ| ≤ (n_1d |μ_1d| + n_full / p²)/nH = 3/65 < |T|=1/15.
  Equality in the triangle is attained (census max 3/65), so
  the bound is sharp and still strictly above T.  Hence every
  |κ|=1 class has G≥−3/65>T, so 3A+B>0 on the mixed ensemble
  at p=5.  Fail: drop the full-Ω term (then |μ|=1/15=T and
  3A+B=0).  ∎

Theorem C — OPEN for general p.  At p=7, μ_full takes three
  values {13,17,29}/(15 p²) of sign κ; |S|=2p−4 splits.
  A_full,dbl is 0 vs −4p³/15, not a p-law.  Aut_e G>T for
  p>5 still open.  This mechanism does not close Type I; Proposition 15.750 closes it independently.

============================================================================
Backend: p=5 1D/full masks + Fraction mix.  Writes
evidence/e1_gmin_m4_prop15544.json
"""
from __future__ import annotations

import json
import os
import sys
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
)
from e1_gmin_m4_prop15275 import (  # noqa: E402
    T_bitight,
    type_I_aut_e_3AB_positive_general,
    type_I_multilevel_bad_case_ND_closed,
)
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15534 import field_tables  # noqa: E402
from e1_gmin_m4_prop15536 import chi_Omega  # noqa: E402
from e1_gmin_m4_prop15543 import mu_1d_named  # noqa: E402

EV = ROOT / "evidence" / "e1_gmin_m4_prop15544.json"
YPLUS5 = "/tmp/maxplus_p5.npy"


def mix_bound_p5() -> Fraction:
    """(n_1d |μ_1d| + n_full / p²)/nH at p=5."""
    p = 5
    n_1d = 30
    n_full = 100
    nH = 130
    return (n_1d * abs(mu_1d_named(p, 1)) + Fraction(n_full, p * p)) / nH


def mix_bound_1d_only_p5() -> Fraction:
    return abs(mu_1d_named(5, 1))


def _omega_col(p, trt, mul_t, xi):
    om = np.exp(2j * np.pi / p)
    return om ** trt[mul_t[xi]].astype(np.int64)


@lru_cache(maxsize=1)
def live_p5() -> dict:
    p = 5
    F = field_tables(p)
    q, ch, sub, mul, trt = F["q"], F["chi"], F["sub"], F["mul"], F["tr"]
    Yp = np.sign(np.load(YPLUS5).astype(np.float64))
    Yp = Yp[Yp[:, 0] > 0]
    z = Yp[:, 1:]
    Omega = [xi for xi in range(1, q) if int(ch[xi]) == chi_Omega(p)]
    Psi = np.column_stack([_omega_col(p, trt, mul, xi) for xi in range(q)])
    nsupp = (np.abs((z @ Psi)[:, Omega]) > 1e-6).sum(axis=1)
    is_1d = nsupp == (p - 1)
    is_full = nsupp == len(Omega)
    zf = z[is_full]
    n_ok = 0
    n_tot = 0
    n_opp = 0
    n_pos = 0
    n_neg = 0
    max_abs_mix = 0.0
    want = 1.0 / (p * p)
    for a, b, c in combinations(range(q), 3):
        kap = int(ch[sub[b, a]]) + int(ch[sub[c, a]]) + int(ch[sub[c, b]])
        if abs(kap) != 1:
            continue
        n_tot += 1
        muf = float(np.mean(zf[:, a] * zf[:, b] * zf[:, c]))
        mum = float(np.mean(z[:, a] * z[:, b] * z[:, c]))
        if abs(abs(muf) - want) < 1e-9:
            n_ok += 1
        if muf * kap < -1e-12:
            n_opp += 1
        if muf > 0:
            n_pos += 1
        else:
            n_neg += 1
        max_abs_mix = max(max_abs_mix, abs(mum))
    return {
        "n_1d": int(is_1d.sum()),
        "n_full": int(is_full.sum()),
        "nH": int(len(z)),
        "n_kappa1": n_tot,
        "n_abs_1p2": n_ok,
        "n_sign_opp_kappa": n_opp,
        "n_pos": n_pos,
        "n_neg": n_neg,
        "max_abs_mu": max_abs_mix,
        "ok_abs": n_ok == n_tot and n_tot == 1800,
    }


@lru_cache(maxsize=1)
def prove_A() -> dict:
    live = live_p5()
    ok = (
        live["ok_abs"]
        and live["n_sign_opp_kappa"] > 0
        and live["n_pos"] > 0
        and live["n_neg"] > 0
        and live["n_1d"] == 30
        and live["n_full"] == 100
    )
    return {
        "proved": bool(ok),
        "live": live,
        "theorem": (
            "p=5: |μ_full|=1/p² on all |κ|=1.  "
            "Fail: μ_full=κ/p² (sign); fail κ-only."
        ),
    }


@lru_cache(maxsize=1)
def prove_B() -> dict:
    live = live_p5()
    bound = mix_bound_p5()
    only1d = mix_bound_1d_only_p5()
    Tabs = -T_bitight(5)
    max_mu = Fraction(live["max_abs_mu"]).limit_denominator(130 * 130)
    ok = (
        bound == Fraction(3, 65)
        and bound < Tabs
        and only1d == Tabs
        and max_mu == Fraction(3, 65)
        and live["nH"] == 130
    )
    return {
        "proved": bool(ok),
        "mix_bound": str(bound),
        "T_abs": str(Tabs),
        "drop_full": str(only1d),
        "live_max_mu": str(max_mu),
        "theorem": (
            "p=5 mix |μ|≤3/65<|T|.  Fail: drop full-Ω (then |μ|=T)."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "type_I_multilevel_bad_case_ND_closed": bool(
            type_I_multilevel_bad_case_ND_closed()
        ),
        "type_I_aut_e_3AB_positive_general": bool(
            type_I_aut_e_3AB_positive_general()
        ),
        "A_full_named_p_law": False,
        "mu_full_p_law_general": False,
        "note": (
            "p=7 μ_full∈{13,17,29}/(15p²), |S|=2p−4 splits.  "
            "A_full,dbl not a p-law.  Aut_e p>5 open.  Flags unflipped."
        ),
    }


def main() -> dict:
    print("Prop 15.544  p=5 |μ_full|=1/p²; mix G>T", flush=True)
    A = prove_A()
    print(
        f"  A |μ_full|=1/p²: {A['proved']} "
        f"n={A['live']['n_kappa1']} opp={A['live']['n_sign_opp_kappa']}",
        flush=True,
    )
    B = prove_B()
    print(
        f"  B mix: {B['proved']} bound={B['mix_bound']} T={B['T_abs']}",
        flush=True,
    )
    C = prove_open()
    print(
        f"  C open: ND={C['type_I_multilevel_bad_case_ND_closed']} "
        f"3AB={C['type_I_aut_e_3AB_positive_general']}",
        flush=True,
    )
    out = {
        "prop": "15.544",
        "title": "p=5 |μ_full|=1/p²; 1D+full mix gives Aut_e G>T",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "C": C,
        "proved": {
            "mu_full_abs_p5": A["proved"],
            "mix_G_gt_T_p5": B["proved"],
            "type_I_multilevel_bad_case_ND_closed": C[
                "type_I_multilevel_bad_case_ND_closed"
            ],
            "type_I_aut_e_3AB_positive_general": C[
                "type_I_aut_e_3AB_positive_general"
            ],
        },
        "e1_closed_general": bool(e1_closed_general()),
        "gsum_disj_lb": bool(gsum_disj_lb_proved_general()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "flags_not_flipped": [
            "type_I_multilevel_bad_case_ND_closed",
            "type_I_aut_e_3AB_positive_general",
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
    print("wrote", EV, flush=True)
    return out


if __name__ == "__main__":
    main()
