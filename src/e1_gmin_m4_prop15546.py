#!/usr/bin/env python3
"""
Prop 15.546 — p=7 mix |μ|≤109/2863<|T|; Aut_e G>T on the ensemble.

Does **not** flip type_I / phi_F / e1 / L / Aut-Schur / Gsum /
pairing / residual_ii / 15.279–15.545 flags.  Does **not**
edit leftover-1 15.522–15.545 or residual-ii 15.528.

============================================================================
Setup.  15.543: μ_1d=κ/(p(p−2)) on |κ|=1.  15.544: p=5 |μ_full|=1/p²
and mix 3/65<|T|.  At p=7 the y_∞=+1 slice splits 1D / full-Ω /
partial (nsupp=18).  T=−(p−2)/(p(2p−1)), |T|=5/91.

============================================================================
Theorem A — PROVED (finite from C; all 14112 |κ|=1 ∞-sets).
  At p=7,
      |μ_full| ≤ (4p+1)/(15 p²) = 29/735,
      |μ_part| ≤ (p−2)/(3 p²) = 5/147,
  with equality attained.  μ_full takes exactly the six values
  ±{13,17,29}/(15 p²) of sign κ.  Fail: |μ_full|=1/p²
  (29/735 ≠ 1/49).  Fail: |μ_part|=1/p² (5/147 ≠ 1/49).
  S(ρ) does not pin the three full magnitudes (|S|=2p−4 splits
  13 vs 29).  ∎

Theorem B — PROVED (p=7 mix; 15.543 A + live counts).
  n_1d=140, n_full=4410, n_part=1176, nH=5726.
      |μ| ≤ (n_1d |μ_1d| + n_full (4p+1)/(15 p²)
              + n_part (p−2)/(3 p²)) / nH
          = 109/2863 < |T|=5/91.
  Equality in the triangle is attained (census max 109/2863).
  Hence every |κ|=1 class has G≥−109/2863>T, so 3A+B>0 on
  the mixed ensemble at p=7.  Fail: drop the 15 in the full
  majorant (use (4p+1)/((p−2)p²)); then the bound is
  283/2863 > |T| and does not prove G>T.  ∎

Theorem C — OPEN for p≥11.  The two majorants are not shown
  to be p-laws (15 is 3(p−2) only at p=7; no partial slice
  at p=5).  A_full still takes three magnitudes.  Aut_e
  3A+B for general p still G>T.  type_I flags stay False.

============================================================================
Backend: p=7 1D/full/part masks + Fraction mix.  Writes
evidence/e1_gmin_m4_prop15546.json
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

EV = ROOT / "evidence" / "e1_gmin_m4_prop15546.json"
YPLUS7 = "/tmp/maxplus_p7.npy"

N1D, NFULL, NPART, NH = 140, 4410, 1176, 5726


def C_full_p7() -> Fraction:
    """(4p+1)/(15 p²) at p=7."""
    p = 7
    return Fraction(4 * p + 1, 15 * p * p)


def C_part_p7() -> Fraction:
    """(p−2)/(3 p²) at p=7."""
    p = 7
    return Fraction(p - 2, 3 * p * p)


def C_full_wrong_drop15() -> Fraction:
    """Fail-when-wrong: drop the 15, use (4p+1)/((p−2)p²)."""
    p = 7
    return Fraction(4 * p + 1, (p - 2) * p * p)


def mix_bound_p7() -> Fraction:
    p = 7
    return (
        N1D * abs(mu_1d_named(p, 1))
        + NFULL * C_full_p7()
        + NPART * C_part_p7()
    ) / NH


def mix_bound_wrong_drop15() -> Fraction:
    p = 7
    return (
        N1D * abs(mu_1d_named(p, 1))
        + NFULL * C_full_wrong_drop15()
        + NPART * C_part_p7()
    ) / NH


def _omega_col(p, trt, mul_t, xi):
    om = np.exp(2j * np.pi / p)
    return om ** trt[mul_t[xi]].astype(np.int64)


@lru_cache(maxsize=1)
def live_p7() -> dict:
    p = 7
    F = field_tables(p)
    q, ch, sub, mul, trt = F["q"], F["chi"], F["sub"], F["mul"], F["tr"]
    Yp = np.sign(np.load(YPLUS7).astype(np.float64))
    Yp = Yp[Yp[:, 0] > 0]
    z = Yp[:, 1:]
    Omega = [xi for xi in range(1, q) if int(ch[xi]) == chi_Omega(p)]
    Psi = np.column_stack([_omega_col(p, trt, mul, xi) for xi in range(q)])
    nsupp = (np.abs((z @ Psi)[:, Omega]) > 1e-6).sum(axis=1)
    is_1d = nsupp == (p - 1)
    is_full = nsupp == len(Omega)
    is_part = (~is_1d) & (~is_full)
    zf = z[is_full]
    zp = z[is_part]
    n_tot = 0
    max_f = max_p = max_m = 0.0
    seen_f: set[str] = set()
    for a, b, c in combinations(range(q), 3):
        kap = int(ch[sub[b, a]]) + int(ch[sub[c, a]]) + int(ch[sub[c, b]])
        if abs(kap) != 1:
            continue
        n_tot += 1
        muf = float(np.mean(zf[:, a] * zf[:, b] * zf[:, c]))
        mup = float(np.mean(zp[:, a] * zp[:, b] * zp[:, c]))
        mum = float(np.mean(z[:, a] * z[:, b] * z[:, c]))
        max_f = max(max_f, abs(muf))
        max_p = max(max_p, abs(mup))
        max_m = max(max_m, abs(mum))
        seen_f.add(str(Fraction(muf).limit_denominator(NFULL * NFULL)))
    return {
        "n_1d": int(is_1d.sum()),
        "n_full": int(is_full.sum()),
        "n_part": int(is_part.sum()),
        "nH": int(len(z)),
        "n_kappa1": n_tot,
        "max_full": str(Fraction(max_f).limit_denominator(NFULL * NFULL)),
        "max_part": str(Fraction(max_p).limit_denominator(NPART * NPART)),
        "max_mix": str(Fraction(max_m).limit_denominator(NH * NH)),
        "n_full_vals": len(seen_f),
    }


@lru_cache(maxsize=1)
def prove_A() -> dict:
    live = live_p7()
    ok = (
        live["n_1d"] == N1D
        and live["n_full"] == NFULL
        and live["n_part"] == NPART
        and live["nH"] == NH
        and live["n_kappa1"] == 14112
        and Fraction(live["max_full"]) == C_full_p7()
        and Fraction(live["max_part"]) == C_part_p7()
        and C_full_p7() != Fraction(1, 49)
        and C_part_p7() != Fraction(1, 49)
        and live["n_full_vals"] == 6
    )
    return {
        "proved": bool(ok),
        "live": live,
        "C_full": str(C_full_p7()),
        "C_part": str(C_part_p7()),
        "theorem": (
            "p=7: |μ_full|≤(4p+1)/(15p²), |μ_part|≤(p−2)/(3p²).  "
            "Fail: either equals 1/p²."
        ),
    }


@lru_cache(maxsize=1)
def prove_B() -> dict:
    live = live_p7()
    bound = mix_bound_p7()
    wrong = mix_bound_wrong_drop15()
    Tabs = -T_bitight(7)
    max_mu = Fraction(live["max_mix"])
    ok = (
        bound == Fraction(109, 2863)
        and bound < Tabs
        and wrong > Tabs
        and max_mu == Fraction(109, 2863)
        and live["nH"] == NH
    )
    return {
        "proved": bool(ok),
        "mix_bound": str(bound),
        "T_abs": str(Tabs),
        "drop_15": str(wrong),
        "live_max_mu": str(max_mu),
        "theorem": (
            "p=7 mix |μ|≤109/2863<|T|.  "
            "Fail: drop 15 in the full majorant (bound > T)."
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
            "Majorants not shown for p≥11.  A_full still three "
            "magnitudes.  Aut_e general still G>T.  Flags unflipped."
        ),
    }


def main() -> dict:
    print("Prop 15.546  p=7 mix |μ|≤109/2863<|T|", flush=True)
    A = prove_A()
    print(
        f"  A majorants: {A['proved']} "
        f"Cfull={A['C_full']} Cpart={A['C_part']}",
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
        "prop": "15.546",
        "title": "p=7 mix |μ|≤109/2863<|T|; Aut_e G>T on the ensemble",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "C": C,
        "proved": {
            "mu_majorants_p7": A["proved"],
            "mix_G_gt_T_p7": B["proved"],
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
