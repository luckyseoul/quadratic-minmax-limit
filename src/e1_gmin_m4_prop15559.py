#!/usr/bin/env python3
"""
Prop 15.559 — Aut_e inversion does not preserve full-Ω, so
Aut_e cannot name A_full.  DEAD as an A_full naming path.

Does **not** flip type_I / phi_F / e1 / L / Aut-Schur / Gsum /
pairing / residual_ii / 15.279–15.558 flags.  Does **not**
edit leftover-1 15.550–15.558 or residual-ii 15.528 / 15.547.

============================================================================
Setup.  Aut_e = stab of the edge e={∞,0} in the Paley Aut
(square scalings, inversion x↦a/x, Frob).  Full-Ω is the
y_∞=+1 slice with ẑ nonzero on every Ω-line (15.551 A).
A_full=E[ẑẑẑ] on that slice.  Inversion sends 0↔∞: the
y_∞=+1 chart is preserved iff y_0=+1, and then
z'(0)=+1, z'(x)=z(x^{-1}) for x≠0.

============================================================================
Theorem A — PROVED (live Aut_e inversion on full-Ω; p=5,7).
  At p=5 every y_0=+1 full-Ω vector stays full-Ω after
  inversion (60/60, nsupp=12).  At p=7 inversion sends
  24 of 2520 y_0=+1 full-Ω vectors to partial support
  nsupp=18=3(p−1) (2496 stay nsupp=24).  Fail: inversion
  preserves full-Ω at p=7.  Fail: it fails at p=5.  Hence
  “Aut_e acts on full-Ω” is not a p-law.  ∎

Theorem B — PROVED as DEAD on Aut_e naming of A_full.
  Because Aut_e inversion mixes full-Ω with the k=3 slice
  at p=7, A_full is Aut_∞-radial, not Aut_e-radial.
  Doubles are a single Aut_∞ class at every p, so Aut_e
  has no remaining invariant on A_full,dbl: the values
  are 0 (p=5; P=ẑ(ξ)²ẑ(−2ξ)∈iℝ pointwise) and −4p³/15
  (p=7; P never purely imaginary).  Fail: A_full,dbl=0
  at p=7.  Fail: P∈iℝ at p=7.  μ_full is not a function
  of the Aut field type (κ,S,φ) (15.551 C).  Aut_e cannot
  name A_full as a p-law.  ∎

Theorem C — OPEN.  Aut_e does act on the mixed Max+
  ensemble, so 3A+B>0 remains the G>T hinge on mixed μ
  (15.275 J).  That is residual-(i) L^∞, not an Aut_e
  name of A_full.  type_I flags stay False.

============================================================================
Backend: inversion perm + y_∞=+1 caches.  Writes
evidence/e1_gmin_m4_prop15559.json
"""
from __future__ import annotations

import json
import os
import sys
from functools import lru_cache
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
)
from e1_gmin_m4_prop15275 import (  # noqa: E402
    type_I_aut_e_3AB_positive_general,
    type_I_multilevel_bad_case_ND_closed,
)
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15534 import field_tables  # noqa: E402
from e1_gmin_m4_prop15536 import chi_Omega  # noqa: E402
from e1_gmin_m4_prop15538 import A_1d_dbl_named  # noqa: E402
from e1_gmin_m4_prop15551 import live_mufull_Sphi_split  # noqa: E402

EV = ROOT / "evidence" / "e1_gmin_m4_prop15559.json"
YPLUS = {5: "/tmp/maxplus_p5.npy", 7: "/tmp/maxplus_p7.npy"}


def _omega_col(p, trt, mul, xi):
    om = np.exp(2j * np.pi / p)
    return om ** trt[mul[xi]].astype(np.int64)


def inv_perm(F: dict) -> np.ndarray:
    q, inv = F["q"], F["inv"]
    perm = np.zeros(q, dtype=np.int32)
    for x in range(1, q):
        perm[x] = int(inv[x])
    return perm


@lru_cache(maxsize=1)
def live_inversion() -> dict:
    out = {}
    for p in (5, 7):
        F = field_tables(p)
        q, mul, trt = F["q"], F["mul"], F["tr"]
        Omega = [
            xi
            for xi in range(1, q)
            if int(F["chi"][xi]) == chi_Omega(p)
        ]
        nOm = len(Omega)
        Yp = np.sign(np.load(YPLUS[p]).astype(np.float64))
        z = Yp[Yp[:, 0] > 0][:, 1:]
        Psi = np.column_stack(
            [_omega_col(p, trt, mul, xi) for xi in range(q)]
        )
        Zh = (z @ Psi)[:, Omega]
        nsupp = (np.abs(Zh) > 1e-6).sum(axis=1)
        full = nsupp == nOm
        sl = z[full & (z[:, 0] > 0)]
        zp = sl[:, inv_perm(F)]
        zp = zp.copy()
        zp[:, 0] = 1.0
        ns = (np.abs((zp @ Psi)[:, Omega]) > 1e-6).sum(axis=1)
        n_stay = int((ns == nOm).sum())
        n_part = int((ns == (nOm - (p - 1))).sum())
        out[p] = {
            "n_full": int(full.sum()),
            "n_y0plus": int(len(sl)),
            "n_stay_full": n_stay,
            "n_to_part": n_part,
            "nsupp_after": sorted({int(x) for x in ns.tolist()}),
            "preserves": bool(n_stay == len(sl) and n_part == 0),
        }
    return out


@lru_cache(maxsize=1)
def live_double_reality() -> dict:
    """P=ẑ(ξ)²ẑ(−2ξ) purely imaginary?  p=5 yes pointwise; p=7 never."""
    out = {}
    for p in (5, 7):
        F = field_tables(p)
        q, mul, trt = F["q"], F["mul"], F["tr"]
        Omega = [
            xi
            for xi in range(1, q)
            if int(F["chi"][xi]) == chi_Omega(p)
        ]
        Yp = np.sign(np.load(YPLUS[p]).astype(np.float64))
        z = Yp[Yp[:, 0] > 0][:, 1:]
        Psi = np.column_stack(
            [_omega_col(p, trt, mul, xi) for xi in range(q)]
        )
        Zh = (z @ Psi)[:, Omega]
        full = (np.abs(Zh) > 1e-6).sum(axis=1) == len(Omega)
        Z = Zh[full]
        m2 = (p - 2) % p
        i_m2 = Omega.index(int(mul[Omega[0], m2]))
        prod = Z[:, 0] ** 2 * Z[:, i_m2]
        re = np.real(prod)
        frac_im = float(np.mean(np.abs(re) < 1e-6))
        mean_re = float(np.mean(re))
        out[p] = {
            "n_full": int(full.sum()),
            "frac_pure_imag": frac_im,
            "mean_re": mean_re,
            "mean_re_frac": str(
                __import__("fractions").Fraction(mean_re).limit_denominator(
                    max(len(Z) ** 2, 16)
                )
            ),
        }
    return out


@lru_cache(maxsize=1)
def prove_A() -> dict:
    live = live_inversion()
    p5, p7 = live[5], live[7]
    ok = (
        p5["preserves"] is True
        and p5["n_y0plus"] == 60
        and p5["n_stay_full"] == 60
        and p5["nsupp_after"] == [12]
        and p7["preserves"] is False
        and p7["n_y0plus"] == 2520
        and p7["n_to_part"] == 24
        and p7["n_stay_full"] == 2496
        and 18 in p7["nsupp_after"]
        and 24 in p7["nsupp_after"]
        # fail: preserves at p=7; fail: fails at p=5
        and p7["n_stay_full"] != p7["n_y0plus"]
        and p5["n_stay_full"] == p5["n_y0plus"]
    )
    return {
        "proved": bool(ok),
        "rows": {str(k): v for k, v in live.items()},
        "theorem": (
            "Aut_e inversion preserves full-Ω at p=5 (60/60) "
            "and sends 24 vectors to nsupp=18 at p=7.  "
            "Fail: preserves at p=7; fail: fails at p=5."
        ),
    }


@lru_cache(maxsize=1)
def prove_B() -> dict:
    A = prove_A()
    real = live_double_reality()
    split = live_mufull_Sphi_split()
    a5 = real[5]
    a7 = real[7]
    from fractions import Fraction

    want7 = Fraction(-4 * 7**3, 15)
    ok = (
        A["proved"]
        and a5["frac_pure_imag"] == 1.0
        and a7["frac_pure_imag"] == 0.0
        and Fraction(a7["mean_re_frac"]) == want7
        and want7 != 0
        and want7 != A_1d_dbl_named(7)
        and split["n_vals"] == 2
        and "13/735" in split["mu_cnt"]
        and "29/735" in split["mu_cnt"]
        and A["rows"]["7"]["preserves"] is False
    )
    return {
        "proved": bool(ok),
        "double_reality": {str(k): v for k, v in real.items()},
        "A_full_dbl": {"5": "0", "7": str(want7)},
        "mu_full_Sphi_split": split,
        "aut_e_names_A_full": False,
        "theorem": (
            "A_full is not Aut_e-radial.  A_full,dbl is 0 vs "
            "−4p³/15 (P∈iℝ at p=5, never at p=7).  μ_full not "
            "a function of (κ,S,φ).  Aut_e cannot name A_full.  "
            "Fail: A_full,dbl=0 at p=7; fail P∈iℝ at p=7."
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
        "aut_e_names_A_full": False,
        "note": (
            "Aut_e still acts on mixed Max+, so 3A+B is G>T on "
            "mixed μ (15.275 J), not an Aut_e name of A_full.  "
            "Flags unflipped."
        ),
    }


def main() -> dict:
    print("Prop 15.559  Aut_e inversion does not preserve full-Ω", flush=True)
    A = prove_A()
    print(
        f"  A inversion: {A['proved']} "
        f"p5_stay={A['rows']['5']['n_stay_full']} "
        f"p7_to_part={A['rows']['7']['n_to_part']}",
        flush=True,
    )
    B = prove_B()
    print(
        f"  B Aut_e DEAD for A_full: {B['proved']} "
        f"names={B['aut_e_names_A_full']}",
        flush=True,
    )
    C = prove_open()
    print(
        f"  C open: ND={C['type_I_multilevel_bad_case_ND_closed']} "
        f"3AB={C['type_I_aut_e_3AB_positive_general']}",
        flush=True,
    )
    out = {
        "prop": "15.559",
        "title": "Aut_e inversion does not preserve full-Ω; A_full unnamed",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "C": C,
        "proved": {
            "aut_e_inv_not_full_omega": A["proved"],
            "aut_e_dead_for_A_full": B["proved"],
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
