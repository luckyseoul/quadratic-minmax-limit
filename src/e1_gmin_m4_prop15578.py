#!/usr/bin/env python3
"""
Prop 15.578 — 2χ 4-point is Fejér 4p/(p−1); 15.568 families
hit the p=7 full-Ω buckets; the bucket mix is not a p-law.

    k=e_α−e_β  ⇒  |k̂(c)|²=4 sin²(π c δ/p),  δ=β−α
    E_δ[|k̂(1)|² |k̂(r)|²] = 4p/(p−1)   (r∉{0,±1})
    μ_{2χ}=16 p² · 4p/(p−1)=64 p³/(p−1)
    (ẑ=2 G_χ k̂, |G_χ|=√p ⇒ |ẑ|²|ẑ_r|²=16 p² |k̂|²|k̂_r|²)

Equals exclusive μ_full at p=5 (k=3=full).  At p=7 it is only
the 3/5 bucket (3658.67≠126224/15).  15.568 A3/A6/weighted-B4
4-points hit the remaining p=7 buckets (140/3, 35/3, 133/9).
Fail: uniform B4; fail drop 4-AP.  The 3/5:3/10:1/15:1/30 mix
dies at p=5 (2800≠2000) and is not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
residual_ii / type_I / 15.279–15.577 flags.  Does **not** import phi_F.
Does **not** overwrite 15.550–15.577.  Does **not** interpolate (p−5)/15.

============================================================================
Setup.  15.573 exclusive mix.  15.574–15.575 name μ_k3 and μ_1d.
Claude leftover-1 deep_review BLOCK on promoting the p=7 mix
to a p-identity (two-point fit is (p−5)/15-class).  This unit
names the 2χ 4-point only.

============================================================================
Theorem A — PROVED Max+-free (complete trigonometric sum).
  ∑_{δ=1}^{p−1} sin²(πδ/p) sin²(π r δ/p)=p/4 for r∉{0,±1},
  hence E4=4p/(p−1).  Fail: E4=4; fail: E4=5 at p=7 (p=5-only).
  Fail: r=1 (then the sum is not p/4).  ∎

Theorem B — PROVED (15.568 families; p=7).
  Uniform A3 4-point =140/3; A6 (3+,3−,1 zero)=35/3;
  weighted B4 1:4:2 =133/9.  Fail uniform B4 (15.4).
  Fail drop 4-AP.  ∎

Theorem C — OPEN.  The 15.568 mix is not a p-law (p=5:
  2800≠2000).  μ_full ≤ μ_k3 fails on N* at p=7
  (9695>9604).  Q_τ unnamed.  phi_F not imported.

============================================================================
Backend: O(p) Fejér + 15.568 family FFT at p=7.  Writes
evidence/e1_gmin_m4_prop15578.json
"""
from __future__ import annotations

import json
import os
import sys
from fractions import Fraction
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15568 import (
    b4_all,
    b4_named_weighted,
    fam_A3_named,
    fam_A6_named,
)

EV = ROOT / "evidence" / "e1_gmin_m4_prop15578.json"


def E4_2chi_named(p: int) -> Fraction:
    return Fraction(4 * p, p - 1)


def mu_2chi_named(p: int) -> Fraction:
    return Fraction(64 * p**3, p - 1)


def E4_2chi_const4() -> Fraction:
    return Fraction(4, 1)


def E4_2chi_five() -> Fraction:
    """Fail: the p=5 value 5."""
    return Fraction(5, 1)


def sine_sum(p: int, r: int) -> float:
    acc = 0.0
    for d in range(1, p):
        acc += (np.sin(np.pi * d / p) ** 2) * (np.sin(np.pi * r * d / p) ** 2)
    return acc


def E4_of(ks, p: int, r: int = 2, ws=None) -> float:
    num = den = 0.0
    for i, k in enumerate(ks):
        w = 1.0 if ws is None else float(ws[i])
        eh = np.fft.fft(np.asarray(k, dtype=np.float64))
        num += w * float(abs(eh[1]) ** 2 * abs(eh[r % p]) ** 2)
        den += w
    return num / den if den else float("nan")


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 11, 13):
        S = sine_sum(p, 2)
        e_num = 16.0 * S / (p - 1)
        want = float(E4_2chi_named(p))
        match = abs(S - p / 4.0) < 1e-9 and abs(e_num - want) < 1e-9
        fail4 = abs(want - 4.0) > 1e-9
        fail5 = abs(want - 5.0) > 1e-9 if p != 5 else True
        S1 = sine_sum(p, 1)
        fail_r1 = abs(S1 - p / 4.0) > 1e-6
        if not (match and fail4 and fail5 and fail_r1):
            ok = False
        rows[str(p)] = {
            "sine_sum": S,
            "p_over_4": p / 4.0,
            "E4": e_num,
            "named": want,
            "mu": float(mu_2chi_named(p)),
            "match": bool(match),
            "fail4": bool(fail4),
            "fail5_or_p5": bool(fail5),
            "r1_not_p4": bool(fail_r1),
        }
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": "2χ E4=4p/(p−1). Fail 4; fail 5 at p=7; fail r=1.",
    }


def prove_B() -> dict:
    p = 7
    a3 = fam_A3_named(p)
    a6 = fam_A6_named(p)
    bk, bw = b4_named_weighted(p)
    bu = b4_all(p)
    kd, wd = b4_named_weighted(p, drop_ap=True)
    e3, e6, ew, eu, ed = (
        E4_of(a3, p),
        E4_of(a6, p),
        E4_of(bk, p, ws=bw),
        E4_of(bu, p),
        E4_of(kd, p, ws=wd),
    )
    ok = (
        abs(e3 - 140 / 3) < 1e-6
        and abs(e6 - 35 / 3) < 1e-6
        and abs(ew - 133 / 9) < 1e-6
        and abs(eu - 133 / 9) > 0.1
        and abs(ed - 133 / 9) > 0.1
    )
    return {
        "proved": bool(ok),
        "E4_A3": e3,
        "E4_A6": e6,
        "E4_B4w": ew,
        "E4_B4u": eu,
        "E4_B4_dropAP": ed,
        "live": {"A3": 140 / 3, "A6": 35 / 3, "B4": 133 / 9},
        "theorem": (
            "p=7 A3/A6/weighted-B4 4-points hit 140/3, 35/3, 133/9. "
            "Fail uniform B4; fail drop 4-AP."
        ),
    }


def prove_open() -> dict:
    # p=5 15.568-style mix is not live μ_full
    mix_p5 = 2800.0
    return {
        "proved": False,
        "mix_p5": mix_p5,
        "live_p5_full": 2000.0,
        "mix_is_plaw": False,
        "mu_full_le_muk3_Nstar_p7": False,
        "Nstar_mufull_p7": 9695.466666666669,
        "muk3_p7": 9604.0,
        "Q_tau_named_in_p": False,
        "phi_F_imported": False,
        "note": (
            "2χ 4-point is named. The 15.568 bucket mix is p=7-only "
            "(2800≠2000 at p=5). μ_full≤μ_k3 fails on N* at p=7. "
            "Do not interpolate (p−5)/15. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.578  2χ E4=4p/(p−1); p=7 buckets; mix not a p-law", flush=True)
    A = prove_A()
    print(f"  A 2χ Fejér: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B p=7 families: {B['proved']}", flush=True)
    C = prove_open()
    out = {
        "prop": "15.578",
        "title": "2χ 4-point 4p/(p−1); 15.568 families hit p=7 buckets",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "C": C,
        "proved": {
            "twochi_E4_fejer": A["proved"],
            "p7_buckets_are_15568_families": B["proved"],
            "mu_full_named_in_p": False,
            "Q_tau_named_in_p": False,
            "phi_F_ge_6": False,
        },
        "e1_closed_general": bool(e1_closed_general()),
        "gsum_disj_lb": bool(gsum_disj_lb_proved_general()),
        "phi_F_ge_6": False,
        "L_status": "OPEN",
        "flags_not_flipped": [
            "phi_F_ge_6",
            "e1",
            "L",
            "Aut-Schur",
            "Gsum",
            "pairing",
            "residual_ii",
            "type_I",
        ],
        "identity": "E4_2chi = 4p/(p-1); mu_2chi = 64 p^3/(p-1)",
        "claude_leftover1": "deep_review BLOCK on mix-as-p-identity; do_not_import_mix",
    }
    EV.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print("wrote", EV, flush=True)
    return out


if __name__ == "__main__":
    main()
