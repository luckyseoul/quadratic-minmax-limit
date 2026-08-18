#!/usr/bin/env python3
"""
Prop 15.564 — F=−2(3p²+2) is a complete character sum for every
odd prime p; hence Q3_02=−4N(2p²+1)/p is a p-law.

    I := Σ_{3-dist} e_θ(u) e_φ(v) = 2
    χ(−1)=1 on F_{p²},  G(χ,1)²=q
    J_u+J_v+J_w = −6 G χ_Ω on every Ω-triple
    F = −2 I + 2κ (J_u+J_v+J_w) = −2(3p²+2)
    Q3_02 = −2 N p + (N/p) F = −4 N (2p²+1)/p

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
residual_ii / type_I / 15.279–15.563 flags.  Does **not** import phi_F.
Does **not** overwrite 15.550–15.563.  Does **not** interpolate (p−5)/15.

============================================================================
Setup.  15.558: S_Ω(δ)=|Ω| (δ=0) else −1/2+(χ_Ω/2)G(χ,δ), and
F:=Σ_{3-dist}[1+2S_Ω(u)+2S_Ω(v)+2S_Ω(v−u)] e_θ(u)e_φ(v)
was certified =−2(3p²+2) at p=5..23, not proved.  Ω={G(χ,ξ)=p}.
On 3-distinct arguments S_Ω never hits 0, so
    1+2S+2S+2S = −2 + 2κ (χ(u)+χ(v)+χ(w)),  κ=χ_Ω G/2.

============================================================================
Theorem A — PROVED (inclusion-exclusion; any θ,φ,θ+φ≠0).
  Σ_{all u,v} e_θ(u)e_φ(v)=0.  Each of {u=0}, {v=0}, {u=v}
  sums to 0.  The three pairwise intersections and the triple
  intersection are the single point (0,0) with value 1.
  Σ_bad=−2, hence I=2.  Fail: I=0.  ∎

Theorem B — PROVED (quadratic χ on F_{p²}; certified p=5,7,11).
  (q−1)/2 is even, so χ(−1)=1.  G(χ,1)²=χ(−1)q=q.
  Fail: χ(−1)=−1.  Fail: G²=−q.  ∎

Theorem C — PROVED (Gauss inclusion-exclusion on an Ω-triple).
  J_u=−G(χ(θ)+χ(θ+φ)), J_v=−G(χ(φ)+χ(θ+φ)),
  J_w=−G(χ(φ)+χ(−1)χ(θ)).  On an Ω-triple θ,φ,−(θ+φ)∈Ω
  one has χ=χ_Ω on all three (since −Ω=Ω), hence
  J_u+J_v+J_w=−6 G χ_Ω.  Fail: drop the 6.  ∎

Theorem D — PROVED (A+B+C; field p=5,7,11).
  F=−4 + 2κ(−6 G χ_Ω)=−4−6 G²=−4−6q=−2(3p²+2).
  Therefore Q3_02=−4N(2p²+1)/p for every odd p (15.558 formula
  is now a p-law, not a p≤23 certificate).  Fail: drop +1.  ∎

Theorem E — OPEN.  Q3_dbl and Q3_n3 still carry A_full for m>3.
  Live ẑ=2 G_χ k̂ reconstructs A_dbl, but the measure on k is
  not uniform (uniform A3/A6/B4/2χ averages vanish).  Do not
  interpolate (p−5)/15.  phi_F not imported.

============================================================================
Backend: field tables O(q²) at p=5,7,11.  Writes
evidence/e1_gmin_m4_prop15564.json
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15279 import _e_tr, _kernel_field
from e1_gmin_m4_prop15550 import omega_set_local
from e1_gmin_m4_prop15558 import LIVE_N, q302_drop_plus1, q302_named

EV = ROOT / "evidence" / "e1_gmin_m4_prop15564.json"
LIVE_Q302 = {5: -5304.0, 7: -323928.0}


def gauss_G(F) -> complex:
    acc = 0j
    for x in range(1, F["q"]):
        acc += F["chi"][x] * _e_tr(F, x)
    return acc


def first_omega_triple(F, Omega):
    add, neg = F["add"], F["neg"]
    omset = set(Omega)
    for a in Omega:
        for b in Omega:
            c = neg[add[a][b]]
            if c in omset:
                return a, b, c
    raise RuntimeError("no Ω-triple")


def three_dist_sums(F, theta: int, phi: int):
    q = F["q"]
    chi, add, neg, mul = F["chi"], F["add"], F["neg"], F["mul"]
    I = Ju = Jv = Jw = 0j
    for u in range(1, q):
        for v in range(1, q):
            if u == v:
                continue
            w = add[v][neg[u]]
            if w == 0:
                continue
            e = _e_tr(F, mul[theta][u]) * _e_tr(F, mul[phi][v])
            I += e
            Ju += chi[u] * e
            Jv += chi[v] * e
            Jw += chi[w] * e
    return I, Ju, Jv, Jw


def F_named(p: int) -> float:
    return -2.0 * (3 * p * p + 2)


def F_wrong_I0(p: int) -> float:
    """Fail: I=0 removes the −4."""
    return -6.0 * p * p


def q302_from_F(p: int, N: int) -> float:
    return -2.0 * N * p + (N / p) * F_named(p)


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 11):
        Fld = _kernel_field(p)
        Omega = omega_set_local(Fld)
        _a, b, c = first_omega_triple(Fld, Omega)
        I, *_ = three_dist_sums(Fld, b, c)
        match = abs(I - 2) < 1e-8
        i0_fails = abs(I) > 1.0
        if not (match and i0_fails):
            ok = False
        rows[str(p)] = {
            "I_re": float(I.real),
            "err": float(abs(I - 2)),
            "match": bool(match),
            "I0_fails": bool(i0_fails),
        }
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": "I=Σ_{3-dist} e_θ(u)e_φ(v)=2. Fail I=0.",
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 11):
        Fld = _kernel_field(p)
        q = Fld["q"]
        chi_m1 = int(Fld["chi"][Fld["neg"][1]])
        G = gauss_G(Fld)
        match_m1 = chi_m1 == 1
        match_G2 = abs(G * G - q) < 1e-6
        g2_neg_fails = abs(G * G + q) > 1.0
        if not (match_m1 and match_G2 and g2_neg_fails):
            ok = False
        rows[str(p)] = {
            "chi_m1": chi_m1,
            "G_re": float(G.real),
            "G2_err": float(abs(G * G - q)),
            "match": bool(match_m1 and match_G2),
            "G2_neg_fails": bool(g2_neg_fails),
        }
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": "χ(−1)=1 on F_{p²} and G²=q. Fail χ(−1)=−1; fail G²=−q.",
    }


def prove_C() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 11):
        Fld = _kernel_field(p)
        Omega = omega_set_local(Fld)
        chi_om = int(Fld["chi"][Omega[0]])
        G = gauss_G(Fld)
        _a, b, c = first_omega_triple(Fld, Omega)
        _I, Ju, Jv, Jw = three_dist_sums(Fld, b, c)
        jsum = Ju + Jv + Jw
        want = -6.0 * G * chi_om
        match = abs(jsum - want) < 1e-6
        drop6 = abs(jsum + G * chi_om) > 1.0
        if not (match and drop6):
            ok = False
        rows[str(p)] = {
            "chi_om": chi_om,
            "err": float(abs(jsum - want)),
            "match": bool(match),
            "drop6_fails": bool(drop6),
        }
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": "J_u+J_v+J_w=−6 G χ_Ω on Ω-triples. Fail drop 6.",
    }


def prove_D() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 11):
        Fld = _kernel_field(p)
        Omega = omega_set_local(Fld)
        chi_om = int(Fld["chi"][Omega[0]])
        G = gauss_G(Fld)
        _a, b, c = first_omega_triple(Fld, Omega)
        I, Ju, Jv, Jw = three_dist_sums(Fld, b, c)
        kappa = chi_om * G / 2.0
        Fgot = -2.0 * I + 2.0 * kappa * (Ju + Jv + Jw)
        want = F_named(p)
        match_F = abs(Fgot - want) < 1e-6
        i0_fails = abs(want - F_wrong_I0(p)) > 1.0
        match_q = True
        live = None
        drop1 = True
        if p in LIVE_N:
            N = LIVE_N[p]
            named = q302_named(p, N)
            from_F = q302_from_F(p, N)
            live = LIVE_Q302[p]
            match_q = abs(named - from_F) < 1e-8 and abs(named - live) < 1e-6
            drop1 = abs(named - q302_drop_plus1(p, N)) > 1.0
        if not (match_F and i0_fails and match_q and drop1):
            ok = False
        rows[str(p)] = {
            "F_re": float(Fgot.real),
            "F_named": want,
            "match_F": bool(match_F),
            "I0_fails": bool(i0_fails),
            "live_Q302": live,
            "match_Q302": bool(match_q),
            "drop_plus1_fails": bool(drop1),
        }
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "F=−2(3p²+2); Q3_02=−4N(2p²+1)/p is a p-law. Fail drop +1."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "Q3_double_named_in_p": False,
        "Q3_n3_named_in_p": False,
        "A_full_dbl_named_in_p": False,
        "NUM_SUM_16pA_general": False,
        "phi_F_imported": False,
        "note": (
            "Q3_02 is a p-law. Q3_dbl / Q3_n3 still carry A_full for m>3. "
            "Uniform k-families give A_dbl=0; live k reconstructs A_dbl "
            "via ẑ=2 G_χ k̂ but the measure on k is unnamed. "
            "Do not interpolate (p−5)/15. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.564  F=−2(3p²+2) names Q3_02 for every odd p", flush=True)
    A = prove_A()
    print(f"  A I=2: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B χ(−1)=1, G²=q: {B['proved']}", flush=True)
    C = prove_C()
    print(f"  C Jsum=−6G χ_Ω: {C['proved']}", flush=True)
    D = prove_D()
    print(f"  D F and Q3_02: {D['proved']}", flush=True)
    E = prove_open()
    out = {
        "prop": "15.564",
        "title": "F=−2(3p²+2); Q3_02=−4N(2p²+1)/p is a p-law",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "C": C,
        "D": D,
        "E": E,
        "proved": {
            "I_equals_2": A["proved"],
            "chi_m1_and_G2": B["proved"],
            "Jsum_named": C["proved"],
            "F_and_Q302_named": D["proved"],
            "Q3_double_named_in_p": False,
            "Q3_n3_named_in_p": False,
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
        "identity": "F=-2(3p^2+2); Q3_02=-4N(2p^2+1)/p",
    }
    EV.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print("wrote", EV, flush=True)
    return out


if __name__ == "__main__":
    main()
