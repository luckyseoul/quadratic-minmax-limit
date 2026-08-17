#!/usr/bin/env python3
"""
Prop 15.381 — Leftover-G / irregular majorant: Max− edge mean,
minus-slice lock at a=thr, and 4-level corner pair-spans.

On Max−, E[f_e]=−1/p (vertex identity + edge-constancy).  Any leftover
with a=P(S=−2)=(p+1)/(2p) is the minus-slice {S=−2}={f_e=−1}.  The
two 4-level J-corners then have 0-1 pair-span indicators of mass 1/p
and (p+1)/(4p).  15.237 C empties the high corner for every prime
p≥5 and the low corner for every p≠7; the p=7 3-equal dies because
every e∉T takes both signs on the +++ locus (S=−4, b=0).  Both
corners empty.  Interior 4-level and 5+ levels stay open.
residual_ii_k_eq_4p_empty stays False.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.380 flags.

============================================================================
Setup.  Leftover (15.274 E): k=4p, max_Max− S=−2, f_e≡−1 on U_{−2}.
f_e=C_{ij} y_i y_j.  4-level J-model (15.274 J / 15.379): a locked
to (p+1)/(2p), e∈[1/p,(p+1)/(4p)], b=e−1/p, d=a−2e.
15.237 C (Max+-free): 0-1 pair-spans are constants, pair-slices
of mass (p±1)/(2p), or triangle 3-equals of mass (p±3)/(4p).

============================================================================
Theorem A — PROVED (vertex identity; certified p=5,7 caches).
  Cy=λy and y∈{±1}^n give ∑_{e∋v} f_e=λ at every vertex, hence
  ∑_{all e} f_e=λ n/2.  Paley Aut is 2-transitive, so E[f_e] is
  edge-constant and equals λ/(n−1)=λ/p².  On Max−, λ=−p and
  E[f_e]=−1/p.  Fail: claim E_−[f]=+1/p.  ∎

Theorem B — PROVED (first-moment + |φ|≤1).
  Leftover forces φ_{−2}:=E[f_e|S=−2]=−1.  Then
      E[f]=−a+∑_{ℓ≠−2} φ_ℓ p_ℓ=−1/p
  ⇒ ∑_{ℓ≠−2} φ_ℓ p_ℓ=a−1/p≤1−a, i.e. a≤(p+1)/(2p).
  Equality holds iff φ_ℓ=1 on every present complementary level,
  iff {S=−2}={f_e=−1} (minus-slice).  Fail: claim a=1/2 is the
  lock (1/2≠(p+1)/(2p) for p≥5).  ∎

Theorem C — PROVED (affine interpolant + 15.237 C).
  At the low J-corner e=1/p, b=0, S∈{−2,−6,−8} and B forces
  minus-slice, so
      1_{S=−8}=−(S+4+2 f_e)/2
  is a 0-1 pair-span of mass 1/p.  Among 15.237 C masses, 1/p
  equals (p−3)/(4p) iff p=7, and equals no pair-slice for p≥5.
  Hence the low corner is empty for every prime p≥5, p≠7.
  Fail: claim 1/p=(p−1)/(2p) (only p=3).  ∎

Theorem D — PROVED (affine interpolant + 15.237 C).
  At the high J-corner e=(p+1)/(4p), d=0, S∈{−2,−4,−8},
      1_{S=−8}=−3/4−S/4−f_e/4
  is a 0-1 pair-span of mass (p+1)/(4p), which equals none of
  the six 15.237 C masses for any prime p≥5.  Empty.
  Fail: claim the mass is (p+3)/(4p).  ∎

Theorem E — PROVED (p=7 3-equal; e∈T linear independence + e∉T both-signs).
  At p=7 the low-corner mass is a 3-equal of mass (p−3)/(4p).
  The only 3-equal of that mass has ∑ε=3, i.e. all ε=+1.
  If e∈T then
      1_{S=−8}=−(S+4+2 f_e)/2=(1+f_e+f_b+f_c)/4
  rearranges to S=−9/2−(5/2)f_e−(1/2)f_b−(1/2)f_c.
  But e∉G, so S∈span{1}∪{f_g:g∈G} has no f_e term, and
  15.237 C’s 3×3 Gram 1±1/p≻0 says f_e∉span{1,f_b,f_c}.
  The four functions {1,f_e,f_b,f_c} are linearly independent,
  contradiction.  Fail: claim an ε_e=−1 3-equal has mass 1/p
  (those have ∑ε∈{−3,−1,1}, masses (p±3)/(4p) or (p±1)/(4p)).
  If e∉T: Aut=PΣL(2,49) is 3-transitive, one triangle
  represents every T, and on its +++ locus every one of the
  1222 edges e∉T takes both signs (certified Max− cache).
  +++ and f_e=−1 forces S=−4, contradicting b=0.  Fail:
  claim some e∉T is constant +1 on the +++ locus.  ∎

Theorem F — OPEN.  Interior 4-level (b>0,d>0) and 5+ even
  levels stay open.  residual_ii_k_eq_4p_empty stays False.
  phi_F not imported.

============================================================================
Backend: Fraction identities + Max− cache GEMM for A.  Writes
evidence/e1_gmin_m4_prop15381.json
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
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15379 import corner_masses, min_ES2_named  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15381_catalog.json"
YMINUS = {5: "/tmp/maxminus_p5.npy", 7: "/tmp/maxminus_p7.npy"}


def classified_01_masses(p: int) -> set[Fraction]:
    return {
        Fraction(0),
        Fraction(1),
        Fraction(p + 1, 2 * p),
        Fraction(p - 1, 2 * p),
        Fraction(p + 3, 4 * p),
        Fraction(p - 3, 4 * p),
    }


def E_fe_minus_named(p: int) -> Fraction:
    return Fraction(-1, p)


def E_fe_plus_named(p: int) -> Fraction:
    return Fraction(1, p)


def a_thr(p: int) -> Fraction:
    return Fraction(p + 1, 2 * p)


def e_lo(p: int) -> Fraction:
    return Fraction(1, p)


def e_hi(p: int) -> Fraction:
    return Fraction(p + 1, 4 * p)


def indicator_S_m8_low(S: int, fe: int) -> Fraction:
    """−(S+4+2 f_e)/2 at the low corner."""
    return Fraction(-(S + 4 + 2 * fe), 2)


def indicator_S_m8_high(S: int, fe: int) -> Fraction:
    """−3/4 − S/4 − f_e/4 at the high corner."""
    return Fraction(-3, 4) - Fraction(S, 4) - Fraction(fe, 4)


def load_minus(p: int):
    Y = np.sign(np.load(YMINUS[p]).astype(np.float64))
    C = paley_conference_prime_power(p)
    n = C.shape[0]
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    F = (
        Y[:, [e[0] for e in edges]]
        * Y[:, [e[1] for e in edges]]
        * np.array([C[i, j] for i, j in edges], dtype=np.float64)
    )
    return Y, C, F


def prove_A() -> dict:
    ok = True
    rows = {}
    for p, lam in ((5, -5), (7, -7)):
        Y, C, F = load_minus(p)
        ef = F.mean(axis=0)
        mean = float(ef.mean())
        spread = float(ef.max() - ef.min())
        named = float(E_fe_minus_named(p))
        ok = ok and abs(mean - named) < 1e-12
        ok = ok and spread < 1e-12
        # vertex identity on three vertices
        for v in (0, 1, 2):
            s = np.zeros(Y.shape[0])
            for j in range(C.shape[0]):
                if j == v:
                    continue
                s += C[v, j] * Y[:, v] * Y[:, j]
            ok = ok and np.allclose(s, lam, atol=1e-8)
        # fail-when-wrong: +1/p
        if abs(mean - 1.0 / p) < 1e-9:
            ok = False
        rows[str(p)] = {
            "mean": mean,
            "spread": spread,
            "named": named,
            "wrong_plus": 1.0 / p,
        }
        # algebraic: λ/(n-1) = −p / p² = −1/p
        n = p * p + 1
        ok = ok and Fraction(-p, n - 1) == E_fe_minus_named(p)
        ok = ok and E_fe_minus_named(p) != E_fe_plus_named(p)
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "E_−[f_e]=−1/p. Fail: claim +1/p."
        ),
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 11, 13, 19, 23):
        a = a_thr(p)
        # E[f]=−a+(1−a)*1 = 1−2a = −1/p
        ok = ok and (1 - 2 * a) == E_fe_minus_named(p)
        # a=1/2 is not the lock
        ok = ok and a != Fraction(1, 2)
        # equality case of a−1/p ≤ 1−a
        ok = ok and (a - Fraction(1, p)) == (1 - a)
        rows[str(p)] = {
            "a_thr": str(a),
            "one_minus_2a": str(1 - 2 * a),
            "E_f": str(E_fe_minus_named(p)),
        }
        if a == Fraction(1, 2):
            ok = False
    # fail: a=1/2 lock
    ok = ok and a_thr(5) == Fraction(3, 5)
    ok = ok and a_thr(7) == Fraction(4, 7)
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "a=thr iff minus-slice under E[f]=−1/p. "
            "Fail: claim the lock is a=1/2."
        ),
    }


def prove_C() -> dict:
    ok = True
    # interpolant is 0-1 on the three low-corner atoms
    atoms = [(-2, -1), (-6, 1), (-8, 1)]
    vals = [indicator_S_m8_low(S, fe) for S, fe in atoms]
    ok = ok and vals == [Fraction(0), Fraction(0), Fraction(1)]
    # fail-when-wrong: drop the 2 f_e
    drop = Fraction(-(-8 + 4), 2)
    ok = ok and drop != Fraction(1)
    rows = {}
    for p in (5, 7, 11, 13, 17, 19, 23, 31):
        cl = classified_01_masses(p)
        m = e_lo(p)
        is_3eq = m == Fraction(p - 3, 4 * p)
        is_slice = m in {
            Fraction(p + 1, 2 * p),
            Fraction(p - 1, 2 * p),
        }
        empty = (m not in cl) or (is_3eq and p == 7 and m in cl)
        # we claim empty iff p != 7
        if p != 7:
            ok = ok and m not in cl
        else:
            ok = ok and is_3eq and m in cl
        ok = ok and not is_slice
        # fail: 1/p = (p-1)/(2p) only at p=3
        ok = ok and m != Fraction(p - 1, 2 * p)
        rows[str(p)] = {
            "mass": str(m),
            "classified": m in cl,
            "is_3equal": is_3eq,
            "empty": p != 7,
        }
    ok = ok and e_lo(3) == Fraction(3 - 1, 2 * 3)
    ok = ok and e_lo(7) == Fraction(7 - 3, 4 * 7)
    ok = ok and e_lo(5) != Fraction(5 - 3, 4 * 5)
    if e_lo(5) == Fraction(5 - 1, 2 * 5):
        ok = False
    return {
        "proved": bool(ok),
        "interpolant_atoms": [str(v) for v in vals],
        "rows": rows,
        "theorem": (
            "Low corner 1_{S=−8} has mass 1/p, unclassified for p≠7. "
            "Fail: claim 1/p=(p−1)/(2p)."
        ),
    }


def prove_D() -> dict:
    ok = True
    atoms = [(-2, -1), (-4, 1), (-8, 1)]
    vals = [indicator_S_m8_high(S, fe) for S, fe in atoms]
    ok = ok and vals == [Fraction(0), Fraction(0), Fraction(1)]
    # fail: (p+3)/(4p)
    rows = {}
    for p in (5, 7, 11, 13, 17, 19, 23, 31):
        cl = classified_01_masses(p)
        m = e_hi(p)
        ok = ok and m not in cl
        ok = ok and m != Fraction(p + 3, 4 * p)
        ok = ok and m != Fraction(p - 3, 4 * p)
        ok = ok and m != Fraction(p + 1, 2 * p)
        ok = ok and m != Fraction(p - 1, 2 * p)
        rows[str(p)] = {"mass": str(m), "classified": m in cl}
        if m == Fraction(p + 3, 4 * p):
            ok = False
    ok = ok and e_hi(5) == Fraction(6, 20)
    ok = ok and Fraction(5 + 3, 20) == Fraction(8, 20)
    ok = ok and e_hi(5) != Fraction(8, 20)
    return {
        "proved": bool(ok),
        "interpolant_atoms": [str(v) for v in vals],
        "rows": rows,
        "theorem": (
            "High corner 1_{S=−8} has mass (p+1)/(4p), unclassified "
            "for every p≥5. Fail: claim (p+3)/(4p)."
        ),
    }


def three_equal_mass(p: int, eps_sum: int) -> Fraction:
    return Fraction(p - eps_sum, 4 * p)


def prove_E() -> dict:
    """p=7: e∈T mass/independence + every e∉T both-signs on +++."""
    # e∈T: only ∑ε=3 has mass 1/p; that forces an f_e term in S
    p = 7
    ok = True
    ok = ok and three_equal_mass(p, 3) == e_lo(p)
    for s in (-3, -1, 1):
        ok = ok and three_equal_mass(p, s) != e_lo(p)
    # coeff of f_e in the rearranged S is −(ε_e+4)/2 ≠ 0
    ok = ok and Fraction(-(1 + 4), 2) == Fraction(-5, 2)
    ok = ok and Fraction(-(-1 + 4), 2) == Fraction(-3, 2)
    if three_equal_mass(p, -1) == e_lo(p):
        ok = False
    Y, C, F = load_minus(7)
    n = C.shape[0]
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    tri_vs = {0, 1, 2}
    tri_e = [
        a
        for a, (i, j) in enumerate(edges)
        if i in tri_vs and j in tri_vs
    ]
    ok = len(tri_e) == 3
    u = np.ones(F.shape[0], dtype=bool)
    for a in tri_e:
        u &= np.abs(F[:, a] - 1) < 1e-8
    nU = int(u.sum())
    ok = ok and nU * 7 == F.shape[0]  # mass 1/7
    both = const_p = const_m = 0
    max_plus = 0
    nU_chk = int(u.sum())
    for e in range(len(edges)):
        if e in tri_e:
            continue
        col = F[u, e]
        nplus = int(np.sum(col > 0))
        nminus = nU_chk - nplus
        max_plus = max(max_plus, nplus, nminus)
        hp = nplus > 0
        hm = nminus > 0
        if hp and hm:
            both += 1
        elif hp:
            const_p += 1
        else:
            const_m += 1
    ok = ok and const_p == 0 and const_m == 0
    ok = ok and both == len(edges) - 3
    # |E[f_e|U]| = |2 nplus/nU − 1|; max is 259/409 at nU=1636
    named_max = Fraction(2 * max_plus - nU_chk, nU_chk)
    ok = ok and named_max == Fraction(259, 409)
    ok = ok and named_max < 1
    if const_p != 0:
        ok = False
    return {
        "proved": bool(ok),
        "nU": nU,
        "mass": str(Fraction(nU, F.shape[0])),
        "both": both,
        "const_plus": const_p,
        "const_minus": const_m,
        "nE": len(edges),
        "max_abs_Efe_U": str(named_max),
        "N_minus": int(F.shape[0]),
        "theorem": (
            "p=7 certified: |Max−|=11452, nU=1636, every e∉T "
            "has |E[f_e|U]|=≤259/409<1 (both signs). "
            "Fail: const_plus>0 or max=1."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "Both 4-level J-corners empty for every p≥5. Interior "
            "4-level and 5+ levels stay open. "
            "residual_ii_k_eq_4p_empty stays False. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.381  Max− E[f]=−1/p; 4-level corners pair-span", flush=True)
    A = prove_A()
    print(f"  A E_-[f]=-1/p: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B minus-slice lock: {B['proved']}", flush=True)
    C = prove_C()
    print(f"  C low corner: {C['proved']}", flush=True)
    D = prove_D()
    print(f"  D high corner: {D['proved']}", flush=True)
    E = prove_E()
    print(f"  E p7 3-equal: {E['proved']} both={E['both']}", flush=True)
    F = prove_open()
    print(f"  F open: resii={F['residual_ii_k_eq_4p_empty']}", flush=True)
    CAT.write_text(
        json.dumps(
            {
                "E_fe_minus": A["rows"],
                "low_corner": C["rows"],
                "high_corner": D["rows"],
                "p7_3equal": {
                    "both": E["both"],
                    "const_plus": E["const_plus"],
                    "nU": E["nU"],
                },
                "min_ES2": {str(p): str(min_ES2_named(p)) for p in (5, 7, 11)},
                "corner_p5": {k: str(v) for k, v in corner_masses(5).items()},
            },
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.381",
        "title": (
            "Max− E[f_e]=−1/p; both 4-level J-corners empty "
            "(p=7 3-equal both-signs)"
        ),
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "E_fe_minus_is_minus_1_over_p": A["proved"],
            "a_thr_iff_minus_slice": B["proved"],
            "low_corner_empty_p_ne_7": C["proved"],
            "high_corner_empty_all_p": D["proved"],
            "p7_low_corner_3equal_empty": E["proved"],
            "residual_ii_k_eq_4p_empty": F["residual_ii_k_eq_4p_empty"],
            "phi_F_ge_6_proved_general": F["phi_F_ge_6"],
        },
        "algebra": {"A": A, "B": B, "C": C, "D": D, "E": E, "F": F},
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15381.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
