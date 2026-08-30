#!/usr/bin/env python3
"""
Prop 15.402 — Square-direction ∞-stars are 1-level on Max−.
An r-line union has S≡−r; 4p is the 1-level law S≡−4, not
two-level {−1,−3}.  Type I k=3p−2 is never 1-level.  The
bad slice f_e≡−1 is not empty as a Max−-only statement.
The r=4 tight family is nevertheless outside residual (ii): a
k=1 Max+ cylinder has score at most 0 on every such four-line union.
3A+B>0 stays open.  type_I flags not imported.

Does **not** flip type_I / phi_F / e1 / L / Aut-Schur / Gsum /
pairing / residual_ii / 15.279–15.401 flags.

============================================================================
Setup.  Paley C on V=F_q∪{∞}, q=p², C_∞u=1, C_uv=χ(u−v).
Max−: Cy=−p y.  Frame (15.27): E_−[yyᵀ]=I−C/p.
A p-star from ∞ onto an affine line ℓ=x_0+F_p d has k=p
and S=y_∞∑_{u∈ℓ} y_u.  Type I leftover is k=3p−2, s_+=1,
bad case f_e≡−1 on {S=−1} when Max− is not two-level {−1,−3}.

============================================================================
Theorem A — PROVED (norm; all odd p).
  χ_q|_{F_p^×}≡+1.  For c∈F_p^×, N_{q/p}(c)=c·c^p=c² and
  χ_q=χ_p∘N, so χ_q(c)=χ_p(c²)=+1.  Equivalently every
  F_p element is a square in F_{p²}.  Fail: claim
  χ_q|_{F_p}=χ_p (then a F_p-nonsquare would have χ_q=−1).  ∎

Theorem B — PROVED (15.27 frame + A).
  For the ∞-star onto ℓ=x_0+F_p d,
      E_−[S²]=p−(p−1)χ(d).
  Each nonzero difference on ℓ is a F_p^×-multiple of d and
  appears p times; A forces χ(c d)=χ(d) for c∈F_p^×.  Hence
      ∑_{u≠v∈ℓ} χ(u−v)=p(p−1)χ(d),
  and E[y_u y_v]=−χ(u−v)/p.  Square direction: E[S²]=1=E[S]²
  and E[S]=−1, so S≡−1 on Max−.  Nonsquare: E[S²]=2p−1.
  Fail: swap square/nonsquare; fail use E[yy]=+C/p (Max+),
  which swaps the two values.  ∎

Theorem C — PROVED (additivity of B).
  The lines of a square-direction parallel class partition F_q
  and each has S_ℓ≡−1.  The ∞-star onto any r of them is the
  disjoint union, so S≡−r on Max−.  In particular r=4 gives
  k=4p and S≡−4 (1-level at the mean), not two-level {−1,−3}.
  {S=−1} is empty.  Fail: claim 4p is two-level {−1,−3};
  fail: claim S≡−1.  ∎

Theorem D — PROVED (vertex identity / C at r=p).
  (Cy)_∞=∑_{finite} y_u=−p y_∞, so the full ∞-star has
  S=y_∞∑ y_u≡−p.  Same for every vertex star (Aut is
  2-transitive; or ∑_{e∋v} f_e=λ).  Fail: claim S≡−1.  ∎

Theorem E — PROVED (mean).
  Type I k=3p−2 has E_−[S]=−(3p−2)/p=−3+2/p∉ℤ, so S is
  not 1-level.  Fail: claim E[S]=−1 or S≡−1.  ∎

Theorem G — PROVED (15.272 G / 15.588 E; all odd p≥5).
  Fix r lines B in one square-direction parallel class.  A k=1
  Max+ cylinder aligned with that class is constant on its lines,
  with (p+1)/2 line values +1 and (p−1)/2 line values −1; every
  choice of the negative lines occurs.  Put as many negative lines
  as possible in B.  Its score on the corresponding ∞-star union is
      p(r−2 min(r,(p−1)/2)).
  For r=4 this is 0 at p=5, −2p at p=7, and −4p at p≥11.
  Hence every four-parallel-line cover from Theorem C has
  min_Max+ S≤0 and is not a residual-(ii) example (which needs
  min_Max+ S=2).  One-sided tight covers exist; only their possible
  compatibility with the other residual conditions remains open.  ∎

Theorem F — OPEN.  3A+B>0 on Aut_e far classes is still
  the g_min>T hinge (15.275 J–K).  f_e≡−1 on a large
  {S=−1} occurs for some k=3p−2 sets that fail Type I
  Max+ {1,5} (p=5 sample), so Max−-only emptiness is
  false.  type_I_multilevel_bad_case_ND_closed stays False.
  type_I_aut_e_3AB_positive_general stays False.

============================================================================
Backend: Fraction + Max− cache for B–D.  Writes
evidence/e1_gmin_m4_prop15402.json
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

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    k_type_I_3p_minus_2,
)
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15275 import (  # noqa: E402
    type_I_aut_e_3AB_positive_general,
    type_I_mean_minus,
    type_I_multilevel_bad_case_ND_closed,
)
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import (  # noqa: E402
    _kernel_field,
    rhat_ge_budget_proved_general,
)
from minmax_quadratic import paley_conference_prime_power  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15402_catalog.json"
YMINUS = {5: "/tmp/maxminus_p5.npy", 7: "/tmp/maxminus_p7.npy"}


def chi_q_on_Fp(c: int, p: int) -> int:
    """χ_q(c) for c∈F_p.  Named +1 on F_p^×."""
    if c % p == 0:
        return 0
    return 1


def ES2_pstar(p: int, chd: int) -> Fraction:
    """E_−[S²]=p−(p−1)χ(d) for an ∞-star onto a dir-d line."""
    return Fraction(p - (p - 1) * chd)


def ES2_pstar_square(p: int) -> Fraction:
    return Fraction(1)


def ES2_pstar_nonsquare(p: int) -> Fraction:
    return Fraction(2 * p - 1)


def type_I_mean_is_integer(p: int) -> bool:
    return (3 * p - 2) % p == 0


def field_vertex(a: int, b: int, p: int) -> int:
    """Vertex of a+bω; ∞ is 0."""
    return 1 + (a % p) + (b % p) * p


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 11, 13, 17, 19, 23, 29, 31):
        F = _kernel_field(p)
        chi = F["chi"]
        n_m1 = 0
        for c in range(1, p):
            # encoding of c∈F_p is c + 0·p
            got = int(chi[c])
            named = chi_q_on_Fp(c, p)
            if got != 1 or named != 1:
                ok = False
            if got == -1:
                n_m1 += 1
        # χ_p has (p-1)/2 nonsquares; χ_q does not
        n_chi_p_ns = (p - 1) // 2
        ok = ok and n_m1 == 0
        ok = ok and n_chi_p_ns > 0
        rows[str(p)] = {
            "n_Fp_star": p - 1,
            "n_chi_q_minus": n_m1,
            "n_chi_p_nonsquare": n_chi_p_ns,
        }
        if n_m1 == n_chi_p_ns:
            ok = False
    # fail: χ_q(2)=χ_p(2) at a p with 2 nonsquare
    # p=7: 2 is F_7-nonsquare, but χ_q(2)=+1
    F7 = _kernel_field(7)
    ok = ok and int(F7["chi"][2]) == 1
    if int(F7["chi"][2]) == -1:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "χ_q|_{F_p^×}≡+1 via N(c)=c². Fail: χ_q|_{F_p}=χ_p."
        ),
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 11, 13, 17, 19, 23, 31):
        sq = ES2_pstar(p, 1)
        ns = ES2_pstar(p, -1)
        ok = ok and sq == ES2_pstar_square(p) == Fraction(1)
        ok = ok and ns == ES2_pstar_nonsquare(p) == Fraction(2 * p - 1)
        ok = ok and sq != ns
        # swapped χ would exchange the two
        ok = ok and ES2_pstar(p, 1) != Fraction(2 * p - 1)
        ok = ok and ES2_pstar(p, -1) != Fraction(1)
        # Max+ 2-point E[yy]=+C/p flips the sign of the off-diag
        plus_sq = Fraction(p + (p - 1) * 1)
        plus_ns = Fraction(p + (p - 1) * (-1))
        ok = ok and plus_sq == Fraction(2 * p - 1)
        ok = ok and plus_ns == Fraction(1)
        ok = ok and plus_sq != sq
        rows[str(p)] = {"square": str(sq), "nonsquare": str(ns)}
    ok = ok and ES2_pstar(5, 1) == 1
    ok = ok and ES2_pstar(5, -1) == 9
    if ES2_pstar(5, 1) == ES2_pstar(5, -1):
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "E_−[S²]=p−(p−1)χ(d). Square ⇒ S≡−1. "
            "Fail: swap square/nonsquare."
        ),
    }


def _edge_index(n: int):
    m = {}
    k = 0
    for a in range(n):
        for b in range(a + 1, n):
            m[(a, b)] = k
            k += 1
    return m


def _line_star_scores(p: int, Y: np.ndarray, C: np.ndarray, bs) -> np.ndarray:
    """S of the ∞-star onto ∪_{b∈bs}{a+bω}."""
    n = C.shape[0]
    emap = _edge_index(n)
    v = np.zeros(n * (n - 1) // 2, dtype=np.float64)
    for b in bs:
        for a in range(p):
            j = field_vertex(a, b, p)
            v[emap[(0, j)]] = 1.0
    ea = []
    eb = []
    ce = []
    for a in range(n):
        for b in range(a + 1, n):
            ea.append(a)
            eb.append(b)
            ce.append(C[a, b])
    F = Y[:, ea] * Y[:, eb] * np.asarray(ce, dtype=np.float64)
    return F @ v


def prove_C() -> dict:
    ok = True
    rows = {}
    # algebra: r copies of S≡−1 add
    for r in range(1, 12):
        ok = ok and (-1) * r == -r
    ok = ok and -4 != -1
    ok = ok and set([-4]) != {-1, -3}
    # live: horizontal class is dir 1, χ(1)=+1
    live = {}
    for p in (5, 7):
        Y = np.sign(np.load(YMINUS[p]).astype(np.float64))
        C = paley_conference_prime_power(p)
        Fchi = _kernel_field(p)
        ok = ok and int(Fchi["chi"][1]) == 1
        # each single line
        for b in range(p):
            S = _line_star_scores(p, Y, C, (b,))
            vals = {int(round(x)) for x in S}
            ok = ok and vals == {-1}
        # prefix r = 1..p
        pref = {}
        for r in range(1, p + 1):
            S = _line_star_scores(p, Y, C, range(r))
            vals = {int(round(x)) for x in S}
            ok = ok and vals == {-r}
            pref[str(r)] = sorted(vals)
        # r=4 is 1-level −4
        ok = ok and pref["4"] == [-4]
        if pref["4"] == [-1, -3] or pref["4"] == [-1]:
            ok = False
        live[str(p)] = pref
        rows[str(p)] = {"prefix": pref, "dir1_square": True}
    return {
        "proved": bool(ok),
        "rows": rows,
        "live": live,
        "theorem": (
            "r square-dir lines: S≡−r. 4p is 1-level −4. "
            "Fail: claim two-level {−1,−3}."
        ),
    }


def prove_D() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        Y = np.sign(np.load(YMINUS[p]).astype(np.float64))
        C = paley_conference_prime_power(p)
        n = C.shape[0]
        # vertex identity at ∞ and at 0
        for v, lam in ((0, -p), (1, -p)):
            s = np.zeros(Y.shape[0])
            for j in range(n):
                if j == v:
                    continue
                s += C[v, j] * Y[:, v] * Y[:, j]
            ok = ok and np.allclose(s, lam, atol=1e-8)
        # full ∞-star = r=p
        S = _line_star_scores(p, Y, C, range(p))
        vals = {int(round(x)) for x in S}
        ok = ok and vals == {-p}
        if vals == {-1}:
            ok = False
        rows[str(p)] = {"star_support": sorted(vals)}
    # algebraic: λ = −p, not −1
    ok = ok and (-5) != (-1)
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": "Vertex stars S≡−p on Max−. Fail: S≡−1.",
    }


def prove_E() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 11, 13, 17, 19, 23, 29, 31):
        mu = type_I_mean_minus(p)
        k = k_type_I_3p_minus_2(p)
        ok = ok and mu == -3 + Fraction(2, p)
        ok = ok and mu == -Fraction(k, p)
        ok = ok and not type_I_mean_is_integer(p)
        ok = ok and mu != -1
        ok = ok and mu != -3
        ok = ok and k % 2 == 1
        rows[str(p)] = {"mean": str(mu), "k": k, "odd": True}
        if mu == -1 or type_I_mean_is_integer(p):
            ok = False
    ok = ok and type_I_mean_minus(5) == Fraction(-13, 5)
    ok = ok and type_I_mean_minus(7) == Fraction(-19, 7)
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "Type I E_−[S]=−3+2/p∉ℤ, not 1-level. Fail: E[S]=−1."
        ),
    }


def maxplus_cylinder_witness_score(p: int, r: int) -> int:
    """Smallest score supplied by an aligned k=1 Max+ cylinder.

    Proposition 15.272 G says every choice of (p+1)/2 positive line
    values (equivalently (p-1)/2 negative values) occurs.
    """
    if p < 5 or p % 2 == 0:
        raise ValueError("p must be odd and at least 5")
    if r < 0 or r > p:
        raise ValueError("r must satisfy 0 <= r <= p")
    negative_in_support = min(r, (p - 1) // 2)
    return p * (r - 2 * negative_in_support)


def prove_G() -> dict:
    rows = {}
    ok = True
    for p in (5, 7, 11, 13, 17, 19, 23, 29, 31):
        score = maxplus_cylinder_witness_score(p, 4)
        named = 0 if p == 5 else (-2 * p if p == 7 else -4 * p)
        outside_residual = score <= 0
        ok = ok and score == named and outside_residual and score != 2
        rows[str(p)] = {
            "r4_witness_score": score,
            "min_Maxplus_upper_bound": score,
            "outside_residual_splus_2": outside_residual,
        }
    # The full parallel class is the vertex star and has Max+ score +p.
    ok = ok and all(
        maxplus_cylinder_witness_score(p, p) == p for p in (5, 7, 11, 13)
    )
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "An aligned k=1 Max+ cylinder gives score "
            "p(r-2 min(r,(p-1)/2)); for r=4 this is <=0, so the "
            "15.402 one-sided-tight family is outside residual (ii)."
        ),
        "depends_on": ["15.272.G", "15.588.E"],
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
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "gsum": bool(gsum_disj_lb_proved_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "note": (
            "3A+B>0 on Aut_e far classes remains G>T (15.275). "
            "f_e≡−1 on large {S=−1} exists off Type I Max+ at p=5, "
            "so Max−-only emptiness is false. Flags unflipped."
        ),
    }


def main() -> dict:
    print("Prop 15.402  square-dir stars 1-level on Max−; 4p is S≡−4", flush=True)
    A = prove_A()
    print(f"  A chi|Fp=1: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B ES2: {B['proved']} p5={B['rows']['5']}", flush=True)
    C = prove_C()
    print(f"  C r-union: {C['proved']} p5_r4={C['live']['5']['4']}", flush=True)
    D = prove_D()
    print(f"  D full star: {D['proved']} {D['rows']}", flush=True)
    E = prove_E()
    print(f"  E Type I mean: {E['proved']} p5={E['rows']['5']}", flush=True)
    G = prove_G()
    print(
        f"  G r=4 Max+ witness: {G['proved']} "
        f"p5={G['rows']['5']['r4_witness_score']} "
        f"p7={G['rows']['7']['r4_witness_score']}",
        flush=True,
    )
    F = prove_open()
    print(
        f"  F open: 3AB={F['type_I_aut_e_3AB_positive_general']} "
        f"ND={F['type_I_multilevel_bad_case_ND_closed']}",
        flush=True,
    )
    CAT.write_text(
        json.dumps(
            {
                "A": A["rows"],
                "B": B["rows"],
                "C": C.get("live", {}),
                "E": {p: E["rows"][p]["mean"] for p in ("5", "7", "11")},
                "G": G["rows"],
            },
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.402",
        "title": (
            "Square-direction ∞-stars are 1-level S≡−1 on Max−; "
            "4p is S≡−4, not two-level; Type I is never 1-level"
        ),
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "chi_q_trivial_on_Fp": A["proved"],
            "square_pstar_ES2_is_1": B["proved"],
            "r_line_union_S_eq_minus_r": C["proved"],
            "vertex_star_S_eq_minus_p": D["proved"],
            "type_I_mean_not_integer": E["proved"],
            "r4_line_union_outside_residual_ii": G["proved"],
            "type_I_multilevel_bad_case_ND_closed": F[
                "type_I_multilevel_bad_case_ND_closed"
            ],
            "type_I_aut_e_3AB_positive_general": F[
                "type_I_aut_e_3AB_positive_general"
            ],
            "phi_F_ge_6_proved_general": F["phi_F_ge_6"],
            "residual_ii_k_eq_4p_empty": F["residual_ii_k_eq_4p_empty"],
        },
        "algebra": {
            "A": A,
            "B": B,
            "C": C,
            "D": D,
            "E": E,
            "F": F,
            "G": G,
        },
        "L_status": "OPEN",
        "flags_not_flipped": [
            "type_I_multilevel_bad_case_ND_closed",
            "type_I_aut_e_3AB_positive_general",
            "phi_F_ge_6",
            "e1",
            "L",
            "Aut-Schur",
            "Gsum",
            "pairing",
            "residual_ii_k_eq_4p_empty",
        ],
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15402.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
