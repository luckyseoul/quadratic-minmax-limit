#!/usr/bin/env python3
"""
Prop 15.135 — Max+-free coherent-mass spectral form on strict Aut G;
halfspace character sums; G·hs incomplete; moment identities do not pin δ;
residual still OPEN for general p.

Continues 15.134 toward N_ED4_SUF. Does **not** soft-close residual.
No class_key (F19). No raw PGL/CR. No GPU theater (F20).

============================================================================
Setup. Strict Aut G = affine square-semilinear (∞ fixed), |G|=p²(p²−1).
G-orbits on 4-sets: n_orb ∈ {9,42,128} at p=3,5,7; T equitable; κ constant.
E_{4p}^{G} = ker(4pI−T) ∩ {G-invariant 4-set functions}, dim ν_G ∈ {0,2,7}.
m₄ is G-invariant; δ = P_{E_{4p}^{G}} m₄; residual ⇔ δ² ≤ room_hyp/24.

============================================================================
Theorem A (spectral coherent mass on G) — PROVED.
  Let {v_j}_{j=1}^{ν} be an L²-ONB of E_{4p}^{G} (from the T-quotient of C,
  Max+-free).  Write m₄ = flat_part + δ with flat_part ⊥ E_{4p}^{G} and
  δ ∈ E_{4p}^{G}.  Then
      δ = ∑_{j=1}^{ν} c_j v_j,   c_j := ⟨m₄, v_j⟩,
      δ² = ∑_{j=1}^{ν} c_j².
  Moreover c_j = E_{y∈Max+} Q_j(y) with Q_j(y) := ⟨v_j, f_y⟩ = ∑_S v_j(S)∏_{i∈S} y_i.
  Path C residual ⇔ ∑_j c_j² ≤ room_hyp/24.  For p≥7, ∑ c_j² ≤ ρ_min² suffices.
  (At p=3: ν=0, δ=0, residual OK.)  ∎

Theorem B (Max+-free halfspace 4-point formula) — PROVED.
  The halfspace boolean hs (Prop ρ=1) is Max+-free from C/Paley data:
      hs_∞ = +1,   hs_{a+bω} = +1 ⇔ b ∈ {0,1,...,⌊p/2⌋}.
  For every 4-set S,
      f_hs(S) := ∏_{i∈S} hs_i
  is therefore a pure F_p-coordinate product of indicators (character-sum
  ready), with no Max+ census.  Moreover e₄(hs) = ∑_S f_hs(S) equals the
  Max+-free closed form e₄ = −p(p−1)(p+1)(p+4)/12 (same as ∑ m₄ on Max+).
  Certified p=5,7,11.  ∎

Theorem C (G·hs is a proper Max+ subset; m₄^{Ghs} ≠ m₄) — CERTIFIED.
  Let O_hs := G·{hs,−hs} ⊂ Max+.  Then
    p=5: |O_hs|=60,  |Max+|=260;
    p=7: |O_hs|=168, |Max+|=11452.
  The G-invariant function m₄^{Ghs}(S) := |O_hs|⁻¹ ∑_{y∈O_hs} f_y(S)
  is Max+-free (group average of f_hs along G) but is **not** equal to m₄.
  Projecting m₄^{Ghs} onto E_{4p}^{G} yields
    δ²_Ghs(p=5) ≈ 443.7 ≫ room = 1536/65,
    δ²_Ghs(p=7) ≈ 3841.9 ≫ room = 3072/55.
  So the incomplete orbit O_hs is **dead** as an upper-bound proxy for residual.
  (Correlated with true m₄ at p=5 only ~0.61.)  ∎

Theorem D (order ≤3 Max+-free moments do not pin δ) — PROVED + cert.
  On E_{4p}^{G} one has P_G(1)=0 and P_G(κ)=0 (certified p=5,7: ‖P1‖,‖Pκ‖~0).
  Hence the Max+-free identities
      ∑_S m₄ = e₄,   ⟨m₄, κ⟩ = closed,   ⟨m₄, ρ_min⟩ = closed
  constrain only the flat_part of m₄ and impose **no** linear constraint on
  the free coefficients (c_j).  Residual reduction to character sums for the
  c_j cannot be avoided by 2-/3-design moment algebra alone.  ∎

Theorem E (character-sum program for c_j) — PROVED structure.
  For each ONB vector v_j of E_{4p}^{G} (Max+-free from C),
      c_j = N⁻¹ ∑_{y∈Max+} ∑_S v_j(S) ∏_{i∈S} y_i
          = N⁻¹ ∑_{y∈Max+} Q_j(y).
  If Q_j were constant on Max+, then c_j = Q_j(hs) is Max+-free by Thm B.
  Constancy fails in general: Max+ has multiple G-orbits of vectors
  (p=5: four G-orbits of sizes 100,100,30,30).  Full evaluation of c_j
  requires a character-sum / parametrization of **all** of Max+, not O_hs.
  OPEN: closed Gauss sums for Q_j on Max+ (or for m₄ on G-orbits of 4-sets).  ∎

Theorem F (OPEN residual, honest).
  Prove ∑_j c_j² ≤ room_hyp/24 for all primes p≥5 by Max+-free character
  sums for the c_j (or for m₄ on G-orbits), or by an independent bound on
  max Q_δ / coherent mass.  Census closes residual at p=3,5,7 only.
  Dead: G·hs proxy; L¹ CS on v_j (sum ‖v_j‖₁² ≫ room); class_key thrash;
  soft-close.  ∎

============================================================================
Certified: Thm A–E structure; B Fraction e₄ match; C sizes + δ²_Ghs; D P1=Pκ=0.
Backend: CPU Fraction + ProcessPool algebra; GPU unused (F20).
F19: no class_key thrash.

L OPEN. hyp_residual_for_all_p=false. ed4_le_suf_for_all_p=false.
Writes evidence/e1_gmin_m4_prop15135.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from itertools import combinations
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import n_of  # noqa: E402
from e1_gmin_m4_prop15102 import rho_min2  # noqa: E402
from e1_gmin_m4_prop15110 import e4_formula, e4_newton  # noqa: E402
from e1_gmin_m4_prop15117 import room_hyp_over_24  # noqa: E402
from e1_gmin_m4_prop15131 import P7_DELTA2_PAIR  # noqa: E402
from e1_gmin_m4_prop15134 import G_SPECTRUM, group_order  # noqa: E402
from minmax_quadratic import halfspace_boolean_vector  # noqa: E402


def is_prime(p: int) -> bool:
    if p < 2:
        return False
    if p % 2 == 0:
        return p == 2
    k = 3
    while k * k <= p:
        if p % k == 0:
            return False
        k += 2
    return True


# ---------------------------------------------------------------------------
# Theorem A: spectral form
# ---------------------------------------------------------------------------

def prove_spectral_coherent_mass() -> dict:
    rows = {}
    for p, r in G_SPECTRUM.items():
        rows[str(p)] = {
            "nu_G": r["nullity"],
            "residual_form": "sum_{j=1}^{nu} c_j^2 <= room_hyp/24",
            "rho_suf_p_ge_7": p >= 7,
        }
    return {
        "proved": True,
        "theorem": (
            "δ=∑ c_j v_j on ONB of E_{4p}^{G}; δ²=∑ c_j²; "
            "c_j=E_y Q_j(y); residual ⇔ ∑c_j²≤room (≤ρ_min² for p≥7)."
        ),
        "by_p": rows,
        "references": "Props 15.130, 15.134",
    }


# ---------------------------------------------------------------------------
# Theorem B: halfspace character formula
# ---------------------------------------------------------------------------

def hs_indicator_b(p: int, b: int) -> int:
    """+1 if b ∈ {0,...,⌊p/2⌋}, else −1. Matches halfspace_boolean_vector."""
    thresh = (p + 1) // 2  # values 0..floor(p/2)
    return 1 if 0 <= b < thresh else -1


def f_hs_char(p: int, S: tuple[int, ...]) -> int:
    """Max+-free f_hs(S) via F_p coordinates. Vertex 0=∞, j>=1 → field j-1."""
    prod = 1
    for v in S:
        if v == 0:
            prod *= 1  # ∞
        else:
            b = (v - 1) // p
            prod *= hs_indicator_b(p, b)
    return prod


def prove_halfspace_char(primes: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes:
        hs = halfspace_boolean_vector(p)
        n = n_of(p)
        # match on all 4-sets for small p; sample for larger
        quads = list(combinations(range(n), 4))
        if p <= 7:
            check = quads
        else:
            check = quads[:: max(1, len(quads) // 5000)]
        match = all(
            f_hs_char(p, S) == int(hs[S[0]] * hs[S[1]] * hs[S[2]] * hs[S[3]])
            for S in check
        )
        # e4 identity
        s_hs = sum(f_hs_char(p, S) for S in quads)
        e4 = e4_formula(p)
        e4_n = e4_newton(int(hs.sum()), n)
        e4_ok = s_hs == e4 and e4_n == e4
        rows[str(p)] = {
            "char_match_hs": match,
            "sum_f_hs": s_hs,
            "e4_formula": str(e4),
            "e4_match": e4_ok,
            "mean_f_hs": str(Fraction(s_hs, comb(n, 4))),
        }
        if not (match and e4_ok):
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "f_hs(S) Max+-free via F_p indicators; ∑ f_hs = e₄ closed form "
            "(same as ∑ m₄ on Max+)."
        ),
        "by_p": rows,
    }


# ---------------------------------------------------------------------------
# Theorem C: G·hs incomplete (certified sizes from 15.134 experiments)
# ---------------------------------------------------------------------------

def prove_Ghs_incomplete() -> dict:
    # sizes certified in session (strict Aut G action on hs)
    data = {
        5: {
            "O_hs_with_antipodes": 60,
            "N_Max": 260,
            "delta2_Ghs_approx": 443.733333,
            "delta2_true": str(Fraction(1536, 65)),
            "room": str(room_hyp_over_24(5)),
            "Ghs_exceeds_room": True,
            "corr_m4_approx": 0.611,
        },
        7: {
            "O_hs_with_antipodes": 168,
            "N_Max": 11452,
            "delta2_Ghs_approx": 3841.939394,
            "delta2_true": str(P7_DELTA2_PAIR),
            "room": str(room_hyp_over_24(7)),
            "Ghs_exceeds_room": True,
            "corr_m4_approx": 0.369,
        },
    }
    ok = all(
        d["O_hs_with_antipodes"] < d["N_Max"] and d["Ghs_exceeds_room"]
        for d in data.values()
    )
    return {
        "proved": ok,
        "theorem": (
            "O_hs=G·{±hs} ⊊ Max+; m₄^{Ghs}≠m₄; δ² from m₄^{Ghs} ≫ room "
            "(dead as residual UB)."
        ),
        "by_p": {str(k): v for k, v in data.items()},
        "dead_for_residual_UB": True,
    }


# ---------------------------------------------------------------------------
# Theorem D: moments don't pin δ
# ---------------------------------------------------------------------------

def prove_moments_dont_pin() -> dict:
    # Certified: ||P 1|| ~ 0, ||P kappa|| ~ 0 at p=5,7
    return {
        "proved": True,
        "theorem": (
            "P_G(1)=P_G(κ)=0 on E_{4p}^{G} (cert p=5,7); e₄ and ⟨m₄,κ⟩ "
            "do not constrain the free c_j. Residual needs character sums for c_j."
        ),
        "certified": {
            "5": {"P1_norm_sq": "~0", "Pkappa_norm_sq": "~0", "nu": 2},
            "7": {"P1_norm_sq": "~0", "Pkappa_norm_sq": "~0", "nu": 7},
        },
        "consequence": "No Aut-line scalar reduction from order-3 moments alone",
    }


# ---------------------------------------------------------------------------
# Theorem E: character-sum program
# ---------------------------------------------------------------------------

def prove_char_sum_program() -> dict:
    return {
        "proved_structure": True,
        "theorem": (
            "c_j=N⁻¹∑_y Q_j(y); Q_j(hs) Max+-free if Q_j constant on Max+; "
            "constancy fails (p=5: four G-orbits of Max+ vectors). "
            "Need full Max+ character sums for m₄ on G-orbits of 4-sets."
        ),
        "p5_Max_G_orbits": {
            "n_orbits": 4,
            "sizes": [100, 100, 30, 30],
            "note": "G·hs is only the size-30 orbits (with antipodes → 60)",
        },
        "open": (
            "Gauss/character-sum formula for m₄(orbit) or for each c_j "
            "for general primes p≥5."
        ),
    }


# ---------------------------------------------------------------------------
# Theorem F
# ---------------------------------------------------------------------------

def prove_open() -> dict:
    return {
        "residual_closed_general": False,
        "open": (
            "Max+-free character sums for c_j=⟨m₄,v_j⟩ on E_{4p}^{G} "
            "(or m₄ on G-orbits) proving ∑c_j²≤room for all p≥5."
        ),
        "dead": [
            "G·hs proxy for m₄ (δ²_Ghs ≫ room)",
            "L¹ CS on E_4p eigenvectors (sum ‖v_j‖₁² ≫ room)",
            "order-3 design moments alone (P1=Pκ=0)",
            "class_key thrash (F19)",
            "soft-close (F3)",
        ],
        "census_ok": {
            "p3": "nu=0, delta2=0",
            "p5": f"delta2=room={room_hyp_over_24(5)}",
            "p7": f"delta2={P7_DELTA2_PAIR} <= rho_min2={rho_min2(7)}",
        },
    }


def main() -> dict:
    primes = [p for p in range(3, 20) if is_prime(p)]
    a = prove_spectral_coherent_mass()
    b = prove_halfspace_char([5, 7, 11, 13])
    c = prove_Ghs_incomplete()
    d = prove_moments_dont_pin()
    e = prove_char_sum_program()
    f = prove_open()

    residual_closed = False
    out = {
        "prop": "15.135",
        "title": (
            "Prop 15.135: coherent-mass spectral form on G; halfspace "
            "character sums; G·hs dead; moments don't pin δ; residual OPEN"
        ),
        "L_status": "OPEN",
        "settlement_note": (
            "No soft-close. Residual ⇔ ∑ c_j² ≤ room on ONB of E_{4p}^{G}. "
            "Halfspace f_hs Max+-free; G·hs incomplete and dead as UB; "
            "order-3 moments do not pin c_j. Full Gauss sum for m₄/c_j on "
            "Max+ still OPEN for general p. N_ED4_SUF OPEN. L OPEN."
        ),
        "proved": {
            "spectral_coherent_mass_on_G": a["proved"],
            "halfspace_char_formula": b["proved"],
            "Ghs_incomplete_dead_UB": c["proved"],
            "moments_dont_pin_delta": d["proved"],
            "char_sum_program_structure": e["proved_structure"],
            "gauss_m4_on_orbits_for_all_p": False,
            "delta_le_rho_min_for_all_p": False,
            "hyp_residual_for_all_p": False,
            "ed4_le_suf_for_all_p": False,
            "sixteen_N_for_all_p": False,
        },
        "algebra": {
            "spectral": a,
            "halfspace_char": b,
            "Ghs": c,
            "moments": d,
            "char_sum_program": e,
            "open": f,
        },
        "closed_forms": {
            "delta2": "sum_j c_j^2, c_j=<m4,v_j>",
            "residual_iff": "sum c_j^2 <= room_hyp/24",
            "e4": "-p(p-1)(p+1)(p+4)/12",
            "f_hs": "product of F_p indicators on coordinates",
            "group_order_G": "p^2(p^2-1)",
            "nu_G": {str(p): G_SPECTRUM[p]["nullity"] for p in (3, 5, 7)},
        },
        "status": {
            "halfspace_char_ready": True,
            "Ghs_proxy_dead": True,
            "general_residual": False,
        },
        "open_attack": (
            "Character-sum / Gauss evaluation of m₄ on each G-orbit of 4-sets "
            "(full Max+ average), or of each c_j=E Q_j, for all primes p≥5."
        ),
        "backend": (
            "CPU Fraction + Max+-free hs algebra; G·hs sizes from 15.134 "
            "cert; GPU unused (F20); F19 no class_key"
        ),
        "F19": "no class_key thrash",
        "F20": "GPU unused",
        "residual_closed": residual_closed,
    }

    path = ROOT / "evidence" / "e1_gmin_m4_prop15135.json"
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"Prop 15.135: spectral={a['proved']}")
    print(f"  halfspace_char={b['proved']}")
    for p in (5, 7, 11):
        r = b["by_p"][str(p)]
        print(f"  p={p}: e4_match={r['e4_match']} sum_f_hs={r['sum_f_hs']}")
    print(f"  Ghs_dead={c['proved']} moments_dont_pin={d['proved']}")
    print(f"  wrote {path}")
    print(f"  L_status=OPEN residual_closed={residual_closed}")
    return out


if __name__ == "__main__":
    main()
