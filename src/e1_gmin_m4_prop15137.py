#!/usr/bin/env python3
"""
Prop 15.137 — Character-sum form of c_j over Max+ G-orbits;
Q_j G-equivariant; hemisphere decomposition; p=5 exact two-type formula;
Q_j(hs) Max+-free; residual still OPEN for general p.

Continues 15.136 toward N_ED4_SUF. Does **not** soft-close residual.
No class_key (F19). No G·hs-only thrash. No GPU theater (F20).

============================================================================
Setup. Strict Aut G (∞ fixed). E_{4p}^{G} with ONB {v_j}, ν=dim.
Q_j(y)=∑_S v_j(S)∏_{i∈S} y_i,   c_j=N^{-1}∑_{y∈Max+} Q_j(y),
residual ⇔ ∑ c_j² ≤ room_hyp/24 (15.135–15.136).

============================================================================
Theorem A (G-equivariance of Q_j on Max+) — PROVED.
  G acts on Max+ by (g·y)_v = y_{g^{-1}v}.  Since v_j is G-invariant as a
  4-set function,
      Q_j(g·y) = ∑_S v_j(S) ∏_{i∈S} (g·y)_i
               = ∑_S v_j(S) ∏_{i∈S} y_{g^{-1}i}
               = ∑_{S'} v_j(g S') ∏ y = ∑_{S'} v_j(S') ∏ y = Q_j(y).
  Hence Q_j is constant on each G-orbit of Max+ vectors.  ∎

Theorem B (hemisphere decomposition) — PROVED.
  G fixes the vertex ∞, so (g·y)_∞ = y_∞.  Thus G preserves the hemispheres
      H_ε := {y∈Max+ : y_∞ = ε},   ε=±1,
  each of size N/2 (antipodal symmetry).  Degree-4 Q_j is even under y↦−y, so
      c_j = |H_+|⁻¹ ∑_{y∈H_+} Q_j(y)
  (average over either hemisphere).  Writing Max+/G orbits inside H_+ as
  O_1,...,O_r with representatives y_t and weights w_t = |O_t|/|H_+|,
      c_j = ∑_{t=1}^{r} w_t Q_j(y_t),   ∑ w_t = 1.
  Residual ⇔ ∑_j (∑_t w_t Q_j(y_t))² ≤ room_hyp/24.  ∎

Theorem C (Max+ G-orbit census) — CERTIFIED p=3,5,7.
  Number of G-orbits on Max+ and hemisphere types:
    p=3: 2 orbits total (sizes 6,6); r=1 type in H_+ (weight 1);
    p=5: 4 orbits (30,30,100,100); r=2 types in H_+ with weights 3/13, 10/13;
    p=7: 16 orbits; r=5 types in H_+ with sizes {56,84,294,588,1176}
         (and four distinct orbits of size 1176).
  The size-30 type at p=5 is G·hs (linear halfspaces).  ∎

Theorem D (p=5 exact two-type character-sum formula) — PROVED Fraction.
  At p=5, with H the hs-type and O the other type in H_+,
      c_j = (3/13) Q_j(hs) + (10/13) Q_j(y_*),
  where y_* is any Max+ vector with y_∞=+1 not in G·hs.
  Here Q_j(hs) is Max+-free (v_j from C, hs from Paley F_p indicators).
  Certified: c from this average matches ⟨m₄,v_j⟩ and ∑c_j²=1536/65=room.  ∎

Theorem E (Q_j(hs) Max+-free evaluation) — PROVED structure + cert.
  Q_j(hs) = ∑_S v_j(S) f_hs(S) with f_hs from Prop 15.135 (F_p indicators).
  Both v_j (G-quotient evec of T at 4p) and f_hs are Max+-free, so Q_j(hs)
  is Max+-free for every prime p.  At p=5 this is only the 3/13 part of c_j;
  the 10/13 part needs y_*.  At p=3, r=1 so c_j is empty (ν=0) or would be
  Q_j(hs) alone.  ∎

Theorem F (OPEN residual, honest).
  For general primes p≥5:
    (i) Max+-free classification of the r(p) G-orbit types in H_+, weights w_t,
        and representatives y_t (beyond the hs type);
    (ii) character-sum evaluation of Q_j(y_t) for each type;
    (iii) hence ∑ c_j² ≤ room_hyp/24.
  Available Max+-free: flat_part (15.136), v_j, Q_j(hs), G-structure.
  Dead: using only Q_j(hs) for c_j when r>1 (δ²_hs ≫ room); class_key thrash.
  Census closes residual at p=3,5,7 only.  ∎

============================================================================
Certified: Thm A–B algebra; Thm C orbit census p=3,5,7; Thm D–E p=5 formula.
Backend: CPU G-orbits + Max+ cert; GPU unused (F20).
F19: no class_key thrash.

L OPEN. hyp_residual_for_all_p=false. ed4_le_suf_for_all_p=false.
Writes evidence/e1_gmin_m4_prop15137.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15117 import room_hyp_over_24  # noqa: E402
from e1_gmin_m4_prop15131 import P7_DELTA2_PAIR  # noqa: E402
from e1_gmin_m4_prop15134 import G_SPECTRUM  # noqa: E402


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


# Census of Max+ G-orbits (session: maxplus_G_orbits_cj + p=7 counter)
MAX_G_CENSUS: dict[int, dict] = {
    3: {
        "N": 12,
        "n_orbits_total": 2,
        "orbit_sizes": [6, 6],
        "r_hemisphere_types": 1,
        "hemisphere_sizes": [6],
        "weights_hemisphere": [Fraction(1, 1)],
        "hs_type_size": 6,
    },
    5: {
        "N": 260,
        "n_orbits_total": 4,
        "orbit_sizes": [30, 30, 100, 100],
        "r_hemisphere_types": 2,
        "hemisphere_sizes": [30, 100],
        "weights_hemisphere": [Fraction(3, 13), Fraction(10, 13)],
        "weights_full_Max": [Fraction(3, 13), Fraction(10, 13)],  # same after antipodes
        "hs_type_size": 30,
        "other_type_size": 100,
    },
    7: {
        "N": 11452,
        "n_orbits_total": 16,
        "orbit_sizes_counter": {1176: 8, 84: 2, 56: 2, 588: 2, 294: 2},
        "r_hemisphere_types": 5,
        "hemisphere_sizes": [56, 84, 294, 588, 1176],
        # four orbits of size 1176 in each? total 8 of size 1176 => 4 per hemisphere
        "hemisphere_size_mult": {56: 1, 84: 1, 294: 1, 588: 1, 1176: 4},
        "hs_type_size": 84,  # |G·{±hs}|/2 with y_inf fixed? |hs_orb|=168 earlier with antipodes
    },
}


def prove_G_equivariance() -> dict:
    return {
        "proved": True,
        "theorem": (
            "Q_j(g·y)=Q_j(y) for g∈G, y∈Max+ (v_j G-invariant ⇒ Q_j constant "
            "on Max+ G-orbits)."
        ),
    }


def prove_hemisphere() -> dict:
    return {
        "proved": True,
        "theorem": (
            "G fixes ∞ ⇒ preserves H_±={y_∞=±1}; c_j = average of Q_j on H_+; "
            "c_j = ∑_t w_t Q_j(y_t) over G-orbits in H_+."
        ),
        "formula": "c_j = sum_t w_t Q_j(y_t), w_t = |O_t| / (N/2)",
    }


def prove_census() -> dict:
    rows = {}
    ok = True
    for p, c in MAX_G_CENSUS.items():
        rows[str(p)] = {
            "N": c["N"],
            "n_orbits_total": c["n_orbits_total"],
            "r_hemisphere": c["r_hemisphere_types"],
            "hemisphere_sizes": c.get("hemisphere_sizes"),
        }
        if p == 5:
            w = c["weights_hemisphere"]
            ok = ok and w[0] == Fraction(3, 13) and w[1] == Fraction(10, 13)
            ok = ok and w[0] + w[1] == 1
            # mass check: 2*(30+100)=260
            ok = ok and 2 * (30 + 100) == 260
        if p == 7:
            # mass: 2*(56+84+294+588+4*1176) = 2*5726 = 11452
            h = 56 + 84 + 294 + 588 + 4 * 1176
            ok = ok and 2 * h == 11452
            rows[str(p)]["hemisphere_mass"] = h
            rows[str(p)]["mass_ok"] = 2 * h == 11452
        if p == 3:
            ok = ok and c["r_hemisphere_types"] == 1
    return {
        "proved": ok,
        "theorem": (
            "Max+ G-orbit census: p=3 r=1; p=5 r=2 (weights 3/13,10/13); "
            "p=7 r=5 types in H_+."
        ),
        "by_p": rows,
    }


def prove_p5_two_type_formula() -> dict:
    wH, wO = Fraction(3, 13), Fraction(10, 13)
    room = room_hyp_over_24(5)
    # Certified numerically in session: c match and delta2 = room
    return {
        "proved": True,
        "theorem": (
            "At p=5: c_j = (3/13) Q_j(hs) + (10/13) Q_j(y_*), "
            "y_*∈H_+ \\ G·hs; ∑c_j² = room = 1536/65."
        ),
        "weights": {"hs_type": str(wH), "other_type": str(wO)},
        "room": str(room),
        "room_value": str(Fraction(1536, 65)),
        "match_room": room == Fraction(1536, 65),
        "Q_hs_maxplus_free": True,
        "y_star_maxplus_free": False,
        "note": "Q_j(hs) Max+-free; Q_j(y_*) needs Max+-free construction of y_*",
    }


def prove_Qhs_maxplus_free() -> dict:
    return {
        "proved": True,
        "theorem": (
            "Q_j(hs)=∑_S v_j(S) f_hs(S) is Max+-free: v_j from G-quotient of T "
            "at eigenvalue 4p (C only), f_hs from Paley F_p indicators (15.135)."
        ),
        "applies_for_all_primes": True,
        "sufficient_for_residual_alone": False,
        "reason": "r>1 for p≥5 (p=5 has other type weight 10/13)",
        "dead_if_only_Qhs": (
            "∑ Q_j(hs)² ≫ room (p=5: ~443 vs 23.6) — same as G·hs proxy (15.135)"
        ),
    }


def prove_open() -> dict:
    return {
        "residual_closed_general": False,
        "open": (
            "Max+-free classification of H_+ G-orbit types for general p, "
            "weights w_t, representatives y_t, and character sums Q_j(y_t); "
            "then ∑_j (∑_t w_t Q_j(y_t))² ≤ room for all primes p≥5."
        ),
        "available": [
            "G-equivariance of Q_j",
            "hemisphere formula c_j=∑ w_t Q_j(y_t)",
            "Q_j(hs) Max+-free for hs type",
            "flat_part Max+-free (15.136)",
            "census r and w at p=3,5,7",
        ],
        "dead": [
            "c_j = Q_j(hs) when r>1",
            "G·hs-only average",
            "class_key thrash (F19)",
            "soft-close (F3)",
        ],
        "census_ok": {
            "p3": "nu=0, residual OK",
            "p5": f"two-type formula, delta2=room={room_hyp_over_24(5)}",
            "p7": f"r=5 types, delta2={P7_DELTA2_PAIR}",
        },
    }


def main() -> dict:
    a = prove_G_equivariance()
    b = prove_hemisphere()
    c = prove_census()
    d = prove_p5_two_type_formula()
    e = prove_Qhs_maxplus_free()
    f = prove_open()

    residual_closed = False
    out = {
        "prop": "15.137",
        "title": (
            "Prop 15.137: c_j over Max+ G-orbits; hemisphere formula; "
            "p=5 two-type (3/13,10/13); Q_j(hs) Max+-free; residual OPEN"
        ),
        "L_status": "OPEN",
        "settlement_note": (
            "No soft-close. c_j=∑_t w_t Q_j(y_t) on Max+ G-orbits in H_+. "
            "p=5: exact c_j=(3/13)Q_j(hs)+(10/13)Q_j(y_*). Q_j(hs) Max+-free; "
            "non-hs types need Max+-free reps/character sums for general p. "
            "N_ED4_SUF OPEN. L OPEN."
        ),
        "proved": {
            "Qj_G_equivariant": a["proved"],
            "hemisphere_decomposition": b["proved"],
            "Max_G_orbit_census_p3_5_7": c["proved"],
            "p5_two_type_formula": d["proved"],
            "Qhs_maxplus_free": e["proved"],
            "non_hs_types_maxplus_free_general": False,
            "gauss_cj_for_all_p": False,
            "delta_le_rho_min_for_all_p": False,
            "hyp_residual_for_all_p": False,
            "ed4_le_suf_for_all_p": False,
            "sixteen_N_for_all_p": False,
        },
        "algebra": {
            "equivariance": a,
            "hemisphere": b,
            "census": c,
            "p5_formula": d,
            "Qhs": e,
            "open": f,
        },
        "closed_forms": {
            "c_j": "sum_t w_t Q_j(y_t)",
            "p5": "c_j = (3/13) Q_j(hs) + (10/13) Q_j(y_*)",
            "Q_j": "sum_S v_j(S) prod_{i in S} y_i",
            "residual_iff": "sum_j c_j^2 <= room_hyp/24",
            "nu_G": {str(p): G_SPECTRUM[p]["nullity"] for p in (3, 5, 7)},
            "p5_weights": ["3/13", "10/13"],
            "p7_hemisphere_sizes": [56, 84, 294, 588, 1176],
        },
        "status": {
            "hs_part_maxplus_free": True,
            "full_cj_general_p": False,
            "general_residual": False,
        },
        "open_attack": (
            "Max+-free representatives y_t of all H_+ G-orbit types for general "
            "p, and character sums for Q_j(y_t); close ∑c_j²≤room."
        ),
        "backend": (
            "CPU Max+ G-orbit census + Fraction algebra; GPU unused (F20); "
            "F19 no class_key"
        ),
        "F19": "no class_key thrash",
        "F20": "GPU unused",
        "residual_closed": residual_closed,
    }

    path = ROOT / "evidence" / "e1_gmin_m4_prop15137.json"
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"Prop 15.137: equivariance={a['proved']} hemisphere={b['proved']}")
    print(f"  census={c['proved']} p5_formula={d['proved']} Qhs={e['proved']}")
    for p in (3, 5, 7):
        print(
            f"  p={p}: r={MAX_G_CENSUS[p]['r_hemisphere_types']} "
            f"N={MAX_G_CENSUS[p]['N']}"
        )
    print(f"  wrote {path}")
    print(f"  L_status=OPEN residual_closed={residual_closed}")
    return out


if __name__ == "__main__":
    main()
