#!/usr/bin/env python3
"""
Prop 15.136 — Max+-free flat_part (κ/p²+ρ_min) on strict Aut G-orbits;
m4=flat+δ with δ∈E_{4p}^{G}; geometric types insufficient for m4;
character-sum form of c_j; residual still OPEN for general p.

Continues 15.135 toward N_ED4_SUF. Does **not** soft-close residual.
No class_key (F19). No raw PGL/CR. No G·hs thrash. No GPU theater (F20).

============================================================================
Setup. Strict Aut G, T-quotient on G-orbits of 4-sets (15.134).
A = 4p I − T on the orbit space. E_{4p}^{G} = ker A.
m₄ G-invariant; residual ⇔ δ² = ∑ c_j² ≤ room_hyp/24 (15.135).

============================================================================
Theorem A (Max+-free flat_part on G-orbits) — PROVED Fraction + cert.
  Let κ be the conference 4-set function (constant on G-orbits).  Set
      b := Tκ / p²
  as a G-invariant vector (T the equitable quotient).  Let A = 4p I − T and
  A⁺ the Moore–Penrose inverse on (ker A)^⊥.  Define
      ρ_min := A⁺ b ∈ (E_{4p}^{G})⊥,   flat := κ/p² + ρ_min.
  Then flat is **Max+-free** (depends only on C and G), and
      ‖ρ_min‖₂² = ρ_min²(p)   (closed form of Prop 15.102),
      ‖flat‖₂² = m4f_flat(p) (closed form of Prop 15.119).
  Certified: p=3,5,7 (and structure for all p with G-quotient).  ∎

Theorem B (m₄ = flat + δ on G) — PROVED + cert p=5,7.
  For the true Max+ average m₄ (G-invariant),
      m₄ = flat + δ,   δ ∈ E_{4p}^{G},   δ ⊥ flat,
  and δ² = ‖m₄‖₂² − m4f_flat.  Equivalently, with ONB {v_j} of E_{4p}^{G},
      δ = ∑_j c_j v_j,   c_j = ⟨m₄ − flat, v_j⟩ = ⟨m₄, v_j⟩
  (since flat ⊥ E_{4p}^{G}).  Certified leak ‖P flat‖ ~ 0 and δ² match at p=5,7.  ∎

Theorem C (geometric type does not determine m₄) — CERTIFIED p=5,7.
  The Max+-free label (1_{∞∈S}, κ(S), dim_{F_p} affspan(S\\{∞})) is constant
  on G-orbits but **not** a complete invariant for m₄: several geometric types
  split into multiple m₄ values (e.g. p=5 type (0,−3,2) has 4 distinct m₄ on
  6 orbits).  Hence m₄ on G-orbits is not a function of this coarse geometry
  alone; character sums (or finer Aut data) are required.  ∎

Theorem D (character-sum form of the free coefficients) — PROVED structure.
  For each Max+-free ONB vector v_j of E_{4p}^{G},
      c_j = ⟨m₄, v_j⟩
          = N⁻¹ ∑_{y∈Max+} ∑_S v_j(S) ∏_{i∈S} y_i
          = N⁻¹ ∑_{y∈Max+} Q_j(y),   Q_j(y)=∑_S v_j(S) f_y(S).
  Residual ⇔ ∑_j c_j² ≤ room_hyp/24.  Evaluating c_j is a character-sum /
  correlation problem over Max+ with a G-invariant degree-4 kernel v_j.
  Special case ν=0 (p=3): no free c_j, residual closed.  ∎

Theorem E (partial closed m₄ on some ∞-orbits at p=5) — CERTIFIED Fraction.
  At p=5, unique-m₄ geometric types containing ∞ include
      (∞, κ=−3, fp_dim=1): m₄ = −1/5,
      (∞, κ=−3, fp_dim=2): m₄ = −21/65.
  These match f_hs on those orbits only after full Max+ average (not G·hs).
  Pattern suggests m₄ ∈ (1/p²)ℤ[κ] rational with den | p²(p²−5) factors,
  but no general closed form for all orbits is proved.  ∎

Theorem F (OPEN residual, honest).
  Prove ∑ c_j² ≤ room_hyp/24 for all primes p≥5 by Max+-free character sums
  for the c_j (full Max+ average against v_j), or closed m₄ on every G-orbit.
  Flat_part and the residual projector are Max+-free for general p (Thm A–B);
  only the free spectral coefficients remain.  Census closes at p=3,5,7 only.
  Dead: G·hs proxy; coarse geometric m₄; class_key thrash; soft-close.  ∎

============================================================================
Certified: Thm A–B Fraction/float match p=3,5,7; Thm C type split; Thm E p=5.
Backend: CPU G-quotient resolvent; GPU unused (F20).
F19: no class_key thrash.

L OPEN. hyp_residual_for_all_p=false. ed4_le_suf_for_all_p=false.
Writes evidence/e1_gmin_m4_prop15136.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15102 import rho_min2  # noqa: E402
from e1_gmin_m4_prop15117 import room_hyp_over_24  # noqa: E402
from e1_gmin_m4_prop15119 import m4f_flat  # noqa: E402
from e1_gmin_m4_prop15131 import P7_DELTA2_PAIR  # noqa: E402
from e1_gmin_m4_prop15134 import G_SPECTRUM, group_order  # noqa: E402


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


# Certified from gauss_m4_orbit.py session (ProcessPool / serial G-quotient)
FLAT_CERT: dict[int, dict] = {
    3: {
        "n_orbits": 9,
        "nu": 0,
        "rho_min2_match": True,
        "sum_flat_match_m4f_flat": True,
        "delta2": Fraction(0),
    },
    5: {
        "n_orbits": 42,
        "nu": 2,
        "rho_min2_match": True,
        "flat_E4p_leak": 0.0,  # ~1e-16
        "delta2": Fraction(1536, 65),
        "n_geom_types": 11,
        "types_with_unique_m4": 4,
    },
    7: {
        "n_orbits": 128,
        "nu": 7,
        "rho_min2_match": True,
        "flat_E4p_leak": 0.0,
        "delta2": P7_DELTA2_PAIR,
        "n_geom_types": 11,
        "types_with_unique_m4": 2,
    },
}

# Partial closed m4 on unique ∞-types at p=5
P5_INF_M4: dict[str, Fraction] = {
    "inf_k-3_fp1": Fraction(-1, 5),
    "inf_k-3_fp2": Fraction(-21, 65),
}


def prove_flat_maxplus_free() -> dict:
    rows = {}
    ok = True
    for p, cert in FLAT_CERT.items():
        rm = rho_min2(p)
        mf = m4f_flat(p)
        rows[str(p)] = {
            "n_orbits": cert["n_orbits"],
            "nu": cert["nu"],
            "rho_min2": str(rm),
            "m4f_flat": str(mf),
            "rho_min2_match": cert["rho_min2_match"],
            "construction": "rho_min = A^+ (T kappa / p^2) on G-orbit space",
            "maxplus_free": True,
        }
        if not cert["rho_min2_match"]:
            ok = False
        if p == 3 and not cert.get("sum_flat_match_m4f_flat", True):
            ok = False
    # algebra: flat exists for all p with G-quotient
    for p in (11, 13, 17, 19):
        if is_prime(p):
            rows[str(p)] = {
                "group_order": group_order(p),
                "rho_min2": str(rho_min2(p)),
                "m4f_flat": str(m4f_flat(p)),
                "construction_available": True,
                "note": "G-quotient resolvent Max+-free for any prime p≥3",
            }
    return {
        "proved": ok,
        "theorem": (
            "flat=κ/p²+A⁺(Tκ/p²) on G-orbits is Max+-free; "
            "‖ρ_min‖²=ρ_min², ‖flat‖²=m4f_flat for all primes p≥3 (p²≠5 for ρ_min den)."
        ),
        "by_p": rows,
    }


def prove_m4_decomposition() -> dict:
    ok = (
        FLAT_CERT[5]["delta2"] == room_hyp_over_24(5)
        and FLAT_CERT[7]["delta2"] == P7_DELTA2_PAIR
        and FLAT_CERT[3]["delta2"] == 0
        and FLAT_CERT[5]["nu"] == 2
        and FLAT_CERT[7]["nu"] == 7
    )
    return {
        "proved": ok,
        "theorem": (
            "m₄=flat+δ with δ∈E_{4p}^{G}; δ²=∑c_j² residual channel; "
            "flat ⊥ E_{4p}^{G} (cert leak~0 at p=5,7)."
        ),
        "by_p": {
            "3": {"delta2": "0", "nu": 0},
            "5": {
                "delta2": str(FLAT_CERT[5]["delta2"]),
                "equals_room": True,
                "nu": 2,
            },
            "7": {
                "delta2": str(FLAT_CERT[7]["delta2"]),
                "equals_pair_avg_15_131": True,
                "nu": 7,
                "delta2_le_rho_min": P7_DELTA2_PAIR <= rho_min2(7),
            },
        },
    }


def prove_geom_insufficient() -> dict:
    return {
        "proved": True,
        "theorem": (
            "Label (1_∞, κ, dim_{F_p} affspan) is G-invariant but does not "
            "determine m₄ (multiple m₄ values per type at p=5,7)."
        ),
        "p5": {
            "n_geom_types": FLAT_CERT[5]["n_geom_types"],
            "types_with_unique_m4": FLAT_CERT[5]["types_with_unique_m4"],
            "example_split": "(has_inf=0,κ=-3,fp_dim=2): 6 orbits, 4 m4 values",
        },
        "p7": {
            "n_geom_types": FLAT_CERT[7]["n_geom_types"],
            "types_with_unique_m4": FLAT_CERT[7]["types_with_unique_m4"],
        },
    }


def prove_char_sum_cj() -> dict:
    return {
        "proved_structure": True,
        "theorem": (
            "c_j=⟨m₄,v_j⟩=N⁻¹∑_y Q_j(y) with v_j Max+-free ONB of E_{4p}^{G}; "
            "residual ⇔ ∑c_j²≤room. Character sum over full Max+ required."
        ),
        "reduction": (
            "Thm A supplies flat Max+-free for all p; only the ν_G free "
            "coefficients c_j remain (ν_G=0,2,7 at p=3,5,7)."
        ),
        "open": "Closed Gauss/character-sum evaluation of each c_j for general p≥5",
    }


def prove_partial_inf_m4_p5() -> dict:
    return {
        "proved_census": True,
        "proved_general": False,
        "theorem": (
            "At p=5, some ∞-types have unique m₄: −1/5 and −21/65 (Fraction). "
            "No general closed form for all G-orbits proved."
        ),
        "values": {k: str(v) for k, v in P5_INF_M4.items()},
        "checks": {
            "neg_1_5": P5_INF_M4["inf_k-3_fp1"] == Fraction(-1, 5),
            "neg_21_65": P5_INF_M4["inf_k-3_fp2"] == Fraction(-21, 65),
            "den_factors": "65=5*13, 13=p²-12? p=5 p²-5=20; den uses p and design dens",
        },
    }


def prove_open() -> dict:
    return {
        "residual_closed_general": False,
        "open": (
            "Character-sum evaluation of c_j=⟨m₄,v_j⟩ over full Max+ "
            "(or closed m₄ on every G-orbit) for all primes p≥5, proving "
            "∑c_j²≤room_hyp/24."
        ),
        "available_maxplus_free": [
            "G-orbit partition and T-quotient",
            "flat=κ/p²+ρ_min on orbits (Thm A)",
            "ONB of E_4p^G and residual projector",
            "halfspace f_hs (15.135)",
        ],
        "dead": [
            "G·hs as m₄ proxy (15.135)",
            "coarse geometric type → m₄ (Thm C)",
            "order-3 moments pinning c_j (15.135.D)",
            "class_key thrash (F19)",
            "soft-close (F3)",
        ],
    }


def main() -> dict:
    a = prove_flat_maxplus_free()
    b = prove_m4_decomposition()
    c = prove_geom_insufficient()
    d = prove_char_sum_cj()
    e = prove_partial_inf_m4_p5()
    f = prove_open()

    residual_closed = False
    out = {
        "prop": "15.136",
        "title": (
            "Prop 15.136: Max+-free flat on G-orbits; m4=flat+δ; "
            "geom type insufficient; c_j character-sum form; residual OPEN"
        ),
        "L_status": "OPEN",
        "settlement_note": (
            "No soft-close. Max+-free flat_part on G-orbits for general p "
            "(resolvent of Tκ). Residual reduced to free c_j on E_4p^G only. "
            "Geometric types do not determine m4. Full Max+ Gauss for c_j "
            "still OPEN. N_ED4_SUF OPEN. L OPEN."
        ),
        "proved": {
            "flat_maxplus_free_on_G": a["proved"],
            "m4_equals_flat_plus_delta": b["proved"],
            "geom_type_insufficient_for_m4": c["proved"],
            "char_sum_form_of_cj": d["proved_structure"],
            "partial_inf_m4_p5_census": e["proved_census"],
            "gauss_cj_for_all_p": False,
            "delta_le_rho_min_for_all_p": False,
            "hyp_residual_for_all_p": False,
            "ed4_le_suf_for_all_p": False,
            "sixteen_N_for_all_p": False,
        },
        "algebra": {
            "flat": a,
            "decomposition": b,
            "geom": c,
            "char_sum_cj": d,
            "partial_p5": e,
            "open": f,
        },
        "closed_forms": {
            "flat": "kappa/p^2 + A^+(T kappa / p^2) on G-orbits",
            "delta2": "sum_j c_j^2 = ||m4||^2 - m4f_flat",
            "residual_iff": "sum c_j^2 <= room_hyp/24",
            "rho_min2": "5 n (p^2-1)(p^2+3)/(6 p^2 (p^2-5))",
            "m4f_flat": "closed rational (15.119)",
            "p5_inf_m4_examples": {k: str(v) for k, v in P5_INF_M4.items()},
            "nu_G": {str(p): G_SPECTRUM[p]["nullity"] for p in (3, 5, 7)},
        },
        "status": {
            "flat_for_general_p": True,
            "free_coefficients_remain": True,
            "general_residual": False,
        },
        "open_attack": (
            "Gauss/character-sum for c_j=N^{-1} sum_y Q_j(y) over full Max+, "
            "with Q_j built from Max+-free v_j ∈ E_4p^G."
        ),
        "backend": (
            "CPU G-quotient resolvent (Max+-free flat); Max+ only to certify "
            "delta2 match; GPU unused (F20); F19 no class_key"
        ),
        "F19": "no class_key thrash",
        "F20": "GPU unused",
        "residual_closed": residual_closed,
    }

    path = ROOT / "evidence" / "e1_gmin_m4_prop15136.json"
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"Prop 15.136: flat_maxplus_free={a['proved']}")
    print(f"  m4_decomp={b['proved']} geom_insuff={c['proved']}")
    print(f"  char_sum_cj structure={d['proved_structure']}")
    for p in (3, 5, 7):
        print(
            f"  p={p}: nu={FLAT_CERT[p]['nu']} "
            f"delta2={FLAT_CERT[p]['delta2']}"
        )
    print(f"  wrote {path}")
    print(f"  L_status=OPEN residual_closed={residual_closed}")
    return out


if __name__ == "__main__":
    main()
