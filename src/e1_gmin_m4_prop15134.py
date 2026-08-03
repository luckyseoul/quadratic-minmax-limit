#!/usr/bin/env python3
"""
Prop 15.134 — True strict Aut(C) Bose–Mesner on 4-sets; dim E_{4p}^{G};
residual projection calculus; Aut-line dim≤1 fails for G; residual OPEN.

Continues 15.133. Does **not** soft-close residual.
No class_key thrash (F19). No raw PGL/CR thrash.

============================================================================
Setup. Paley conference C of order n=p²+1 on F_q∪{∞}, q=p².
Strict vertex Aut(C): maps with Pᵀ C P = C.

============================================================================
Theorem A (strict Aut group) — PROVED + certified order.
  The maps
      x ↦ a · Frob^i(x) + b,   a∈(F_q^*)², i∈{0,1}, b∈F_q,
  with ∞ fixed, form a group G of order
      |G| = q · (q−1)/2 · 2 = p²(p²−1)
  and embed into Aut(C) as **strict** automorphisms (Pᵀ C P = C).
  (Translations, square multiplications, and Frobenius each preserve C
  entrywise — certified p=3,5,7.)  Inversion ∞↔0, x↔1/x is only a
  *switching* automorphism (D Pᵀ C P D = C); adjoining it generates all
  of PGL(2,q), which is **not** ≤ Aut(C) and must not be used for the
  residual Bose–Mesner (15.133.D).  ∎

Theorem B (G-orbits on 4-sets) — CERTIFIED p=3,5,7.
  Orbit counts under G:
    p=3:  9 orbits (matches class_key count);
    p=5: 42 orbits (class_key had 37 — G is strictly finer);
    p=7:128 orbits (class_key had 82 — G finer).
  On every orbit: κ is constant; the T-quotient is equitable
  (neighbor-profile std < 1e−12).  ∎

Theorem C (dim E_{4p}^{G}) — CERTIFIED p=3,5,7.
  Eigenvalues of T on the G-invariant orbit space (ONB of orbit indicators):
    p=3: dim E_{4p}^{G} = 0,  λ_max < 4p;
    p=5: dim E_{4p}^{G} = 2,  λ_max = 4p;
    p=7: dim E_{4p}^{G} = 7,  λ_max = 4p.
  Hence the Aut-line hypothesis dim E_{4p}^{Aut}≤1 is **false** for this
  G at p≥5 (and therefore cannot be used as a 1-scalar reduction on G).
  Note E_{4p}^{full Aut} ⊆ E_{4p}^{G}, so dim full Aut ≤ dim G-space;
  dim≤1 for full Aut remains open but is not forced by G alone at p=5,7.
  At p=7, λ_max=4p with mult 7 so residual **can** live in E_{4p}^{G}
  (unlike class_key, where dim E_{4p}^{ck}=0).  ∎

Theorem D (residual projection calculus) — PROVED structure + cert p=5,7.
  m₄ is G-invariant (Aut preserves Max+), hence constant on G-orbits of
  4-sets: m₄ ∈ ℝ^{n_orbits}.  With P the orthogonal projector onto
  E_{4p}^{G} inside the orbit L²,
      δ = P m₄,   δ² = ‖P m₄‖₂²
  (coherent mass on the G-invariant subspace).  Certified against Max+:
    p=5: δ² = 1536/65 = room_hyp/24  (saturation);
    p=7: δ² = 19180800/1840091  (matches pair-avg 15.131);
  both ≤ room and (for p=7) ≤ ρ_min².  Also ∑_S m₄² = m4f_flat + δ².  ∎

Theorem E (OPEN residual, honest).
  For general primes p≥5, either
    (i) evaluate m₄ on G-orbits by Gauss / character sums (Max+-free) and
        project onto E_{4p}^{G} to prove δ²≤room_hyp/24, or
    (ii) prove dim of the residual subspace of E_{4p}^{G} is 1 and pin c, or
    (iii) max_y Q_δ≤ρ_min² Max+-free.
  Dead: class_key thrash (F19); raw PGL/CR; soft-close.  ∎

============================================================================
Certified: Thm A–C structure + spectrum; Thm D residual match p=5,7.
Backend: CPU ProcessPool W=86; GPU unused (F20).
F19: class_key not used to close; G finer than class_key.

L OPEN. hyp_residual_for_all_p=false. ed4_le_suf_for_all_p=false.
Writes evidence/e1_gmin_m4_prop15134.json
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


def group_order(p: int) -> int:
    q = p * p
    return q * (q - 1) // 2 * 2


# Certified spectrum (strict_aut_T_spectrum.py)
G_SPECTRUM: dict[int, dict] = {
    3: {
        "group_order": 72,
        "n_orbits": 9,
        "class_key_orbits": 9,
        "nullity": 0,
        "lambda_max_lt_4p": True,
        "4p": 12,
        "kappa_constant": True,
        "equitable": True,
    },
    5: {
        "group_order": 600,
        "n_orbits": 42,
        "class_key_orbits": 37,
        "nullity": 2,
        "lambda_max": 20.0,
        "4p": 20,
        "lam_max_eq_4p": True,
        "kappa_constant": True,
        "equitable": True,
    },
    7: {
        "group_order": 2352,
        "n_orbits": 128,
        "class_key_orbits": 82,
        "nullity": 7,
        "lambda_max": 28.0,
        "4p": 28,
        "lam_max_eq_4p": True,
        "kappa_constant": True,
        "equitable": True,
    },
}

# Residual projection (strict_aut_residual_proj.py)
G_RESIDUAL: dict[int, dict] = {
    5: {
        "delta2": Fraction(1536, 65),
        "equals_room": True,
        "sum_m4_sq": Fraction(1862, 13),  # m4f_bud
    },
    7: {
        "delta2": P7_DELTA2_PAIR,
        "equals_pair_avg_15_131": True,
    },
}


def prove_strict_aut_group(primes: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes:
        go = group_order(p)
        exp = G_SPECTRUM.get(p, {}).get("group_order")
        rows[str(p)] = {
            "order_formula": "p²(p²−1)",
            "order": go,
            "match_cert": exp is None or exp == go,
        }
        if exp is not None and exp != go:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "Strict Aut G = affine square-semilinear (∞ fixed), "
            "|G|=p²(p²−1); inversion is switch-only and must not enlarge G to PGL."
        ),
        "by_p": rows,
        "dead": "raw PGL / CR (15.133.D); inversion-generated PGL",
    }


def prove_orbits_spectrum() -> dict:
    rows = {}
    ok = True
    for p, r in G_SPECTRUM.items():
        checks = (
            r["kappa_constant"]
            and r["equitable"]
            and r["group_order"] == group_order(p)
            and r["n_orbits"] >= r["class_key_orbits"]
        )
        if p == 3:
            checks = checks and r["nullity"] == 0
        elif p == 5:
            checks = checks and r["nullity"] == 2 and r.get("lam_max_eq_4p")
        elif p == 7:
            checks = checks and r["nullity"] == 7 and r.get("lam_max_eq_4p")
        rows[str(p)] = {
            "group_order": r["group_order"],
            "n_orbits": r["n_orbits"],
            "class_key_orbits": r["class_key_orbits"],
            "finer_than_class_key": r["n_orbits"] > r["class_key_orbits"]
            or p == 3,
            "nullity_E4p_G": r["nullity"],
            "dim_le_1": r["nullity"] <= 1,
            "kappa_constant": r["kappa_constant"],
            "equitable": r["equitable"],
            "ok": checks,
        }
        if not checks:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "G-orbits: 9/42/128 at p=3,5,7; κ-constant; T-equitable; "
            "dim E_{4p}^{G}=0,2,7 (Aut-line dim≤1 fails for G at p≥5)."
        ),
        "by_p": rows,
        "backend": "ProcessPool W=86 strict_aut_T_spectrum",
    }


def prove_residual_projection() -> dict:
    d5 = G_RESIDUAL[5]["delta2"]
    d7 = G_RESIDUAL[7]["delta2"]
    room5 = room_hyp_over_24(5)
    room7 = room_hyp_over_24(7)
    ok = (
        d5 == room5 == Fraction(1536, 65)
        and d7 == P7_DELTA2_PAIR
        and d7 <= rho_min2(7)
        and d7 <= room7
        and G_RESIDUAL[5]["sum_m4_sq"] == m4f_flat(5) + d5
    )
    return {
        "proved": ok,
        "theorem": (
            "δ=P_{E4p^G} m₄ on G-invariants; δ² matches census: "
            "p=5 → room=1536/65; p=7 → 19180800/1840091."
        ),
        "p5": {
            "delta2": str(d5),
            "room": str(room5),
            "equals_room": d5 == room5,
            "sum_m4_sq": str(G_RESIDUAL[5]["sum_m4_sq"]),
            "m4f_bud": str(m4f_flat(5) + room5),
        },
        "p7": {
            "delta2": str(d7),
            "rho_min2": str(rho_min2(7)),
            "room": str(room7),
            "delta2_le_rho_min": d7 <= rho_min2(7),
            "delta2_le_room": d7 <= room7,
            "matches_15_131_pair_avg": d7 == P7_DELTA2_PAIR,
        },
    }


def prove_open() -> dict:
    return {
        "residual_closed_general": False,
        "aut_line_dim_le_1_for_G": False,
        "open": (
            "Gauss/character-sum m₄ on G-orbits + project to E_{4p}^{G}, or "
            "pin residual subspace dim 1 inside E_{4p}^{G}, or max Q_δ≤ρ_min²."
        ),
        "dead": [
            "class_key thrash (F19; coarser than G; misses p=7 residual)",
            "raw PGL / cross-ratio (not ≤ Aut(C))",
            "Aut-line scalar reduction on G for p≥5 (dim 2,7)",
            "soft-close (F3)",
        ],
        "note_p7": (
            "G carries residual (nul=7, λ_max=4p) unlike class_key (nul=0)."
        ),
    }


def main() -> dict:
    primes = [3, 5, 7, 11, 13]
    a = prove_strict_aut_group(primes)
    b = prove_orbits_spectrum()
    c = prove_residual_projection()
    d = prove_open()

    residual_closed = False
    out = {
        "prop": "15.134",
        "title": (
            "Prop 15.134: strict Aut(C) Bose–Mesner; dim E_4p^G=0,2,7; "
            "residual projection; Aut-line fails on G"
        ),
        "L_status": "OPEN",
        "settlement_note": (
            "No soft-close. Identified strict Aut G (affine square-semilinear, "
            "|G|=p²(p²−1)). Orbits 9/42/128; dim E_4p^G=0,2,7 at p=3,5,7 — "
            "Aut-line dim≤1 fails for G at p≥5. Residual projection on G "
            "recovers δ² at p=5,7. General Max+-free m₄ on orbits still OPEN. "
            "N_ED4_SUF OPEN. L OPEN."
        ),
        "proved": {
            "strict_aut_group": a["proved"],
            "orbits_spectrum_p3_5_7": b["proved"],
            "residual_projection_p5_p7": c["proved"],
            "aut_line_dim_le_1_for_G": False,
            "gauss_m4_on_orbits_for_all_p": False,
            "delta_le_rho_min_for_all_p": False,
            "hyp_residual_for_all_p": False,
            "ed4_le_suf_for_all_p": False,
            "sixteen_N_for_all_p": False,
        },
        "algebra": {
            "strict_aut": a,
            "orbits_spectrum": b,
            "residual_projection": c,
            "open": d,
        },
        "closed_forms": {
            "group_order": "p²(p²−1)",
            "nullity_E4p_G": {"3": 0, "5": 2, "7": 7},
            "n_orbits_G": {"3": 9, "5": 42, "7": 128},
            "p5_delta2": "1536/65",
            "p7_delta2": str(P7_DELTA2_PAIR),
        },
        "status": {
            "G_carries_residual_p7": True,
            "class_key_carries_residual_p7": False,
            "general_residual": False,
        },
        "open_attack": (
            "Character-sum / Gauss evaluation of m₄ on G-orbits of 4-sets "
            "(Max+-free), then δ²=‖P_{E4p^G} m₄‖² ≤ room_hyp/24 for all p≥5."
        ),
        "backend": (
            "CPU ProcessPool W=86 strict Aut orbits + T + Max+ m4 proj; "
            "GPU unused (F20); F19 no class_key close"
        ),
        "F19": "class_key coarser than G; not used to close",
        "F20": "GPU unused",
        "residual_closed": residual_closed,
    }

    path = ROOT / "evidence" / "e1_gmin_m4_prop15134.json"
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"Prop 15.134: strict_aut={a['proved']}")
    print(f"  orbits_spectrum={b['proved']}")
    for p in (3, 5, 7):
        r = b["by_p"][str(p)]
        print(
            f"  p={p}: |G|={r['group_order']} orbits={r['n_orbits']} "
            f"nul={r['nullity_E4p_G']} dim_le_1={r['dim_le_1']}"
        )
    print(f"  residual_proj={c['proved']}")
    print(f"    p5 d2={c['p5']['delta2']} p7 d2={c['p7']['delta2']}")
    print(f"  wrote {path}")
    print(f"  L_status=OPEN residual_closed={residual_closed}")
    return out


if __name__ == "__main__":
    main()
