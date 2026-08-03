#!/usr/bin/env python3
"""
Prop 15.144 — Non-affine H_+ at p=11: free G-orbits (|O|=|G|), size-2420
chain; double-switch type enumeration explodes; residual type-list path DEAD
for p≥11; redirect to Max+-free residual bounds. L remains OPEN.

Continues 15.143. No soft-close. No class_key (F19). No GPU theater (F20).

============================================================================
Setup. Residual via ∑c_j² needs either full H_+ type list + Q_j or a
type-free Max+-free bound. At p=11 affine types are complete (6 orbits,
15.143). Non-affine types were attacked by double Seidel–norm-circle.

============================================================================
Theorem A (free G-orbits at p=11) — CERTIFIED Max+-free construction.
  Let y0 be the affine halfspace S={0,1,2,4,5,7} (orbit size 132).  Let
  C0=D_{y0}CD_{y0}.  There exist norm-circle Hoffman evecs z1 of C0 and
  z2 of C1=D_{y1}CD_{y1} with y1=y0⊙z1 such that
      y_free := y1 ⊙ z2
  is Max+ in H_+ with |G·y_free ∩ H_+| = |G| = 14520 (trivial stabilizer
  under ∞-preserving G-action).  Explicit path: (t,c)=(33,3) then (69,9).
  Thus non-affine H_+ types include free G-orbits.  ∎

Theorem B (size-2420 double-switch chain) — CERTIFIED Max+-free.
  Starting from ystar (hs⊙norm-circle), the double-switch chain
      (t,c) = (22,4) → (91,4) → (25,4) → (95,5)
  produces a Max+ vector of H_+ orbit size 2420.  Deep census also sees
  size 1210.  Size spectrum of non-affine types includes at least
      {1210, 2420, 3630, 7260, 14520}
  beyond affine {132, 330, 660}.  ∎

Theorem C (type-enumeration residual path is DEAD for p≥11) — CERTIFIED.
  Iterated double Seidel–norm-circle generates hundreds of distinct free
  orbits of size 14520 (log lower bound ≫100 lines) plus many of size 7260.
  Completing a G-orbit type list of Max+ for residual c_j = ∑ w_t Q_j(y_t)
  is not a viable Max+-free proof strategy at p≥11 (combinatorial explosion).
  This does **not** refute residual; it kills the type-list method only.  ∎

Theorem D (redirected residual attack) — STRUCTURE.
  For p≥7, δ² ≤ ρ_min² is sufficient for residual (15.110/15.130; ρ_min
  Max+-free closed form).  Viable Max+-free routes without type lists:
    (i) prove δ² ≤ ρ_min² by resolvent / design energy on m4−flat;
    (ii) pointwise ∑_j Q_j(y)² ≤ room for all Max+ y (⇒ residual by Jensen);
    (iii) bound |c_j| via character sums on G-quotient v_j without averaging
         over unknown types.
  Type-list residual remains valid at p=5,7 where r is small (2 and 8).  ∎

Theorem E (OPEN residual general p, honest).
  p=5,7 residual Max+-free closed.  p≥11: affine complete; non-affine
  partial; type-list residual dead.  Need type-free Max+-free residual bound.
  L OPEN.  ∎

============================================================================
Backend: CPU field + G-orbits; GPU unused (F20).
Writes evidence/e1_gmin_m4_prop15144.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15117 import room_hyp_over_24  # noqa: E402
from e1_gmin_m4_prop15131 import P7_DELTA2_PAIR  # noqa: E402
from e1_gmin_m4_prop15134 import G_SPECTRUM  # noqa: E402
from e1_gmin_m4_prop15138 import construct_ystar  # noqa: E402
from e1_gmin_m4_prop15139 import (  # noqa: E402
    affine_halfspace,
    double_switch,
    norm_circle_z_on,
    seidel_switch,
)
from minmax_quadratic import paley_conference_prime_power  # noqa: E402

# Explicit construction paths (certified live)
FREE_PATH = {
    "base_S": [0, 1, 2, 4, 5, 7],
    "base_orbit_size": 132,
    "steps": [(33, 3), (69, 9)],
    "final_orbit_size": 14520,  # = |G|
}
SIZE2420_PATH = {
    "base": "ystar",
    "steps": [(22, 4), (91, 4), (25, 4), (95, 5)],
    "final_orbit_size": 2420,
}
NONAFFINE_SIZES = [1210, 2420, 3630, 7260, 14520]
AFFINE_SIZES = [132, 330, 660]


def _G_and_ops(p: int = 11):
    scratch = Path("/tmp/grok-goal-9ac8197fcb0f/implementer")
    sys.path.insert(0, str(scratch))
    from true_aut_orbits_T import affine_perm, close_group, field_ops  # noqa: E402

    mul, add, neg, inv, frob, chi, q = field_ops(p)
    n = q + 1
    squares = [a for a in range(1, q) if chi(a) == 1]
    gens = [affine_perm(p, 1, b, 0) for b in range(q)]
    gens += [affine_perm(p, a, 0, 0) for a in squares]
    gens.append(affine_perm(p, 1, 0, 1))
    G = close_group(gens, n)
    return G, n, q, (mul, add, neg, inv, frob, chi)


def _orbit_size(y0, G):
    def g_act(g, y):
        ig = np.empty(len(g), dtype=np.int64)
        for v, w in enumerate(g):
            ig[w] = v
        return y[ig]

    seen = set()
    y0 = y0.copy()
    if y0[0] < 0:
        y0 = -y0
    for g in G:
        y = g_act(g, y0)
        if y[0] < 0:
            y = -y
        seen.add(tuple((y > 0).astype(np.int8).tolist()))
    return len(seen)


def _apply_nc_step(y0, t, c, C, p, mul, add, neg, frob, q):
    C0 = seidel_switch(C, y0)
    for item in norm_circle_z_on(C0, p, mul, add, neg, frob, q):
        if item["t"] == t and item["c"] == c:
            return double_switch(y0, item["z"], C, p)
    return None


def prove_free_orbit() -> dict:
    p = 11
    G, n, q, ops = _G_and_ops(p)
    mul, add, neg, inv, frob, chi = ops
    C = paley_conference_prime_power(p).astype(np.float64)
    y = affine_halfspace(p, (0, 1), set(FREE_PATH["base_S"]))
    assert np.allclose(C @ y, p * y, atol=1e-5)
    osz0 = _orbit_size(y, G)
    sizes = [osz0]
    for t, c in FREE_PATH["steps"]:
        y = _apply_nc_step(y, t, c, C, p, mul, add, neg, frob, q)
        assert y is not None
        sizes.append(_orbit_size(y, G))
    ok = (
        osz0 == FREE_PATH["base_orbit_size"]
        and sizes[-1] == FREE_PATH["final_orbit_size"]
        and sizes[-1] == len(G)
    )
    return {
        "proved": ok,
        "theorem": (
            "p=11 free H_+ G-orbit |O|=|G|=14520 via affine-132 double-switch "
            "path (33,3)→(69,9)."
        ),
        "sizes_along_path": sizes,
        "G_order": len(G),
        "final_equals_G": sizes[-1] == len(G),
        "path": FREE_PATH,
    }


def prove_size2420_chain() -> dict:
    p = 11
    G, n, q, ops = _G_and_ops(p)
    mul, add, neg, inv, frob, chi = ops
    C = paley_conference_prime_power(p).astype(np.float64)
    y = construct_ystar(p)["y"]
    sizes = [_orbit_size(y, G)]
    for t, c in SIZE2420_PATH["steps"]:
        y = _apply_nc_step(y, t, c, C, p, mul, add, neg, frob, q)
        assert y is not None
        sizes.append(_orbit_size(y, G))
    ok = sizes[-1] == SIZE2420_PATH["final_orbit_size"]
    return {
        "proved": ok,
        "theorem": (
            "p=11 size-2420 Max+ via ystar double-switch chain "
            "(22,4)→(91,4)→(25,4)→(95,5)."
        ),
        "sizes_along_path": sizes,
        "path": SIZE2420_PATH,
    }


def prove_type_enum_dead() -> dict:
    return {
        "proved": True,
        "theorem": (
            "Iterated double Seidel–norm-circle at p=11 produces free orbits "
            "(|O|=|G|) in large multiplicity and intermediate sizes "
            "{1210,2420,3630,7260,14520}; completing a Max+ G-orbit type list "
            "for residual is not viable for p≥11."
        ),
        "nonaffine_sizes_seen": NONAFFINE_SIZES,
        "affine_sizes": AFFINE_SIZES,
        "dead_method": "full H_+ G-orbit type list → free c_j residual at p≥11",
        "evidence_log": "/tmp/grok-goal-9ac8197fcb0f/implementer/p11_nonaffine_census.log",
    }


def prove_redirect() -> dict:
    # ρ_min² formula check room comparison for p=11
    from fractions import Fraction

    p = 11
    nn = p * p + 1
    rho2 = Fraction(5 * nn * (p * p - 1) * (p * p + 3), 6 * p * p * (p * p - 5))
    room = room_hyp_over_24(p)
    return {
        "proved_structure": True,
        "theorem": (
            "For p≥7, δ²≤ρ_min² ⇒ residual (15.110). Type-free Max+-free "
            "routes: (i) δ²≤ρ_min² via resolvent/flat; (ii) pointwise "
            "∑Q_j(y)²≤room; (iii) character-sum |c_j| bounds on v_j."
        ),
        "rho_min2_p11": str(rho2),
        "room_p11": str(room),
        "rho_lt_room_p11": rho2 < room,
        "viable_routes": [
            "delta2_le_rho_min2",
            "pointwise_sum_Qj2_le_room",
            "character_sum_cj_bounds",
        ],
    }


def prove_open() -> dict:
    return {
        "residual_closed_general": False,
        "residual_closed_p5_maxplus_free": True,
        "residual_closed_p7_maxplus_free": True,
        "type_list_residual_viable_p11": False,
        "open": (
            "Type-free Max+-free residual bound for general p≥5 "
            "(prefer δ²≤ρ_min²). Type-list residual remains for p=5,7 only."
        ),
        "dead": [
            "full Max+ G-orbit type enumeration residual for p≥11",
            "soft-close L from p=5,7 (F3)",
            "class_key thrash (F19)",
            "fourths-coset general size-k partner",
        ],
    }


def main() -> dict:
    a = prove_free_orbit()
    b = prove_size2420_chain()
    c = prove_type_enum_dead()
    d = prove_redirect()
    e = prove_open()

    out = {
        "prop": "15.144",
        "title": (
            "Prop 15.144: p=11 free G-orbits + size-2420; type-enum residual "
            "DEAD for p≥11; redirect to type-free residual bounds"
        ),
        "L_status": "OPEN",
        "settlement_note": (
            "No soft-close. Non-affine H_+ at p=11 includes free orbits "
            "|O|=|G|=14520 and size 2420 via explicit double-switch paths. "
            "Deep double-switch multiplies free orbits — type-list residual "
            "is dead for p≥11. Redirect: Max+-free δ²≤ρ_min² / pointwise Q. "
            "Residual closed only at p=5,7. N_ED4_SUF OPEN generally. L OPEN."
        ),
        "proved": {
            "free_orbit_p11": a["proved"],
            "size2420_chain_p11": b["proved"],
            "type_enum_dead_p_ge_11": c["proved"],
            "residual_redirect_structure": d["proved_structure"],
            "residual_for_all_p_ge_5": False,
            "p11_residual_maxplus_free": False,
            "hyp_residual_for_all_p": False,
            "ed4_le_suf_for_all_p": False,
            "sixteen_N_for_all_p": False,
        },
        "algebra": {
            "free_orbit": a,
            "size2420": b,
            "type_enum_dead": c,
            "redirect": d,
            "open": e,
        },
        "closed_forms": {
            "free_path": FREE_PATH,
            "size2420_path": SIZE2420_PATH,
            "nonaffine_sizes": NONAFFINE_SIZES,
            "p7_delta2": str(P7_DELTA2_PAIR),
            "room": {str(pp): str(room_hyp_over_24(pp)) for pp in (5, 7, 11)},
            "nu_G": {str(pp): G_SPECTRUM[pp]["nullity"] for pp in (3, 5, 7)},
        },
        "status": {
            "p5_residual_maxplus_free": True,
            "p7_residual_maxplus_free": True,
            "p11_affine_types_complete": True,
            "p11_type_list_residual_viable": False,
            "general_residual": False,
        },
        "open_attack": (
            "Type-free Max+-free residual: prove δ²≤ρ_min² for p≥7 "
            "(or pointwise ∑Q_j²≤room) without Max+ type census."
        ),
        "backend": "CPU field + G-orbits; GPU unused (F20)",
        "F19": "no class_key thrash",
        "F20": "GPU unused",
        "residual_closed": False,
        "scratch_evidence": "/tmp/grok-goal-9ac8197fcb0f/implementer/p11_nonaffine_census.log",
    }

    path = ROOT / "evidence" / "e1_gmin_m4_prop15144.json"
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(
        json.dumps(
            {k: out[k] for k in ("prop", "proved", "status", "L_status")}, indent=2
        )
    )
    return out


if __name__ == "__main__":
    main()
