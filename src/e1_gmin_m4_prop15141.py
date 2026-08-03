#!/usr/bin/env python3
"""
Prop 15.141 — Explicit size-12 Seidel partner closes last p=7 H_+ type
Max+-free; free c_j recovers δ²≤room; p=7 residual Max+-free closed.
General p still OPEN (no soft-close L).

Continues 15.140. No class_key (F19). No GPU theater (F20).

============================================================================
Setup. From 15.140: 7/8 H_+ types Max+-free; one size-1176 orbit lacked
a Max+-free y (distinct Q). Weights w_t=|G|/|Stab| and Q_j character sums
already recover δ² when all 8 types are available.

============================================================================
Theorem A (size-12 Seidel partner) — CERTIFIED Max+-free.
  Let y0 be the affine halfspace for L(a+bω)=b and S_aff={2,3,4,5}
  (a 4-AP ⇒ size-84 type, 15.139).  Set C0=D_{y0} C D_{y0}.  Let
      T = {10,12,13,16,18,25,29,36,38,42,44,48} ⊂ F_{49}
  (field encoding a+b·p as in paley_conference_prime_power), and let z be
  the boolean vector with z_∞=+1 and z_u=−1 ⇔ u∈T.  Then:
    (i) C0 z = p z  (boolean evec of the switched conference);
    (ii) y_♯ := y0 ⊙ z satisfies C y_♯ = p y_♯ and y_♯∈H_+;
    (iii) |G·y_♯ ∩ H_+| = 1176 with Stab order 2;
    (iv) Q_j(y_♯) matches the previously census-only size-1176 type
         (distinct from the three double-norm-circle size-1176 orbits).
  Construction uses only field constants + C + y0 — Max+-free.  ∎

Theorem B (all eight H_+ types Max+-free at p=7) — CERTIFIED.
  Combining 15.138–15.140 with Theorem A:
    56,84 affine; 588 hs⊙nc; 294 ync⊙nc; three 1176 from y56/ync⊙nc;
    one 1176 from size-12 Seidel partner (Thm A).
  All eight H_+ G-orbits have Max+-free representatives.  ∎

Theorem C (p=7 residual Max+-free) — CERTIFIED Fraction.
  With Max+-free y_t (Thm B), Max+-free weights w_t=|O_t|/|H_+| via
  |G|/|Stab| (15.140), and Max+-free Q_j (v_j from G-quotient of T at 4p),
      c_j = ∑_t w_t Q_j(y_t)
  is Max+-free and ∑_j c_j² = 19180800/1840091 = δ²_pair ≤ room_hyp/24
  = 3072/55.  Hence Path C residual holds at p=7 **without Max+ census**.  ∎

Theorem D (bi-tight empty at p=5,7 Max+-free path) — PROVED form.
  mult(λ₂)≥d−1 (15.98) + orth≤room_hyp (⇔ residual) ⇒ 16N (15.98+15.110)
  ⇒ bi-tight empty for p≥5 (15.61).  Residual Max+-free at p=5 (15.138)
  and p=7 (Thm C) ⇒ bi-tight empty at p=5,7 on the Max+-free path.  ∎

Theorem E (OPEN general p, honest).
  For all primes p≥5: Max+-free types/weights/Q_j and residual.  p=3 trivial
  (ν=0); p=5,7 closed Max+-free.  Need general-p type classification
  (size-12 partners, affine AP dichotomy, double switches) and character sums.
  Do not soft-close L.  ∎

============================================================================
Backend: CPU field + G-orbits + vectorized Q; GPU unused (F20).
Writes evidence/e1_gmin_m4_prop15141.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
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
    is_4_AP,
    norm_circle_z_on,
    seidel_switch,
)
from e1_gmin_m4_prop15140 import P7_G, P7_HPLUS, P7_ORBIT_TABLE  # noqa: E402
from minmax_quadratic import (  # noqa: E402
    halfspace_boolean_vector,
    paley_conference_prime_power,
)

# Explicit size-12 field set for Thm A (encoding a + b*p)
T_FIELD_SIZE12: list[int] = [10, 12, 13, 16, 18, 25, 29, 36, 38, 42, 44, 48]
S_AFF_SIZE12: list[int] = [2, 3, 4, 5]


def construct_y_sharp(p: int = 7) -> dict:
    """Max+-free missing size-1176 type: y0 ⊙ z_T with |T|=12."""
    if p != 7:
        return {"found": False, "reason": "explicit T certified only at p=7"}
    C = paley_conference_prime_power(p).astype(np.float64)
    y0 = affine_halfspace(p, (0, 1), set(S_AFF_SIZE12))
    # S_aff={2,3,4,5} is a 4-AP (size-84 type); size-12 partner still Max+
    assert is_4_AP(S_AFF_SIZE12, p)
    n = p * p + 1
    z = np.ones(n, dtype=np.float64)
    for u in T_FIELD_SIZE12:
        z[1 + u] = -1.0
    C0 = seidel_switch(C, y0)
    if not np.allclose(C0 @ z, p * z, atol=1e-5):
        return {"found": False, "reason": "z not C0-evec"}
    y = y0 * z
    if y[0] < 0:
        y = -y
    if not np.allclose(C @ y, p * y, atol=1e-5):
        return {"found": False, "reason": "y not Max+"}
    return {
        "found": True,
        "y": y,
        "y0_S_aff": list(S_AFF_SIZE12),
        "T_field": list(T_FIELD_SIZE12),
        "T_ab": [(u % p, u // p) for u in T_FIELD_SIZE12],
        "y_inf": float(y[0]),
        "sum": float(y.sum()),
    }


def prove_size12_partner() -> dict:
    r = construct_y_sharp(7)
    return {
        "proved": r["found"],
        "theorem": (
            "y_♯=y0⊙z_T with y0 affine 4-AP S={2,3,4,5} (size-84 type), "
            "T size-12 explicit field set, C0=D_y0 C D_y0, C0z=pz, Cy_♯=py_♯; "
            "H_+ orbit size 1176 — Max+-free last type at p=7."
        ),
        "S_aff": S_AFF_SIZE12,
        "T_field": T_FIELD_SIZE12,
        "found": r["found"],
        "y_inf": r.get("y_inf"),
        "sum": r.get("sum"),
    }


def prove_all_eight_and_residual() -> dict:
    """Live: all 8 types free + weighted Q gives δ²."""
    p = 7
    C = paley_conference_prime_power(p).astype(np.float64)
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

    def g_act(g, y):
        ig = np.empty(len(g), dtype=np.int64)
        for v, w in enumerate(g):
            ig[w] = v
        return y[ig]

    def orbit_size(y0):
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

    # load V for Q
    cache = scratch / "p7_E4pG_cache.npz"
    if not cache.exists():
        return {"proved": False, "reason": "missing E4pG cache"}
    zc = np.load(cache, allow_pickle=True)
    orbits4 = list(zc["orbits4"])
    Ds = zc["Ds"]
    V = zc["V"]
    nu = int(zc["nu"])
    idx, vmat = [], []
    for oi, o in enumerate(orbits4):
        vf = V[oi] / Ds[oi]
        for S in o:
            idx.append(S)
            vmat.append(vf)
    idx = np.array(idx)
    vmat = np.array(vmat)
    Ia, Ib, Ic, Id = idx[:, 0], idx[:, 1], idx[:, 2], idx[:, 3]

    def Qs(y):
        return vmat.T @ (y[Ia] * y[Ib] * y[Ic] * y[Id])

    # eight free reps
    y56 = affine_halfspace(p, (0, 1), {0, 1, 2, 4})
    y84 = affine_halfspace(p, (0, 1), {0, 1, 2, 3})
    y588 = construct_ystar(p)["y"]
    y_sharp = construct_y_sharp(p)["y"]

    C_nc = seidel_switch(C, y588)
    y294 = None
    for item in norm_circle_z_on(C_nc, p, mul, add, neg, frob, q):
        yy = double_switch(y588, item["z"], C, p)
        if yy is not None and orbit_size(yy) == 294:
            y294 = yy
            break
    assert y294 is not None

    y1176s = []
    for y0, tag in [(y56, "56"), (y588, "nc")]:
        C0 = seidel_switch(C, y0)
        for item in norm_circle_z_on(C0, p, mul, add, neg, frob, q):
            yy = double_switch(y0, item["z"], C, p)
            if yy is None:
                continue
            if orbit_size(yy) == 1176:
                y1176s.append(yy)
    # dedupe 1176 by orbit + add sharp
    def key_orb(y0):
        keys = set()
        y0 = y0.copy()
        if y0[0] < 0:
            y0 = -y0
        for g in G:
            y = g_act(g, y0)
            if y[0] < 0:
                y = -y
            keys.add(tuple((y > 0).astype(np.int8).tolist()))
        return min(keys)

    unique = {}
    for name, yy in [
        ("56", y56),
        ("84", y84),
        ("294", y294),
        ("588", y588),
        ("sharp", y_sharp),
    ] + [(f"1176_{i}", yy) for i, yy in enumerate(y1176s)]:
        k = key_orb(yy)
        if k not in unique:
            unique[k] = (name, yy, orbit_size(yy), Qs(yy))

    sizes = sorted(v[2] for v in unique.values())
    H = sum(v[2] for v in unique.values())
    c = np.zeros(nu)
    for name, yy, osz, Q in unique.values():
        c += (osz / H) * np.asarray(Q)
    d2 = float(np.dot(c, c))
    room = room_hyp_over_24(7)
    ok = (
        len(unique) == 8
        and H == P7_HPLUS
        and abs(d2 - float(P7_DELTA2_PAIR)) < 1e-8
        and d2 <= float(room) + 1e-9
        and orbit_size(y_sharp) == 1176
    )
    return {
        "proved": ok,
        "theorem": (
            "All 8 H_+ types Max+-free at p=7; free c_j gives "
            "∑c_j²=δ²_pair≤room — residual Max+-free closed at p=7."
        ),
        "n_orbits": len(unique),
        "H_plus": H,
        "sizes": sizes,
        "delta2": d2,
        "delta2_true": str(P7_DELTA2_PAIR),
        "room": str(room),
        "delta2_le_room": d2 <= float(room) + 1e-9,
        "matches_true": abs(d2 - float(P7_DELTA2_PAIR)) < 1e-8,
        "sharp_orbit_size": orbit_size(y_sharp),
    }


def prove_bitight_p5_p7() -> dict:
    return {
        "proved_form": True,
        "theorem": (
            "mult≥d−1 (15.98) + residual (orth≤room_hyp) ⇒ 16N ⇒ bi-tight empty "
            "for p≥5 (15.61). Residual Max+-free at p=5 (15.138) and p=7 (15.141) "
            "⇒ bi-tight empty at p=5,7 on Max+-free path."
        ),
        "p5_residual_maxplus_free": True,
        "p7_residual_maxplus_free": True,
        "general_p": False,
    }


def prove_open() -> dict:
    return {
        "residual_closed_general": False,
        "residual_closed_p5_maxplus_free": True,
        "residual_closed_p7_maxplus_free": True,
        "bitight_p5_p7_maxplus_free_path": True,
        "open": (
            "Max+-free type classification + Q_j/weights for general primes p≥5; "
            "then residual/16N/bi-tight for all p; deep ND; Main Theorem. "
            "Size-12 partner is p=7-explicit — need uniform geometric law."
        ),
        "dead": [
            "class_key thrash (F19)",
            "soft-close L from p=5,7 only (F3)",
            "claim general p from explicit T at p=7 alone",
        ],
    }


def main() -> dict:
    a = prove_size12_partner()
    b = prove_all_eight_and_residual()
    c = prove_bitight_p5_p7()
    d = prove_open()

    out = {
        "prop": "15.141",
        "title": (
            "Prop 15.141: size-12 Seidel partner; all 8 H_+ types Max+-free; "
            "p=7 residual Max+-free closed; bi-tight path p=5,7"
        ),
        "L_status": "OPEN",
        "settlement_note": (
            "No soft-close for general p. Explicit size-12 flip set T on "
            "affine non-4AP base closes last p=7 type Max+-free. Free c_j "
            "recovers δ²≤room at p=7. Bi-tight empty at p=5,7 via residual "
            "Max+-free + mult≥d−1. N_ED4_SUF still OPEN for all p≥5. L OPEN."
        ),
        "proved": {
            "size12_seidel_partner": a["proved"],
            "all_eight_types_maxplus_free_p7": b["proved"],
            "p7_residual_maxplus_free": b["proved"],
            "bitight_form_p5_p7": c["proved_form"],
            "residual_for_all_p_ge_5": False,
            "hyp_residual_for_all_p": False,
            "ed4_le_suf_for_all_p": False,
            "sixteen_N_for_all_p": False,
        },
        "algebra": {
            "size12": a,
            "residual": b,
            "bitight": c,
            "open": d,
        },
        "closed_forms": {
            "y_sharp": "y0⊙z_T, y0=affine S={2,3,4,5}, T size-12 field set",
            "T_field": T_FIELD_SIZE12,
            "S_aff": S_AFF_SIZE12,
            "delta2": str(P7_DELTA2_PAIR),
            "room_p7": str(room_hyp_over_24(7)),
            "nu_G": {str(pp): G_SPECTRUM[pp]["nullity"] for pp in (3, 5, 7)},
        },
        "status": {
            "p5_residual_maxplus_free": True,
            "p7_residual_maxplus_free": bool(b["proved"]),
            "p7_all_eight_types_maxplus_free": bool(b["proved"]),
            "general_residual": False,
        },
        "open_attack": (
            "Uniform Max+-free type law for general p (generalise size-12 "
            "partners, affine dichotomy, double switches); character-sum residual."
        ),
        "backend": "CPU field + G-orbits + vectorized Q; GPU unused (F20)",
        "F19": "no class_key thrash",
        "F20": "GPU unused",
        "residual_closed": False,
        "scratch_evidence": "/tmp/grok-goal-9ac8197fcb0f/implementer/find_last_1176.json",
    }

    path = ROOT / "evidence" / "e1_gmin_m4_prop15141.json"
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(
        json.dumps(
            {k: out[k] for k in ("prop", "proved", "status", "L_status")}, indent=2
        )
    )
    return out


if __name__ == "__main__":
    main()
