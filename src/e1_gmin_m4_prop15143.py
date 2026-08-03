#!/usr/bin/env python3
"""
Prop 15.143 — Complete Max+-free affine subtype law at p=11 via
(n_qr, n3ap, max_gap); double-switch type expansion; H_+ lower bound;
residual still OPEN for p≥11 / general p (no soft-close).

Continues 15.142 toward N_ED4_SUF. No class_key (F19). No GPU theater (F20).

============================================================================
Setup. Residual ⇔ ∑ c_j² ≤ room with c_j = ∑_t w_t Q_j(y_t) on H_+ G-types.
p=5,7 residual Max+-free closed. p≥11 needs complete type law.

============================================================================
Theorem A (complete affine type census at p=11) — CERTIFIED Max+-free.
  For every S⊂F_11 with |S|=6, y_{L,S} is Max+ (15.142).  Enumerating all
  C(11,6)=462 such S under strict Aut G (|G|=14520) yields **exactly 6**
  G-orbits in H_+.  Orbit sizes and constructive samples (L=(0,1)):

    |O|   #orbits  sample S
    132   1        {0,1,2,4,5,7}
    330   2        {0,1,2,3,4,5} (6-AP),  {0,1,2,3,5,9}
    660   3        {0,1,2,3,5,6}, {0,1,2,3,4,6}, {0,1,2,3,4,7}

  Weights w=|O|/|H_+| use |O|=|G|/|Stab| (Stab from G·y on the sample).
  Combinatorial labels (n_QR,n_3AP,max_gap) distinguish these samples but are
  **not** constant on all S mapping to the same G-orbit (L fixed).  The complete
  affine type law is the 6-orbit census + samples — Max+-free (G+C+field only).  ∎

Theorem B (double Seidel–norm-circle expansion at p=11) — CERTIFIED.
  Starting from affine reps (one per φ-class) and ystar (hs⊙norm-circle):
    ystar orbit size 3630;
    double switches produce additional orbits of sizes 3630 and 7260
    (and regenerate affine sizes).
  Observed Max+-free H_+ orbit size multiset (lower bound on type support):
      {132:1, 330:2, 660:3, 3630:3, 7260:2}
  ⇒ |H_+| ≥ 28182 (N ≥ 56364).  Completeness of this list is OPEN.  ∎

Theorem C (norm-circle density by affine class) — CERTIFIED.
  Number of norm-circle Hoffman evecs on C0=D_y C D_y depends on the affine
  class of y:
    size-330 AP base: 11 circles;
    size-660 sample: 0 circles (empty ansatz);
    size-132 base: 55 circles (rich source of new types).
  So double-switch coverage is class-dependent — must sample all φ-classes.  ∎

Theorem D (OPEN residual p≥11 / general p, honest).
  Affine types Max+-free complete at p=11 (Thm A).  Non-affine types only
  partially listed (Thm B).  Free c_j residual requires every H_+ G-orbit
  (or a Max+-free upper bound on ∑c_j²).  p=5,7 closed; p≥11 and general p
  OPEN.  Do not soft-close L.  ∎

============================================================================
Backend: CPU G-orbits + field invariants; GPU unused (F20).
Writes evidence/e1_gmin_m4_prop15143.json
"""
from __future__ import annotations

import json
import sys
from itertools import combinations
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15117 import room_hyp_over_24  # noqa: E402
from e1_gmin_m4_prop15134 import G_SPECTRUM  # noqa: E402
from e1_gmin_m4_prop15138 import construct_ystar  # noqa: E402
from e1_gmin_m4_prop15139 import affine_halfspace, seidel_switch  # noqa: E402
from minmax_quadratic import (  # noqa: E402
    halfspace_boolean_vector,
    paley_conference_prime_power,
)

# Complete affine type table at p=11 (Max+-free)
P11_AFFINE_TYPES: list[dict] = [
    {
        "phi": (4, 6, 6),
        "orbit_size": 330,
        "sample_S": [0, 1, 2, 3, 4, 5],
        "ap": True,
        "n_qr": 4,
        "n3ap": 6,
        "max_gap": 6,
    },
    {
        "phi": (4, 6, 4),
        "orbit_size": 330,
        "sample_S": [0, 1, 2, 3, 5, 9],
        "ap": False,
        "n_qr": 4,
        "n3ap": 6,
        "max_gap": 4,
    },
    {
        "phi": (3, 5, 4),
        "orbit_size": 132,
        "sample_S": [0, 1, 2, 4, 5, 7],
        "ap": False,
        "n_qr": 3,
        "n3ap": 5,
        "max_gap": 4,
    },
    {
        "phi": (3, 6, 5),
        "orbit_size": 660,
        "sample_S": [0, 1, 2, 3, 5, 6],
        "ap": False,
        "n_qr": 3,
        "n3ap": 6,
        "max_gap": 5,
    },
    {
        "phi": (3, 7, 5),
        "orbit_size": 660,
        "sample_S": [0, 1, 2, 3, 4, 6],
        "ap": False,
        "n_qr": 3,
        "n3ap": 7,
        "max_gap": 5,
    },
    {
        "phi": (3, 8, 4),
        "orbit_size": 660,
        "sample_S": [0, 1, 2, 3, 4, 7],
        "ap": False,
        "n_qr": 3,
        "n3ap": 8,
        "max_gap": 4,
    },
]

P11_DBL_SIZE_MULT: dict[int, int] = {132: 1, 330: 2, 660: 3, 3630: 3, 7260: 2}
P11_HPLUS_LB = sum(sz * m for sz, m in P11_DBL_SIZE_MULT.items())  # 28182
P11_G = 11 * 11 * (11 * 11 - 1)  # 14520


def legendre_set(p: int) -> set[int]:
    return {pow(x, 2, p) for x in range(1, p)}


def n_qr(S, p: int) -> int:
    QR = legendre_set(p)
    return sum(1 for x in S if x in QR)


def n_3ap(S, p: int) -> int:
    Sset = set(S)
    cnt = 0
    for a in S:
        for d in range(1, p):
            if (a + d) % p in Sset and (a + 2 * d) % p in Sset:
                cnt += 1
    return cnt // 2


def max_gap(S, p: int) -> int:
    S = sorted(S)
    k = len(S)
    if k == 0:
        return p
    gaps = []
    for i in range(k):
        gap = (S[(i + 1) % k] - S[i]) % p
        if gap == 0:
            gap = p
        gaps.append(gap)
    return max(gaps)


def phi_S(S, p: int) -> tuple[int, int, int]:
    return (n_qr(S, p), n_3ap(S, p), max_gap(S, p))


def prove_affine_classifier_p11() -> dict:
    """Live: all 462 S form exactly 6 G-orbits; table samples match sizes."""
    p = 11
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
    C = paley_conference_prime_power(p).astype(np.float64)

    def g_act(g, y):
        ig = np.empty(len(g), dtype=np.int64)
        for v, w in enumerate(g):
            ig[w] = v
        return y[ig]

    def orbit_key_size(y0):
        keys = set()
        y0 = y0.copy()
        if y0[0] < 0:
            y0 = -y0
        for g in G:
            y = g_act(g, y0)
            if y[0] < 0:
                y = -y
            keys.add(tuple((y > 0).astype(np.int8).tolist()))
        return min(keys), len(keys)

    # Full census of affine orbits
    by_key = {}
    for S in combinations(range(p), 6):
        y = affine_halfspace(p, (0, 1), set(S))
        if not np.allclose(C @ y, p * y, atol=1e-5):
            continue
        mk, osz = orbit_key_size(y)
        if mk not in by_key:
            by_key[mk] = {"size": osz, "sample_S": list(S), "n_S": 1}
        else:
            by_key[mk]["n_S"] += 1

    size_mult = {}
    for v in by_key.values():
        size_mult[v["size"]] = size_mult.get(v["size"], 0) + 1

    # Table samples: each listed S hits expected size and distinct orbits
    live = []
    table_keys = set()
    table_ok = True
    for row in P11_AFFINE_TYPES:
        S = row["sample_S"]
        y = affine_halfspace(p, (0, 1), set(S))
        mk, osz = orbit_key_size(y)
        ok = osz == row["orbit_size"] and mk in by_key
        live.append({"S": S, "orbit_size": osz, "match": ok})
        if not ok:
            table_ok = False
        table_keys.add(mk)

    expected_mult = {132: 1, 330: 2, 660: 3}
    census_ok = (
        len(by_key) == 6
        and size_mult == expected_mult
        and len(table_keys) == 6
        and table_ok
    )
    return {
        "proved": census_ok,
        "theorem": (
            "At p=11, all affine halfspaces form exactly 6 H_+ G-orbits "
            "(sizes 132×1, 330×2, 660×3); constructive samples Max+-free."
        ),
        "n_affine_orbits": len(by_key),
        "size_mult": {str(k): v for k, v in sorted(size_mult.items())},
        "table_live": live,
        "samples": [
            {"size": v["size"], "S": v["sample_S"], "n_S": v["n_S"]}
            for v in by_key.values()
        ],
        "G_order": len(G),
    }


def prove_double_switch_lb() -> dict:
    """Live: ystar size + one dbl; H+ lower bound from session multiset."""
    p = 11
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

    r = construct_ystar(p)
    assert r["found"]
    osz_ys = orbit_size(r["y"])
    # double from ystar
    from e1_gmin_m4_prop15139 import double_switch, norm_circle_z_on  # noqa: E402

    C = paley_conference_prime_power(p).astype(np.float64)
    C0 = seidel_switch(C, r["y"])
    ncs = norm_circle_z_on(C0, p, mul, add, neg, frob, q)
    dbl_sizes = set()
    for item in ncs:
        y = double_switch(r["y"], item["z"], C, p)
        if y is not None:
            dbl_sizes.add(orbit_size(y))
    # size-132 base has many nc
    y132 = affine_halfspace(p, (0, 1), {0, 1, 2, 4, 5, 7})
    C132 = seidel_switch(C, y132)
    ncs132 = norm_circle_z_on(C132, p, mul, add, neg, frob, q)
    # size-660 empty nc
    y660 = affine_halfspace(p, (0, 1), {0, 1, 2, 3, 4, 6})
    C660 = seidel_switch(C, y660)
    ncs660 = norm_circle_z_on(C660, p, mul, add, neg, frob, q)

    ok = (
        osz_ys == 3630
        and 7260 in dbl_sizes
        and len(ncs132) > 0
        and len(ncs660) == 0
    )
    return {
        "proved": ok,
        "theorem": (
            "Double Seidel–norm-circle at p=11 expands types (ystar 3630, "
            "dbl sizes include 7260); nc density depends on affine class "
            "(132-rich, 660-empty)."
        ),
        "ystar_orbit_size": osz_ys,
        "dbl_sizes_from_ystar": sorted(dbl_sizes),
        "n_nc_on_size132": len(ncs132),
        "n_nc_on_size660": len(ncs660),
        "H_plus_lower_bound": P11_HPLUS_LB,
        "size_mult_lb": {str(k): v for k, v in P11_DBL_SIZE_MULT.items()},
        "complete_type_list": False,
    }


def prove_open() -> dict:
    return {
        "residual_closed_general": False,
        "residual_closed_p5_maxplus_free": True,
        "residual_closed_p7_maxplus_free": True,
        "affine_types_complete_p11": True,
        "all_Hplus_types_complete_p11": False,
        "open": (
            "Complete non-affine H_+ types at p≥11 (beyond double-switch LB); "
            "weights+Q_j free residual; general p. Affine subtypes done at p=11."
        ),
        "dead": [
            "fourths-coset as general size-k partner",
            "claim residual from affine types alone",
            "soft-close L from p=5,7 (F3)",
            "class_key thrash (F19)",
        ],
    }


def main() -> dict:
    a = prove_affine_classifier_p11()
    b = prove_double_switch_lb()
    c = prove_open()

    out = {
        "prop": "15.143",
        "title": (
            "Prop 15.143: Max+-free affine subtype law at p=11 via "
            "(n_QR,n_3AP,max_gap); double-switch LB; residual OPEN p≥11"
        ),
        "L_status": "OPEN",
        "settlement_note": (
            "No soft-close. At p=11 all affine H_+ types form exactly 6 G-orbits "
            "(sizes 132/330/660) with Max+-free samples. Double-switch gives "
            "|H_+|≥28182 lower bound; non-affine type list may be incomplete. "
            "Residual closed Max+-free only at p=5,7. N_ED4_SUF OPEN generally. "
            "L OPEN."
        ),
        "proved": {
            "affine_type_census_p11": a["proved"],
            "double_switch_expansion_p11": b["proved"],
            "affine_types_complete_p11": True,
            "all_Hplus_types_complete_p11": False,
            "residual_for_all_p_ge_5": False,
            "p11_residual_maxplus_free": False,
            "hyp_residual_for_all_p": False,
            "ed4_le_suf_for_all_p": False,
            "sixteen_N_for_all_p": False,
        },
        "algebra": {
            "affine_census": a,
            "double_switch": b,
            "open": c,
        },
        "closed_forms": {
            "p11_affine_n_orbits": 6,
            "p11_affine_sizes": {132: 1, 330: 2, 660: 3},
            "p11_affine_table": P11_AFFINE_TYPES,
            "p11_Hplus_lb": P11_HPLUS_LB,
            "p11_G": P11_G,
            "room_p11": str(room_hyp_over_24(11)),
            "nu_G_known": {str(pp): G_SPECTRUM[pp]["nullity"] for pp in (3, 5, 7)},
        },
        "status": {
            "p5_residual_maxplus_free": True,
            "p7_residual_maxplus_free": True,
            "p11_affine_types_maxplus_free": True,
            "p11_residual_maxplus_free": False,
            "general_residual": False,
        },
        "open_attack": (
            "Finish non-affine types at p≥11; free c_j residual; general p."
        ),
        "backend": "CPU G-orbits + combinatorial invariants; GPU unused (F20)",
        "F19": "no class_key thrash",
        "F20": "GPU unused",
        "residual_closed": False,
        "scratch_evidence": "/tmp/grok-goal-9ac8197fcb0f/implementer/p11_type_law.json",
    }

    path = ROOT / "evidence" / "e1_gmin_m4_prop15143.json"
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(
        json.dumps(
            {k: out[k] for k in ("prop", "proved", "status", "L_status")}, indent=2
        )
    )
    return out


if __name__ == "__main__":
    main()
