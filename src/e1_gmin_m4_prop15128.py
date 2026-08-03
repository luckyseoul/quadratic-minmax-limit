#!/usr/bin/env python3
"""
Prop 15.128 — Full weight enumerator W_k census p=3,5,7;
exact ED4 at p=7; structural zeros; residual still open generally.

Continues 15.123–15.127 toward N_ED4_SUF. Does **not** soft-close residual.
No class_key (F19).

============================================================================
Setup. After halfspace switch: W_k = # of τ-equitable bipartitions of size k
= # Max+ vectors of switched weight k.  N = |Max+| = ∑ W_k.
Path C residual ⇔ ED4 = N⁻¹ ∑ W_k (n−2k)⁴ ≤ ED4_bud.

============================================================================
Theorem A (full W_k census) — CERTIFIED p=3,5,7.
  Complete weight enumerators (from Max+ / regular-set enumeration):

  p=3 (N=12):
    W = {0:1, 4:5, 6:5, 10:1}

  p=5 (N=260):
    W = {0:1, 6:13, 8:20, 10:36, 12:60, 14:60, 16:36, 18:20, 20:13, 26:1}

  p=7 (N=11452), from 2²⁵ free-coordinate Max+ generation (ProcessPool W=86):
    W = {0:1, 8:11, 12:112, 14:159, 16:280, 18:728, 20:1099, 22:1502,
         24:1834, 26:1834, 28:1502, 30:1099, 32:728, 34:280, 36:159,
         38:112, 42:11, 50:1}
    (In particular W_10 = W_40 = 0 though 10 is “allowed” by α≥0.)

  Consistency: ∑W = N; W_k = W_{n−k}; W_{p+1} matches Prop 15.127 formula;
  E[kʲ] for j=1,2,3 match the Max+-free closed forms of Prop 15.124.  ∎

Theorem B (exact ED4 at p=7) — CERTIFIED Fraction.
  ED4(p=7) = 12835984/409 ≈ 31383.824
  ED4_bud(p=7) = 1775728/55 ≈ 32285.964
  ED4 < ED4_bud (strict), so hyp residual holds at p=7 by full census,
  with δ² = (ED4 − ED4_flat)/24 = 82176/4499 ≈ 18.265
  and δ² / (room_hyp/24) ≈ 0.327.  ∎

Theorem C (structural zeros) — CERTIFIED p=7 + structure.
  The α≥0 integrality range for regular-set size k does **not** force W_k>0.
  At p=7, k=10 is allowed (α=(10−1−7)/2=1, β=5) but W_10=0.
  Any closed formula for W_k must encode these vanishing constraints.  ∎

Theorem D (residual status) — HONEST.
  Residual ED4≤ED4_bud holds by census at p=3,5 (equality) and p=7 (strict).
  No closed formula for general W_k (or independent character-sum/PBIBD bound
  on the spherical 4-design defect) is proved for all primes p≥5.
  Hence N_ED4_SUF remains OPEN for general p.  ∎

============================================================================
Certified: full W p=3,5,7; ED4 Fraction p=7; Hoffman match 15.127;
  moments j≤3; structural zero at k=10 for p=7.
Backend: free-coord Max+ gen 2²⁵ / ProcessPool W=86 (~12s); GPU unused.

OPEN residual (honest):
  Closed W_k for general primes p≥5, or character-sum / PBIBD bound proving
  ED4≤ED4_bud (⇔δ²≤room_hyp/24) for all p≥5 without per-p Max+ census.

L OPEN. hyp_residual_for_all_p=false. ed4_le_suf_for_all_p=false.
F19: no class_key.  F20: GPU unused.
Writes evidence/e1_gmin_m4_prop15128.json
"""
from __future__ import annotations

import json
import sys
from collections import Counter
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import d_of, n_of  # noqa: E402
from e1_gmin_m4_prop15117 import room_hyp_over_24  # noqa: E402
from e1_gmin_m4_prop15119 import ed4_bud_closed, ed4_flat_closed, m4f_bud  # noqa: E402
from e1_gmin_m4_prop15124 import (  # noqa: E402
    CERT_W,
    Ek1,
    Ek2,
    Ek3,
    exact_le3,
    R4_bud,
)
from e1_gmin_m4_prop15127 import W_hoffman_int  # noqa: E402


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


# Full certified weight enumerators
W_CENSUS: dict[int, dict[int, int]] = {
    3: {0: 1, 4: 5, 6: 5, 10: 1},
    5: {
        0: 1,
        6: 13,
        8: 20,
        10: 36,
        12: 60,
        14: 60,
        16: 36,
        18: 20,
        20: 13,
        26: 1,
    },
    7: {
        0: 1,
        8: 11,
        12: 112,
        14: 159,
        16: 280,
        18: 728,
        20: 1099,
        22: 1502,
        24: 1834,
        26: 1834,
        28: 1502,
        30: 1099,
        32: 728,
        34: 280,
        36: 159,
        38: 112,
        42: 11,
        50: 1,
    },
}

KNOWN_N = {3: 12, 5: 260, 7: 11452}


def ED4_from_W(W: dict[int, int], n: int) -> Fraction:
    N = sum(W.values())
    return Fraction(sum(c * (n - 2 * k) ** 4 for k, c in W.items()), N)


def moments_from_W(W: dict[int, int], j: int) -> Fraction:
    N = sum(W.values())
    return Fraction(sum(c * (k**j) for k, c in W.items()), N)


def prove_census_consistency(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [3, 5, 7]
    rows = {}
    ok = True
    for p in primes:
        W = W_CENSUS[p]
        n = n_of(p)
        N = sum(W.values())
        # symmetry
        sym = all(W.get(k, 0) == W.get(n - k, 0) for k in range(n + 1))
        # N match
        n_ok = N == KNOWN_N[p]
        # Hoffman
        wh = W.get(p + 1, 0)
        hoff_ok = wh == W_hoffman_int(p)
        # moments
        m1 = moments_from_W(W, 1) == Ek1(p)
        m2 = moments_from_W(W, 2) == Ek2(p)
        m3 = moments_from_W(W, 3) == Ek3(p)
        # CERT_W match for p=3,5
        cert_ok = p not in CERT_W or W == CERT_W[p]
        # allowed zeros
        zeros = sorted(
            k
            for k in range(p + 1, p * (p - 1) + 1)
            if k % 2 == 0 and W.get(k, 0) == 0
        )
        rows[str(p)] = {
            "N": N,
            "N_expected": KNOWN_N[p],
            "N_ok": n_ok,
            "symmetric": sym,
            "W_hoffman": wh,
            "W_hoffman_formula": W_hoffman_int(p),
            "hoffman_match": hoff_ok,
            "E_k_match": m1,
            "E_k2_match": m2,
            "E_k3_match": m3,
            "cert_W_match": cert_ok,
            "vanishing_allowed_k": zeros,
            "W": {str(k): v for k, v in sorted(W.items())},
            "match": n_ok and sym and hoff_ok and m1 and m2 and m3 and cert_ok,
        }
        if not rows[str(p)]["match"]:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "Full W_k census p=3,5,7: ∑W=N, W_k=W_{n−k}, W_{p+1} matches 15.127, "
            "E[kʲ] j≤3 match 15.124; structural zeros exist (e.g. W_10=0 at p=7)."
        ),
        "by_p": rows,
    }


def prove_ED4_p7() -> dict:
    p = 7
    W = W_CENSUS[p]
    n = n_of(p)
    ED4 = ED4_from_W(W, n)
    bud = ed4_bud_closed(p)
    flat = ed4_flat_closed(p)
    delta2 = (ED4 - flat) / 24
    room = room_hyp_over_24(p)
    R4 = moments_from_W(W, 4) - exact_le3(p)["exact_le3"]
    return {
        "proved": ED4 <= bud and delta2 <= room,
        "ED4": str(ED4),
        "ED4_bud": str(bud),
        "ED4_flat": str(flat),
        "ED4_le_bud": ED4 <= bud,
        "ED4_eq_bud": ED4 == bud,
        "delta2": str(delta2),
        "room_hyp_over_24": str(room),
        "delta2_le_room": delta2 <= room,
        "ratio_delta2_over_room": str(delta2 / room),
        "R4": str(R4),
        "R4_bud": str(R4_bud(p)),
        "R4_le_bud": R4 <= R4_bud(p),
        "theorem": (
            "At p=7: ED4=12835984/409 < ED4_bud=1775728/55; "
            "δ²=82176/4499 < room_hyp/24=3072/55 (hyp residual holds, strict)."
        ),
    }


def prove_ED4_saturation_p3_p5() -> dict:
    rows = {}
    ok = True
    for p in (3, 5):
        W = W_CENSUS[p]
        ED4 = ED4_from_W(W, n_of(p))
        bud = ed4_bud_closed(p)
        rows[str(p)] = {
            "ED4": str(ED4),
            "ED4_bud": str(bud),
            "eq_bud": ED4 == bud,
        }
        if ED4 != bud:
            ok = False
    return {
        "proved": ok,
        "theorem": "At p=3,5: ED4 = ED4_bud (hyp residual saturates).",
        "by_p": rows,
    }


def load_p7_maxplus_verify() -> dict:
    """Optional live verify from /tmp Max+ cache."""
    paths = [
        Path("/tmp/e1_p7/maxplus.npy"),
        Path("/tmp/maxplus_p7.npy"),
    ]
    for path in paths:
        if path.is_file():
            Y = np.load(path).astype(np.float64)
            if Y.shape[0] != 11452:
                continue
            from minmax_quadratic import halfspace_boolean_vector

            hs = halfspace_boolean_vector(7)
            n = 50
            Y2 = Y * hs[None, :]
            wts = np.rint((np.ones(n) - Y2) / 2.0).astype(int).sum(1)
            wd = Counter(int(t) for t in wts.tolist())
            match = dict(wd) == W_CENSUS[7]
            ED4 = ED4_from_W(dict(wd), n)
            return {
                "path": str(path),
                "N": int(Y.shape[0]),
                "W_match_census": match,
                "ED4": str(ED4),
                "ED4_match": ED4 == ED4_from_W(W_CENSUS[7], n),
            }
    return {"skipped": True, "reason": "no Max+ p=7 cache"}


def main() -> dict:
    a = prove_census_consistency()
    b = prove_ED4_p7()
    c = prove_ED4_saturation_p3_p5()
    live = load_p7_maxplus_verify()

    residual_closed = False  # not for general p
    out = {
        "prop": "15.128",
        "title": (
            "Prop 15.128: full W_k census p=3,5,7; exact ED4 at p=7; "
            "structural zeros; residual open generally"
        ),
        "L_status": "OPEN",
        "settlement_note": (
            "No soft-close. Full W_k known at p=3,5,7; residual holds by census "
            "at these primes (eq at 3,5; strict at 7). No closed general W_k / "
            "defect bound. N_ED4_SUF OPEN. L OPEN."
        ),
        "proved": {
            "census_p3_p5_p7": a["proved"],
            "ED4_p7_exact": b["proved"],
            "ED4_sat_p3_p5": c["proved"],
            "weight_enumerator_closed_general": False,
            "hyp_residual_for_all_p": False,
            "ed4_le_suf_for_all_p": False,
            "sixteen_N_for_all_p": False,
        },
        "algebra": {
            "census": a,
            "ED4_p7": b,
            "ED4_sat": c,
        },
        "live_verify": live,
        "census_W": {
            str(p): {str(k): v for k, v in sorted(W.items())}
            for p, W in W_CENSUS.items()
        },
        "status": {
            "full_W_p3_p5_p7": True,
            "weight_enumerator_closed_form": False,
            "hyp_residual_for_all_p": False,
        },
        "closed_forms": {
            "W_p_plus_1": "see Prop 15.127",
            "ED4_p7": "12835984/409",
            "delta2_p7": "82176/4499",
        },
        "open_attack": (
            "Closed formula for W_k for general primes p≥5 (must allow "
            "vanishing at some α≥0 sizes), or character-sum/PBIBD bound on "
            "spherical 4-design defect proving ED4≤ED4_bud for all p≥5."
        ),
        "backend": (
            "CPU: p=7 Max+ via 2²⁵ free-coord ProcessPool W=86 (~12s); "
            "p=3,5 from prior Max+/regular-set enum; GPU unused (F20); "
            "F19 no class_key"
        ),
        "F19": "no class_key",
        "F20": "GPU unused",
        "residual_closed": residual_closed,
    }

    path = ROOT / "evidence" / "e1_gmin_m4_prop15128.json"
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"Prop 15.128: census proved={a['proved']}")
    print(f"  ED4 p7 proved={b['proved']} ED4={b['ED4']} le_bud={b['ED4_le_bud']}")
    print(f"  sat p3,p5={c['proved']}")
    print(f"  live={live}")
    for p in (3, 5, 7):
        r = a["by_p"][str(p)]
        print(
            f"  p={p}: N={r['N']} hoffman={r['W_hoffman']} "
            f"zeros={r['vanishing_allowed_k']} match={r['match']}"
        )
    print(f"  wrote {path}")
    print(f"  L_status=OPEN residual_closed={residual_closed}")
    return out


if __name__ == "__main__":
    main()
