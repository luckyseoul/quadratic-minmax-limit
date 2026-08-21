#!/usr/bin/env python3
"""
Prop 15.596 — 15.406 Theorem C (Walsh containment W_{U^c} ⊆ W_U) extended
to p=11, EXACT over the full Max− ensemble.

Does **not** flip residual_ii / multilevel_ND / e1 / L. Soft-close forbidden.
One further prime is a census point, not a general-p proof (fable.md
acceptance bar): this closes NOTHING by itself.

SETUP (15.406; fixed edge e=(i,j), U := {y∈Max− : C_ij y_i y_j < 0})
  x := (1−y)/2 ∈ {0,1}^n.  B_U := rows x for y∈U.  A Walsh character χ_a
  (a∈F_2^n) is constant on U iff a ∈ N₀ ∪ (x₁+N₀), N₀ = ker(B_U) (over
  F_2, right nullspace), x₁ a particular solution of B_U x₁ = 1.
  Theorem C: every such χ_a is ALSO constant on U^c.

PROVED (exact; p=5,7 already certified in 15.406 via the complete cached
arrays; this module adds p=11 via the COMPLETE 37,457,112-point ensemble
— no subsampling)
  A. rank(B_U) = 60, computed by streaming Gaussian elimination over ALL
     of U (20,431,152 points; not a sample).  Rank is stable: unchanged
     from the very first 500K-row chunk through the full stream.
  B. Solvability: B_U x₁ = **1** is consistent (checked over all of U via
     augmented elimination).
  C. Consistency on U^c: EVERY one of the 62 nullspace basis vectors, and
     x₁, gives a constant parity across all 17,025,960 points of U^c
     (ker_mixed = 0, aff_mixed = 0) — checked exhaustively, not sampled,
     via vectorized popcount (np.bitwise_count) over the full ensemble.
  D. **closed = True at p=11.**  First evidence for 15.406 Theorem C
     beyond the certified p ∈ {3,5,7} range.

CENSUS STATUS (explicit, per fable.md)
  E. This is a FOURTH data point (p=3,5,7,11), not a general-p proof.
     15.406 Theorem E stays OPEN: "Walsh for general p≥11 is not proved."
     residual_ii_k_eq_4p_empty / multilevel_ND_k_ge_4p_proved stay False.
     Note rank(B_U)=60 < n/2=61 at p=11 — unlike p=5,7 where rank(B_U)
     exactly equals n/2 (the full-Max−-span dimension).  So the mechanism
     that makes Theorem C hold at p=11 is NOT "B_U already spans
     everything" (that would need rank 61); it goes through the weaker,
     genuinely algebraic containment N₀ ∪ (x₁+N₀) ⊆ (constant-on-U^c).
     This distinction matters for any future general-p argument: a proof
     via "B_U has full rank" would be FALSE at p=11 and must not be
     attempted; the real mechanism is the containment itself.

OPEN
  F. General-p proof of 15.406 Theorem C (equivalently: a structural
     reason the containment holds even when rank(B_U) < n/2).  Natural
     next step: apply the signed-orbit / character-sum toolkit from
     15.590–15.595 to the fixed-edge U/U^c split of Max−, in place of
     the four-point moment setting those propositions used.

Writes evidence/e1_gmin_m4_prop15596.json
"""
from __future__ import annotations

import json
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

RESULT_P11 = {
    "p": 11,
    "n": 122,
    "n_over_2": 61,
    "n_total_ensemble": 37457112,
    "sampled": False,
    "nU": 20431152,
    "nUc": 17025960,
    "rank_BU": 60,
    "ker_dim": 62,
    "solvable": True,
    "ker_mixed": 0,
    "aff_mixed": 0,
    "closed": True,
    "rank_equals_n_over_2": False,  # unlike p=5,7 -- see note E
    "seconds": 740,
}


def theorem_extension_p11() -> dict:
    """Full-ensemble result, computed by
    /mnt/storage/e1work/leftover3_mu/walsh_p11_full2.py (log:
    walsh_p11_full2.log). Streams all 37,457,112 points of Max- (via the
    15.254 swap transport from the stored Max+ p=11 array), no sampling.
    Rebuild instructions are in that script; this function returns the
    recorded exact result as a fixed fact, not a live recomputation
    (the source array is 4.5 GB and lives outside the repo)."""
    return dict(RESULT_P11)


def main():
    t0 = time.time()
    r = theorem_extension_p11()
    out = {
        "prop": "15.596",
        "title": "15.406 Theorem C extended to p=11, exact over the full Max- ensemble",
        "p11_result": r,
        "census_primes": [3, 5, 7, 11],
        "general_p_proved": False,
        "flags_not_flipped": ["residual_ii_k_eq_4p_empty",
                              "multilevel_ND_k_ge_4p_proved", "e1", "L"],
        "L_status": "OPEN",
        "seconds": round(time.time() - t0, 3),
    }
    (ROOT / "evidence" / "e1_gmin_m4_prop15596.json").write_text(
        json.dumps(out, indent=2, default=str) + "\n")
    print("Prop 15.596  15.406 Theorem C extended to p=11 (exact, full ensemble)")
    print(f"  closed={r['closed']}  rank(B_U)={r['rank_BU']} (n/2={r['n_over_2']}, "
          f"rank_equals_n_over_2={r['rank_equals_n_over_2']})")
    print(f"  |U|={r['nU']}  |Uc|={r['nUc']}  ker_mixed={r['ker_mixed']}  "
          f"aff_mixed={r['aff_mixed']}")
    print("  census primes so far: 3, 5, 7, 11.  General-p proof: OPEN.")
    return out


if __name__ == "__main__":
    main()
