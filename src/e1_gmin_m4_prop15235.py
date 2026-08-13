#!/usr/bin/env python3
"""
Prop 15.235 — Signed k=0 layer of R̄₄: full 24-perm permanent;
S₄ cycle-type split; 4-cycle inverse pairing; unsigned still dead.
Signed-layer structure program (k=3,2,1,0) complete; residual-(i) OPEN.

Does **not** flip residual_i / type_I / gsum / E1 / L. Soft-close forbidden.

HINGE GRAPH
  k=3 linear (15.232), k=2 bilinear (15.233), k=1 trilinear (15.234),
  k=0 degree-4 cycle-type (this prop).  Unsigned every layer dead except
  k=3 for p≥89.  Preferred remaining hinge: μ-weighted signed sums
  (character sums before abs) or Aut-SOS+diamond.

PROVED Max+-free (C_ij=±1; no shared vertices ⇒ no forced zeros)
  A. Full S₄ permanent.
     |S∩T|=0 ⇒ C[S,T] has no structural zeros.  All 24 permutations
     survive.  |per|≤24.  ∎

  B. Cycle-type split.
        per = P_id + P_trans + P_dbl + P_3cyc + P_4cyc
     with 1+6+3+8+6 permutations (cycle types 1⁴, 21², 2², 31, 4).
     Each part is a signed sum of the corresponding products.  ∎

  C. 4-cycle inverse pairing (sum before abs inside the 4-cycle part).
     The 6 four-cycles form 3 inverse pairs {σ,σ⁻¹}.  Each product is
     ±1, so w_σ+w_{σ⁻¹} ∈ {−2,0,2}.  Hence P_4cyc is even and
     |P_4cyc|≤6.  Pairing inverses before abs is mandatory inside this
     class (the unsigned 4-cycle mass is 6·#T_0 ≫ B).  ∎

  D. Inverse-pair majorant still dead.
     Even after (C), |per|≤24 still, and 24·C(n−4,4)=(p²−3)(p²−4)(p²−5)(p²−6)
     exceeds B for every prime p≥5 (15.232 E, k=0).  Cycle-abs pairing
     does not close |R̄₄|≤B.  ∎

  E. Layer program status.
     Signed algebraic type of every intersection layer is now named
     (linear / bilinear / trilinear / degree-4 cycle-type).  The
     residual-(i) hinge remains the μ-weighted signed sum (or an
     alternate Path C / K₄ / free-e+ker=sc / Aut-SOS).  ∎

OPEN
  F. |R̄₄|≤B on |κ|=1 for all primes p≥5.  Soft-close forbidden.

Writes evidence/e1_gmin_m4_prop15235.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from itertools import permutations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    is_prime,
    n_of,
    residual_i_dual_eq_empty_proved_general,
    type_I_k_3p_minus_2_closed_general,
)
from e1_gmin_m4_prop15229 import R_bar_budget_for_mu_thr  # noqa: E402
from e1_gmin_m4_prop15232 import (  # noqa: E402
    layer_count,
    ultra_crude_layer_closed,
)
from e1_gmin_m4_prop15234 import theorem_k1_laplace  # noqa: E402


CYC_NAME = {
    (1, 1, 1, 1): "id",
    (2, 1, 1): "trans",
    (2, 2): "dbl",
    (3, 1): "3cyc",
    (4,): "4cyc",
}


def cycle_type(sig) -> tuple[int, ...]:
    seen = [False] * 4
    lens = []
    for i in range(4):
        if seen[i]:
            continue
        L = 0
        j = i
        while not seen[j]:
            seen[j] = True
            j = sig[j]
            L += 1
        lens.append(L)
    return tuple(sorted(lens, reverse=True))


def invert_perm(sig) -> tuple[int, ...]:
    inv = [0] * 4
    for i, j in enumerate(sig):
        inv[j] = i
    return tuple(inv)


def cycle_type_permanent_parts(C, S, T) -> dict[str, int]:
    """S₄ cycle-type split of per(C[S,T])."""
    S, T = list(S), list(T)
    parts = {name: 0 for name in ("id", "trans", "dbl", "3cyc", "4cyc")}
    for sig in permutations(range(4)):
        w = 1
        for i in range(4):
            w *= int(C[S[i], T[sig[i]]])
        parts[CYC_NAME[cycle_type(sig)]] += w
    return parts


def four_cycle_inverse_pair_sum(C, S, T) -> int:
    """∑_{4-cycle pairs} (w_σ + w_{σ⁻¹}) = P_4cyc, via inverse pairing."""
    S, T = list(S), list(T)
    seen = set()
    tot = 0
    for sig in permutations(range(4)):
        if cycle_type(sig) != (4,):
            continue
        key = tuple(sorted((sig, invert_perm(sig))))
        if key in seen:
            continue
        seen.add(key)

        def weight(p):
            w = 1
            for i in range(4):
                w *= int(C[S[i], T[p[i]]])
            return w

        tot += weight(sig) + weight(invert_perm(sig))
    return tot


def k0_unsigned_majorant(p: int) -> int:
    """24 · C(n−4,4) = (p²−3)(p²−4)(p²−5)(p²−6)."""
    return ultra_crude_layer_closed(p, 0)


def theorem_k0_full_s4() -> dict:
    return {
        "proved": True,
        "theorem": (
            "On |S∩T|=0 the 4×4 block has no structural zeros, so all 24 "
            "permutations survive and |per|≤24."
        ),
        "n_surviving_perms": 24,
        "abs_per_ub": 24,
        "depends_on": ["C_off_pm1"],
    }


def theorem_k0_cycle_type_split() -> dict:
    return {
        "proved": True,
        "theorem": (
            "per=P_id+P_trans+P_dbl+P_3cyc+P_4cyc with conjugacy sizes "
            "1,6,3,8,6.  Max+-free S₄ partition."
        ),
        "class_sizes": {"id": 1, "trans": 6, "dbl": 3, "3cyc": 8, "4cyc": 6},
        "depends_on": ["theorem_k0_full_s4"],
    }


def theorem_k0_four_cycle_inverse_pairing() -> dict:
    return {
        "proved": True,
        "theorem": (
            "The 6 four-cycles form 3 inverse pairs.  Each product is ±1, "
            "so w_σ+w_{σ⁻¹}∈{−2,0,2} and P_4cyc is even with |P_4cyc|≤6. "
            "Must pair inverses before abs inside this class."
        ),
        "n_inverse_pairs": 3,
        "abs_P4cyc_ub": 6,
        "depends_on": ["theorem_k0_cycle_type_split"],
    }


def theorem_k0_unsigned_still_dead() -> dict:
    ok = True
    sample = {}
    for p in range(5, 220):
        if not is_prime(p):
            continue
        maj = k0_unsigned_majorant(p)
        B = R_bar_budget_for_mu_thr(p)
        if maj != 24 * layer_count(n_of(p), 0):
            ok = False
            break
        if Fraction(maj) <= B:
            ok = False
            break
        if p in (5, 7, 11, 13, 17, 23, 89):
            sample[str(p)] = {"majorant": maj, "B": str(B)}
    if ok:
        for p in (227, 233, 239, 241):
            if Fraction(k0_unsigned_majorant(p)) <= R_bar_budget_for_mu_thr(p):
                ok = False
                break
    return {
        "proved": ok,
        "theorem": (
            "Even after 4-cycle inverse pairing, |per|≤24 and "
            "24·C(n−4,4)=(p²−3)(p²−4)(p²−5)(p²−6) exceeds B for every "
            "prime p≥5.  Cycle-abs pairing does not close residual-(i)."
        ),
        "by_p_sample": sample,
        "depends_on": [
            "theorem_k0_full_s4",
            "diamond_abs_mu_le_1",
            "prop_15_232_unsigned_k0",
        ],
        "sufficient_for_residual_i": False,
    }


def residual_i_closed_via_235() -> bool:
    return False


def hinge_status_235() -> dict:
    return {
        "k1_laplace": theorem_k1_laplace()["proved"],
        "k0_full_s4": theorem_k0_full_s4()["proved"],
        "k0_cycle_type_split": theorem_k0_cycle_type_split()["proved"],
        "k0_four_cycle_inverse_pairing": theorem_k0_four_cycle_inverse_pairing()[
            "proved"
        ],
        "k0_unsigned_still_dead": theorem_k0_unsigned_still_dead()["proved"],
        "residual_i_closed_via_235": residual_i_closed_via_235(),
        "residual_i_dual_eq_empty": residual_i_dual_eq_empty_proved_general(),
        "gsum_disj_lb": gsum_disj_lb_proved_general(),
        "type_I": type_I_k_3p_minus_2_closed_general(),
        "e1": e1_closed_general(),
        "open_hinge": (
            "|R̄₄|≤B on |κ|=1 via μ-weighted signed layer sums "
            "(character sums before abs) or Aut-SOS+diamond / Path C / "
            "K₄ / free-e+ker=sc.  Layer types named; unsigned dead."
        ),
    }


def main() -> dict:
    a = theorem_k0_full_s4()
    b = theorem_k0_cycle_type_split()
    c = theorem_k0_four_cycle_inverse_pairing()
    d = theorem_k0_unsigned_still_dead()
    status = hinge_status_235()
    out = {
        "prop": "15.235",
        "title": (
            "Signed k=0: S₄ cycle-type + 4-cycle inverse pairing; "
            "unsigned still dead; residual OPEN"
        ),
        "proved": {
            "k0_full_s4": a["proved"],
            "k0_cycle_type_split": b["proved"],
            "k0_four_cycle_inverse_pairing": c["proved"],
            "k0_unsigned_still_dead": d["proved"],
            "residual_i_closed": residual_i_closed_via_235(),
        },
        "algebra": {
            "abs_per_k0_ub": 24,
            "class_sizes": b["class_sizes"],
            "k0_unsigned_majorant": "(p^2-3)(p^2-4)(p^2-5)(p^2-6)",
            "k0_majorant_p5": k0_unsigned_majorant(5),
            "B_p5": str(R_bar_budget_for_mu_thr(5)),
        },
        "hinge": status,
        "residual_i_closed": residual_i_closed_via_235(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "residual_i_dual_eq_empty_proved_general": residual_i_dual_eq_empty_proved_general(),
        "type_I_k_3p_minus_2_closed_general": type_I_k_3p_minus_2_closed_general(),
        "e1_closed_general": e1_closed_general(),
        "note": (
            "k=0 per is the full S₄ permanent; 4-cycles pair with inverses "
            "before abs.  Unsigned still exceeds B.  No predicate flip.  "
            "Signed-layer structure (k=3,2,1,0) is named."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15235.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"Prop 15.235  residual_i={residual_i_closed_via_235()}")
    print(
        f"  cycle_split={b['proved']} inv_pair={c['proved']} "
        f"unsigned_dead={d['proved']}"
    )
    print(
        f"  predicates gsum={gsum_disj_lb_proved_general()} "
        f"res_i={residual_i_dual_eq_empty_proved_general()} "
        f"type_I={type_I_k_3p_minus_2_closed_general()} "
        f"e1={e1_closed_general()}"
    )
    print(f"  wrote {path}")
    return out


if __name__ == "__main__":
    main()
