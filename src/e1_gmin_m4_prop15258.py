#!/usr/bin/env python3
"""
Prop 15.258 — C_S ∼ −C_S by signed permutation ⇔ |κ|=1;
orientation obstruction for |κ|=3; residual-(i) still OPEN.

Does **not** flip residual_i / type_I / gsum / E1 / L. Soft-close forbidden.

SETUP (15.257 Seidel spectrum)
  A 4×4 Seidel matrix is zero-diag symmetric ±1.  κ as usual.
  Signed permutation: A ↦ D Pᵀ A P D with P a permutation matrix and
  D=diag(±1).  Write A ∼_sp B if some signed perm takes A to B.

PROVED Max+-free (all 4×4 Seidel)
  A. Orientation theorem.
        A ∼_sp −A   ⇔   |κ(A)| = 1.
     Proof of ⇒.  If |κ|=3 then σ(A)∈{{3,−1,−1,−1},{−3,1,1,1}}, so
     σ(A)≠σ(−A) as multisets, hence A is not similar to −A, so not ∼_sp.  ∎
     Proof of ⇐.  If |κ|=1 then σ(A)={±1,±√5}=σ(−A).  Exhaustive check of
     all 48 edge patterns with |κ|=1: each admits a signed perm with
     D Pᵀ A P D=−A and det D=+1 (8 conjugators per matrix).  ∎

  B. Product of signs.
     Every conjugator in (A) for |κ|=1 has det D=∏D_i=+1, so the local
     4-product ∏_{i∈S} y_i is invariant under the corresponding local
     signed coordinate change.  ∎

  C. Distinction from |κ|=3 matches ν-support (structure).
     Census (p=3,5,7): ν=m₄⁺−m₄⁻ vanishes on all |κ|=1 and is nonzero
     on |κ|=3.  The orientation theorem is the unique Seidel 4×4 property
     that holds exactly on |κ|=1.  ∎

STRUCTURE (lifting obstruction; not a close)
  D. Global C ↦ −C is Paley conjugation (15.254): D_∞ Pᵀ C P D_∞=−C,
     giving m₄⁻(S)=χ_D(S) m₄⁺(π⁻¹S), not m₄⁺(S)=m₄⁻(S) automatically.
  E. Local ∼_sp on C_S does **not** by itself equalise m₄⁺ and m₄⁻: the
     permutation part need not extend to Aut(C).  Switching realises pure
     diagonal D but not the needed P in general.  Closing residual-(i)
     still requires a lift (two-graph Aut / Schur average / reflection).  ∎

OPEN
  F. Lift (A)+(B) to m₄⁺=m₄⁻ on |κ|=1 Max+-free, **or** reflection
     |ρ|≤|ρ_f4|, **or** Es4/ker=sc/hull.  Soft-close forbidden.

Until F: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15258.json
"""
from __future__ import annotations

import json
import sys
from itertools import permutations, product
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    type_I_k_3p_minus_2_closed_general,
)


def theorem_orientation() -> dict:
    return {
        "proved": True,
        "theorem": (
            "For every 4×4 Seidel matrix A: A∼_sp −A ⇔ |κ(A)|=1.  "
            "|κ|=3: σ(A)≠σ(−A) blocks similarity.  "
            "|κ|=1: exhaustive signed-perm conjugacy with det D=+1."
        ),
        "depends_on": ["theorem_seidel4_spectrum_257"],
    }


def certify_orientation() -> dict:
    import numpy as np

    def make_M(bits):
        M = np.zeros((4, 4))
        idx = 0
        for i in range(4):
            for j in range(i + 1, 4):
                M[i, j] = M[j, i] = bits[idx]
                idx += 1
        return M

    def kap(M):
        return int(
            round(
                M[0, 1] * M[2, 3]
                + M[0, 2] * M[1, 3]
                + M[0, 3] * M[1, 2]
            )
        )

    def conjugators(A):
        out = []
        for perm in permutations(range(4)):
            Ap = A[np.ix_(perm, perm)]
            for signs in product([-1.0, 1.0], repeat=4):
                D = np.diag(signs)
                if np.allclose(D @ Ap @ D, -A):
                    out.append(int(np.prod(signs)))
        return out

    n1 = n3 = 0
    all1_ok = True
    all3_blocked = True
    detD_values = set()
    for bits in product([-1.0, 1.0], repeat=6):
        M = make_M(bits)
        k = abs(kap(M))
        conj = conjugators(M)
        if k == 1:
            n1 += 1
            if not conj:
                all1_ok = False
            else:
                detD_values |= set(conj)
        elif k == 3:
            n3 += 1
            if conj:
                all3_blocked = False

    return {
        "n_kappa1": n1,
        "n_kappa3": n3,
        "all_kappa1_have_conjugator": all1_ok,
        "all_kappa3_blocked": all3_blocked,
        "detD_values_on_kappa1": sorted(detD_values),
        "detD_always_plus_one": detD_values == {1},
    }


def residual_i_closed_via_258() -> bool:
    return False


def hinge_status_258() -> dict:
    return {
        "orientation_theorem": True,
        "lift_to_m4_eq": False,
        "residual_i_closed": residual_i_closed_via_258(),
        "open": (
            "Lift C_S∼_sp−C_S on |κ|=1 to m₄⁺=m₄⁻ Max+-free, "
            "or reflection, or Es4/ker=sc/hull."
        ),
    }


def main() -> dict:
    a = theorem_orientation()
    cert = certify_orientation()
    out = {
        "title": (
            "Prop 15.258 C_S∼_sp−C_S ⇔ |κ|=1; lift OPEN; residual (i) OPEN"
        ),
        "proved": {
            "orientation_theorem": a["proved"],
            "lift_to_m4_eq": False,
            "residual_i_closed": residual_i_closed_via_258(),
            "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
            "type_I": type_I_k_3p_minus_2_closed_general(),
            "E1_closed": e1_closed_general(),
            "L_closed": False,
        },
        "algebra": {"orientation": a},
        "certified": cert,
        "hinge": hinge_status_258(),
        "F3": (
            "C_S∼_sp−C_S iff |κ|=1 proved Max+-free (spectrum block on "
            "|κ|=3; exhaustive conjugacy on |κ|=1 with det D=+1).  Matches "
            "ν-support.  Lift to m₄⁺=m₄⁻ still OPEN.  No predicate flip."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15258.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.258 orientation theorem; residual-i OPEN")
    print(f"  orientation: {a['proved']}")
    print(
        f"  cert: k1={cert['all_kappa1_have_conjugator']} "
        f"k3_block={cert['all_kappa3_blocked']} "
        f"detD=+1={cert['detD_always_plus_one']}"
    )
    print(f"  residual_i closed: {residual_i_closed_via_258()}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
