#!/usr/bin/env python3
"""
Prop 15.628 — Eligible GQR circles close W1, W2, and Walsh.

This is a genuine general-p close of Walsh 15.406 E.  It does not close
the unrelated 5+-level / even-k>4p or other E1 leftovers.

============================================================================
Theorem A — PROVED (all odd primes).
  Nonsquare Miquelian circles meeting e={infinity,0} evenly span
      H0 intersect ker(e_infinity+e_0).
  Tangent pencils reduce the quotient to connectivity of a bipartite Cayley
  graph on M=(F_(p^2)^*)^2 with connection set
      {(1+r sigma)^(-2): r in F_p}.
  Failure to generate M would make a nontrivial even multiplicative
  character constant on an affine F_p-line.  A norm-one-torus
  parameterisation turns the line sum into Katz's t=-2 Soto-Andrade sum.
  Its two exceptional character pairs induce the trivial extension-field
  character, so |sum|<=2sqrt(p)<p for p>=5.  The p=3 group is direct.

Theorem B — PROVED (all odd primes).
  If T subset F_p has |T|=(p+1)/2 and L:F_(p^2)->F_p has square kernel,
  then h_T(infinity)=1 and h_T(u)=2*1_T(Lu)-1 satisfies C h_T=p h_T.
  At a finite point, the zero L-fibre contributes (p-1)s_T(b), every
  nonzero affine square fibre has character sum -1, and sum s_T=1:
      (C h_T)_u = 1+(p-1)s_T(b)-sum_(t!=0)s_T(b+t) = p s_T(b).
  Nonsquare dilation plus the infinity sign switch gives a Max-minus
  completion of the standard nonsquare circle whenever 0 is in T.

Theorem C — PROVED (p>=5; p=3 direct).
  For an outside pair sigma*u,sigma*v, with transverse coordinates b,d,
  the U condition is s_T(b)s_T(d)=chi(u-v).  If b=d, chi(u-v)=1.  If
  b!=d, prescribe equal membership for target +1 and opposite membership
  for target -1, then extend T to size (p+1)/2 while retaining 0 in T.
  Thus every eligible circle is an actual difference of two U-points.
  Theorem A gives
      dir(U)=H0 intersect ker(e_infinity+e_0).
  At p=3 the exact rank is 4 and W2 is vacuous.  Hence W1, W2, and Walsh
  15.406 E hold for every odd prime.

Exact certificates: scripts/w2_affine_circle_close.py.
Long proof: evidence/NOTE_2026-08-24_w2_gqr_circle_route.md.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))


def theorem_A_eligible_span() -> dict:
    return {
        "proved": True,
        "all_odd_primes": True,
        "dimension": "(p^2-1)/2",
        "theorem": (
            "Eligible nonsquare Miquelian circles span "
            "H0 intersect ker(e_infinity+e_0)."
        ),
        "proof": "tangent pencils + connected Cayley tangency graph + Katz",
    }


def theorem_B_affine_completions() -> dict:
    return {
        "proved": True,
        "all_odd_primes": True,
        "theorem": (
            "Every size-(p+1)/2 affine halfspace on a square-kernel "
            "functional is a +p eigenvector; nonsquare transport gives "
            "standard-circle Max-minus completions."
        ),
        "finite_row_identity": (
            "1+(p-1)s(b)-sum_(t!=0)s(b+t)=p*s(b)"
        ),
    }


def theorem_C_walsh_close() -> dict:
    return {
        "proved": True,
        "all_odd_primes": True,
        "W1": True,
        "W2": True,
        "Walsh_15_406_E": True,
        "p3_direct_rank": 4,
        "theorem": (
            "Every eligible circle is a U-difference, hence "
            "dir(U)=H0 intersect ker(e_infinity+e_0)."
        ),
        "not_closed": [
            "5+-level and even-k>4p leftovers",
            "other E1 leftovers",
            "L",
        ],
    }


def main() -> dict:
    from w2_affine_circle_close import run

    checks = run([5, 7, 11, 19])
    out = {
        "prop": "15.628",
        "title": "Eligible GQR circles + affine completions close W1/W2/Walsh",
        "proved": {
            "eligible_circle_span_all_odd_p": True,
            "affine_completion_all_odd_p": True,
            "every_eligible_circle_U_difference_p_ge_5": True,
            "W1_all_odd_p": True,
            "W2_all_odd_p": True,
            "Walsh_15_406_E_all_odd_p": True,
        },
        "algebra": {
            "A": theorem_A_eligible_span(),
            "B": theorem_B_affine_completions(),
            "C": theorem_C_walsh_close(),
        },
        "checks": checks,
        "L_status": "OPEN",
        "residual_ii_k_eq_4p_empty": False,
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15628.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print(f"wrote {dest}", flush=True)
    return out


if __name__ == "__main__":
    main()
