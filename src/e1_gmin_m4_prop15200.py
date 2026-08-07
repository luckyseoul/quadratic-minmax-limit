#!/usr/bin/env python3
"""
Prop 15.200 — C−2/n ∈ ker(Gsum) Max+-free; dual-eq ker-box restated;
residual (i) still OPEN.

Does **not** flip residual_i / gsum / E1 / L. Soft-close forbidden.

SETUP (15.182, 15.189–192)
  Dual equality ⇔ ∃κ∈ker(Gsum) with κ_e=2−α, 1ᵀκ=0, and
  −α≤κ_g≤1−α for all g≠e (box on freeness-edge free).
  Equivalent: max{κ_e : Gsum κ=0, 1ᵀκ=0, box on g≠e} < 2−α
  blocks dual-eq (Aut_e averaging: 15.192).
  Scheme-ker κ_ij=f_i+f_j (∑f=0) ⊆ ker (15.190); max scheme κ_e
  =3(n−2)/(pn)<2−α. Census full max κ_e <2−α at p=5,7; > at p=3.

PROVED Max+-free
  A. Edge feature identity. For y∈Max± and edge vector κ,
        f(y)ᵀκ = ∑_{i<j} κ_ij C_ij y_i y_j = yᵀ A y
     with A_ij=κ_ij C_ij (i≠j), A_ii=0.  ∎

  B. Ker via vanishing quadratics.  Gsum=∑_{Max±} w_y f(y)f(y)ᵀ with w>0
     ⇒ Gsum κ=0 ⇔ f(y)ᵀκ=0 for every y∈Max±
     ⇔ yᵀ A y=0 on Max±.  ∎

  C. C−2/n is in ker(Gsum).  Let κ_ij = C_ij − 2/n (i≠j as edges).
     Then for y∈Max+:
        f(y)ᵀκ = ∑_{i<j} y_i y_j − (2/n) ∑_{i<j} C_ij y_i y_j
                 = ((1ᵀy)² − n)/2 − (2/n)·(yᵀ C y)/2
                 = ((1ᵀy)² − n)/2 − (1/n)·(p n)   (Cy=py)
                 = ((1ᵀy)² − n)/2 − p.
     By 15.189: 1ᵀy=(p+1)y_∞ ⇒ (1ᵀy)²=(p+1)², so
        ((p+1)² − (p²+1))/2 − p = (2p)/2 − p = 0.
     For y∈Max−: Cy=−py, 1ᵀy=−(p−1)y_∞, (1ᵀy)²=(p−1)²,
        fᵀκ = ((1ᵀy)²−n)/2 − (2/n)·(−p n)/2 = ((1ᵀy)²−n)/2 + p
             = ((p−1)² − p² − 1)/2 + p = −p + p = 0.
     Hence f(y)ᵀκ=0 on Max± ⇒ κ∈ker(Gsum).  Uses only 15.189 sum
     identity + Cy=±py (no Max+ enumeration).  ∎

  D. Dual-eq restated as free-e ker-box.  (15.182 algebra recalled)
     Dual-eq feasible iff max{κ_e: Gsumκ=0, 1ᵀκ=0, −α≤κ_g≤1−α (g≠e)}
     ≥ 2−α.  Box does **not** constrain κ_e.  ∎

CERTIFIED
  E. C−2/n ∈ ker at p=3,5 (‖Gsum κ‖≈0); scheme-ker still in ker.
  F. Free-e ker-box LP: max κ_e ≈2.30 (p=3) ≥ need 1.8 (feasible);
     max κ_e ≈0.805 (p=5) < need 1.95 (blocked) — matches 15.190 scale.

OPEN (blocks residual i / predicate flip)
  G. Prove for all primes p≥5 that
        max{κ_e : Gsumκ=0, 1ᵀκ=0, −α≤κ_g≤1−α (g≠e)} < 2−α
     e.g. via Aut_e orbit algebra + closed Gsum, or max κ_e ≤ (3/2)·scheme-max
     (15.192: (3/2)·scheme < need for p≥5), or |μ|≤2/n / K₄≤Wick_hi.
  Scheme + span{C−2/n} alone is too small (misses cross-ker that lifts max κ_e).

Until G: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15200.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    is_prime,
    n_of,
)
from e1_gmin_m4_prop15190 import (  # noqa: E402
    alpha_particular,
    scheme_ker_max_kappa_e,
)


def kappa_C_minus_two_over_n(p: int, edges: list[tuple[int, int]], C) -> "object":
    """Edge vector κ_ij = C_ij − 2/n."""
    import numpy as np

    n = n_of(p)
    return np.array([C[i, j] - 2.0 / n for i, j in edges], dtype=float)


def theorem_C_minus_2n_in_ker() -> dict:
    return {
        "proved": True,
        "theorem": (
            "κ_ij=C_ij−2/n lies in ker(Gsum) for Paley Max±. "
            "Proof: f(y)ᵀκ=((1ᵀy)²−n)/2 − (sign)·p with sign=+ on Max+ / − on Max− "
            "from yᵀCy=±pn; plug 1ᵀy=±(p±1)y_∞ (15.189) ⇒ fᵀκ=0 on Max± "
            "⇒ Gsum κ=0. Max+-free (no enumeration)."
        ),
        "depends_on": ["prop15189_sum_identity", "Cy_pm_p_y"],
    }


def theorem_dual_eq_free_e_box() -> dict:
    return {
        "proved": True,
        "theorem": (
            "Dual equality feasible ⇔ "
            "max{κ_e: Gsumκ=0, 1ᵀκ=0, −α≤κ_g≤1−α (g≠e)} ≥ 2−α. "
            "Box does not constrain the freeness edge e. (15.182.)"
        ),
    }


def theorem_three_halves_scheme_blocks_if_max_le() -> dict:
    """If max κ_e ≤ (3/2) scheme_max then dual-eq blocked for p≥5."""
    primes = [p for p in range(5, 80) if is_prime(p)]
    rows = {}
    ok = True
    for p in primes:
        sch = scheme_ker_max_kappa_e(p)
        thr = Fraction(3, 2) * sch
        need = 2 - alpha_particular(p)
        rows[str(p)] = {
            "scheme_max": str(sch),
            "three_halves": str(thr),
            "need_2_minus_alpha": str(need),
            "three_halves_lt_need": thr < need,
        }
        if not (thr < need):
            ok = False
    return {
        "proved_implication": ok,
        "max_bound_proved": False,
        "theorem": (
            "For all primes p≥5: (3/2)·scheme_max < 2−α (15.192 cubic). "
            "Hence IF max κ_e ≤ (3/2)·scheme_max on the free-e ker-box, "
            "THEN dual-eq is empty. The IF is OPEN (census ratios ~1.46 at p=5,7)."
        ),
        "by_p_sample": {k: rows[k] for k in ("5", "7", "11") if k in rows},
    }


def certify_C_minus_2n_ker(p: int) -> dict:
    import numpy as np
    from e1_bitight_infeas import boolean_evecs
    from minmax_quadratic import paley_conference_prime_power

    if p > 5:
        return {"p": p, "skipped": True}
    C = paley_conference_prime_power(p).astype(float)
    n = C.shape[0]
    Mp = boolean_evecs(C, float(p), 0).astype(np.float64)
    Mm = boolean_evecs(C, float(-p), 0).astype(np.float64)
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    kappa = kappa_C_minus_two_over_n(p, edges, C)

    def feats(M):
        return np.stack([C[i, j] * M[:, i] * M[:, j] for i, j in edges], 1)

    G = feats(Mp).T @ feats(Mp) / len(Mp) + feats(Mm).T @ feats(Mm) / len(Mm)
    gnorm = float(np.linalg.norm(G @ kappa))
    # also check f(y)^T kappa = 0 pointwise on a sample
    Fp, Fm = feats(Mp), feats(Mm)
    max_f = max(
        float(np.max(np.abs(Fp @ kappa))),
        float(np.max(np.abs(Fm @ kappa))),
    )
    return {
        "p": p,
        "||Gsum kappa||": gnorm,
        "max_|f^T kappa|_Maxpm": max_f,
        "in_ker": gnorm < 1e-8 and max_f < 1e-8,
        "sum_kappa": float(kappa.sum()),
    }


def certify_free_e_kerbox(p: int) -> dict:
    """LP: max κ_e free, box on g≠e, Gsumκ=0, 1ᵀκ=0."""
    import numpy as np
    from e1_bitight_infeas import boolean_evecs
    from minmax_quadratic import paley_conference_prime_power
    from scipy.optimize import linprog

    if p > 5:
        return {"p": p, "skipped": True}
    C = paley_conference_prime_power(p).astype(float)
    n = C.shape[0]
    al = float(alpha_particular(p))
    need = float(2 - alpha_particular(p))
    Mp = boolean_evecs(C, float(p), 0).astype(np.float64)
    Mm = boolean_evecs(C, float(-p), 0).astype(np.float64)
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    E = len(edges)
    e0 = edges.index((0, 1))

    def feats(M):
        return np.stack([C[i, j] * M[:, i] * M[:, j] for i, j in edges], 1)

    G = feats(Mp).T @ feats(Mp) / len(Mp) + feats(Mm).T @ feats(Mm) / len(Mm)
    ev, vec = np.linalg.eigh(G)
    K = vec[:, ev < 1e-8]
    ones = np.ones(E)
    sn = ones @ K
    ns = np.linalg.svd([sn], full_matrices=True)[2][1:].T
    N = K @ ns
    free = N.shape[1]
    rows, bs = [], []
    for i in range(E):
        if i == e0:
            continue
        rows += [N[i], -N[i]]
        bs += [1 - al, al]
    res = linprog(
        -N[e0],
        A_ub=np.array(rows),
        b_ub=np.array(bs),
        method="highs",
        bounds=[(None, None)] * free,
    )
    if not res.success:
        return {"p": p, "success": False, "message": res.message}
    vmax = float(N[e0] @ res.x)
    return {
        "p": p,
        "success": True,
        "max_kappa_e": vmax,
        "need_2_minus_alpha": need,
        "blocks_dual_eq": vmax < need - 1e-9,
        "proved_general": False,
    }


def K4_bound_proved_general() -> bool:
    return False


def kerbox_empty_proved_general() -> bool:
    return False


def residual_i_closed_via_kerbox() -> bool:
    return kerbox_empty_proved_general()


def hinge_status_200() -> dict:
    return {
        "C_minus_2n_in_ker": True,
        "dual_eq_free_e_box_form": True,
        "three_halves_scheme_implication": True,
        "kerbox_empty_general": False,
        "residual_i_closed": residual_i_closed_via_kerbox(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "open": (
            "Prove free-e ker-box max κ_e < 2−α for all primes p≥5 "
            "(e.g. ≤(3/2)·scheme-max), or |μ|≤2/n / K₄≤Wick_hi."
        ),
    }


def main() -> dict:
    a = theorem_C_minus_2n_in_ker()
    b = theorem_dual_eq_free_e_box()
    c = theorem_three_halves_scheme_blocks_if_max_le()
    cert_ker = {str(p): certify_C_minus_2n_ker(p) for p in (3, 5)}
    cert_box = {str(p): certify_free_e_kerbox(p) for p in (3, 5)}
    out = {
        "title": (
            "Prop 15.200 C−2/n∈ker(Gsum) Max+-free; free-e dual-eq box; "
            "ker-box empty OPEN"
        ),
        "proved": {
            "C_minus_2n_in_ker": a["proved"],
            "dual_eq_free_e_box": b["proved"],
            "three_halves_implication": c["proved_implication"],
            "kerbox_empty_general": False,
            "residual_i_closed": residual_i_closed_via_kerbox(),
            "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
            "E1_closed": e1_closed_general(),
            "L_closed": False,
        },
        "algebra": {"C_ker": a, "dual_eq_box": b, "three_halves": c},
        "certified": {"C_ker": cert_ker, "free_e_box": cert_box},
        "hinge": hinge_status_200(),
        "F3": (
            "C−2/n∈ker is Max+-free. Do not flip residual/E1/L without "
            "general free-e ker-box emptiness or |μ|/K₄ bound."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15200.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.200 C-2/n in ker(Gsum); dual-eq free-e box")
    print(f"  C-2/n in ker proved: {a['proved']}")
    print(f"  free-e box form: {b['proved']}")
    print(f"  (3/2)scheme implication: {c['proved_implication']}")
    for pk, v in cert_ker.items():
        if v.get("skipped"):
            continue
        print(f"  p={pk} C-ker cert in_ker={v['in_ker']}")
    for pk, v in cert_box.items():
        if v.get("skipped") or not v.get("success"):
            continue
        print(
            f"  p={pk} max_ke={v['max_kappa_e']:.4f} need={v['need_2_minus_alpha']:.4f} "
            f"blocks={v['blocks_dual_eq']}"
        )
    print(f"  residual_i closed: {residual_i_closed_via_kerbox()}")
    print(f"  open: {out['hinge']['open']}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
