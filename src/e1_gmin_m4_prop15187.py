#!/usr/bin/env python3
"""
Prop 15.187 — Global T²κ = −24φ + 48κ (any conference, all 4-sets);
T²φ = 4(p²+3)(φ−2κ) for Paley; residual-(i) still needs |μ_actual|.

Does **not** flip gsum_disj_lb_proved_general. Soft-close forbidden.

PROVED Max+-free
  A. Global T²κ identity (any symmetric conference C, all 4-sets S):
        T²κ(S) = −24 φ(S) + 48 κ(S) = −24(φ(S) − 2 κ(S)).
     Proof. Same C² reduction as 15.184: Tκ = −6 star, expand
     T(Tκ)(S) = −6 ∑_v ∑_r C_vr star(S_{v→r}); leading terms give −24φ;
     K4 remainder equals +48κ by exhaustion of **all 64** edge-labelings
     (both |κ|=1 and |κ|=3). ∎
     (Strengthens 15.184 which stated the identity only on |κ|=1.)

  B. T²φ identity for Paley order n=p²+1 (certified p=3,5 full enum;
     sample p=7): 
        T²φ = 4(p²+3)(φ − 2κ) = 4(n+2)(φ − 2κ).
     Open as pure symbolic C² reduction (φ depends on n external vertices).

  C. Consequences for the master (16p²I − T²)μ₄ = 16κ:
     span{κ, φ} is T²-invariant under (A)+(B). Unique solution in this span:
        μ_part = [(p²−1)κ − 2φ] / (p²(p²−5))   (p²≠5)
     with |μ_part| ≤ (p²+4p−9)/(p²(p²−5)) ≤ 1/(2p) on |κ|=1 when
     |φ|≤2(p−2) (15.186). ∎

  D. Actual Max± μ₄ = μ_part + δ with δ ∈ E_{±4p}(T) when 4p ∈ σ(T).
     At p=5: Aut-invariant δ ≠ 0 on |κ|=1 (actual=(4κ−φ)/(pn) ≠ μ_part).
     At p=5,7: |μ_actual| ≤ 1/(2p) by census (15.176–177 evidence).

OPEN (blocks residual i)
  E. Prove |μ_actual| ≤ (p²+4p−9)/(p²(p²−5)) (or ≤1/(2p)) on |κ|=1
     for all primes p≥5 — i.e. control Aut-invariant δ ∈ E_{±4p}.
     Then Gsum ≥ −1/p, residual-(i) Farkas, flip predicates.

Until E: gsum_disj_lb_proved_general False; residual_i False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15187.json
"""
from __future__ import annotations

import json
import sys
from itertools import combinations, product
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    is_prime,
)
from e1_gmin_m4_prop15186 import (  # noqa: E402
    majorant_mu_part,
    theorem_particular_majorant_beats_thr,
    theorem_phi_bound_general,
)


def theorem_global_T2_kappa() -> dict:
    """T²κ = −24φ + 48κ on every 4-set, any conference (64-labeling)."""
    import sympy as sp

    E = {
        (i, j): sp.Symbol(f"e{i}{j}", commutative=True)
        for i, j in combinations(range(4), 2)
    }

    def EE(i, j):
        if i == j:
            return 0
        return E[tuple(sorted((i, j)))]

    def sum_r_Cvr_star_nonphi(v: int):
        others = [x for x in range(4) if x != v]
        b, c, d = others
        sum_rb = -EE(v, c) * EE(c, b) - EE(v, d) * EE(d, b)
        term2 = EE(b, c) * EE(b, d) * sum_rb
        sum_rc = -EE(v, b) * EE(b, c) - EE(v, d) * EE(d, c)
        term3 = EE(b, c) * EE(c, d) * sum_rc
        sum_rd = -EE(v, b) * EE(b, d) - EE(v, c) * EE(c, d)
        term4 = EE(b, d) * EE(c, d) * sum_rd
        return sp.expand(term2 + term3 + term4)

    total_nonphi = sum(sum_r_Cvr_star_nonphi(v) for v in range(4))
    rest = sp.expand(-6 * total_nonphi)
    kap = EE(0, 1) * EE(2, 3) + EE(0, 2) * EE(1, 3) + EE(0, 3) * EE(1, 2)
    edges = list(E.values())
    ok_all = True
    n1 = n3 = 0
    for bits in product([-1, 1], repeat=6):
        subs = dict(zip(edges, bits))
        kv = int(kap.subs(subs))
        rv = int(rest.subs(subs))
        if rv != 48 * kv:
            ok_all = False
        if abs(kv) == 1:
            n1 += 1
        elif abs(kv) == 3:
            n3 += 1
    return {
        "proved": ok_all and n1 == 48 and n3 == 16,
        "theorem": (
            "Any symmetric conference: T²κ = −24φ + 48κ on every 4-set. "
            "C² reduction of T(−6 star); K4 remainder = 48κ on all 64 "
            "edge-labelings (|κ|=1 and |κ|=3)."
        ),
        "n_kappa1_labelings": n1,
        "n_kappa3_labelings": n3,
        "rest_equals_48_kappa_all": ok_all,
    }


def certify_T2_phi_formula(primes: list[int] | None = None) -> dict:
    """Census T²φ = 4(p²+3)(φ−2κ) for Paley."""
    import numpy as np
    from minmax_quadratic import paley_conference_prime_power
    from scipy import sparse

    if primes is None:
        primes = [3, 5]
    rows = {}
    all_ok = True
    for p in primes:
        C = paley_conference_prime_power(p).astype(np.int8)
        n = C.shape[0]
        sets = list(combinations(range(n), 4))
        N = len(sets)
        idx = {S: i for i, S in enumerate(sets)}
        kap = np.zeros(N)
        phi = np.zeros(N)
        for S, i in idx.items():
            a, b, c, d = S
            kap[i] = (
                C[a, b] * C[c, d]
                + C[a, c] * C[b, d]
                + C[a, d] * C[b, c]
            )
            phi[i] = sum(
                int(np.prod([C[r, j] for j in S]))
                for r in range(n)
                if r not in S
            )
        rows_i = []
        cols_i = []
        data = []
        for S, i in idx.items():
            Sl = list(S)
            for pos, v in enumerate(Sl):
                for r in range(n):
                    if r in S:
                        continue
                    Sp = Sl.copy()
                    Sp[pos] = r
                    j = idx[tuple(sorted(Sp))]
                    rows_i.append(i)
                    cols_i.append(j)
                    data.append(float(C[v, r]))
        T = sparse.csr_matrix(
            (data, (rows_i, cols_i)), shape=(N, N), dtype=np.float64
        )
        T2phi = T @ (T @ phi)
        expect = 4 * (p * p + 3) * (phi - 2 * kap)
        err = float(np.max(np.abs(T2phi - expect)))
        ok = err < 1e-8
        if not ok:
            all_ok = False
        rows[str(p)] = {"max_err": err, "ok": ok, "formula": "4(p^2+3)(phi-2kappa)"}
    return {
        "certified": all_ok,
        "proved_general": False,
        "by_p": rows,
        "note": "Full enum p=3,5. Symbolic C² proof of T²φ still open.",
    }


def hinge_status_187() -> dict:
    return {
        "residual_ii_closed": True,
        "residual_i_closed": False,
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "global_T2_kappa_proved": True,
        "phi_bound_all_p": theorem_phi_bound_general()["proved"],
        "particular_majorant_le_thr": theorem_particular_majorant_beats_thr()[
            "proved"
        ],
        "actual_mu4_bound": False,
        "open": (
            "Control δ=μ_actual−μ_part on |κ|=1 for all p≥5 "
            "(or direct |μ_actual|≤1/(2p)); then flip residual_i."
        ),
        "majorant_p5": str(majorant_mu_part(5)),
        "majorant_p7": str(majorant_mu_part(7)),
    }


def main() -> dict:
    t2k = theorem_global_T2_kappa()
    t2phi = certify_T2_phi_formula([3, 5])
    maj = theorem_particular_majorant_beats_thr(
        [p for p in range(5, 40) if is_prime(p)]
    )
    phi = theorem_phi_bound_general()
    out = {
        "title": (
            "Prop 15.187 global T²κ; T²φ census; "
            "μ_part majorant; residual (i) OPEN (need |μ_actual|)"
        ),
        "proved": {
            "global_T2_kappa": t2k["proved"],
            "T2_phi_paley_general": False,
            "T2_phi_certified_p3_p5": t2phi["certified"],
            "phi_bound_all_p": phi["proved"],
            "particular_majorant_le_1_over_2p": maj["proved"],
            "actual_mu4_le_1_over_2p": False,
            "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
            "residual_i_closed": False,
            "E1_closed": e1_closed_general(),
            "L_closed": False,
        },
        "algebra": {"global_T2_kappa": t2k, "particular_majorant": maj, "phi": phi},
        "certified": {"T2_phi": t2phi},
        "hinge": hinge_status_187(),
        "F3": "Do not flip residual/E1/L without |μ_actual|≤1/(2p).",
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15187.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.187 global T²κ; residual (i) OPEN")
    print(f"  global T2κ proved={t2k['proved']}")
    print(f"  T2φ census p3,5={t2phi['certified']}")
    print(f"  phi bound={phi['proved']}; majorant={maj['proved']}")
    print(f"  gsum={gsum_disj_lb_proved_general()}; e1={e1_closed_general()}")
    print(f"  open: {out['hinge']['open']}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
