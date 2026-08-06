#!/usr/bin/env python3
"""
Prop 15.188 — T²φ identity (Paley census); μ_part solves master; p=7 Max±
census; corrected residual-(i) target |μ|≤2/n; hinge still OPEN.

Does **not** flip gsum_disj_lb_proved_general. Soft-close forbidden.

PROVED Max+-free
  A. Global T²κ = −24φ + 48κ (15.187; any conference, all 4-sets).

  B. T²φ = 4(p²+3)(φ − 2κ) for Paley n=p²+1: certified by full enum at
     p=3,5 and full sparse T² at p=7 (max abs residual 0). Symbolic C²
     reduction still open as pure algebra; identity used as certified for
     p∈{3,5,7} and the master algebra below.

  C. Consequently on span{κ,φ}, the unique master solution is
        μ_part = [(p²−1)κ − 2φ] / (p²(p²−5))   (p²≠5 as integer always)
     and (16p² I − T²) μ_part = 16κ exactly whenever (A)+(B) hold.
     Certified residual of master on μ_part: max|res| < 1e-12 at p=3,5,7.

  D. Particular majorant (15.186): on |κ|=1 with |φ|≤2(p−2),
        |μ_part| ≤ (p²+4p−9)/(p²(p²−5)) ≤ 1/(2p) for all primes p≥5.

  E. Target comparison 2/n ≤ 1/(2p):
        2/(p²+1) ≤ 1/(2p)  ⇔  4p ≤ p²+1  ⇔  p²−4p+1 ≥ 0,
     true for all primes p≥5 (value 6 at p=5). ∎

CERTIFIED (Max± census; not general)
  F. p=5: μ_actual = (4κ−φ)/(p n) exactly on |κ|=1; max|μ|=3/65 < 1/10;
     and 3/65 < 2/n = 1/13.
  G. p=7: |Max±|=11452 each; on |κ|=1, max|μ_actual|=109/2863 ≈ 0.03807
     < 1/14 and < 2/n = 1/25. μ is **not** constant on every (κ,φ) class
     (std > 0 on ( ±1, ∓10 )); mean close to μ_part but **not equal**.
  H. **Correction:** |μ_actual| is **not** ≤ |μ_part| pointwise at p=7
     (≈29400 of 176400 |κ|=1 sets have |act| > |part|). The open claim in
     15.186 that |μ_actual|≤ majorant_of_μ_part is **false** as a general
     inequality. The viable sufficient target is
        |μ_actual| ≤ 2/n = 2/(p²+1) on |κ|=1   (⇒ ≤1/(2p) for p≥5),
     which holds at p=5,7 by census and fails correctly at p=3 (1/3 > 1/5).

OPEN (blocks residual i / predicate flip)
  I. Prove Max+-free |μ_actual| ≤ 2/(p²+1) (or ≤1/(2p)) on every |κ|=1
     four-set for all primes p≥5 — without Max± enumeration.
     Then Gsum ≥ −1/p (15.176) ⇒ residual-(i) Farkas ⇒ flip predicates.

Until I: gsum_disj_lb_proved_general False; residual_i False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15188.json
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
from e1_gmin_m4_prop15176 import theorem_minus_1_over_p_suffices  # noqa: E402
from e1_gmin_m4_prop15177 import mu4_threshold  # noqa: E402
from e1_gmin_m4_prop15186 import (  # noqa: E402
    majorant_mu_part,
    mu_part_formula,
    theorem_particular_majorant_beats_thr,
    theorem_phi_bound_general,
)
from e1_gmin_m4_prop15187 import (  # noqa: E402
    certify_T2_phi_formula,
    theorem_global_T2_kappa,
)


def two_over_n(p: int) -> Fraction:
    """2/n = 2/(p²+1)."""
    return Fraction(2, n_of(p))


def theorem_two_over_n_le_mu4_thr(primes: list[int] | None = None) -> dict:
    """Max+-free Fraction: 2/(p²+1) ≤ 1/(2p) for all primes p≥5."""
    if primes is None:
        primes = [p for p in range(5, 80) if is_prime(p)]
    rows = {}
    ok = True
    for p in primes:
        a = two_over_n(p)
        thr = mu4_threshold(p)
        beats = a <= thr
        # cubic form: p² − 4p + 1 ≥ 0
        disc = p * p - 4 * p + 1
        rows[str(p)] = {
            "two_over_n": str(a),
            "mu4_thr": str(thr),
            "two_over_n_le_thr": beats,
            "p2_minus_4p_plus_1": disc,
        }
        if not beats or disc < 0:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "2/(p²+1) ≤ 1/(2p) for all primes p≥5 "
            "⇔ p²−4p+1 ≥ 0 (value 6 at p=5; increasing for p≥3)."
        ),
        "by_p_sample": {k: rows[k] for k in list(rows)[:6]},
        "n_primes": len(primes),
    }


def theorem_two_over_n_implies_gsum_and_farkas() -> dict:
    """|μ|≤2/n on |κ|=1 ⇒ |μ|≤1/(2p) ⇒ Gsum≥−1/p ⇒ residual-i Farkas."""
    return {
        "proved_implication": True,
        "bound_proved_general": False,
        "chain": (
            "|μ|≤2/n ≤1/(2p) (p≥5) ⇒ Gsum≥−2·(1/(2p))=−1/p "
            "⇒ need 6/p−4 < k·(−1/p) for residual-(i) Farkas (15.176)."
        ),
        "minus_1_over_p_suffices": theorem_minus_1_over_p_suffices(),
    }


def certify_mu_part_solves_master(primes: list[int] | None = None) -> dict:
    """Census: (16p²I−T²)μ_part = 16κ when T²φ identity holds."""
    import numpy as np
    from itertools import combinations

    from minmax_quadratic import paley_conference_prime_power
    from scipy import sparse

    if primes is None:
        primes = [3, 5, 7]
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
        for i, S in enumerate(sets):
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
        rows_i, cols_i, data = [], [], []
        for i, S in enumerate(sets):
            for pos, v in enumerate(S):
                for r in range(n):
                    if r in S:
                        continue
                    Sp = list(S)
                    Sp[pos] = r
                    j = idx[tuple(sorted(Sp))]
                    rows_i.append(i)
                    cols_i.append(j)
                    data.append(float(C[v, r]))
        T = sparse.csr_matrix(
            (data, (rows_i, cols_i)), shape=(N, N), dtype=np.float64
        )
        den = float(p * p * (p * p - 5))
        mu_part = ((p * p - 1) * kap - 2 * phi) / den
        res = 16.0 * p * p * mu_part - (T @ (T @ mu_part)) - 16.0 * kap
        err = float(np.max(np.abs(res)))
        ok = err < 1e-10
        if not ok:
            all_ok = False
        rows[str(p)] = {"max_master_residual": err, "ok": ok}
    return {
        "certified": all_ok,
        "proved_general": False,
        "note": "Holds when T²φ=4(p²+3)(φ−2κ); certified p=3,5,7.",
        "by_p": rows,
    }


def certified_p7_mu_census_constants() -> dict:
    """
    Hard-coded certified constants from full Max± enum at p=7
    (|M±|=11452; full |κ|=1 pass). Re-runnable offline; not a general proof.
    """
    return {
        "p": 7,
        "n": 50,
        "abs_Max_plus": 11452,
        "abs_Max_minus": 11452,
        "n_kappa1": 176400,
        "max_abs_mu_kappa1": str(Fraction(109, 2863)),
        "max_abs_mu_float": float(Fraction(109, 2863)),
        "mu4_thr": str(mu4_threshold(7)),
        "two_over_n": str(two_over_n(7)),
        "majorant_mu_part": str(majorant_mu_part(7)),
        "beats_mu4_thr": Fraction(109, 2863) <= mu4_threshold(7),
        "beats_two_over_n": Fraction(109, 2863) <= two_over_n(7),
        "not_equal_mu_part": True,
        "not_bounded_by_mu_part_pointwise": True,
        "n_sets_act_gt_part_abs": 29400,
        "varies_within_kappa_phi_class": True,
        "not_general_proof": True,
    }


def certified_p5_mu_on_kappa1() -> dict:
    """p=5: actual = (4κ−φ)/(pn) on |κ|=1 (full census)."""
    return {
        "p": 5,
        "formula_on_kappa1": "(4*kappa - phi)/(p*n)",
        "max_abs_mu": str(Fraction(3, 65)),
        "mu4_thr": str(mu4_threshold(5)),
        "two_over_n": str(two_over_n(5)),
        "beats_mu4_thr": Fraction(3, 65) <= mu4_threshold(5),
        "beats_two_over_n": Fraction(3, 65) <= two_over_n(5),
        "not_general_proof": True,
    }


def hinge_status_188() -> dict:
    return {
        "residual_ii_closed": False,  # full open (15.193); affine branch closed by 15.179
        "residual_i_closed": False,
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "phi_bound_all_p": theorem_phi_bound_general()["proved"],
        "particular_majorant_le_thr": theorem_particular_majorant_beats_thr()[
            "proved"
        ],
        "two_over_n_le_thr": theorem_two_over_n_le_mu4_thr()["proved"],
        "actual_mu4_le_two_over_n_general": False,
        "actual_mu4_le_1_over_2p_general": False,
        "open": (
            "Prove |μ_actual|≤2/(p²+1) (or ≤1/(2p)) on |κ|=1 for all primes "
            "p≥5 Max+-free; then flip residual_i / gsum / E1 / L."
        ),
        "corrected": (
            "|μ_actual|≰|μ_part| pointwise (false at p=7); use 2/n target."
        ),
    }


def main() -> dict:
    t2k = theorem_global_T2_kappa()
    t2phi = certify_T2_phi_formula([3, 5])
    part_master = certify_mu_part_solves_master([3, 5])  # p=7 heavy; optional
    two_n = theorem_two_over_n_le_mu4_thr()
    maj = theorem_particular_majorant_beats_thr(
        [p for p in range(5, 40) if is_prime(p)]
    )
    phi = theorem_phi_bound_general()
    p5 = certified_p5_mu_on_kappa1()
    p7 = certified_p7_mu_census_constants()
    out = {
        "title": (
            "Prop 15.188 T²φ census; μ_part master; p=7 Max±; "
            "target |μ|≤2/n; residual (i) OPEN"
        ),
        "proved": {
            "global_T2_kappa": t2k["proved"],
            "T2_phi_certified_p3_p5": t2phi["certified"],
            "mu_part_solves_master_certified": part_master["certified"],
            "phi_bound_all_p": phi["proved"],
            "particular_majorant_le_1_over_2p": maj["proved"],
            "two_over_n_le_1_over_2p": two_n["proved"],
            "actual_mu4_le_two_over_n_general": False,
            "actual_mu4_le_1_over_2p_general": False,
            "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
            "residual_i_closed": False,
            "E1_closed": e1_closed_general(),
            "L_closed": False,
        },
        "algebra": {
            "global_T2_kappa": t2k,
            "two_over_n_le_thr": two_n,
            "particular_majorant": maj,
            "phi": phi,
            "implication": theorem_two_over_n_implies_gsum_and_farkas(),
        },
        "certified": {
            "T2_phi": t2phi,
            "mu_part_master": part_master,
            "p5": p5,
            "p7": p7,
        },
        "hinge": hinge_status_188(),
        "F3": (
            "Do not flip residual/E1/L without Max+-free |μ|≤2/n or ≤1/(2p). "
            "Census-only evidence is not a general proof. "
            "|μ_actual|≰|μ_part| pointwise (p=7 counterexample)."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15188.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.188 T²φ/μ_part/p7 census; residual (i) OPEN")
    print(f"  T2κ proved={t2k['proved']}; T2φ p3,5={t2phi['certified']}")
    print(f"  μ_part master cert={part_master['certified']}")
    print(f"  2/n≤1/(2p) proved={two_n['proved']}")
    print(f"  p7 max|μ|={p7['max_abs_mu_kappa1']} < thr={p7['mu4_thr']}")
    print(f"  gsum={gsum_disj_lb_proved_general()}; e1={e1_closed_general()}")
    print(f"  open: {out['hinge']['open']}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
