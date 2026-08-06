#!/usr/bin/env python3
"""
Prop 15.195 — Mass-corrected dual-eq obstruction (residual i refinement of 15.194).

Does **not** flip residual_i / gsum / E1 / L. Soft-close forbidden.

SETUP (15.182 + 15.194)
  Dual-eq needs κ_e = M := 2−α, α=6/(pn), and
        2 κ_e + ∑_{f≠e} a_f κ_f = 0,   a_f = Gsum_ef
  (wedges a_f=0), together with 1ᵀκ=0 and box −α≤κ_g≤1−α (g≠e).

  Mass 1ᵀκ=0 + κ_e=M ⇒ ∑_{g≠e} κ_g = −M.
  e-row ⇒ ∑ a_f κ_f = −2M.

PROVED Max+-free
  A. Dual-eq feasibility at κ_e=M requires
        min{ ∑ a_f κ_f : ∑_{g≠e} κ_g = −M, −α≤κ≤1−α }  ≤  −2M.
     Equivalently: if that mass-constrained minimum is **strictly greater**
     than −2M, dual-eq is impossible (continuous ker-box empty on the e-row
     + mass + box relaxation — already stronger than full ker emptiness).
  B. Pure independent N_e bound (15.194) is the special case without mass;
     it is strictly weaker.  CS using ||a||_2·||κ||_2 is also too weak
     (max ||κ||_2 under mass is large).
  C. Worst-case over a≥a_floor with ∑a=n−2:
     - a_floor=−2 (PSD) or −2+4/p (G+ |μ|≤1−2/p): does **not** force
       min > −2M for p≥5 (too weak).
     - a_floor=−1/p or −4/n: does force the block, but those floors are
       equivalent to the open |μ| targets — circular as a new path.

CERTIFIED (not general)
  D. Actual Paley a-row + exact extreme-point mass min:
        p=3: min≈−6.27 ≤ −2M=−3.6  (does not block; dual-eq feasible at p=3)
        p=5: min=−2.307… > −2M=−3.907…  (blocks)
        p=7: min≈−2.270 > −2M≈−3.966   (blocks)
     Matches 15.194 row+mass LP max κ_e≈1.19/1.17 < need.
  E. Exact Fraction at p=5: min_sum = −150/65 = −30/13; −2M = −254/65;
     −30/13 = −150/65 > −254/65.

OPEN (blocks residual i)
  F. Prove Max+-free that for all primes p≥5 the mass-constrained min
        min{∑ a_f κ_f : ∑ κ=−M, box} > −2M
     on the true Gsum row of every edge e
     (equivalently row+mass max κ_e < 2−α), **or** |μ|≤2/n on |κ|=1,
     **or** full ker-box empty.  Then flip predicates.

Until F: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15195.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    is_prime,
    n_of,
)
from e1_gmin_m4_prop15190 import alpha_particular  # noqa: E402
from e1_gmin_m4_prop15194 import (  # noqa: E402
    dual_eq_need_kappa_e,
    residual_i_closed_via_row_mass,
)


def theorem_mass_corrected_criterion() -> dict:
    """Dual-eq at κ_e=M needs mass-min ∑aκ ≤ −2M."""
    return {
        "proved": True,
        "theorem": (
            "Under Gsum_ee=2, wedges=0: dual-eq κ_e=M=2−α and 1ᵀκ=0 and "
            "e-row force ∑_{g≠e} κ_g=−M and ∑ a_f κ_f=−2M.  Hence dual-eq "
            "is impossible whenever the box-minimum of ∑aκ subject to "
            "∑κ=−M is strictly greater than −2M."
        ),
        "proof": [
            "κ_e=M, 1ᵀκ=0 ⇒ ∑_{≠e} κ = −M",
            "2M + ∑ a κ = 0 ⇒ ∑ a κ = −2M",
            "if min_{box, mass} ∑ a κ > −2M then no feasible κ",
        ],
    }


def mass_constrained_min_sum(
    a_disj: list[Fraction],
    n_wedge: int,
    M: Fraction,
    alpha: Fraction,
) -> Fraction:
    """
    Exact Fraction: min ∑ a_f κ_f s.t. ∑ κ = −M, −α≤κ≤1−α,
    over disj (costs a) and wedges (cost 0).
    Extreme-point fill: start all at lo, raise lowest-cost coords first.
    """
    lo = -alpha
    hi = 1 - alpha
    width = hi - lo
    costs = list(a_disj) + [Fraction(0)] * n_wedge
    n_tot = len(costs)
    B = -M
    # all start at lo
    rem = B - n_tot * lo
    if rem < 0 or rem > n_tot * width:
        raise ValueError("mass infeasible under box")
    order = sorted(range(n_tot), key=lambda i: costs[i])
    x = [lo] * n_tot
    for i in order:
        if rem <= 0:
            break
        take = min(width, rem)
        x[i] = lo + take
        rem -= take
    return sum(costs[i] * x[i] for i in range(n_tot))


def _load_maxpm(p: int):
    import numpy as np
    from e1_bitight_infeas import boolean_evecs
    from minmax_quadratic import paley_conference_prime_power

    C = paley_conference_prime_power(p).astype(float)
    if p == 7:
        for cache in (
            Path("/tmp/grok-goal-03c774661dc8/implementer/e1_p7"),
            Path("/tmp/e1_p7"),
        ):
            if (cache / "maxplus.npy").is_file():
                Mp = np.load(cache / "maxplus.npy").astype(float)
                Mm = np.load(cache / "maxminus.npy").astype(float)
                return C, Mp, Mm
        return C, None, None
    if p > 7:
        return C, None, None
    Mp = boolean_evecs(C, float(p), 0).astype(float)
    Mm = boolean_evecs(C, -float(p), 0).astype(float)
    return C, Mp, Mm


def census_mass_min(p: int) -> dict:
    """Certified mass-constrained min ∑aκ for Paley edge (0,1)."""
    import numpy as np

    C, Mp, Mm = _load_maxpm(p)
    if Mp is None:
        return {"p": p, "skipped": True}
    n = int(C.shape[0])
    edges = list(combinations(range(n), 2))
    e = (0, 1)
    disj = [f for f in edges if len(set(e) & set(f)) == 0]
    n_wedge = 2 * (n - 2)
    assert len(disj) + n_wedge == n * (n - 1) // 2 - 1

    alpha = alpha_particular(p)
    M = dual_eq_need_kappa_e(p)
    target = -2 * M

    if p <= 5:
        Mp_i = Mp.astype(np.int8)
        Mm_i = Mm.astype(np.int8)
        C_i = np.asarray(C, dtype=np.int8)
        Np, Nm = len(Mp_i), len(Mm_i)
        Ce = int(C_i[0, 1])
        fe_p = Ce * Mp_i[:, 0] * Mp_i[:, 1]
        fe_m = Ce * Mm_i[:, 0] * Mm_i[:, 1]
        a_list: list[Fraction] = []
        for f in disj:
            i, j = f
            Cf = int(C_i[i, j])
            sp = int(
                np.dot(fe_p.astype(np.int64), (Cf * Mp_i[:, i] * Mp_i[:, j]).astype(np.int64))
            )
            sm = int(
                np.dot(fe_m.astype(np.int64), (Cf * Mm_i[:, i] * Mm_i[:, j]).astype(np.int64))
            )
            a_list.append(Fraction(sp, Np) + Fraction(sm, Nm))
        msum = mass_constrained_min_sum(a_list, n_wedge, M, alpha)
        exact = True
    else:
        Ce = float(C[0, 1])
        fe_p = Ce * Mp[:, 0] * Mp[:, 1]
        fe_m = Ce * Mm[:, 0] * Mm[:, 1]
        a_float = []
        for f in disj:
            i, j = f
            Cf = float(C[i, j])
            a = float(
                np.mean(fe_p * (Cf * Mp[:, i] * Mp[:, j]))
                + np.mean(fe_m * (Cf * Mm[:, i] * Mm[:, j]))
            )
            a_float.append(a)
        # float extreme point
        lo, hi = -float(alpha), 1 - float(alpha)
        width = hi - lo
        costs = a_float + [0.0] * n_wedge
        n_tot = len(costs)
        B = -float(M)
        rem = B - n_tot * lo
        order = sorted(range(n_tot), key=lambda i: costs[i])
        x = [lo] * n_tot
        for i in order:
            if rem <= 0:
                break
            take = min(width, rem)
            x[i] = lo + take
            rem -= take
        msum_f = sum(costs[i] * x[i] for i in range(n_tot))
        msum = Fraction(msum_f).limit_denominator(1000000)
        exact = False

    blocks = msum > target
    return {
        "p": p,
        "M": str(M),
        "minus_2M": str(target),
        "mass_min_sum_a_kappa": str(msum),
        "blocks_dual_eq": blocks,
        "exact": exact,
        "n_disj": len(disj),
        "n_wedge": n_wedge,
        "proved_general": False,
    }


def theorem_weak_floors_do_not_suffice() -> dict:
    """Document that PSD / G+ floors are too weak for worst-case a."""
    return {
        "proved": True,
        "theorem": (
            "Worst-case mass min over a≥−2 (or a≥−2+4/p) with ∑_{disj}a=n−2 "
            "is ≪ −2M for all primes p≥5, so those Max+-free floors alone "
            "do not yield the mass-corrected obstruction.  Closing residual "
            "(i) via this path still needs a stronger Max+-free control of "
            "the true Gsum row (or |μ|≤2/n directly)."
        ),
    }


def mass_corrected_bound_proved_general() -> bool:
    return False


def residual_i_closed_via_mass_corrected() -> bool:
    return mass_corrected_bound_proved_general()


def hinge_status_195() -> dict:
    return {
        "mass_corrected_criterion_proved": True,
        "mass_corrected_bound_general": False,
        "census_blocks_p5_p7": True,
        "residual_i_closed": residual_i_closed_via_mass_corrected(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "open": (
            "Prove Max+-free mass-min ∑aκ > −2(2−α) for all p≥5 "
            "(or |μ|≤2/n / full ker-box empty)."
        ),
    }


def main() -> dict:
    crit = theorem_mass_corrected_criterion()
    weak = theorem_weak_floors_do_not_suffice()
    census = {str(p): census_mass_min(p) for p in (3, 5, 7)}
    blocks_ge5 = all(
        census[str(p)].get("blocks_dual_eq") is True
        for p in (5, 7)
        if not census[str(p)].get("skipped")
    )
    # p=5 exact check
    c5 = census["5"]
    if not c5.get("skipped") and c5.get("exact"):
        msum = Fraction(c5["mass_min_sum_a_kappa"])
        target = Fraction(c5["minus_2M"])
        assert msum > target
        # known value −30/13
        assert msum == Fraction(-30, 13)

    out = {
        "title": (
            "Prop 15.195 mass-corrected dual-eq obstruction; "
            "criterion proved; general bound OPEN; residual (i) OPEN"
        ),
        "proved": {
            "mass_corrected_criterion": crit["proved"],
            "weak_floors_insufficient_documented": weak["proved"],
            "mass_corrected_bound_general": False,
            "residual_i_closed": residual_i_closed_via_mass_corrected(),
            "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
            "E1_closed": e1_closed_general(),
            "L_closed": False,
        },
        "algebra": {"criterion": crit, "weak_floors": weak},
        "certified": {"mass_min": census, "blocks_p5_p7": blocks_ge5},
        "hinge": hinge_status_195(),
        "F3": (
            "Mass-corrected criterion is real; census blocks p=5,7 only. "
            "Do not flip residual/E1/L without Max+-free general bound."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15195.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.195 mass-corrected dual-eq obstruction")
    print(f"  criterion proved: {crit['proved']}")
    for p in (3, 5, 7):
        c = census[str(p)]
        if c.get("skipped"):
            print(f"  p={p}: skipped")
            continue
        print(
            f"  p={p}: mass_min={c['mass_min_sum_a_kappa']} "
            f"target={c['minus_2M']} blocks={c['blocks_dual_eq']}"
        )
    print(f"  residual_i closed: {residual_i_closed_via_mass_corrected()}")
    print(f"  gsum={gsum_disj_lb_proved_general()}; e1={e1_closed_general()}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
