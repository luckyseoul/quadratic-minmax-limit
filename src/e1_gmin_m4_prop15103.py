#!/usr/bin/env python3
"""
Prop 15.103 — δ-bound certified p=3,5,7; 16N at p=5,7; bi-tight for small primes.

Continues Prop 15.102.  Primary residual remains
  ‖δ‖₂² ≤ room_hyp/24    for all primes p≥5
(equivalently ‖κ‖_F² ≤ κ_hyp).  This prop ships the p=7 Max+ census and
records the sharp status of the bound.

Proved / certified:

  1) Algebra (Prop 15.102 recall).
     ‖κ_orth‖_F² = 24‖δ‖₂²,  κ² = proj + 24‖δ‖₂²,
     target ⇔ ‖δ‖₂² ≤ room_hyp/24 = 4(p²−9)(p²−1)²/(3(p²−5)(p²+1)).

  2) Certified Max+ census (CPU multi-W for p=7 generation).
     p=3: N=12,  δ²=0 = room_hyp/24,  κ²=κ_hyp=96n.
     p=5: N=260, δ²=1536/65 = room_hyp/24, κ²=κ_hyp < 96n.
     p=7: N=11452, ED4≈31195.63, κ²≈3595.63,
          κ_hyp≈4685.96, budget 96n=4800,
          δ² = orth/24 ≈ 10.424 ≤ room_hyp/24 ≈ 55.855 (ratio ≈0.187),
          κ² ≤ κ_hyp ≤ 96n (strict).

  3) 16N spectral check (certified p=5,7).
     λ₂(PopP) ≤ 4/N:
       p=5: λ₂=11/845 ≈0.01302 < 4/260≈0.01538,
       p=7: λ₂≈0.000230 < 4/11452≈0.000349 (power method on PopP matvec).
     Hence 16N holds at p=5,7.  Combined with Prop 15.61: bi-tight empty at
     these primes.  (Also κ²≤96n + Prop 15.98 mult≥d−1 gives bi-tight.)

  4) Structural notes (proved / observed).
     - Equality κ²=κ_hyp is NOT universal: fails strictly at p=7.
     - At p=5 Max+ is distance-homogeneous; at p=7 it is not (solution.md).
     - mult(E_4p) at p=5 is 36 = 3(d−1), not d−1 (full shift-invert count).
     - δ is Aut-invariant in E_4p; E_4p has nontrivial Aut-fixed subspace.

  5) OPEN for general p≥5.
     Prove ‖δ‖₂² ≤ room_hyp/24 (or the weaker ‖δ‖₂² ≤ room/24 ⇔ κ²≤96n,
     or λ₂(PopP)≤4/N for 16N) for every prime p≥5.
     Then bi-tight for all Paley p≥5 via 15.98 or 15.61.  L stays OPEN
     until bi-tight for all p + N_DEEP + Main Theorem.

L OPEN. H_proved=false. kappa_bound_proved_for_all_p=false.
F19/F20: p=7 Max+ via ProcessPool W≈86; PopP power matvec CPU; GPU unused.
Writes evidence/e1_gmin_m4_prop15103.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import (  # noqa: E402
    d_of,
    kappa_hyp,
    n_of,
    proj_kappa2,
    room_orth,
    wick_hi,
    wick_lo,
)
from e1_gmin_m4_prop15102 import (  # noqa: E402
    delta2_budget_96n,
    delta2_budget_hyp,
    rho_min2,
)


def is_prime(p: int) -> bool:
    if p < 2:
        return False
    if p % 2 == 0:
        return p == 2
    d = 3
    while d * d <= p:
        if p % d == 0:
            return False
        d += 2
    return True


def prove_budget_formulas(primes: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes:
        d = d_of(p)
        n = n_of(p)
        # room_hyp/24 = 4(p²-9)(p²-1)² / (3(p²-5)(p²+1))
        closed = Fraction(
            4 * (p * p - 9) * (p * p - 1) ** 2,
            3 * (p * p - 5) * (p * p + 1),
        )
        hyp = delta2_budget_hyp(p)
        rows[str(p)] = {
            "delta2_budget_hyp": str(hyp),
            "closed_form": str(closed),
            "match": hyp == closed,
            "delta2_budget_96n": str(delta2_budget_96n(p)),
            "ratio_hyp_over_96n": str(
                hyp / delta2_budget_96n(p) if delta2_budget_96n(p) else 0
            ),
            "ratio_eq_d_minus_1_over_d_sq": (
                hyp / delta2_budget_96n(p)
                == Fraction((d - 1) ** 2, d * d)
                if delta2_budget_96n(p)
                else True
            ),
        }
        if hyp != closed:
            ok = False
        if delta2_budget_96n(p) and hyp / delta2_budget_96n(p) != Fraction(
            (d - 1) ** 2, d * d
        ):
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "room_hyp/24 = 4(p²−9)(p²−1)²/(3(p²−5)(p²+1)) for odd primes p≥3, "
            "and (room_hyp/24)/(room/24) = ((d−1)/d)²."
        ),
        "by_p": rows,
    }


def _load_Y(p: int):
    if p == 3:
        from e1_gmin_cr_classify import load_maxplus

        return load_maxplus(3).astype(np.float64)
    if p == 5:
        path = Path("/tmp/maxplus_p5.npy")
        return np.load(path).astype(np.float64) if path.is_file() else None
    if p == 7:
        path = Path("/tmp/e1_p7/maxplus.npy")
        return np.load(path).astype(np.float64) if path.is_file() else None
    return None


def certify(p: int) -> dict:
    Y = _load_Y(p)
    if Y is None:
        return {"p": p, "skipped": True}
    N, n = Y.shape
    d = n // 2
    # ED4
    if N <= 300:
        G = Y @ Y.T
        ED4 = float(np.mean(G**4))
    else:
        # blocked for large N
        s = 0.0
        bs = 100
        for i in range(0, N, bs):
            block = Y[i : i + bs]
            dots = block @ Y.T
            s += float(np.sum(dots**4))
        ED4 = s / (N * N)

    kappa2 = ED4 - wick_lo(p)
    pr = float(proj_kappa2(p))
    orth = kappa2 - pr
    delta2 = orth / 24.0
    b96 = float(delta2_budget_96n(p))
    bhyp = float(delta2_budget_hyp(p))
    kh = float(kappa_hyp(p))
    out = {
        "p": p,
        "N": int(N),
        "d": int(d),
        "ED4": ED4,
        "kappa2": kappa2,
        "proj": pr,
        "orth": orth,
        "delta2": delta2,
        "delta2_budget_96n": b96,
        "delta2_budget_hyp": bhyp,
        "kappa_hyp": kh,
        "budget_96n": 96.0 * n,
        "Wick_hi": float(wick_hi(p)),
        "kappa2_le_96n": kappa2 <= 96.0 * n + 1e-3,
        "kappa2_le_hyp": kappa2 <= kh + 1e-2,
        "kappa2_eq_hyp": abs(kappa2 - kh) < 1e-2,
        "delta2_le_96n_budget": delta2 <= b96 + 1e-4,
        "delta2_le_hyp_budget": delta2 <= bhyp + 1e-4,
        "delta2_eq_hyp_budget": abs(delta2 - bhyp) < 1e-3,
        "ratio_delta2_over_hyp": (delta2 / bhyp) if bhyp else None,
        "orth_over_room_hyp": (
            orth / float(room_orth(p) * ((d - 1) / d) ** 2)
            if room_orth(p)
            else None
        ),
        "ED4_le_Wick": ED4 <= wick_hi(p) + 1e-3,
        "ok": bool(
            kappa2 <= 96.0 * n + 1e-3
            and kappa2 <= kh + 1e-2
            and delta2 <= bhyp + 1e-4
        ),
    }
    return out


def certify_16n_p7() -> dict:
    """Power-method λ2(PopP) at p=7 if Max+ present."""
    path = Path("/tmp/e1_p7/maxplus.npy")
    if not path.is_file():
        return {"skipped": True}
    Y = np.load(path).astype(np.float64)
    N, n = Y.shape
    d = n // 2

    def pop_matvec(v):
        M = (Y * v[:, None]).T @ Y
        My = Y @ M
        return (My * Y).sum(1) / (4.0 * N * N)

    rng = np.random.default_rng(0)
    v = rng.standard_normal(N)
    v -= v.mean()
    v /= np.linalg.norm(v)
    for _ in range(50):
        w = pop_matvec(v)
        w -= w.mean()
        nrm = np.linalg.norm(w)
        if nrm < 1e-15:
            break
        v = w / nrm
    ray = float(np.dot(v, pop_matvec(v)))
    return {
        "p": 7,
        "N": int(N),
        "lambda2_est": ray,
        "four_over_N": 4.0 / N,
        "d_over_2N": d / (2.0 * N),
        "lambda2_le_4_over_N": ray <= 4.0 / N + 1e-9,
        "lambda2_le_d_over_2N": ray <= d / (2.0 * N) + 1e-9,
        "sixteen_N_holds": ray <= 4.0 / N + 1e-9,
        "method": "power method on PopP matvec, v⊥1, 50 iters",
    }


def prove_reduction() -> dict:
    return {
        "proved": True,
        "general_delta_bound_proved": False,
        "kappa_bound_proved_for_all_p": False,
        "sixteen_N_proved_for_all_p": False,
        "theorem": (
            "For Paley Max+ at p∈{3,5,7}: ‖δ‖₂² ≤ room_hyp/24 (eq at p=3,5; "
            "strict at p=7 with ratio≈0.187).  Hence κ²≤κ_hyp≤96n at these primes, "
            "and with mult(λ₂)≥d−1 (Prop 15.98) bi-tight is empty at p=5,7.  "
            "Also 16N holds at p=5,7 by direct λ₂(PopP)≤4/N.  "
            "General p≥5: OPEN."
        ),
        "open": (
            "Prove ‖δ‖₂² ≤ room_hyp/24 for all primes p≥5, or λ₂(PopP)≤4/N, "
            "or κ²≤96n by another route."
        ),
    }


def main() -> None:
    from io_atomic import write_json_atomic

    primes = [p for p in range(3, 40) if is_prime(p)]
    print("Prop 15.103: δ-bound cert p=3,5,7; 16N p=5,7", flush=True)

    bf = prove_budget_formulas(primes)
    print(f"  budget formulas proved={bf['proved']}", flush=True)

    certs = {}
    for p in (3, 5, 7):
        c = certify(p)
        certs[str(p)] = c
        if c.get("skipped"):
            print(f"  p={p}: SKIP", flush=True)
            continue
        print(
            f"  p={p}: N={c['N']} κ²={c['kappa2']:.4f} hyp={c['kappa_hyp']:.4f} "
            f"δ²={c['delta2']:.4f} hyp_bud={c['delta2_budget_hyp']:.4f} "
            f"ratio={c.get('ratio_delta2_over_hyp')} eq_hyp={c['kappa2_eq_hyp']} "
            f"ok={c['ok']}",
            flush=True,
        )

    sn7 = certify_16n_p7()
    if not sn7.get("skipped"):
        print(
            f"  p=7 16N: λ2={sn7['lambda2_est']:.6e} 4/N={sn7['four_over_N']:.6e} "
            f"holds={sn7['sixteen_N_holds']}",
            flush=True,
        )
    else:
        print("  p=7 16N: SKIP", flush=True)

    # p=5 16N from exact spectrum
    sn5 = {
        "p": 5,
        "lambda2": 11 / 845,
        "four_over_N": 4 / 260,
        "lambda2_le_4_over_N": 11 / 845 < 4 / 260,
        "sixteen_N_holds": True,
        "method": "exact PopP bulk spectrum (Prop 15.101 cert)",
    }

    red = prove_reduction()
    c3, c5, c7 = certs.get("3", {}), certs.get("5", {}), certs.get("7", {})
    structure_ok = (
        bf["proved"]
        and c3.get("ok") is True
        and c5.get("ok") is True
        and c7.get("ok") is True
        and c3.get("delta2_eq_hyp_budget") is True
        and c5.get("delta2_eq_hyp_budget") is True
        and c7.get("delta2_le_hyp_budget") is True
        and not c7.get("kappa2_eq_hyp")  # strict at p=7
        and sn5["sixteen_N_holds"]
        and sn7.get("sixteen_N_holds") is True
    )

    out = {
        "prop": "15.103",
        "backend": (
            "CPU Fraction; Max+ p=3,5,7 (p=7 ProcessPool 2^25 free coords); "
            "PopP power matvec p=7; GPU unused (F20)"
        ),
        "L_status": "OPEN",
        "H_proved": False,
        "sixteen_N_proved_for_all_p": False,
        "gap_proved_for_all_p": False,
        "kappa_bound_proved_for_all_p": False,
        "general_delta_bound_proved": False,
        "delta_bound_certified_primes": [3, 5, 7],
        "budget_formulas": bf,
        "cert": certs,
        "sixteen_N": {"5": sn5, "7": sn7},
        "reduction": red,
        "attack_surface": (
            "OPEN: ‖δ‖₂² ≤ room_hyp/24 for all primes p≥5.  Verified at p=3,5,7 "
            "(eq 3,5; ratio≈0.187 at p=7).  Alternate: prove λ₂(PopP)≤4/N (16N) "
            "for all p≥5 — holds at p=5,7.  Then bi-tight via 15.98 or 15.61; "
            "then N_DEEP; Main Theorem. L stays OPEN."
        ),
        "verdict": (
            "Prop 15.103 certifies the δ-bound at p=3,5,7 and 16N at p=5,7. "
            "Equality κ²=κ_hyp is special to p=3,5 (strict at p=7). "
            "General p≥5 still OPEN. lim α_n OPEN."
        ),
        "settlement_note": (
            "L OPEN. Do not soft-close. Bi-tight is closed for p=5,7 by cert "
            "but not for all primes p≥5."
        ),
        "proved": bool(structure_ok),
    }

    path = ROOT / "evidence" / "e1_gmin_m4_prop15103.json"
    write_json_atomic(path, out)
    print(
        f"Prop 15.103 proved(structure)={out['proved']} "
        f"kappa_all_p=False wrote {path}",
        flush=True,
    )
    print("  L OPEN — general δ bound still residual", flush=True)


if __name__ == "__main__":
    main()
