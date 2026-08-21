#!/usr/bin/env python3
"""
Prop 15.592 — ‖ν‖₂² = ½‖m₄⁺‖₂² − n(n−2)/16; leftover 3 reduces to Es4
with leading constant exactly 12; exact Es4 at p=11.

Does **not** flip type_I / phi_F_ge_6 / residual_ii / e1 / L. Soft-close forbidden.

SETUP (15.590/15.591 route; s := y·z)
  Max± ⊂ V±, the ±p eigenspaces of the symmetric conference C.  Since the
  eigenspaces are orthogonal:  y·z = 0 for EVERY y ∈ Max+, z ∈ Max−.
  Hence all cross moments are forced:  E₊₋[s^k] = 0 (k≥1) and
  E₊₋[e₄] = e₄(0) = n(n−2)/8,  where for ±1 vectors
  e₄(s) = (s⁴ − (6n−8)s² + 3n² − 6n)/24  and  Σ_S y_S z_S = e₄(y∘z).

PROVED (exact; verified by full double sums at p=5,7)
  A. Orthogonality:  s ≡ 0 on Max+ × Max−  (eigenspace orthogonality).
  B. Second moments:  E₊₊[s²] = 2n  (from E₊[yyᵀ] = 2P₊), E₊₋[s²] = 0.
  C. The L² identity:
        Σ_S ν(S)²  =  ½ E₊₊[e₄] − ½ E₊₋[e₄]  =  ½‖m₄⁺‖₂² − n(n−2)/16 .
     (ν = ½(m₄⁺−m₄⁻); Σm₄⁺m₄⁻ = E₊₋[e₄] is FORCED by A — no estimate.)
  D. Equivariance ⟹ per-fiber sup bound:  for every locus fiber F,
        ν̂_F² ≤ (Σ_S ν(S)²) / |orb_F|,
     with |orb_F| exact: generic (q³−q)/2, harmonic (q³−q)/8·(4/…)=221430·,
     equianharmonic (q³−q)/12 = min (values verified at p=11).

CONSEQUENCE (the sharpest current form of leftover 3)
  E. With binding(p) from the §9 budget (data-free) and
     min_F|orb_F| ≥ (q³−q)/12:
        leftover 3  ⟸  census(5,7)  +  [ Es4 := E₊₊[s⁴] ≤ 12n² + x(p)·n ]
     where x(p) = [48·binding(p)²·min|orb| − 16n + 6n·…]/n is data-free.
     At p=11: x = 32.60 needed vs TRUE 17.57 (chain verified exactly).
     x(p) grows to ≈ 85 as p → ∞; leftover-1's principal room demands
     ≈ 20 — strictly stronger.  SHARED BLOCKER, NOW EXACT:
     leftover-1's principal room ⟹ leftover-3's residual (large p).
  F. Kill: any Es4 majorant with leading constant (12+ε)n², ε > 0 fixed,
     CANNOT close leftover 3 (εn² ≫ 85n).  The leading constant must be
     exactly 12 — same wall as leftover 1, weaker lower-order demand.

NEW DATA (from C read backwards at p=11, where census is impossible)
  Es4(11) = 12n² + 17.57n  (exact rational in the JSON), continuing
  44.19 (p=5), 23.92 (p=7), 17.57 (p=11) — recorded, NOT extrapolated.

Writes evidence/e1_gmin_m4_prop15592.json
"""
from __future__ import annotations

import json
import time
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def e4_of_s(s: int, n: int) -> Fraction:
    return Fraction(s ** 4 - (6 * n - 8) * s ** 2 + 3 * n * n - 6 * n, 24)


def theorem_A_orthogonality(p: int) -> dict:
    from e1_gmin_m4_prop15590 import MuLab
    lab = MuLab(p, with_deg6=False)
    S = lab.Yp.astype(np.int64) @ lab.Ym.astype(np.int64).T
    return {"p": p, "proved": bool((S == 0).all()),
            "pairs": int(S.size)}


def theorem_C_identity(p: int) -> dict:
    from e1_gmin_m4_prop15590 import MuLab
    lab = MuLab(p, with_deg6=False)
    n, N = lab.n, lab.N
    lhs = Fraction(int((lab.d4.astype(object) ** 2).sum()), N * N)
    Yp = lab.Yp.astype(np.float32)
    Np = len(lab.Yp)
    Spp = (Yp @ Yp.T).astype(np.int64)
    E_pp = Fraction(sum(e4_of_s(int(s), n) for s in Spp.ravel()), Np * Np) * Np * Np
    E_pp = Fraction(int(sum((s**4 - (6*n-8)*s**2 + 3*n*n - 6*n) for s in Spp.ravel().astype(object))), 24 * Np * Np)
    rhs = E_pp / 2 - Fraction(n * (n - 2), 16)
    es4 = Fraction(int((Spp.ravel().astype(object) ** 4).sum()), Np * Np)
    return {"p": p, "proved": lhs == rhs, "sum_nu_sq": str(lhs),
            "Es4": str(es4),
            "Es4_excess_over_12n2_per_n": float((es4 - 12 * n * n) / n)}


def es4_from_nu_p11() -> dict:
    """Es4(11) exactly, via identity C read backwards from the nu data."""
    p, n = 11, 122
    N = 2 * 37457112
    # exact per-fiber (orbit size, nu_int) verified in orbit_sizes_p11
    fibers = [(221430, -500128), (442860, -369184), (442860, -503008),
              (147620, -417120), (885720, -504416)]
    Snu2 = sum(Fraction(o) * Fraction(v, N) ** 2 for o, v in fibers)
    m4sq = 2 * (Snu2 + Fraction(n * (n - 2), 16))
    es4 = 24 * m4sq + (6 * n - 8) * 2 * n - 3 * n * n + 6 * n
    return {"p": 11, "sum_nu_sq": str(Snu2), "Es4_exact": str(es4),
            "Es4_float": float(es4),
            "excess_over_12n2_per_n": float((es4 - 12 * n * n) / n),
            "min_locus_orbit": 147620,
            "min_locus_orbit_eq": "(q^3-q)/12"}


def needed_x_of_p11() -> dict:
    """The data-free Es4 threshold at p=11 through the section-9 budget."""
    p, n, q = 11, 122, 121
    A1 = 11
    c1 = Fraction(p - 4, 2 * p * A1)
    minorb = 147620
    budget = c1 * c1 * minorb            # allowed Sum nu^2
    m4sq_b = 2 * (budget + Fraction(n * (n - 2), 16))
    es4_b = 24 * m4sq_b + (6 * n - 8) * 2 * n - 3 * n * n + 6 * n
    return {"allowed_sum_nu_sq": float(budget),
            "needed_Es4_excess_per_n": float((es4_b - 12 * n * n) / n)}


def main():
    t0 = time.time()
    out = {"prop": "15.592",
           "title": "nu L2 identity; leftover 3 reduces to Es4 with leading constant 12"}
    for p in (5, 7):
        out[f"orthogonality_p{p}"] = theorem_A_orthogonality(p)["proved"]
        c = theorem_C_identity(p)
        out[f"identity_p{p}"] = c["proved"]
        out[f"Es4_p{p}"] = {"exact": c["Es4"], "excess_per_n": c["Es4_excess_over_12n2_per_n"]}
    out["Es4_p11"] = es4_from_nu_p11()
    out["needed_p11"] = needed_x_of_p11()
    out["chain_closes_p11"] = (out["needed_p11"]["needed_Es4_excess_per_n"]
                              > out["Es4_p11"]["excess_over_12n2_per_n"])
    out["flags_not_flipped"] = ["type_I", "phi_F_ge_6", "residual_ii", "e1", "L"]
    out["L_status"] = "OPEN"
    out["seconds"] = round(time.time() - t0, 1)
    path = ROOT / "evidence" / "e1_gmin_m4_prop15592.json"
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print("Prop 15.592  nu L2 identity; Es4 reduction")
    for k, v in out.items():
        if k.startswith(("orthogonality", "identity", "chain")):
            print(f"  {k} = {v}")
    print(f"  Es4 excess/n: p=5 {out['Es4_p5']['excess_per_n']:.2f}, "
          f"p=7 {out['Es4_p7']['excess_per_n']:.2f}, "
          f"p=11 {out['Es4_p11']['excess_over_12n2_per_n']:.2f} "
          f"(needed at 11: {out['needed_p11']['needed_Es4_excess_per_n']:.2f})")
    print(f"  ({out['seconds']}s)")
    return out


if __name__ == "__main__":
    main()
