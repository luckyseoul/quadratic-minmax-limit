#!/usr/bin/env python3
"""
Prop 15.106 — Rest-average-8; mult≥d ⇒ 16-threshold; kurtosis form of residual.

Closes structural gaps around the Path C 16N / κ residual (builds on 15.105).

============================================================================
Setup. Paley conference Max+, n=p²+1, d=n/2, m=dim Z=d(d−3)/2,
Φ|Z eigenvalues λ_1≥⋯≥λ_m≥0, T=∑λ=n(n−2), μ̄=T/m=8(n−2)/(n−6).

From Prop 15.105: ∑(λ_α−μ̄)² = ‖κ_orth‖_F²  (variance = Fickus residual).

============================================================================
Theorem A (rest-average-8 threshold) — PROVED Fraction algebra.
  Assume mult(λ_max) ≥ d.  Write L=λ_max and
      R_avg := (T − d·L) / (m − d)   (mean of the complementary spectrum,
      using the lower bound mult≥d for an upper bound on L).
  Then the following are equivalent:
    (i)  L ≤ 16
    (ii) R_avg ≥ 8
    (iii) d(L−μ̄)² + (m−d)(8−μ̄)² ≤ ∑(λ−μ̄)²   is automatic at L=16, rest≡8;
         and at the variance budget orth=room one has L≤16 (Prop 15.105.E).
  Proof of (i)⇔(ii).
    m−d = d(d−5)/2 for d>5 (and 0 at d=5 where m=d and L=μ̄=16).
    T − 16d = n(n−2)−16d = 2d(2d−2)−16d = 4d(d−1)−16d = 4d(d−5).
    (T − 16d)/(m−d) = 4d(d−5) / (d(d−5)/2) = 8.
    Hence L≤16 ⇔ T−d L ≥ T−16d ⇔ (T−d L)/(m−d) ≥ 8.  ∎

  Interpretation.  The Wick baseline for the Rayleigh of unit B∈Z is 8 in the
  Q/N = 6+2 ray normalisation (ray_Wick=1 gives 8).  Extremal 16N with mult=d
  forces the complementary spectrum to average exactly the Wick value 8.

Theorem B (kurtosis form of κ≤96n) — PROVED as equivalence.
  Let D=y·z for y,z independent uniform on Max+.  E[D]=0, E[D²]=2n (2-design).
  Kurtosis κ₄(D) := E[D⁴] / (E[D²])² = ED4 / (4 n²).
  Then
    ‖κ‖_F² ≤ 96n   ⇔   ED4 ≤ 12n²+48n = wick_hi
                    ⇔   κ₄(D) ≤ 3 + 12/n.
  Proof.  ED4 = wick_lo + ‖κ‖_F² = 12n²−48n + ‖κ‖_F², so
  ‖κ‖_F²≤96n ⇔ ED4≤12n²+48n.  Divide by 4n²:
  ED4/(4n²) ≤ 3 + 12/n.  ∎
  (Gaussian D would have κ₄=3; the design allows excess 12/n.)

Theorem C (mult ≥ d for Paley — structure + certification).
  By Prop 15.98, mult(λ₂(P⊙P)) ≥ d−1 (PSL min nontrivial irrep dim (q−1)/2).
  Since λ₂(P⊙P)=λ_max(Φ|Z)/(4N), same multiplicity.
  Character table of PSL(2,q), q=p² odd: irrep dimensions include both
  (q−1)/2 = d−1 and (q+1)/2 = d.
  Certified: mult(λ_max(Φ))=d exactly at p=3,5,7 (top isotypic = one copy of
  the d-dimensional irrep).  The general identification mult=d for all primes
  p≥5 is the standard consequence of the top zonal spherical function of the
  Aut scheme living in the (q+1)/2-isotypic (Weil / discrete-series type);
  we package it as:
    Conjecture C* (mult=d).  mult(λ_max(Φ|Z))=d for every prime p≥5.
  Proven substitute used below: mult ≥ d−1 (Prop 15.98, unconditional).

Theorem D (16N criterion with mult≥d−1 and rest≥8) — PROVED algebra.
  If mult(λ_max)≥d−1 and every eigenvalue other than the top is ≥8, then
    L ≤ (T − 8(m − (d−1))) / (d−1).
  This upper bound is NOT ≤16 for all p (e.g. p=5: (624−8·53)/12=16.666>16).
  So rest≥8 with only mult≥d−1 is insufficient; one needs mult≥d (Theorem A)
  or a stronger floor than 8 on a large part of the spectrum.

============================================================================
OPEN residual (honest, single target for Path C bi-tight):
  Prove for all primes p≥5 one of:
    (1) κ₄(D) ≤ 3+12/n   (⇔ ‖κ‖_F²≤96n ⇔ δ²≤room/24)
    (2) λ_max(Φ|Z) ≤ 16 directly (e.g. Norton scheme eigenvalues, or Gu_disj≤5N)
    (3) Conjecture C* (mult=d) AND rest-average ≥8 (then Theorem A ⇒ 16N).

  With (1) + Prop 15.98: bi-tight empty for all p≥5.
  With (2) + Prop 15.61: bi-tight empty for all p≥5.
  Certified (1)(2) at p=3,5,7.  L OPEN.

F19/F20: Fraction algebra + Max+ cert; GPU unused.
Writes evidence/e1_gmin_m4_prop15106.json
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
    ED4_flat,
    d_of,
    kappa_hyp,
    n_of,
    proj_kappa2,
    room_orth,
    wick_hi,
    wick_lo,
)
from e1_gmin_m4_prop15105 import (  # noqa: E402
    avg_Phi_Z,
    certify_variance_on_maxplus,
    dim_Z,
    tr_Phi_Z,
)


def is_prime(p: int) -> bool:
    if p < 2:
        return False
    if p % 2 == 0:
        return p == 2
    k = 3
    while k * k <= p:
        if p % k == 0:
            return False
        k += 2
    return True


def prove_rest_average_8(primes: list[int]) -> dict:
    """Theorem A: mult≥d ⇒ (L≤16 ⇔ rest avg≥8); at L=16 rest avg=8 exactly."""
    rows = {}
    ok = True
    for p in primes:
        n = n_of(p)
        d = d_of(p)
        m = dim_Z(p)
        T = tr_Phi_Z(p)
        if d <= 5:
            # p=3: m=d=5, spectrum flat at 16
            rows[str(p)] = {
                "d": d,
                "m": m,
                "flat_spectrum": True,
                "L_eq_16": True,
                "rest_avg_at_L16": None,
            }
            continue
        rest_dim = m - d
        T_minus_16d = T - 16 * d
        rest_avg_at_16 = Fraction(T_minus_16d, rest_dim)
        # algebra: =8
        closed = Fraction(4 * d * (d - 5), rest_dim)  # should be 8
        rows[str(p)] = {
            "d": d,
            "m": m,
            "rest_dim": rest_dim,
            "T_minus_16d": T_minus_16d,
            "rest_avg_at_L16": str(rest_avg_at_16),
            "equals_8": rest_avg_at_16 == 8,
            "closed_4d_dm5_over_rest": str(closed),
            "mu": str(avg_Phi_Z(p)),
        }
        if rest_avg_at_16 != 8:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "If mult(λ_max)≥d then L≤16 ⇔ mean(rest)≥8. "
            "At L=16 the rest mean equals 8 exactly "
            "(T−16d)/(m−d)=8 for conference d=n/2>5)."
        ),
        "by_p": rows,
    }


def prove_kurtosis_equivalence(primes: list[int]) -> dict:
    """Theorem B: κ≤96n ⇔ kurtosis(D)≤3+12/n."""
    rows = {}
    ok = True
    for p in primes:
        n = n_of(p)
        wh = wick_hi(p)
        thr_kurt = 3 + Fraction(12, n)
        # wick_hi / (4n²) = (12n²+48n)/(4n²) = 3+12/n
        check = Fraction(wh, 4 * n * n)
        rows[str(p)] = {
            "wick_hi": wh,
            "kurtosis_threshold": str(thr_kurt),
            "wick_hi_over_4n2": str(check),
            "match": check == thr_kurt,
            "ED2": 2 * n,
            "equivalence": "‖κ‖²≤96n ⇔ ED4≤wick_hi ⇔ κ₄(D)≤3+12/n",
        }
        if check != thr_kurt:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "With E[D]=0, E[D²]=2n on Max×Max: "
            "‖κ‖_F²≤96n ⇔ ED4≤12n²+48n ⇔ kurtosis κ₄(D)≤3+12/n."
        ),
        "by_p": rows,
    }


def prove_mult_d_status() -> dict:
    return {
        "mult_ge_d_minus_1_proved": True,
        "mult_ge_d_minus_1_ref": "Prop 15.98 PSL min irrep (q−1)/2",
        "mult_eq_d_certified_primes": [3, 5, 7],
        "mult_eq_d_proved_for_all_p": False,
        "conjecture_C_star": (
            "mult(λ_max(Φ|Z))=d for every prime p≥5 "
            "(top isotypic = Weil/discrete-series irrep of dim (q+1)/2)."
        ),
        "note": (
            "Prop 15.105.E needs mult≥d for the exact L≤16 criterion with "
            "orth≤room. Prop 15.98 bi-tight only needs mult≥d−1 + κ≤96n."
        ),
    }


def certify_rest_and_kurtosis(p: int) -> dict:
    """On Max+: check kurtosis bound, mult, λ_max≤16, rest average."""
    base = certify_variance_on_maxplus(p)
    if base.get("skipped"):
        return base
    n = n_of(p)
    d = d_of(p)
    m = dim_Z(p)
    T = tr_Phi_Z(p)
    ED4 = base["ED4"]
    kurt = ED4 / (4.0 * n * n)
    thr = 3 + 12.0 / n
    L = base["lambda_max_Phi"]
    mult = base["mult_top_observed"]
    # rest average using observed mult if full spectrum, else using mult=d
    m1 = mult if mult >= d - 1 else d
    if m > m1 and L is not None:
        rest_avg = (T - m1 * L) / (m - m1)
    else:
        rest_avg = None
    return {
        **base,
        "kurtosis": kurt,
        "kurtosis_threshold": thr,
        "kurtosis_le_thr": kurt <= thr + 1e-9,
        "kappa2_le_96n": base["kappa2"] <= 96 * n + 1e-6,
        "rest_avg_using_obs_mult": rest_avg,
        "rest_avg_ge_8": (rest_avg is not None and rest_avg >= 8 - 1e-6),
        "theorem_A_check": (
            (L <= 16 + 1e-8) == (rest_avg is None or rest_avg >= 8 - 1e-6)
            if rest_avg is not None
            else None
        ),
    }


def prove_general_status() -> dict:
    return {
        "rest_average_8_proved": True,
        "kurtosis_equivalence_proved": True,
        "mult_ge_d_minus_1_proved": True,
        "mult_eq_d_proved_for_all_p": False,
        "sixteen_N_proved_for_all_p": False,
        "kappa_bound_proved_for_all_p": False,
        "H_proved": False,
        "open_residual": (
            "Prove κ₄(D)≤3+12/n (⇔‖κ‖²≤96n) OR λ_max(Φ)≤16 OR "
            "(mult=d and rest-average≥8) for all primes p≥5."
        ),
    }


def main() -> dict:
    primes = [p for p in range(3, 40) if is_prime(p)]
    print("Prop 15.106: rest-average-8 + kurtosis residual form", flush=True)

    ra = prove_rest_average_8(primes)
    print(f"  rest-average-8 proved={ra['proved']}", flush=True)

    ku = prove_kurtosis_equivalence(primes)
    print(f"  kurtosis equivalence proved={ku['proved']}", flush=True)

    md = prove_mult_d_status()
    print(f"  mult≥d−1 proved; mult=d all p={md['mult_eq_d_proved_for_all_p']}", flush=True)

    cert = {}
    for p in (3, 5, 7):
        c = certify_rest_and_kurtosis(p)
        cert[str(p)] = c
        if c.get("skipped"):
            print(f"  cert p={p}: skipped", flush=True)
        else:
            print(
                f"  cert p={p}: κ₄={c['kurtosis']:.6f}≤{c['kurtosis_threshold']:.6f}? "
                f"{c['kurtosis_le_thr']} L={c['lambda_max_Phi']:.4f} "
                f"mult={c['mult_top_observed']} rest_avg={c['rest_avg_using_obs_mult']}",
                flush=True,
            )

    gen = prove_general_status()
    out = {
        "prop": "15.106",
        "backend": "CPU Fraction + Max+ cert p=3,5,7; GPU unused (F20)",
        "L_status": "OPEN",
        "H_proved": False,
        "sixteen_N_proved_for_all_p": False,
        "kappa_bound_proved_for_all_p": False,
        "rest_average_8": ra,
        "kurtosis_equivalence": ku,
        "mult_status": md,
        "cert": cert,
        "general_status": gen,
        "attack_surface": gen["open_residual"],
        "verdict": (
            "Proved: rest-average-8 threshold (mult≥d ⇒ L≤16 ⇔ rest≥8; "
            "at L=16 rest mean=8=Wick). Kurtosis form κ₄≤3+12/n ⇔ κ≤96n. "
            "Certified kurtosis+16N at p=3,5,7. "
            "OPEN: κ₄≤3+12/n or λ_max≤16 for general p≥5. L OPEN."
        ),
        "settlement_note": "No soft-close. Path C residual remains for general p≥5.",
        "proved": {
            "rest_average_8": ra["proved"],
            "kurtosis_equivalence": ku["proved"],
            "sixteen_N_for_all_p": False,
            "kappa_bound_for_all_p": False,
        },
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15106.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str))
    print(f"Wrote {path}", flush=True)
    return out


if __name__ == "__main__":
    main()
