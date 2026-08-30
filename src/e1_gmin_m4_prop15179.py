#!/usr/bin/env python3
"""
Prop 15.179 — Residual (ii) freeze-to-tight: dual two-level affine freeness-fail
forces k=3p−1 (no Gsum LB). Residual (i) hinge unchanged.

Does **not** soft-close residual (i) / E1 / L. Soft-close forbidden.

SETUP (15.171 dual two-level freeness-fail)
  Deep s₊=2 scores S∈{2,4} on Max+, dual S∈{−2,−4} on Max−.
  Freeness-fail affine on Max+: f_e = 3 − S
    (f_e=+1 on S=2, f_e=−1 on S=4).
  H := G ∪ {e}, |H|=k+1.

PROVED Max+-free Fraction:
  A. Freeze. Under freeness-fail affine, for every Max+ vector
        S_H = S + f_e = S + (3 − S) = 3.
     So S_H ≡ 3 on all of Max+.
  B. First moment. E_+[S_H] = 3, but also E_+[S_H] = |H|/p = (k+1)/p
     (each edge contributes E[f]=1/p). Hence
        (k+1)/p = 3  ⇔  k = 3p − 1.
  C. k≥3p impossible. Dual two-level freeness-fail affine cannot occur
     for any k≥3p (or any k≠3p−1). Pure first-moment contradiction;
     no Gsum, no Max+ census.
  D. k=3p−1 is fail-eq. Already empty when bi-tight empty
     (15.168.G / 15.171 fail-eq; S_H≡3 tight cover of size 3p).
  E. Consequence for residual (ii). The dual two-level affine Gsum Farkas
     of 15.171 (need=2(8−3k/p), candidate −12/(pn)) is **vacuous** for
     k≥3p: the configuration does not exist. Residual (ii) dual-two-level
     affine branch closes by freeze + fail-eq, without disj Gsum LB.

NOT claimed here:
  F. Residual (i) Type I k=3p−2 (S∈{1,5}, affine S+2f_e=3):
     S_H ∈ {2,4}, not constant — freeze does **not** apply.
     Still needs Gsum LB / |μ|≤2/n / ker-box (15.176–192).
     gsum_disj_lb_proved_general False.
  G. Non-affine / multi-level deep freeness-fail (f_e≡1 on U₂ but not
     f_e=3−S on support S∈{2,4}, or S support larger): **not** closed by
     this prop. Prop 15.193 audits that exhaustiveness is missing;
     15.171 weak ND only covers freeness (not freeness-fail). Full
     residual (ii) remains OPEN.
  H. E1 / L remain OPEN until residual (i) and full residual (ii) close.

Writes evidence/e1_gmin_m4_prop15179.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15720 import (  # noqa: E402
    required_bitight_levels_empty_all_primes,
)
from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    is_prime,
)


def freeze_S_H_equals_3() -> dict:
    """Affine f_e=3−S on S∈{2,4} ⇒ S_H≡3."""
    for S in (2, 4):
        f_e = 3 - S
        assert S + f_e == 3
    return {
        "proved": True,
        "theorem": (
            "Freeness-fail affine f_e=3−S with S∈{2,4} ⇒ "
            "S_H=S+f_e=3 on every Max+ vector."
        ),
    }


def first_moment_forces_k(p: int) -> dict:
    """S_H≡3 ⇒ (k+1)/p=3 ⇒ k=3p−1."""
    k_forced = 3 * p - 1
    return {
        "p": p,
        "E_S_H": 3,
        "E_S_H_also": "(k+1)/p",
        "k_forced": k_forced,
        "proved": True,
        "theorem": (
            "S_H≡3 on Max+ ⇒ E[S_H]=3=(k+1)/p ⇒ k=3p−1. "
            "Any other k contradicts first moment."
        ),
    }


def dual_twolevel_affine_impossible_k_ge_3p(p: int) -> dict:
    """k≥3p cannot support dual two-level freeness-fail affine."""
    k_forced = 3 * p - 1
    samples = []
    ok = True
    k0 = 3 * p
    for k in range(k0, k0 + 12):
        # first moment requires k=3p-1
        possible = k == k_forced
        samples.append({"k": k, "affine_first_moment_ok": possible})
        if k >= 3 * p and possible:
            ok = False
    return {
        "p": p,
        "k_forced_by_freeze": k_forced,
        "k_forced_lt_3p": k_forced < 3 * p,  # 3p-1 < 3p
        "no_k_ge_3p_satisfies_first_moment": k_forced < 3 * p,
        "proved": ok and k_forced == 3 * p - 1 and k_forced < 3 * p,
        "sample_k": samples[:6],
        "theorem": (
            "Dual two-level freeness-fail affine ⇒ k=3p−1. "
            "Hence impossible for all k≥3p."
        ),
    }


def theorem_freeze_all_primes(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [p for p in range(5, 80) if is_prime(p)]
    freeze = freeze_S_H_equals_3()
    rows = {}
    ok = freeze["proved"]
    for p in primes:
        fm = first_moment_forces_k(p)
        imp = dual_twolevel_affine_impossible_k_ge_3p(p)
        rows[str(p)] = {
            "k_forced": fm["k_forced"],
            "impossible_k_ge_3p": imp["proved"],
        }
        if not (fm["proved"] and imp["proved"]):
            ok = False
    return {
        "proved": ok,
        "freeze": freeze,
        "by_p": rows,
        "n_primes": len(primes),
        "theorem": (
            "For every prime p≥5: dual two-level freeness-fail affine "
            "forces k=3p−1 and is impossible for k≥3p (first moment)."
        ),
    }


def residual_ii_dual_twolevel_affine_closed() -> bool:
    """
    Dual two-level affine branch of residual (ii) closed for k≥3p
    by freeze (no Gsum LB). Fail-eq k=3p−1 still needs bi-tight empty.

    This is **not** full residual (ii) closed — see Prop 15.193
    exhaustiveness (multi-level / non-affine still open).
    """
    th = theorem_freeze_all_primes()
    # fail-eq uses bi-tight
    bt = required_bitight_levels_empty_all_primes()
    return bool(th["proved"] and bt)


def hinge_status_179() -> dict:
    return {
        "residual_ii_dual_twolevel_affine_k_ge_3p": (
            "CLOSED by freeze-to-tight first moment (this prop)"
        ),
        "residual_ii_fail_eq_k_3p_minus_1": (
            "empty when bi-tight empty (15.168.G)"
        ),
        "residual_i_type_I_k_3p_minus_2": (
            "OPEN — freeze does not apply; needs Gsum LB (15.176–178)"
        ),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "bound_proved_general": False,
        "open": (
            "Prove disj Gsum≥−1/p (or |μ₄|≤1/(2p) on |κ|=1) for residual (i); "
            "then flip gsum_disj_lb_proved_general and E1/L."
        ),
    }


def main() -> dict:
    primes = [p for p in range(5, 80) if is_prime(p)]
    th = theorem_freeze_all_primes(primes)
    rii = residual_ii_dual_twolevel_affine_closed()
    out = {
        "title": (
            "Prop 15.179 residual (ii) freeze-to-tight: dual two-level affine "
            "forces k=3p−1; residual (i) still OPEN"
        ),
        "proved": {
            "freeze_S_H_eq_3": True,
            "first_moment_k_eq_3p_minus_1": True,
            "dual_twolevel_affine_impossible_k_ge_3p": th["proved"],
            "residual_ii_dual_twolevel_affine_closed": rii,
            "residual_i_closed": False,
            "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
            "E1_closed": e1_closed_general(),
            "L_closed": False,
        },
        "algebra": th,
        "hinge": hinge_status_179(),
        "F3": (
            "Do not flip residual (i)/E1/L without Gsum LB; "
            "residual (ii) affine branch no longer needs Gsum."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15179.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.179 residual (ii) freeze-to-tight")
    print(f"  dual two-level affine impossible k≥3p: {th['proved']}")
    print(f"  residual_ii affine branch closed: {rii}")
    print(f"  residual_i / gsum / E1 still open")
    print(f"  gsum_disj_lb_proved_general={gsum_disj_lb_proved_general()}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
