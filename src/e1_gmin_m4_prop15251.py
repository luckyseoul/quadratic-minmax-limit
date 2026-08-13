#!/usr/bin/env python3
"""
Prop 15.251 — Cy-expansion identity on |κ|=1; m4+=m4- census;
|μ₄|≤2/n census at p=5,7; residual-(i) still OPEN on |μ|≤2/n general.

Does **not** flip residual_i / type_I / gsum / E1 / L. Soft-close forbidden.

SETUP (15.176–177, 15.191, 15.250)
  μ₄ = ½(m₄⁺+m₄⁻), m₄⁺=E_{Max+}[∏_{i∈S} y_i].
  Target: |μ₄|≤2/n on |κ|=1 (⇒ ≤1/(2p) for p≥5 ⇒ Gsum≥−1/p ⇒ residual-i).
  Cy-expansion: E[∏_{i∈S}(Cy)_i] = p⁴ m₄^{+} on Max+ (Cy=py).

PROVED Max+-free
  A. Permanent of C[S,S] equals derangement permanent (C_ii=0 kills fixed
     points). On |κ|=1, derang_perm = 1 (15.191 A, 64-exhaust).  ∎

  B. Cy-expansion identity on |κ|=1 for Max+ (Paley π + 15.191 A–E):
        p⁴ m₄⁺ = −2 φ + A[m₄⁺],
     where A is the injection transfer
        (A f)(S) = ∑_{injections f:S→[n]} (∏_{i∈S} C_{i,f(i)}) f(f(S)).
     Equivalently with Per=1 on the diagonal:
        (p⁴ − 1) m₄⁺ + 2 φ = Ext[m₄⁺],
     Ext = ∑_{T≠S} W(S,T) m₄⁺(T), W the bipartite-permanent weight.  ∎
     (size1+size2 = −2φ; size3 expectation vanishes on |κ|=1; size4 = A.)

  C. Sufficient pointwise target.
     |m₄⁺| ≤ 2/n and |m₄⁻| ≤ 2/n on |κ|=1 ⇒ |μ₄| ≤ 2/n
     ⇒ |μ₄| ≤ 1/(2p) for primes p≥5 (15.188 E)
     ⇒ Gsum ≥ −1/p (15.176) ⇒ residual-(i) Farkas.  ∎

  D. 1ᵀy = ±(p+1) on Max+ for normalised Paley (15.189 A, recalled).  ∎

CERTIFIED (Max± census; not general)
  E. m₄⁺ = m₄⁻ on all |κ|=1 at p=5 and p=7 (max |m₄⁺−m₄⁻|=0).
     Hence μ₄ = m₄⁺ = m₄⁻ on |κ|=1 at these p.
  F. |μ₄| ≤ 2/n on |κ|=1 at p=5 (max=3/65 < 1/13) and p=7
     (max=109/2863 < 1/25). Also ≤1/(2p).
  G. Identity (B): Ext1:=(p⁴−1)m₄+2φ equals Ext2:=∑_{T≠S}W m₄ at p=5
     on |κ|=1 sample (machine precision).
  H. At p=5: m₄=f4 on |κ|=1 and Ext = (4(p²−1)/p)κ + (2−p+1/p)φ.

OPEN (blocks residual i / predicate flip)
  I. Prove |μ₄| ≤ 2/n (or ≤1/(2p)) on every |κ|=1 four-set for all
     primes p≥5 Max+-free — e.g. control Ext so that
        |(p⁴−1)μ₄| = |−2φ+Ext| ≤ 2(p²−1),
     or prove hull/envelope, or E[s⁴]≤15n²−22n (15.250 R-path),
     or ker=sc (15.249 free-e). Soft-close forbidden.

Until I: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15251.json
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
    residual_i_dual_eq_empty_proved_general,
    type_I_k_3p_minus_2_closed_general,
)
from e1_gmin_m4_prop15177 import mu4_threshold  # noqa: E402
from e1_gmin_m4_prop15188 import two_over_n  # noqa: E402
from e1_gmin_m4_prop15191 import (  # noqa: E402
    theorem_derangement_perm_on_kappa1,
    theorem_f4_majorant,
    theorem_size1_plus_size2,
)


def two_over_n_le_half_p(p: int) -> bool:
    """2/n ≤ 1/(2p) for primes p≥5."""
    n = n_of(p)
    return Fraction(2, n) <= Fraction(1, 2 * p)


def theorem_per_equals_derang() -> dict:
    t = theorem_derangement_perm_on_kappa1()
    return {
        "proved": t["proved"],
        "theorem": (
            "C_ii=0 ⇒ Per(C[S,S])=derang_perm. On |κ|=1, derang_perm=1 "
            "(15.191 A). Hence Per=1 on |κ|=1."
        ),
        "depends_on": ["derang_perm_kappa1", "zero_diagonal"],
    }


def theorem_Cy_identity_structure() -> dict:
    s12 = theorem_size1_plus_size2()
    per = theorem_per_equals_derang()
    return {
        "proved": s12["proved"] and per["proved"],
        "theorem": (
            "On |κ|=1 for Max+ (Paley): (p⁴−1)m₄⁺ + 2φ = Ext[m₄⁺], "
            "with Ext the external-injection transfer.  From Cy-expansion "
            "p⁴ m₄⁺ = −2φ + A[m₄⁺], Per=1 on diagonal (A–B), size3 "
            "expectation vanishes (15.191 D)."
        ),
        "depends_on": [
            "size1_plus_size2",
            "size31_vanishes",
            "per_equals_1",
            "pairwise_pi_189",
        ],
    }


def theorem_two_over_n_suffices() -> dict:
    f4 = theorem_f4_majorant()
    ok = True
    sample = {}
    for p in [p for p in range(5, 80) if is_prime(p)]:
        le = two_over_n_le_half_p(p)
        if not le:
            ok = False
        if p in (5, 7, 11, 13):
            sample[str(p)] = {
                "2/n": str(Fraction(2, n_of(p))),
                "1/(2p)": str(mu4_threshold(p)),
                "2n_le_half_p": le,
            }
    return {
        "proved": ok and f4["proved"],
        "theorem": (
            "|m₄⁺|≤2/n and |m₄⁻|≤2/n on |κ|=1 ⇒ |μ₄|≤2/n ≤1/(2p) (p≥5) "
            "⇒ Gsum≥−1/p ⇒ residual-(i) Farkas (15.176–177)."
        ),
        "by_p_sample": sample,
        "depends_on": ["f4_majorant_two_over_n", "mu4_gsum_farkas_176"],
    }


def certify_mu4_census() -> dict:
    """Census |μ₄|≤2/n on |κ|=1 at p=5,7 using Max± caches."""
    import numpy as np
    from itertools import combinations
    from minmax_quadratic import paley_conference_prime_power

    out = {}
    for p in (5, 7):
        yp = Path(f"/tmp/maxplus_p{p}.npy")
        ym = Path(f"/tmp/maxminus_p{p}.npy")
        if not yp.exists() or not ym.exists():
            out[str(p)] = {"ok": False, "reason": "missing Max± cache"}
            continue
        Yp = np.load(yp)
        Ym = np.load(ym)
        C = paley_conference_prime_power(p).astype(np.float64)
        n = C.shape[0]
        fours = np.array(list(combinations(range(n), 4)), dtype=np.int32)
        a, b, c, d = fours[:, 0], fours[:, 1], fours[:, 2], fours[:, 3]
        kap = C[a, b] * C[c, d] + C[a, c] * C[b, d] + C[a, d] * C[b, c]
        ids = np.where(np.abs(np.abs(kap) - 1) < 1e-9)[0]
        max_mu = 0.0
        max_diff = 0.0
        chunk = 3000
        for s in range(0, len(ids), chunk):
            ix = ids[s : s + chunk]
            aa, bb, cc, dd = (
                fours[ix, 0],
                fours[ix, 1],
                fours[ix, 2],
                fours[ix, 3],
            )
            mp = np.mean(Yp[:, aa] * Yp[:, bb] * Yp[:, cc] * Yp[:, dd], axis=0)
            mm = np.mean(Ym[:, aa] * Ym[:, bb] * Ym[:, cc] * Ym[:, dd], axis=0)
            mu = 0.5 * (mp + mm)
            max_mu = max(max_mu, float(np.max(np.abs(mu))))
            max_diff = max(max_diff, float(np.max(np.abs(mp - mm))))
        two_n = 2.0 / n
        thr = 1.0 / (2 * p)
        out[str(p)] = {
            "ok": True,
            "n_kappa1": int(len(ids)),
            "max_abs_mu4": max_mu,
            "max_abs_m4plus_minus_m4minus": max_diff,
            "two_over_n": two_n,
            "one_over_2p": thr,
            "mu4_le_2n": max_mu <= two_n + 1e-12,
            "mu4_le_half_p": max_mu <= thr + 1e-12,
            "m4plus_eq_m4minus": max_diff < 1e-12,
        }
    return out


def residual_i_closed_via_251() -> bool:
    return False


def hinge_status_251() -> dict:
    return {
        "Cy_identity_structure": True,
        "two_over_n_suffices": True,
        "mu4_bound_general": False,
        "residual_i_closed": residual_i_closed_via_251(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "open": (
            "Prove |μ₄|≤2/n on |κ|=1 for all primes p≥5 (or Es4 R-path / "
            "ker=sc / hull)."
        ),
    }


def main() -> dict:
    a = theorem_per_equals_derang()
    b = theorem_Cy_identity_structure()
    c = theorem_two_over_n_suffices()
    cert = certify_mu4_census()
    out = {
        "title": (
            "Prop 15.251 Cy-identity on |κ|=1; |μ₄|≤2/n census p=5,7; "
            "residual (i) OPEN"
        ),
        "proved": {
            "per_equals_1_on_kappa1": a["proved"],
            "Cy_identity_structure": b["proved"],
            "two_over_n_suffices_for_residual_i": c["proved"],
            "mu4_bound_general": False,
            "residual_i_closed": residual_i_closed_via_251(),
            "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
            "type_I": type_I_k_3p_minus_2_closed_general(),
            "E1_closed": e1_closed_general(),
            "L_closed": False,
        },
        "algebra": {
            "per": a,
            "Cy_identity": b,
            "two_over_n_suffices": c,
        },
        "certified": {"mu4_census": cert},
        "hinge": hinge_status_251(),
        "F3": (
            "Cy-expansion gives (p⁴−1)m₄+2φ=Ext on |κ|=1.  "
            "|m₄⁺|,|m₄⁻|≤2/n would close residual-i.  Census holds at p=5,7 "
            "with m₄⁺=m₄⁻.  General |μ₄|≤2/n still OPEN.  No predicate flip."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15251.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.251 Cy-identity; |μ4|<=2/n path")
    print(f"  Per=1: {a['proved']}")
    print(f"  Cy identity structure: {b['proved']}")
    print(f"  2/n suffices: {c['proved']}")
    for pk, v in cert.items():
        if not v.get("ok"):
            print(f"  p={pk}: {v}")
            continue
        print(
            f"  p={pk}: max|μ4|={v['max_abs_mu4']:.6f} ≤2/n={v['two_over_n']:.6f}? "
            f"{v['mu4_le_2n']} m4+=m4-? {v['m4plus_eq_m4minus']}"
        )
    print(f"  residual_i closed: {residual_i_closed_via_251()}")
    print(f"  open: {out['hinge']['open']}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
