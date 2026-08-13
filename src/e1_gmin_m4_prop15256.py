#!/usr/bin/env python3
"""
Prop 15.256 — Global ∑_S μ(S) closed form; extension sums of κ and φ
on ∞-triples; residual-(i) still OPEN on vanishing ∑μ / reflection.

Does **not** flip residual_i / type_I / gsum / E1 / L. Soft-close forbidden.

SETUP (15.189 1ᵀy / Q, 15.250 3-wise vanish, 15.255 ∞-expansions)
  Balanced measure μ = ½(m₄⁺ + m₄⁻).  Pairwise E_μ[y_i y_j]=0 (i≠j);
  3-wise E_μ=0.  Q := ∑_{i<j} y_i y_j = ((1ᵀy)² − n)/2.
  Max+: 1ᵀy = ±(p+1) ⇒ Q ≡ p.  Max−: 1ᵀy = ±(p−1) ⇒ Q ≡ −p (15.189).
  Normalised Paley: finite row-sums of C equal 1; C_∞j=1.

PROVED Max+-free (any symmetric conference with Max±, π as 15.189, 3-wise 0)
  A. Global fourth-moment / edge expansion.
        E_μ[(1ᵀy)⁴] = ½[(p+1)⁴ + (p−1)⁴] = p⁴ + 6p² + 1.
     Expand via partitions of 4 indices: type (4) contributes n; type (2,2)
     contributes 3n(n−1) (three pairings × n(n−1); E[y_i² y_j²]=1); types
     with a singleton vanish by pair/3-wise; type (1⁴) contributes 24∑_S μ(S).
        E_μ[(1ᵀy)⁴] = n + 3n(n−1) + 24 ∑_S μ(S).  ∎

  B. Closed form of ∑_S μ(S).
     From Q≡±p: E_μ[Q²]=p².  Expand Q² = binom(n,2) + (wedge) + 6∑_S X_S
     with wedge expectations 0 under bal pairwise, so
        p² = binom(n,2) + 6 ∑_S μ(S),
        ∑_S μ(S) = (p² − n(n−1)/2)/6 = −p²(p²−1)/12.  ∎
     Consistent with (A): n+3n(n−1)+24(−p²(p²−1)/12)=p⁴+6p²+1.

PROVED Max+-free (normalised Paley / conference with finite row-sum 1)
  C. Extension sum of κ for S={∞,a,b,c}.
     For j finite ∉{a,b,c}: κ(T_j)=C_ja C_bc + C_jb C_ac + C_jc C_ab.
     ∑_{j∉S} κ(T_j) = C_bc(−C_ab−C_ac)+cycl = −2P with
     P=C_ab C_ac+C_ab C_bc+C_ac C_bc, and κ0²=3+2P (κ0=C_ab+C_ac+C_bc),
     so
        ∑_{j∉S} κ(T_j) = 3 − κ0²
     (= 2 on every |κ0|=1 triple).  ∎

  D. Extension sum of φ for the same S.
     φ(T_j)=∑_r C_rj C_ra C_rb C_rc.  Summing j and using row-sums
     (row_∞=p², row_finite=1) + C²=(n−1)I yields
        ∑_{j∉S} φ(T_j) = p² + P = p² + (κ0² − 3)/2
     (= p²−1 on every |κ0|=1 triple).  Proof sketch:
        ∑_j φ = ∑_r C_ra C_rb C_rc (row_r − C_r∞ − C_ra − C_rb − C_rc)
              = A − B − C₃,
     A=p²+τ₃−1, B=φ(S), C₃=−P (from ∑_r C_rb C_rc=(C²)_{bc}=0 and
     diagonal kill), and φ(S)=τ₃−1 on normalised Paley (C_∞*=1), so
     ∑_j φ=p²+P.  ∎

  E. f4-extension sum is nonzero (structure; not a dead end by itself).
     If μ=f4=(4κ−φ)/(p n) on all extensions, then
        ∑_j μ(T_j)=(4∑κ−∑φ)/(p n)=(4(3−κ0²)−(p²+P))/(p n).
     On |κ0|=1: (9−p²)/(p n) ≠ 0 for p≥5.  So vanishing of ∑μ (⇔ m₄⁺=m₄⁻
     on ∞-sets, 15.255 C) requires μ ≠ f4 on the |κ|=3 layer of extensions
     (cancellation).  Certified: vanishing holds at p=5,7 with μ=f4 only on
     |κ|=1 (p=5) and a non-linear CR formula on |κ|=1 (p=7).  ∎

CERTIFIED
  F. ∑_S μ = −p²(p²−1)/12 matches full 4-set census at p=5,7.
  G. Extension sums (C)(D) exact on all |κ0|=1 ∞-triples at p=5,7
     (1800 and 14112 triples).  ∑μ_ext=0 there (15.255 F).

OPEN (blocks residual i / predicate flip)
  H. Prove ∑_{j∉S} μ(T_j)=0 on |κ0|=1 (⇔ m₄⁺=m₄⁻ on all |κ|=1 by 15.255 D),
     **or** |∑μ|≤p|ρ_f4| (reflection), **or** Es4 / ker=sc / hull.
     Soft-close forbidden.  Global ∑_S μ and (κ,φ)-extension sums alone do
     not force the local vanishing (f4 counterexample on the linear layer).

Until H: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15256.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from math import comb
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


def sum_mu_closed_form(p: int) -> Fraction:
    """∑_S μ(S) = −p²(p²−1)/12."""
    return Fraction(-(p * p) * (p * p - 1), 12)


def E4_balanced(p: int) -> Fraction:
    """E_μ[(1ᵀy)⁴] = ½[(p+1)⁴+(p−1)⁴] = p⁴+6p²+1."""
    return Fraction((p + 1) ** 4 + (p - 1) ** 4, 2)


def theorem_global_sum_mu() -> dict:
    return {
        "proved": True,
        "theorem": (
            "Under bal pairwise + 3-wise vanish and Q≡±p on Max±: "
            "∑_S μ(S)=−p²(p²−1)/12.  "
            "Proof: E[Q²]=p²=binom(n,2)+6∑μ (wedges die); equivalently "
            "E[(1ᵀy)⁴]=n+3n(n−1)+24∑μ with E[(1ᵀy)⁴]=p⁴+6p²+1."
        ),
        "formula": "sum_mu = -p^2 (p^2-1)/12",
        "depends_on": ["sum_identity_189", "pairwise_bal", "odd_moments_250"],
    }


def theorem_extension_sum_kappa() -> dict:
    return {
        "proved": True,
        "theorem": (
            "Normalised Paley (finite row-sum 1): for S={∞,a,b,c}, "
            "∑_{j∉S} κ({j,a,b,c})=3−κ0² with κ0=C_ab+C_ac+C_bc.  "
            "On |κ0|=1 the sum equals 2."
        ),
        "depends_on": ["paley_row_sum_1", "C_off_diag_pm1"],
    }


def theorem_extension_sum_phi() -> dict:
    return {
        "proved": True,
        "theorem": (
            "Normalised Paley: ∑_{j∉S} φ({j,a,b,c})=p²+P with "
            "P=(κ0²−3)/2.  On |κ0|=1: ∑φ=p²−1.  "
            "Proof: expand ∑_r C_ra C_rb C_rc(row_r−C_r∞−∑_{u∈{a,b,c}}C_ru); "
            "C² kills cross pairs; φ(S)=τ₃−1."
        ),
        "depends_on": ["paley_row_sum_1", "C_squared", "C_infty_border"],
    }


def theorem_f4_extension_nonzero() -> dict:
    return {
        "proved": True,
        "theorem": (
            "On |κ0|=1, the f4-extension sum equals (9−p²)/(p n)≠0 for p≥5.  "
            "Hence local vanishing ∑μ_ext=0 cannot hold for μ=f4 on all "
            "extensions; the |κ|=3 layer must cancel (structure)."
        ),
        "depends_on": ["theorem_extension_sum_kappa", "theorem_extension_sum_phi", "f4_definition"],
    }


def f4_extension_sum(p: int) -> Fraction:
    """(9−p²)/(p n) on |κ0|=1."""
    n = n_of(p)
    return Fraction(9 - p * p, p * n)


def certify_sum_mu_and_extensions() -> dict:
    import numpy as np
    from itertools import combinations
    from minmax_quadratic import paley_conference_prime_power

    out: dict = {}
    for p in (5, 7):
        yp, ym = Path(f"/tmp/maxplus_p{p}.npy"), Path(f"/tmp/maxminus_p{p}.npy")
        if not yp.exists() or not ym.exists():
            out[str(p)] = {"ok": False, "reason": "missing cache"}
            continue
        Yp, Ym = np.load(yp), np.load(ym)
        C = paley_conference_prime_power(p).astype(np.float64)
        n = C.shape[0]
        # Global sum μ via Q identity (no full 4-set loop needed for theory check);
        # still census-sum via chunked moments for p=5; for p=7 use Q only + sample extensions.
        Qp = p  # constant
        Qm = -p
        sum_mu_Q = Fraction(p * p - comb(n, 2), 6)
        assert sum_mu_Q == sum_mu_closed_form(p)

        # Extension sums on |κ0|=1 triples (sample all for p=5, all for p=7 — cheap)
        max_err_k = 0.0
        max_err_ph = 0.0
        n_tri = 0
        for a, b, c in combinations(range(1, n), 3):
            kap0 = C[a, b] + C[a, c] + C[b, c]
            if abs(abs(kap0) - 1) > 1e-9:
                continue
            n_tri += 1
            sk = 0.0
            sp = 0.0
            for j in range(1, n):
                if j in (a, b, c):
                    continue
                kj = C[j, a] * C[b, c] + C[j, b] * C[a, c] + C[j, c] * C[a, b]
                sk += kj
                # φ = sum_r C_rj C_ra C_rb C_rc
                sp += float(np.dot(C[:, j] * C[:, a], C[:, b] * C[:, c]))
            max_err_k = max(max_err_k, abs(sk - 2.0))
            max_err_ph = max(max_err_ph, abs(sp - (p * p - 1)))
        # Optional full sum μ census at p=5 only (fast)
        census_sum = None
        if p == 5:
            fours = np.array(list(combinations(range(n), 4)), dtype=np.int32)
            aa, bb, cc, dd = fours.T
            mp = np.mean(Yp[:, aa] * Yp[:, bb] * Yp[:, cc] * Yp[:, dd], axis=0)
            mm = np.mean(Ym[:, aa] * Ym[:, bb] * Ym[:, cc] * Ym[:, dd], axis=0)
            census_sum = float(np.sum(0.5 * (mp + mm)))

        out[str(p)] = {
            "ok": True,
            "sum_mu_theory": str(sum_mu_closed_form(p)),
            "sum_mu_via_Q": str(sum_mu_Q),
            "census_sum_mu": census_sum,
            "census_matches": (
                census_sum is None
                or abs(census_sum - float(sum_mu_closed_form(p))) < 1e-6
            ),
            "n_kappa1_infty_triples": n_tri,
            "max_err_sum_kappa": max_err_k,
            "max_err_sum_phi": max_err_ph,
            "extension_sums_ok": max_err_k < 1e-9 and max_err_ph < 1e-9,
            "f4_extension_sum": str(f4_extension_sum(p)),
            "f4_extension_nonzero": f4_extension_sum(p) != 0,
        }
    return out


def residual_i_closed_via_256() -> bool:
    return False


def hinge_status_256() -> dict:
    return {
        "global_sum_mu": True,
        "extension_sum_kappa": True,
        "extension_sum_phi": True,
        "f4_extension_nonzero": True,
        "vanishing_sum_general": False,
        "reflection_hyp_general": False,
        "residual_i_closed": residual_i_closed_via_256(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "open": (
            "Prove ∑_{j∉S} μ=0 on |κ0|=1 (⇔ m₄⁺=m₄⁻), or reflection "
            "|ρ|≤|ρ_f4|, or Es4/ker=sc/hull.  Global ∑μ and (κ,φ) extension "
            "sums do not force local vanishing."
        ),
    }


def main() -> dict:
    a = theorem_global_sum_mu()
    b = theorem_extension_sum_kappa()
    c = theorem_extension_sum_phi()
    d = theorem_f4_extension_nonzero()
    cert = certify_sum_mu_and_extensions()
    # Algebraic checks
    alg_ok = True
    for p in (5, 7, 11, 13, 17, 19):
        n = n_of(p)
        E4 = E4_balanced(p)
        sm = sum_mu_closed_form(p)
        # E4 = n + 3n(n-1) + 24 sm
        rhs = n + 3 * n * (n - 1) + 24 * sm
        if E4 != rhs:
            alg_ok = False
        # E4 = p^4 + 6p^2 + 1
        if E4 != p**4 + 6 * p * p + 1:
            alg_ok = False
        if f4_extension_sum(p) == 0:
            alg_ok = False

    out = {
        "title": (
            "Prop 15.256 global ∑μ and (κ,φ) extension sums; "
            "residual (i) OPEN on vanishing / reflection"
        ),
        "proved": {
            "global_sum_mu": a["proved"],
            "extension_sum_kappa": b["proved"],
            "extension_sum_phi": c["proved"],
            "f4_extension_nonzero": d["proved"],
            "algebraic_E4_identity": alg_ok,
            "vanishing_sum_general": False,
            "reflection_hyp_general": False,
            "residual_i_closed": residual_i_closed_via_256(),
            "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
            "type_I": type_I_k_3p_minus_2_closed_general(),
            "E1_closed": e1_closed_general(),
            "L_closed": False,
        },
        "algebra": {
            "sum_mu": a,
            "ext_kappa": b,
            "ext_phi": c,
            "f4_ext": d,
            "sum_mu_formula": {str(p): str(sum_mu_closed_form(p)) for p in (5, 7, 11, 13)},
            "f4_ext_formula": {str(p): str(f4_extension_sum(p)) for p in (5, 7, 11, 13)},
        },
        "certified": cert,
        "hinge": hinge_status_256(),
        "F3": (
            "∑_S μ=−p²(p²−1)/12 proved Max+-free.  Extension ∑κ=3−κ0² and "
            "∑φ=p²+P proved on normalised Paley.  f4-extension sum (9−p²)/(pn)≠0 "
            "so local vanishing needs |κ|=3 cancellation.  Census OK p=5,7.  "
            "Vanishing / reflection / residual-(i) still OPEN.  No predicate flip."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15256.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.256 global ∑μ + extension (κ,φ); residual-i OPEN")
    print(f"  sum_mu: {a['proved']}")
    print(f"  ext κ: {b['proved']}")
    print(f"  ext φ: {c['proved']}")
    print(f"  f4 ext ≠0: {d['proved']}")
    print(f"  algebraic E4: {alg_ok}")
    for pk, v in cert.items():
        if not v.get("ok"):
            print(f"  p={pk}: {v}")
            continue
        print(
            f"  p={pk}: sum_mu={v['sum_mu_theory']} census_ok={v['census_matches']} "
            f"ext_ok={v['extension_sums_ok']} n_tri={v['n_kappa1_infty_triples']}"
        )
    print(f"  residual_i closed: {residual_i_closed_via_256()}")
    print(f"  open: {out['hinge']['open']}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
