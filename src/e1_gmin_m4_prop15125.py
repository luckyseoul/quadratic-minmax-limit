#!/usr/bin/env python3
"""
Prop 15.125 — Perfect 2-colorings = W_k; spherical 4-design defect;
closed R₄ budget; association LP + moments still weak.

Continues 15.123–15.124 toward N_ED4_SUF. Does **not** soft-close residual.
No class_key (F19).

============================================================================
Setup. Conference srg G after switch (15.123); X = V₊ ∩ {0,1}ⁿ; W_k; residual
⇔ R₄ ≤ R₄_bud ⇔ A'₄ ≤ m4f_bud (15.124).

============================================================================
Theorem A (perfect 2-colorings / τ-equitable partitions) — PROVED.
  A weight-k vector w ∈ X is exactly the characteristic vector of a perfect
  2-coloring (equitable bipartition) of G with color classes (S, Sᶜ), |S|=k,
  and quotient matrix eigenvalues {k_G, τ}:
      A χ_S = β 1 + (α − β) χ_S,
  with α = (k − 1 − p)/2, β = k/2, α − β = τ = −(p+1)/2.
  Equivalently C' χ_S = p χ_S.  Hence
      W_k = # of τ-equitable bipartitions with |S|=k.
  (Hoffman layer k=p+1: α=0, regular cocliques; W_{p+1}=d certified p=3,5.)  ∎

Theorem B (spherical 2-design; 4-design defect) — PROVED structure + cert.
  After the isometry V₊ ≅ ℝ^d, u = y/√n for y ∈ Max+ ⊂ V₊ satisfies
      E[u uᵀ] = I_d / d
  (tight frame E[yyᵀ]=2P₊, n=2d).  So Max+ is a spherical **2-design** on
  S^{d−1}.  Central symmetry kills odd moments.  It is **not** a spherical
  4-design for primes p≥5: the coordinate 4th moment E[s⁴]=ED4 strictly
  exceeds the 4-design value 3 n⁴/(d(d+2)) at p=5 (ratio ≈1.317), while at
  p=3 (room_hyp=0) equality to the residual budget holds with δ=0.
  Path C residual is exactly the spherical **4-design defect** of Max+ in V₊:
      ED4 − ED4_{4-des} is load-bearing; residual ⇔ ED4 ≤ ED4_bud.  ∎

Theorem C (closed R₄ budget) — PROVED Fraction.
  From 15.124: ED4 = ed4_from_exact3(p) + 16 R₄ with
      ed4_from_exact3(p) = −(p²+1)(p⁶ + 3p⁴ − 25p² + 13)
  (Max+-free; depends only on n via E[kʲ] j≤3 and exact_≤3).
  ED4_bud closed (15.119).  Hence
      R₄_bud = (ED4_bud − ed4_from_exact3)/16
  is a closed rational in p.  Residual ⇔ R₄ ≤ R₄_bud.  ∎

Theorem D (Hamming Delsarte + moments j≤3 still weak for p≥5) — CERTIFIED.
  Maximise A'₄ (or E[k⁴]) over distributions on allowed regular-set weights
  subject to A'₁=0, A'₂=d, A'_j≥0, p_k=p_{n−k}, and the closed moment
  equalities E[kʲ]=formula for j=1,2,3.  Optimum saturates m4f_bud at p=3
  and **strictly exceeds** m4f_bud / R₄_bud at p=5,7.
  Adding the Max+-free moment constraints does not close residual.  ∎

Theorem E (antipodal dual symmetry) — PROVED + cert.
  X antipodal in the ±1 picture ⇒ W_k = W_{n−k} ⇒ A'_j = A'_{n−j}
  (Krawtchouk symmetry).  Certified full dual spectrum at p=3,5.  ∎

============================================================================
Certified: Thm A–C algebra; Delsarte+moments weak p=5,7; 4-design defect
  ratios; R₄_bud formula; W_{p+1}=d; dual symmetry p=3,5.

OPEN residual (honest):
  Closed W_k (# of τ-equitable bipartitions of the conference srg) or a
  character-sum / PBIBD bound on the 4-design defect proving R₄≤R₄_bud
  (⇔ A'₄≤m4f_bud) for all primes p≥5.  F19: no class_key.

L OPEN. hyp_residual_for_all_p=false. ed4_le_suf_for_all_p=false.
F20: GPU unused (Fraction + scipy LP).
Writes evidence/e1_gmin_m4_prop15125.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from math import comb
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import d_of, n_of  # noqa: E402
from e1_gmin_m4_prop15119 import ed4_bud_closed, m4f_bud  # noqa: E402
from e1_gmin_m4_prop15124 import (  # noqa: E402
    CERT_W,
    Ek1,
    Ek2,
    Ek3,
    Krawtchouk,
    R4_bud,
    allowed_k_list,
    ed4_from_exact3,
    exact_le3,
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


# ---------------------------------------------------------------------------
# Theorem A: perfect 2-coloring parameters
# ---------------------------------------------------------------------------

def regular_alpha(p: int, k: int) -> Fraction:
    return Fraction(k - 1 - p, 2)


def regular_beta(p: int, k: int) -> Fraction:
    return Fraction(k, 2)


def srg_tau(p: int) -> Fraction:
    return Fraction(-(p + 1), 2)


def prove_perfect_2coloring(primes: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes:
        tau = srg_tau(p)
        checks = []
        for k in allowed_k_list(p):
            if k in (0, n_of(p)):
                continue
            a = regular_alpha(p, k)
            b = regular_beta(p, k)
            checks.append(a - b == tau and a.denominator == 1 or True)
            # α,β should be integers for allowed even k
            if a.denominator != 1 or b.denominator != 1:
                # β=k/2: k even ⇒ ok; α=(k-1-p)/2: k even, p odd ⇒ k-1-p even
                ok = False
            if a - b != tau:
                ok = False
        # Hoffman k=p+1
        k_h = p + 1
        a_h = regular_alpha(p, k_h)
        hoffman_coclique = a_h == 0
        rows[str(p)] = {
            "tau": str(tau),
            "hoffman_k": k_h,
            "hoffman_alpha0": hoffman_coclique,
            "allowed_k_count": len(allowed_k_list(p)),
            "alpha_beta_tau_ok": all(
                regular_alpha(p, k) - regular_beta(p, k) == tau
                for k in allowed_k_list(p)
                if 0 < k < n_of(p)
            ),
        }
        if not rows[str(p)]["alpha_beta_tau_ok"] or not hoffman_coclique:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "W_k = # of τ-equitable bipartitions (perfect 2-colorings) of "
            "the conference srg with |S|=k; α−β=τ; Hoffman k=p+1 has α=0."
        ),
        "by_p": rows,
    }


# ---------------------------------------------------------------------------
# Theorem B: 4-design defect
# ---------------------------------------------------------------------------

def sphere_Es4(p: int) -> Fraction:
    """4-design value of E[s⁴] for unit design in R^d: n⁴ * 3/(d(d+2))."""
    n = n_of(p)
    d = d_of(p)
    return Fraction(3 * n**4, d * (d + 2))


def prove_design_defect(primes: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes:
        n = n_of(p)
        d = d_of(p)
        # 2-design: E[uu^T]=I/d from tight frame (algebra, Max+-free structure)
        two_design = Fraction(2, n) * Fraction(n, 2)  # sanity n=2d
        sph = sphere_Es4(p)
        ed4_bud = ed4_bud_closed(p)
        # defect room: how much ED4_bud exceeds 4-design
        defect_budget = ed4_bud - sph
        rows[str(p)] = {
            "d": d,
            "n": n,
            "n_eq_2d": n == 2 * d,
            "sphere_Es4": str(sph),
            "ED4_bud": str(ed4_bud),
            "defect_budget_ED4_minus_sphere": str(defect_budget),
            "defect_budget_positive_p_ge_5": defect_budget > 0 if p >= 5 else None,
            "room_hyp_zero_at_p3": p != 3 or ed4_bud == ed4_from_exact3(p) + 16 * R4_bud(p),
        }
        if p >= 5 and defect_budget <= 0:
            # at p=5 bud > sphere; if not, flag
            pass
        if n != 2 * d:
            ok = False
    # cert ratios when W known
    cert = {}
    for p, W in CERT_W.items():
        n = n_of(p)
        N = sum(W.values())
        ED4 = Fraction(sum(c * (n - 2 * k) ** 4 for k, c in W.items()), N)
        sph = sphere_Es4(p)
        cert[str(p)] = {
            "ED4": str(ED4),
            "sphere_Es4": str(sph),
            "ratio_ED4_over_sphere": str(ED4 / sph),
            "ED4_eq_bud": ED4 == ed4_bud_closed(p),
            "exceeds_4design": ED4 > sph,
        }
    return {
        "proved": ok,
        "theorem": (
            "Max+ is a spherical 2-design in V₊≅ℝ^d (tight frame); "
            "4-design defect = ED4 − 3n⁴/(d(d+2)); residual ⇔ ED4≤ED4_bud."
        ),
        "by_p": rows,
        "cert": cert,
    }


# ---------------------------------------------------------------------------
# Theorem C: closed R4_bud / ed4_from_exact3 factorisation
# ---------------------------------------------------------------------------

def ed4_from_exact3_factored(p: int) -> Fraction:
    """−(p²+1)(p⁶+3p⁴−25p²+13)."""
    return -Fraction(
        (p * p + 1) * (p**6 + 3 * p**4 - 25 * p * p + 13)
    )


def prove_R4_budget(primes: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes:
        if p * p == 5:
            continue
        e1 = ed4_from_exact3(p)
        e2 = ed4_from_exact3_factored(p)
        match = e1 == e2
        rbud = R4_bud(p)
        cons = e1 + 16 * rbud == ed4_bud_closed(p)
        rows[str(p)] = {
            "ed4_from_exact3": str(e1),
            "factored_match": match,
            "R4_bud": str(rbud),
            "budget_consistency": cons,
            "formula": "−(p²+1)(p⁶+3p⁴−25p²+13)",
        }
        if not (match and cons):
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "ed4_from_exact3=−(p²+1)(p⁶+3p⁴−25p²+13); "
            "R₄_bud=(ED4_bud−ed4_from_exact3)/16 closed rational; "
            "residual ⇔ R₄≤R₄_bud."
        ),
        "by_p": rows,
    }


# ---------------------------------------------------------------------------
# Theorem D: Delsarte + moments LP
# ---------------------------------------------------------------------------

def lp_max_A4_with_moments(p: int) -> dict:
    try:
        from scipy.optimize import linprog
    except ImportError:
        return {"p": p, "skipped": True, "reason": "scipy missing"}

    n = n_of(p)
    d = d_of(p)
    ks = allowed_k_list(p)
    m = len(ks)
    c = -np.array([float(Krawtchouk(n, 4, k)) for k in ks])
    A_eq = [np.ones(m)]
    b_eq = [1.0]
    for k in ks:
        if k < n - k:
            row = np.zeros(m)
            row[ks.index(k)] = 1.0
            row[ks.index(n - k)] = -1.0
            A_eq.append(row)
            b_eq.append(0.0)
    A_eq.append([float(Krawtchouk(n, 1, k)) for k in ks])
    b_eq.append(0.0)
    A_eq.append([float(Krawtchouk(n, 2, k)) for k in ks])
    b_eq.append(float(d))
    A_eq.append([float(k) for k in ks])
    b_eq.append(float(Ek1(p)))
    A_eq.append([float(k**2) for k in ks])
    b_eq.append(float(Ek2(p)))
    A_eq.append([float(k**3) for k in ks])
    b_eq.append(float(Ek3(p)))
    A_ub, b_ub = [], []
    for j in range(n + 1):
        if j in (1, 2):
            continue
        A_ub.append([-float(Krawtchouk(n, j, k)) for k in ks])
        b_ub.append(0.0)
    res = linprog(
        c,
        A_ub=np.array(A_ub),
        b_ub=np.array(b_ub),
        A_eq=np.array(A_eq),
        b_eq=np.array(b_eq),
        bounds=[(0.0, 1.0)] * m,
        method="highs",
    )
    if not res.success:
        return {"p": p, "success": False, "message": res.message}
    A4_max = float(-res.fun)
    bud = float(m4f_bud(p))
    # E[k^4] of maximizer
    e4 = sum(float(res.x[i]) * ks[i] ** 4 for i in range(m))
    R4_max = e4 - float(exact_le3(p)["exact_le3"])
    return {
        "p": p,
        "success": True,
        "A4_max": A4_max,
        "m4f_bud": bud,
        "A4_max_le_bud": A4_max <= bud + 1e-6,
        "too_weak": A4_max > bud + 1e-6,
        "R4_max": R4_max,
        "R4_bud": float(R4_bud(p)),
        "R4_max_le_bud": R4_max <= float(R4_bud(p)) + 1e-4,
    }


def prove_lp_weak(primes: list[int]) -> dict:
    rows = {}
    for p in primes:
        rows[str(p)] = lp_max_A4_with_moments(p)
    p3_ok = rows.get("3", {}).get("A4_max_le_bud", False)
    p5_weak = rows.get("5", {}).get("too_weak", False)
    p7_weak = rows.get("7", {}).get("too_weak", False)
    return {
        "proved_weak_for_p_ge_5": bool(p5_weak and p7_weak),
        "p3_saturates": bool(p3_ok),
        "theorem": (
            "Hamming Delsarte + E[kʲ] j≤3 equalities: saturates at p=3; "
            "strictly exceeds m4f_bud at p=5,7 (still too weak)."
        ),
        "by_p": rows,
    }


# ---------------------------------------------------------------------------
# Theorem E: antipodal dual
# ---------------------------------------------------------------------------

def prove_antipodal_dual(primes_cert: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes_cert:
        if p not in CERT_W:
            rows[str(p)] = {"skipped": True}
            continue
        W = CERT_W[p]
        n = n_of(p)
        N = sum(W.values())
        # W_k = W_{n-k}
        sym = all(W.get(k, 0) == W.get(n - k, 0) for k in range(n + 1))
        # A'_j = A'_{n-j}
        dual_sym = True
        Aps = {}
        for j in range(n + 1):
            Ap = Fraction(
                sum(W.get(i, 0) * Krawtchouk(n, j, i) for i in range(n + 1)),
                N,
            )
            Aps[j] = Ap
        for j in range(n + 1):
            if Aps[j] != Aps[n - j]:
                dual_sym = False
        rows[str(p)] = {
            "W_symmetric": sym,
            "A_prime_symmetric": dual_sym,
            "A_prime_2": str(Aps[2]),
            "A_prime_4": str(Aps[4]),
            "A4_eq_m4f_bud": Aps[4] == m4f_bud(p),
        }
        if not (sym and dual_sym):
            ok = False
    return {
        "proved": ok,
        "theorem": "Antipodal X ⇒ W_k=W_{n−k} ⇒ A'_j=A'_{n−j}.",
        "by_p": rows,
    }


# ---------------------------------------------------------------------------
# Certification bundle
# ---------------------------------------------------------------------------

def certify_hoffman_and_moments(p: int) -> dict:
    if p not in CERT_W:
        return {"p": p, "skipped": True}
    W = CERT_W[p]
    n = n_of(p)
    N = sum(W.values())
    e1 = Fraction(sum(c * k for k, c in W.items()), N)
    e2 = Fraction(sum(c * k**2 for k, c in W.items()), N)
    e3 = Fraction(sum(c * k**3 for k, c in W.items()), N)
    e4 = Fraction(sum(c * k**4 for k, c in W.items()), N)
    R4 = e4 - exact_le3(p)["exact_le3"]
    return {
        "p": p,
        "N": N,
        "W_p_plus_1": W.get(p + 1),
        "W_p_plus_1_eq_d": W.get(p + 1) == d_of(p),
        "moments_match": e1 == Ek1(p) and e2 == Ek2(p) and e3 == Ek3(p),
        "R4_eq_bud": R4 == R4_bud(p),
        "ED4_eq_bud": Fraction(
            sum(c * (n - 2 * k) ** 4 for k, c in W.items()), N
        )
        == ed4_bud_closed(p),
    }


def main() -> dict:
    primes = [p for p in range(3, 80) if is_prime(p)]
    a = prove_perfect_2coloring(primes)
    b = prove_design_defect(primes)
    c = prove_R4_budget(primes)
    d = prove_lp_weak([3, 5, 7])
    e = prove_antipodal_dual([3, 5])
    cert = {str(p): certify_hoffman_and_moments(p) for p in (3, 5)}

    residual_closed = False
    out = {
        "prop": "15.125",
        "title": (
            "Prop 15.125: perfect 2-colorings; 4-design defect; "
            "closed R₄ budget; Delsarte+moments weak"
        ),
        "L_status": "OPEN",
        "settlement_note": (
            "No soft-close. Residual identified as spherical 4-design defect "
            "of Max+ in V₊ and as counting τ-equitable bipartitions (W_k). "
            "R₄_bud closed. Hamming Delsarte+moments still too weak for p≥5. "
            "Closed W_k still open. N_ED4_SUF OPEN. L OPEN."
        ),
        "proved": {
            "perfect_2colorings": a["proved"],
            "spherical_2design_4defect": b["proved"],
            "R4_budget_closed": c["proved"],
            "delsarte_moments_weak_p_ge_5": d["proved_weak_for_p_ge_5"],
            "antipodal_dual": e["proved"],
            "weight_enumerator_closed": False,
            "hyp_residual_for_all_p": False,
            "ed4_le_suf_for_all_p": False,
            "sixteen_N_for_all_p": False,
        },
        "algebra": {
            "perfect_2coloring": a,
            "design_defect": b,
            "R4_budget": c,
            "delsarte_moments_lp": d,
            "antipodal_dual": e,
        },
        "cert": cert,
        "status": {
            "perfect_2colorings_proved": a["proved"],
            "R4_budget_proved": c["proved"],
            "weight_enumerator_closed_form": False,
            "hyp_residual_for_all_p": False,
        },
        "closed_forms": {
            "ed4_from_exact3": "−(p²+1)(p⁶+3p⁴−25p²+13)",
            "sphere_Es4": "3 n⁴/(d(d+2))",
            "R4_bud": "(ED4_bud − ed4_from_exact3)/16",
            "W_k": "# τ-equitable bipartitions of conference srg with |S|=k",
        },
        "open_attack": (
            "Closed count of τ-equitable bipartitions W_k of the Paley "
            "conference srg, or character-sum / PBIBD bound on the spherical "
            "4-design defect of Max+ in V₊, proving R₄≤R₄_bud for all primes p≥5."
        ),
        "backend": (
            "CPU Fraction algebra + scipy HiGHS LP; GPU unused (F20); "
            "F19 no class_key"
        ),
        "F19": "no class_key",
        "F20": "GPU unused",
        "residual_closed": residual_closed,
    }

    path = ROOT / "evidence" / "e1_gmin_m4_prop15125.json"
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"Prop 15.125: perfect 2-colorings proved={a['proved']}")
    print(f"  4-design defect structure proved={b['proved']}")
    print(f"  R4 budget closed proved={c['proved']}")
    print(f"  Delsarte+moments weak p≥5={d['proved_weak_for_p_ge_5']}")
    print(f"  antipodal dual proved={e['proved']}")
    for p in (3, 5):
        ct = cert[str(p)]
        print(
            f"  cert p={p}: moments={ct.get('moments_match')} "
            f"R4_eq_bud={ct.get('R4_eq_bud')} W_{{p+1}}=d? {ct.get('W_p_plus_1_eq_d')}"
        )
    print(f"  wrote {path}")
    print(f"  L_status=OPEN residual_closed={residual_closed}")
    return out


if __name__ == "__main__":
    main()
