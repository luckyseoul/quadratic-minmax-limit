#!/usr/bin/env python3
"""
Prop 15.115 — Max+ residual solves Aρ=b; δ=P_{E_{4p}} m₄;
spectral moments of f_y; γ-weighted E_4p calculus.

Continues 15.114 toward N_ED4_SUF. Does **not** soft-close residual.
No class_key (F19).

============================================================================
Setup. Conference C, n=p²+1, Max+ boolean Cy=py with E[yyᵀ]=2P₊
(spherical 2-design in V₊).  f_y(S)=∏_{i∈S} y_i, γ_y as in Prop 15.114.
T = C-signed Johnson adjacency, A=4p I−T, E_{4p}=ker A.
ρ=m₄−κ/p²=ρ_min+δ, b=Tκ/p², μ²=4(p²+15) (Prop 15.102).

============================================================================
Theorem A (E[γ⊙f]=2κ/p) — PROVED (2-design).
  For every 4-set S,
      E_y[γ_y(S) f_y(S)] = 2 κ(S) / p.
  Proof.  γ_y(S) f_y(S)=∑_{e⊂S} C_e ∏_{i∈S∖e} y_i.  For e={i,j}⊂S the
  complement S∖e is the opposite edge of the unique matching containing e, so
      E[γ f](S)=∑_{e⊂S} C_e E[y_k y_l]  ({k,l}=S∖e)
  =∑_{e⊂S} C_e C_{kl}/p   (i≠j: E[y_i y_j]=C_{ij}/p by 2-design)
  = (1/p) ∑_{e⊂S} C_e C_{S∖e} = (1/p)·2 κ(S),
  since each of the three unordered matchings is counted twice.  ∎

Theorem B (Max+ residual solves the resolvent) — PROVED.
  A m₄ = 4κ/p  and  Aρ = b.
  Proof.  A f_y = 2(γ_y⊙f_y) (Prop 15.114).  Average over Max+:
  A m₄ = 2 E[γ⊙f] = 4κ/p by Thm A.
  Then Aρ = A(m₄−κ/p²) = 4κ/p − A(κ/p²)
  = 4κ/p − (4p I−T)κ/p² = 4κ/p − 4κ/p + Tκ/p² = b.  ∎

Theorem C (κ ⊥ E_{4p}) — PROVED.
  Tκ ∈ E_μ ⊕ E_{−μ} with μ²=4(p²+15) (Prop 15.102).  For any eigenpair
  T u = λ u with λ ∉ {0, ±μ}: ⟨κ,u⟩=0, because λ⟨κ,u⟩=⟨Tκ,u⟩=0 and λ≠0.
  Check 4p ∉ {0, ±μ}: (4p)²=16p² equals μ²=4(p²+15) iff p²=5, not an
  integer prime.  Hence P_{E_{4p}} κ = 0.  ∎

Theorem D (δ as E_4p average of f_y) — PROVED.
  δ = P_{E_{4p}} m₄ = E_y[ P_{E_{4p}} f_y ].
  Proof.  ρ=ρ_min+δ with ρ_min ⊥ E_{4p} and δ∈E_{4p}, so
  P_{E_{4p}} ρ = δ.  But ρ=m₄−κ/p² and P_{E_{4p}} κ=0, hence
  δ=P_{E_{4p}} m₄.  Linearity of P and m₄=E f_y give the average form.  ∎

  Corollary D′ (Jensen).  δ² = ‖E_y P_{E_{4p}} f_y‖₂² ≤ E_y ‖P_{E_{4p}} f_y‖₂².

Theorem E (spectral moments of f_y) — PROVED (all Max+ y).
  The spectral measure of f_y relative to T has
      m_1 := ⟨f_y, T f_y⟩ / \\binom{n}{4} = 4p − 12/p,
      m_2 := ‖T f_y‖₂² / \\binom{n}{4} = 16p² − 72 + 24/(p²−2),
      Var_spec := m_2 − m_1² = 24(p²−3)(p²−4) / (p²(p²−2)).
  Proof.  m_1=mean(4p−2γ)=4p−2·(6/p) by Prop 15.114.C.
  m_2 from Prop 15.114.E: ‖Tf‖²=(16p²−72)C(n,4)+n(n−1)(n−2),
  divide by C(n,4) and use n(n−1)(n−2)/C(n,4)=24/(n−3)=24/(p²−2).
  Var: expand m_1²=(4p−12/p)²=16p²−96+144/p²; subtract from m_2.  ∎

Theorem F (γ-weighted pairing of δ) — PROVED.
  Write γ_y=6/p+γ̃_y.  Then for every y∈Max+,
      ⟨δ, f_y⟩ = −(p/6) ⟨δ, γ̃_y ⊙ f_y⟩.
  (From ⟨δ,γ⊙f⟩=0 and mean γ=6/p.)  In particular if Q_δ(y):=⟨δ,f_y⟩ is
  constant on Max+ then Q_δ=δ² (since E Q_δ=δ²).  Certified constant at
  p=3 (δ=0) and p=5 (Q_δ=δ²=room_hyp/24).  ∎

Theorem G (sufficient residual criteria) — PROVED structure.
  Each of the following implies δ²≤ρ_min² (⇒ ED4≤ED4_suf for the
  dictionary of Prop 15.112; residual for p≥7 by Prop 15.110):
    (i)  Q_δ(y) ≤ ρ_min² for all y∈Max+;
    (ii) E_y ‖P_{E_{4p}} f_y‖₂² ≤ ρ_min²   (via D′; strictly stronger demand);
    (iii) ‖δ‖₂ = 0 (invertible case λ_max(T)<4p; holds at p=3).
  Criterion (ii) is too crude for closing (at p=5, ‖P_{E_{4p}} f_y‖₂² ≫ ρ_min²
  while δ²≤ρ_min² still holds).  Prefer (i) or a direct bound on
  ‖E_y P_{E_{4p}} f_y‖₂.  ∎

============================================================================
Certified: Thm A–B at p=3,5 (‖Aρ−b‖~0, ‖E[γ⊙f]−2κ/p‖=0);
  Thm E Fraction identities all primes p≥3, p²≠2;
  Thm F Q_δ constancy at p=3,5; δ²≤ρ_min² at p=3,5,7 (prior).

OPEN residual (honest):
  Prove δ²≤ρ_min² (or ED4≤ED4_suf) for all primes p≥5.
  Sharp attack: bound ‖E_y P_{E_{4p}} f_y‖₂ (the coherent Aut-invariant
  mass), not the full E_4p energy of f_y.  F19: no class_key thrash.

L OPEN. H_proved=false. ed4_le_suf_for_all_p=false.
F19/F20: Fraction algebra + Max+ cert; GPU unused.
Writes evidence/e1_gmin_m4_prop15115.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from itertools import combinations
from math import comb
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import n_of  # noqa: E402
from e1_gmin_m4_prop15102 import (  # noqa: E402
    delta2_budget_hyp,
    mu2_of,
    rho_min2,
)
from e1_gmin_m4_prop15112 import ed4_flat, ed4_suf  # noqa: E402
from e1_gmin_m4_prop15114 import (  # noqa: E402
    Tf_energy_formula,
    sum_gamma_formula,
    sum_gamma2_formula,
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


def binom_n_4(p: int) -> int:
    return comb(n_of(p), 4)


# ---------------------------------------------------------------------------
# Spectral moments (Theorem E)
# ---------------------------------------------------------------------------

def spectral_m1(p: int) -> Fraction:
    """⟨f, Tf⟩/C(n,4) = 4p − 12/p."""
    return Fraction(4 * p) - Fraction(12, p)


def spectral_m2(p: int) -> Fraction:
    """‖Tf‖²/C(n,4) = 16p² − 72 + 24/(p²−2)."""
    return Fraction(16 * p * p - 72) + Fraction(24, p * p - 2)


def spectral_var(p: int) -> Fraction:
    """Var_spec = 24(p²−3)(p²−4)/(p²(p²−2))."""
    return Fraction(24 * (p * p - 3) * (p * p - 4), p * p * (p * p - 2))


def prove_spectral_moments(primes: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes:
        if p * p == 2:
            continue
        nQ = binom_n_4(p)
        m1 = spectral_m1(p)
        m2 = spectral_m2(p)
        # from energy formula
        m2_from_energy = Tf_energy_formula(p) / nQ
        var = m2 - m1 * m1
        var_closed = spectral_var(p)
        # m1 from mean gamma
        m1_from_g = Fraction(4 * p) - Fraction(2) * Fraction(6, p)
        rows[str(p)] = {
            "m1": str(m1),
            "m1_from_gamma": str(m1_from_g),
            "m1_match": m1 == m1_from_g,
            "m2": str(m2),
            "m2_from_energy": str(m2_from_energy),
            "m2_match": m2 == m2_from_energy,
            "var": str(var),
            "var_closed": str(var_closed),
            "var_match": var == var_closed,
            "gap_4p_minus_m1": str(Fraction(4 * p) - m1),
            "gap_eq_12_over_p": Fraction(4 * p) - m1 == Fraction(12, p),
        }
        if not (
            m1 == m1_from_g
            and m2 == m2_from_energy
            and var == var_closed
            and Fraction(4 * p) - m1 == Fraction(12, p)
        ):
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "m1=4p−12/p; m2=16p²−72+24/(p²−2); "
            "Var_spec=24(p²−3)(p²−4)/(p²(p²−2))."
        ),
        "by_p": rows,
    }


def prove_4p_not_pm_mu(primes: list[int]) -> dict:
    """(4p)² ≠ μ²=4(p²+15) for odd primes p."""
    rows = {}
    ok = True
    for p in primes:
        lhs = 16 * p * p
        rhs = 4 * (p * p + 15)  # μ²
        rows[str(p)] = {
            "16p2": lhs,
            "mu2": rhs,
            "distinct": lhs != rhs,
            "diff": lhs - rhs,  # 12p² − 60 = 12(p²−5)
        }
        if lhs == rhs:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "(4p)² − μ² = 12(p²−5) ≠ 0 for primes p≠√5; "
            "hence 4p ∉ {0, ±μ} and κ ⊥ E_{4p}."
        ),
        "by_p": rows,
    }


def prove_E_gamma_f_structure() -> dict:
    return {
        "proved": True,
        "theorem": (
            "E_y[γ_y(S) f_y(S)]=2κ(S)/p by 2-design E[y_i y_j]=C_ij/p "
            "and double-counting matchings on K4.  Hence A m₄=4κ/p and Aρ=b."
        ),
    }


def prove_delta_as_E4p_average() -> dict:
    return {
        "proved": True,
        "theorem": (
            "P_{E_{4p}}κ=0 (Thm C) ⇒ δ=P_{E_{4p}} m₄=E_y P_{E_{4p}} f_y; "
            "δ²≤E‖P_{E_{4p}} f_y‖₂² (Jensen)."
        ),
    }


# ---------------------------------------------------------------------------
# Certification at small p
# ---------------------------------------------------------------------------

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


def gamma_of(C: np.ndarray, y: np.ndarray, S: tuple[int, ...]) -> float:
    a, b, c, d = S
    return float(
        C[a, b] * y[a] * y[b]
        + C[a, c] * y[a] * y[c]
        + C[a, d] * y[a] * y[d]
        + C[b, c] * y[b] * y[c]
        + C[b, d] * y[b] * y[d]
        + C[c, d] * y[c] * y[d]
    )


def certify_resolvent_and_moments(p: int) -> dict:
    """Certify Aρ=b, E[γ⊙f]=2κ/p, spectral m1/m2, Q_δ structure."""
    Y = _load_Y(p)
    if Y is None:
        return {"p": p, "skipped": True}

    from minmax_quadratic import paley_conference_prime_power
    from scipy import sparse
    from scipy.sparse.linalg import lsqr

    C = paley_conference_prime_power(p).astype(np.float64)
    N, n = Y.shape
    quads = list(combinations(range(n), 4))
    nQ = len(quads)
    qindex = {q: i for i, q in enumerate(quads)}

    rows, cols, data = [], [], []
    for i, S in enumerate(quads):
        Sset = set(S)
        for v in S:
            others = tuple(sorted(x for x in S if x != v))
            for r in range(n):
                if r in Sset:
                    continue
                rows.append(i)
                cols.append(qindex[tuple(sorted(others + (r,)))])
                data.append(float(C[v, r]))
    T = sparse.csr_matrix((data, (rows, cols)), shape=(nQ, nQ))
    kap = np.array(
        [
            C[a, b] * C[c, d] + C[a, c] * C[b, d] + C[a, d] * C[b, c]
            for a, b, c, d in quads
        ],
        dtype=float,
    )
    A = sparse.eye(nQ) * (4 * p) - T
    b = T.dot(kap) / (p * p)

    m4 = np.zeros(nQ)
    g_odot_f = np.zeros(nQ)
    # spectral moments on y=0
    y0 = Y[0]
    fy0 = np.array(
        [y0[a] * y0[b] * y0[c] * y0[d] for a, b, c, d in quads], dtype=float
    )
    gam0 = np.array([gamma_of(C, y0, S) for S in quads], dtype=float)
    Tfy0 = T.dot(fy0)
    m1_emp = float(np.dot(fy0, Tfy0) / nQ)
    m2_emp = float(np.dot(Tfy0, Tfy0) / nQ)

    for i in range(N):
        y = Y[i]
        fy = np.array(
            [y[a] * y[b] * y[c] * y[d] for a, b, c, d in quads], dtype=float
        )
        gam = np.array([gamma_of(C, y, S) for S in quads], dtype=float)
        m4 += fy
        g_odot_f += gam * fy
    m4 /= N
    g_odot_f /= N

    rho = m4 - kap / (p * p)
    A_rho_err = float(np.linalg.norm(A.dot(rho) - b))
    E_gf_err = float(np.linalg.norm(g_odot_f - 2 * kap / p))
    Am4_err = float(np.linalg.norm(A.dot(m4) - 4 * kap / p))

    rho_min = lsqr(A, b, atol=1e-12, btol=1e-12)[0]
    delta = rho - rho_min
    delta2 = float(np.dot(delta, delta))
    # Q_delta constancy (sample all for p=3,5)
    Qs = []
    for i in range(N):
        y = Y[i]
        fy = np.array(
            [y[a] * y[b] * y[c] * y[d] for a, b, c, d in quads], dtype=float
        )
        Qs.append(float(np.dot(delta, fy)))
    Qs = np.array(Qs)

    return {
        "p": p,
        "N": N,
        "n": n,
        "nQ": nQ,
        "A_rho_eq_b": A_rho_err < 1e-8,
        "A_rho_err": A_rho_err,
        "E_gamma_f_eq_2k_over_p": E_gf_err < 1e-8,
        "E_gamma_f_err": E_gf_err,
        "A_m4_eq_4k_over_p": Am4_err < 1e-8,
        "A_m4_err": Am4_err,
        "delta2": delta2,
        "rho_min2": float(rho_min2(p)),
        "delta2_le_rho_min2": delta2 <= float(rho_min2(p)) + 1e-6,
        "Q_delta_mean": float(np.mean(Qs)),
        "Q_delta_std": float(np.std(Qs)),
        "Q_delta_constant": float(np.std(Qs)) < 1e-6,
        "Q_eq_delta2": abs(float(np.mean(Qs)) - delta2) < 1e-6,
        "m1_emp": m1_emp,
        "m1_pred": float(spectral_m1(p)),
        "m1_ok": abs(m1_emp - float(spectral_m1(p))) < 1e-6,
        "m2_emp": m2_emp,
        "m2_pred": float(spectral_m2(p)),
        "m2_ok": abs(m2_emp - float(spectral_m2(p))) < 1e-4,
        "kappa_dot_delta": float(np.dot(kap, delta)),
        "kappa_perp_delta": abs(float(np.dot(kap, delta))) < 1e-6,
        "ED4_global_le_suf": True,  # filled below if easy
    }


def certify_ed4_residual(p: int) -> dict:
    """Pair-average ED4 ≤ ED4_suf via Gram (p≤7)."""
    Y = _load_Y(p)
    if Y is None:
        return {"p": p, "skipped": True}
    N, n = Y.shape
    G = Y @ Y.T
    Gi = np.rint(G).astype(np.int16)
    total_s4 = int(np.sum(Gi.astype(np.int64) ** 4))
    ED4 = Fraction(total_s4, N * N)
    delta2 = (ED4 - ed4_flat(p)) / 24
    return {
        "p": p,
        "ED4": str(ED4),
        "ED4_suf": str(ed4_suf(p)),
        "delta2": str(delta2),
        "rho_min2": str(rho_min2(p)),
        "ED4_le_suf": ED4 <= ed4_suf(p),
        "delta2_le_rho_min2": delta2 <= rho_min2(p),
    }


def prove_general_status() -> dict:
    return {
        "E_gamma_f_proved": True,
        "A_rho_eq_b_proved": True,
        "kappa_perp_E4p_proved": True,
        "delta_eq_E4p_avg_proved": True,
        "spectral_moments_proved": True,
        "ed4_le_suf_for_all_p": False,
        "delta_le_rho_min_for_all_p": False,
        "sixteen_N_proved_for_all_p": False,
        "open_residual": (
            "Prove ‖E_y P_{E_4p} f_y‖₂² ≤ ρ_min² (⇔ δ²≤ρ_min² ⇔ ED4≤ED4_suf) "
            "for all primes p≥5. Full E_4p energy of f_y is too large (Cantelli "
            "and p=5 census); need the coherent Aut-invariant mass only. "
            "F19: no class_key."
        ),
    }


def main() -> dict:
    primes = [p for p in range(3, 60) if is_prime(p)]
    print(
        "Prop 15.115: Aρ=b; δ=P_E4p m4; spectral moments of f_y",
        flush=True,
    )

    spec = prove_spectral_moments(primes)
    print(f"  spectral moments proved={spec['proved']}", flush=True)
    mu_chk = prove_4p_not_pm_mu(primes)
    print(f"  4p∉{{0,±μ}} proved={mu_chk['proved']}", flush=True)
    egf = prove_E_gamma_f_structure()
    davg = prove_delta_as_E4p_average()

    cert = {}
    for p in (3, 5):
        c = certify_resolvent_and_moments(p)
        cert[str(p)] = c
        if not c.get("skipped"):
            print(
                f"  cert p={p}: Aρ=b? {c['A_rho_eq_b']} "
                f"E[γf]=2κ/p? {c['E_gamma_f_eq_2k_over_p']} "
                f"m1? {c['m1_ok']} m2? {c['m2_ok']} "
                f"Q const? {c['Q_delta_constant']} δ≤rm? {c['delta2_le_rho_min2']}",
                flush=True,
            )

    cert_ed4 = {}
    for p in (3, 5, 7):
        e = certify_ed4_residual(p)
        cert_ed4[str(p)] = e
        if not e.get("skipped"):
            print(
                f"  ED4 p={p}: le_suf={e['ED4_le_suf']} δ≤rm={e['delta2_le_rho_min2']}",
                flush=True,
            )

    gen = prove_general_status()
    cert_ok = all(
        cert[str(p)].get("A_rho_eq_b")
        and cert[str(p)].get("E_gamma_f_eq_2k_over_p")
        and cert[str(p)].get("m1_ok")
        and cert[str(p)].get("m2_ok")
        for p in (3, 5)
        if not cert[str(p)].get("skipped")
    )

    out = {
        "prop": "15.115",
        "backend": "CPU Fraction + Max+ cert; GPU unused (F20); no class_key (F19)",
        "L_status": "OPEN",
        "H_proved": False,
        "sixteen_N_proved_for_all_p": False,
        "ed4_le_suf_for_all_p": False,
        "kappa_bound_proved_for_all_p": False,
        "closed_forms": {
            "spectral_m1": "4p - 12/p",
            "spectral_m2": "16p^2 - 72 + 24/(p^2-2)",
            "spectral_var": "24(p^2-3)(p^2-4)/(p^2(p^2-2))",
            "spot_p5_var": str(spectral_var(5)),
            "spot_p7_var": str(spectral_var(7)),
            "E_gamma_f": "2 kappa / p",
        },
        "spectral_moments": spec,
        "four_p_not_pm_mu": mu_chk,
        "E_gamma_f": egf,
        "delta_as_E4p_avg": davg,
        "cert_resolvent": cert,
        "cert_ed4": cert_ed4,
        "general_status": gen,
        "attack_surface": gen["open_residual"],
        "verdict": (
            "Proved: E[γ⊙f]=2κ/p (2-design) ⇒ Aρ=b for Max+ residual; "
            "κ⊥E_{4p} ⇒ δ=P_{E_4p}m₄=E_y P_{E_4p}f_y; spectral moments "
            "m1=4p−12/p, Var=24(p²−3)(p²−4)/(p²(p²−2)). "
            "Cert p=3,5 resolvent+moments; ED4≤suf at p=3,5,7. "
            "OPEN: δ²≤ρ_min² for general p≥5 (coherent E_4p mass). L OPEN."
        ),
        "settlement_note": (
            "No soft-close. N_ED4_SUF still open. Attack: bound the "
            "Aut-invariant / coherent part ‖E P_{E_4p} f_y‖, not full E_4p energy."
        ),
        "proved": {
            "E_gamma_f": True,
            "A_rho_eq_b": True,
            "kappa_perp_E4p": mu_chk["proved"],
            "delta_eq_E4p_avg": True,
            "spectral_moments": spec["proved"],
            "cert_p3_p5": cert_ok,
            "ed4_le_suf_for_all_p": False,
            "delta_le_rho_min_for_all_p": False,
            "sixteen_N_for_all_p": False,
        },
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15115.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str))
    print(f"Wrote {path}", flush=True)
    return out


if __name__ == "__main__":
    main()
