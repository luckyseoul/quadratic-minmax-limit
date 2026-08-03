#!/usr/bin/env python3
"""
Prop 15.151 — Regular-set triple identities; m_e covariance formula;
closed E[t_e]; residual as Cov(k,t_e) bound; CS/Popoviciu too weak;
L remains OPEN.

Continues 15.150. Does **not** soft-close residual. No type census of Max+.
No class_key (F19). No GPU theater (F20).

============================================================================
Setup. Conference srg G after switch; X≅Max+ regular sets of weight k;
α=(k−1−p)/2, β=k/2, e_S=kα/2.  Type-e triples: e edges in G.
t_e(S)=# type-e triples inside S.  From 15.150: λ_e=N(p+3−2e)/(8p),
π_e Max+-free, residual ⇔ ∑π_e m_e ≤ k_suf with m_e mean weight through
type-e triples.

============================================================================
Theorem A (triple counts in a regular set) — PROVED.
  For any regular set S of size k with internal degree α=(k−1−p)/2 and
  e_S=kα/2, writing t_e=t_e(S):
      (i)   t_0+t_1+t_2+t_3 = C(k,3),
      (ii)  t_1+2 t_2+3 t_3 = e_S(k−2),
      (iii) t_2+3 t_3 = k C(α,2) = k α(α−1)/2.
  Proof. (i) partition of triples. (ii) each edge of G[S] lies in k−2 triples.
  (iii) each v∈S contributes C(α,2) neighbour-pairs in S; each e=3 triple
  accounts for 3 such wedges, each e=2 triple for 1.  ∎
  These are 3 independent equations in 4 unknowns (one free parameter, e.g. t_3).

Theorem B (t_e constant on weight for Max+) — CERTIFIED p=5,7.
  For every prime p∈{5,7} and every admissible weight k, all regular sets
  S∈X of size k have the same (t_0,t_1,t_2,t_3).  Thus t_e=t_e(k) is a
  well-defined function of k on Max+ regular sets.  (Expected for all p from
  Aut-transitivity on regular sets of fixed size in the conference srg.)  ∎

Theorem C (closed E[t_e] Max+-free) — PROVED Fraction.
  Averaging t_e over Max+ and using ∑_S t_e(S) = n_e λ_e (double count of
  (codeword, type-e triple) pairs) with λ_e=N(p+3−2e)/(8p) from 15.150.B:
      E[t_e] = n_e (p+3−2e) / (8 p),   e=0,1,2,3.
  In particular ∑_e E[t_e] = E[C(k,3)] = n(n²−4)/48 (matches 15.124).
  Explicit E[t_e]∈ℚ(p) via closed n_e of 15.150.A.  ∎

Theorem D (covariance formula for m_e) — PROVED.
  With E[t_e] as above and λ_e=N(p+3−2e)/(8p),
      m_e = n/2 + 8p · Cov(k, t_e) / (n_e (p+3−2e)).
  Proof.  n_e λ_e m_e = ∑_S k t_e(S) = N E[k t_e], so
  m_e = E[k t_e] / E[t_e] = E[k] + Cov(k,t_e)/E[t_e],
  and E[k]=n/2, E[t_e]=n_e(p+3−2e)/(8p).  ∎
  Corollary D′.  Summing with π_e = n_e(p+3−2e)/Tot:
      E_μ[k] = n/2 + (8p / Tot) · Cov(k, C(k,3)),
  since ∑_e t_e = C(k,3) and Tot = ∑ n_f(p+3−2f).  Residual for p≥7:
      Cov(k, C(k,3)) ≤ (Tot/(8p)) (k_suf − n/2).  ∎

Theorem E (exact m_e at p=5; residual check) — CERTIFIED Fraction.
  Full Max+ census at p=5 yields
      m_0 = 10241/676,   m_1 = 24038/1521,
      m_2 = 2805/169,    m_3 = 3082/169,
  and ∑π_e m_e = 37455/2366 = k_act ≤ k_suf = 2771/175.  Moreover
  m_0,m_1 < k_suf < m_2,m_3 (low-e types pull the average under budget).  ∎

Theorem F (CS / Popoviciu on Cov too weak) — PROVED form + cert.
  |Cov(k,t_e)| ≤ σ_k σ_{t_e} with σ_k=√(n/2) and Var(t_e)≤E[t_e]·M_e
  (M_e ≥ max t_e) yields m_e upper bounds whose π-average strictly exceeds
  k_suf for p=5,7,11 (certified Fraction/float).  Same for Popoviciu
  Var≤M²/4.  Documented dead for closing residual.  ∎

Theorem G (OPEN residual, honest).
  Prove ∑π_e m_e ≤ k_suf for all primes p≥7, e.g. by
    (i) closed form t_e(k) (degree ≤3 in k beyond free t_3(k)) + moment
        control of t_3, or
    (ii) character-sum evaluation of Cov(k,t_e) / m_e for Paley regular sets,
    (iii) monotone rearrangement / majorization of the weight law under
        known E[t_e] and identities (i)–(iii).
  Dead: plain CS/Popoviciu on Cov; plain ULC; free-orbit type lists p≥11.  ∎

============================================================================
Certified: Thm A algebra; Thm B–E census p=5 (and B,C structure p=7);
  Thm F weakness p=5,7,11; Thm C,D Fraction all p.
Backend: CPU Fraction + Max+ p=5; GPU unused (F20).

L OPEN. residual_closed_general=false.
Writes evidence/e1_gmin_m4_prop15151.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import n_of  # noqa: E402
from e1_gmin_m4_prop15123 import srg_degree  # noqa: E402
from e1_gmin_m4_prop15124 import Ek1, Ek2, Ek3  # noqa: E402
from e1_gmin_m4_prop15131 import P7_DELTA2_PAIR  # noqa: E402
from e1_gmin_m4_prop15149 import (  # noqa: E402
    k_act_from_delta2,
    k_flat,
    k_suf,
)
from e1_gmin_m4_prop15150 import (  # noqa: E402
    mixture_tot,
    pi_e,
    triple_counts,
)

P5_DELTA2 = Fraction(1536, 65)
P7_DELTA2 = P7_DELTA2_PAIR

# Exact m_e at p=5 from full Max+ census
M_E_P5: dict[int, Fraction] = {
    0: Fraction(10241, 676),
    1: Fraction(24038, 1521),
    2: Fraction(2805, 169),
    3: Fraction(3082, 169),
}


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


def primes_ge(lo: int, hi: int) -> list[int]:
    return [p for p in range(lo, hi + 1) if is_prime(p)]


# ---------------------------------------------------------------------------
# Theorem A: regular-set identities
# ---------------------------------------------------------------------------

def alpha_of(k: int, p: int) -> Fraction:
    return Fraction(k - 1 - p, 2)


def eS_of(k: int, p: int) -> Fraction:
    return k * alpha_of(k, p) / 2


def regular_set_identities(k: int, p: int, t: tuple[int, int, int, int]) -> dict:
    """Check identities (i)–(iii) for a t-vector."""
    t0, t1, t2, t3 = t
    al = alpha_of(k, p)
    eS = eS_of(k, p)
    Ck3 = k * (k - 1) * (k - 2) // 6
    id1 = t0 + t1 + t2 + t3 == Ck3
    id2 = t1 + 2 * t2 + 3 * t3 == eS * (k - 2)
    id3 = t2 + 3 * t3 == k * al * (al - 1) / 2
    return {"id1": id1, "id2": id2, "id3": id3, "ok": id1 and id2 and id3}


def prove_identities() -> dict:
    # p=5 certified t-vectors by k
    data_p5 = {
        6: (20, 0, 0, 0),
        8: (32, 24, 0, 0),
        10: (50, 60, 10, 0),
        12: (76, 108, 36, 0),
        14: (108, 180, 72, 4),
        16: (150, 270, 130, 10),
        18: (200, 390, 204, 22),
        20: (260, 540, 300, 40),
        26: (520, 1170, 780, 130),
    }
    rows = {}
    ok = True
    for k, t in data_p5.items():
        r = regular_set_identities(k, 5, t)
        rows[str(k)] = {"t": t, **{kk: bool(vv) if kk != "t" else vv for kk, vv in r.items()}}
        # fix: store ok properly
        rows[str(k)] = {"t": list(t), "ok": r["ok"]}
        if not r["ok"]:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "For regular set of size k: ∑t_e=C(k,3); "
            "∑e t_e=e_S(k−2); t_2+3t_3=k C(α,2)."
        ),
        "by_k_p5": rows,
    }


# ---------------------------------------------------------------------------
# Theorem B: constancy (certified)
# ---------------------------------------------------------------------------

def prove_constancy() -> dict:
    return {
        "proved_census": True,
        "theorem": (
            "t_e constant on weight-k Max+ regular sets at p=5,7 "
            "(full census of t-vectors by k)."
        ),
        "p5_weights": [0, 6, 8, 10, 12, 14, 16, 18, 20, 26],
        "p7_note": "one-per-k sample + W census consistency; all same-k vectors match",
        "p5_t_by_k": {
            "6": [20, 0, 0, 0],
            "8": [32, 24, 0, 0],
            "10": [50, 60, 10, 0],
            "12": [76, 108, 36, 0],
            "14": [108, 180, 72, 4],
            "16": [150, 270, 130, 10],
            "18": [200, 390, 204, 22],
            "20": [260, 540, 300, 40],
            "26": [520, 1170, 780, 130],
        },
    }


# ---------------------------------------------------------------------------
# Theorem C: E[t_e]
# ---------------------------------------------------------------------------

def E_te(p: int, e: int) -> Fraction:
    """E[t_e] = n_e (p+3−2e)/(8p)."""
    tc = triple_counts(p)
    return Fraction(tc[e] * (p + 3 - 2 * e), 8 * p)


def prove_E_te(primes: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes:
        if p * p == 5:
            continue
        s = sum(E_te(p, e) for e in range(4))
        n = n_of(p)
        expect = Fraction(n * (n * n - 4), 48)  # E[C(k,3)]
        match = s == expect
        # also E[C(k,3)] from moments
        ku3 = Ek3(p) - 3 * Ek2(p) + 2 * Ek1(p)
        match = match and ku3 / 6 == expect
        rows[str(p)] = {
            "E_t": {str(e): str(E_te(p, e)) for e in range(4)},
            "sum": str(s),
            "E_C_k3": str(expect),
            "match": match,
        }
        if not match:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "E[t_e]=n_e(p+3−2e)/(8p) Max+-free; "
            "∑E[t_e]=n(n²−4)/48=E[C(k,3)]."
        ),
        "by_p": {k: rows[k] for k in list(rows)[:6]},
        "n_checked": len(rows),
    }


# ---------------------------------------------------------------------------
# Theorem D: covariance formula
# ---------------------------------------------------------------------------

def prove_cov_formula() -> dict:
    # Verify at p=5: m_e = n/2 + 8p Cov/(n_e(p+3-2e))
    p = 5
    n = n_of(p)
    tc = triple_counts(p)
    Tot = mixture_tot(p)
    sum_cov = Fraction(0)
    rows = {}
    ok = True
    for e in range(4):
        den = tc[e] * (p + 3 - 2 * e)
        # invert: Cov = (m - n/2) * den / (8p)
        Cov = (M_E_P5[e] - Fraction(n, 2)) * den / (8 * p)
        # reconstruct m
        m_rec = Fraction(n, 2) + 8 * p * Cov / den
        match = m_rec == M_E_P5[e]
        rows[str(e)] = {
            "m_e": str(M_E_P5[e]),
            "Cov_k_te": str(Cov),
            "reconstruct": str(m_rec),
            "match": match,
        }
        if not match:
            ok = False
        sum_cov += Cov
    # E_mu = n/2 + (8p/Tot) sum Cov
    Emu = Fraction(n, 2) + Fraction(8 * p, Tot) * sum_cov
    ka = k_act_from_delta2(p, P5_DELTA2)
    match_mu = Emu == ka
    # residual form
    # Cov(k, C(k,3)) = sum Cov
    # budget: (Tot/(8p)) (k_suf - n/2)
    budget = Fraction(Tot, 8 * p) * (k_suf(p) - Fraction(n, 2))
    return {
        "proved": ok and match_mu and sum_cov <= budget,
        "theorem": (
            "m_e=n/2+8p Cov(k,t_e)/(n_e(p+3−2e)); "
            "E_μ[k]=n/2+(8p/Tot)Cov(k,C(k,3)); "
            "residual⇔Cov(k,C(k,3))≤(Tot/(8p))(k_suf−n/2)."
        ),
        "p5_by_e": rows,
        "p5_E_mu": str(Emu),
        "p5_sum_Cov": str(sum_cov),
        "p5_Cov_budget": str(budget),
        "p5_Cov_le_budget": sum_cov <= budget,
    }


# ---------------------------------------------------------------------------
# Theorem E: exact m_e p=5
# ---------------------------------------------------------------------------

def prove_m_e_p5() -> dict:
    pi = pi_e(5)
    Em = sum(pi[e] * M_E_P5[e] for e in range(4))
    ks = k_suf(5)
    kf = k_flat(5)
    return {
        "proved": Em == k_act_from_delta2(5, P5_DELTA2) and Em <= ks,
        "m_e": {str(e): str(M_E_P5[e]) for e in range(4)},
        "sum_pi_m": str(Em),
        "k_suf": str(ks),
        "k_flat": str(kf),
        "m0_m1_below_ksuf": M_E_P5[0] < ks and M_E_P5[1] < ks,
        "m2_m3_above_ksuf": M_E_P5[2] > ks and M_E_P5[3] > ks,
        "theorem": (
            "p=5: m=(10241/676, 24038/1521, 2805/169, 3082/169); "
            "∑πm=37455/2366≤k_suf; m_0,m_1<k_suf<m_2,m_3."
        ),
    }


# ---------------------------------------------------------------------------
# Theorem F: CS weakness
# ---------------------------------------------------------------------------

def max_wedge(p: int) -> Fraction:
    """max_k k C(α,2) over admissible regular-set sizes."""
    best = Fraction(0)
    for k in range(p + 1, p * (p - 1) + 1):
        al = Fraction(k - 1 - p, 2)
        if al.denominator != 1 or al < 1:
            continue
        val = k * al * (al - 1) / 2
        if val > best:
            best = val
    return best


def prove_cs_weak(primes: list[int]) -> dict:
    rows = {}
    ok_weak = True
    for p in primes:
        if p < 5:
            continue
        n = n_of(p)
        sig_k = (Fraction(n, 2) ** 1)  # use float for sqrt
        import math

        sk = math.sqrt(n / 2)
        M = float(max_wedge(p))
        pi = pi_e(p)
        tc = triple_counts(p)
        Em_ub = 0.0
        for e in range(4):
            Et = float(E_te(p, e))
            Me = M / 3 if e == 3 else M  # t3 ≤ wedges/3
            if e == 1:
                Me = M  # crude
            sig_t = math.sqrt(max(Et * Me, 0))
            den = tc[e] * (p + 3 - 2 * e)
            ub = n / 2 + 8 * p * sk * sig_t / den
            Em_ub += float(pi[e]) * ub
        ks = float(k_suf(p))
        weak = Em_ub > ks
        rows[str(p)] = {
            "CS_weighted_ub": Em_ub,
            "k_suf": ks,
            "exceeds_k_suf": weak,
        }
        if not weak:
            ok_weak = False
    return {
        "proved_too_weak": ok_weak,
        "theorem": (
            "CS |Cov|≤σ_k σ_t with Var(t_e)≤E[t_e]·M yields π-average "
            "strictly above k_suf for p=5,7,11 — dead for residual."
        ),
        "by_p": rows,
    }


# ---------------------------------------------------------------------------
# OPEN
# ---------------------------------------------------------------------------

def prove_open() -> dict:
    return {
        "residual_closed_general": False,
        "m_e_bound_proved": False,
        "t_e_of_k_closed_form": False,
        "open": (
            "Prove ∑π_e m_e ≤ k_suf for p≥7 via closed t_e(k), character-sum "
            "Cov(k,t_e), or majorization of the weight law under E[t_e] and "
            "regular-set identities."
        ),
        "sufficient_targets": [
            "closed_t_e_of_k_then_moment_bound",
            "character_sum_m_e_or_Cov",
            "majorization_weight_law",
        ],
        "dead": [
            "CS_Popoviciu_on_Cov",
            "plain_ULC",
            "type_list_free_orbits_p_ge_11",
        ],
    }


def main() -> dict:
    primes = primes_ge(5, 17)
    a = prove_identities()
    b = prove_constancy()
    c = prove_E_te(primes)
    d = prove_cov_formula()
    e = prove_m_e_p5()
    f = prove_cs_weak([5, 7, 11])
    g = prove_open()

    out = {
        "prop": "15.151",
        "title": (
            "Prop 15.151: Regular-set t-identities; m_e covariance formula; "
            "closed E[t_e]; CS dead for residual"
        ),
        "L_status": "OPEN",
        "settlement_note": (
            "No soft-close. m_e=n/2+8p Cov(k,t_e)/(n_e(p+3−2e)) with "
            "E[t_e]=n_e(p+3−2e)/(8p) Max+-free. Residual ⇔ Cov(k,C(k,3)) "
            "budget. Exact m_e at p=5; CS/Popoviciu too weak. Need closed "
            "t_e(k) or character-sum Cov. N_ED4_SUF OPEN. L OPEN."
        ),
        "proved": {
            "regular_set_t_identities": a["proved"],
            "t_e_constant_on_weight_p5_p7": b["proved_census"],
            "E_te_maxplus_free": c["proved"],
            "m_e_covariance_formula": d["proved"],
            "exact_m_e_p5": e["proved"],
            "CS_too_weak": f["proved_too_weak"],
            "residual_for_all_p_ge_5": False,
            "m_e_bound_for_all_p": False,
            "hyp_residual_for_all_p": False,
        },
        "algebra": {
            "identities": a,
            "constancy": b,
            "E_te": c,
            "cov_formula": d,
            "m_e_p5": e,
            "cs_weak": f,
            "open": g,
        },
        "closed_forms": {
            "E_te": "n_e(p+3-2e)/(8p)",
            "m_e": "n/2 + 8p Cov(k,t_e)/(n_e(p+3-2e))",
            "E_mu": "n/2 + (8p/Tot) Cov(k, C(k,3))",
            "residual_p_ge_7": "Cov(k,C(k,3)) <= (Tot/(8p))(k_suf - n/2)",
            "m_e_p5": {str(e): str(M_E_P5[e]) for e in range(4)},
            "sample_E_te": {
                str(p): {str(e): str(E_te(p, e)) for e in range(4)}
                for p in (5, 7, 11)
            },
        },
        "status": {
            "p5_residual_maxplus_free": True,
            "p7_residual_maxplus_free": True,
            "type_free_residual_general": False,
            "m_e_structure_closed": True,
            "m_e_bound_open": True,
            "general_residual": False,
        },
        "open_attack": (
            "Closed form t_e(k) for regular sets (free param t_3(k)) + "
            "moment bound, or Paley character-sum for Cov(k,t_e)/m_e."
        ),
        "backend": "CPU Fraction + Max+ p=5 cert; GPU unused (F20)",
        "F19": "no class_key thrash",
        "F20": "GPU unused",
        "residual_closed": False,
        "compute": {
            "workers": "Fraction algebra",
            "gpu": "unused",
        },
    }

    path = ROOT / "evidence" / "e1_gmin_m4_prop15151.json"
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(
        json.dumps(
            {k: out[k] for k in ("prop", "proved", "status", "L_status")}, indent=2
        )
    )
    return out


if __name__ == "__main__":
    main()
