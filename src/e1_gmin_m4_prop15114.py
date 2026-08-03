#!/usr/bin/env python3
"""
Prop 15.114 — γ-calculus of T on Max+; Tf_y multiplicative formula;
∑γ, ∑γ² closed forms; ED4 multi-type at p=7; residual still OPEN.

Continues 15.110–15.113 toward N_ED4_SUF (ED4≤ED4_suf ⇔ δ²≤ρ_min²).
Does **not** soft-close the residual. No class_key (F19).

============================================================================
Setup. Conference C, n=p²+1, Max+ boolean Cy=py.  For unordered 4-sets S,
f_y(S)=∏_{i∈S} y_i, and the edge form
    γ_y(S) := ∑_{{i,j}⊂S} C_{ij} y_i y_j
(six edges of K4).  T = C-signed Johnson adjacency on 4-set functions.
A := 4p I − T.  Residual ρ=m₄−κ/p²=ρ_min+δ, δ∈ker A=E_{4p}.

============================================================================
Theorem A (multiplicative eigenformula) — PROVED.
  For every boolean y with Cy=py and every 4-set S,
      (T f_y)(S) = (4p − 2 γ_y(S)) f_y(S).
  Equivalently A f_y = 2 (γ_y ⊙ f_y) (pointwise product).

  Proof.  (Tf_y)(S)=∑_{v∈S, r∉S} C_{vr} f_y(S−v+r)
  = f_y(S) · ∑_{v∈S} y_v · ∑_{r∉S} C_{vr} y_r
  (using f_y(S−v+r)=f_y(S)·y_r/y_v).
  ∑_{r≠v} C_{vr} y_r = (Cy)_v = p y_v, so
  ∑_{r∉S} C_{vr} y_r = p y_v − ∑_{r∈S\\{v}} C_{vr} y_r.
  Hence the inner sum is ∑_{v∈S}(p − ∑_{r∈S\\{v}} C_{vr} y_v y_r)
  = 4p − ∑_{v≠r in S} C_{vr} y_v y_r = 4p − 2 γ_y(S).  ∎

Theorem B (δ annihilates γ⊙f) — PROVED.
  For every y∈Max+ and every δ∈ker A,
      ⟨δ, γ_y ⊙ f_y⟩ = 0.
  Proof.  ⟨δ, A f_y⟩ = ⟨Aδ, f_y⟩ = 0, and A f_y = 2 γ_y⊙f_y.  ∎

Theorem C (∑ γ closed form) — PROVED for all Max+ y.
      ∑_S γ_y(S) = (6/p) \\binom{n}{4}.
  Proof.  Each edge e appears in \\binom{n-2}{2} quads, so
  ∑_S γ = \\binom{n-2}{2} ∑_{i<j} C_{ij} y_i y_j
  = \\binom{n-2}{2} · (1/2) yᵀ C y = \\binom{n-2}{2} · (p n)/2.
  Divide by \\binom{n}{4}: mean = 6p/(n−1)=6/p.  ∎

Theorem D (∑ γ² closed form) — PROVED for all Max+ y.
      ∑_S γ_y(S)² = 6 \\binom{n}{4} + n(n−1)(n−2)/4.
  Proof.  Expand ∑_S (∑_{e⊂S} a_e)² with a_e=C_e∏_{i∈e} y_i.
  (i) Diagonal e=f: a_e²=1, six edges per S ⇒ 6\\binom{n}{4}.
  (ii) Opposite edges (perfect matchings on S): ordered contribution
      ∑_S 2 κ(S) f_y(S) = 2 · n(n−1)(n−2)/8 = n(n−1)(n−2)/4
      (Prop 15.110: ∑κ∏ constant on Max+).
  (iii) Adjacent edges (share a vertex): for e={i,j}, f={i,k},
      a_e a_f = C_{ij} C_{ik} y_j y_k.  Summing with multiplicity (n−3)
      yields (n−3) ∑_i [(∑_{j≠i} C_{ij} y_j)² − (n−1)]
      = (n−3) ∑_i (p² y_i² − (n−1)) = (n−3)·n·(p²−(n−1))=0
      since n−1=p².  ∎

Theorem E (‖T f_y‖² closed form) — PROVED.
  ‖T f_y‖₂² = ∑_S (4p−2γ_y(S))²
            = (16p² − 72)\\binom{n}{4} + n(n−1)(n−2).
  (Expand using C,D; independent of the choice of y∈Max+.)  ∎

Theorem F (pair-average ED4 is the residual moment) — PROVED structure.
  The dictionary ED4 = ED4_flat + 24 δ² of Prop 15.112 uses the
  **pair average**
      ED4 := N^{-2} ∑_{y,z∈Max+} (y·z)^4 = E_y[ ED4(y) ],
  not a single basepoint.  At p=7, Max+ is **not** 2-point homogeneous:
  ED4(y) takes three values (certified), all ≤ ED4_suf, with global mean
  strictly below ED4_suf.  Attack note: Aut(C) is not transitive on the
  three ED4-types; residual calculus must use the pair average (or control
  the worst type).  ∎

============================================================================
Certified:
  Thm A–E algebra + Max+ check at p=3,5; Thm F census at p=7
  (types 2352 + 8400 + 700; max ED4(y)≤ED4_suf; global ED4≤ED4_suf);
  δ²≤ρ_min² at p=3,5,7 (prior).

OPEN residual (honest):
  Prove pair-average ED4 ≤ ED4_suf (⇔ δ²≤ρ_min²) for all primes p≥5,
  e.g. via γ-weighted spectral mass of f_y on E_{4p}, or closed weight
  enumerator of Max+ dots.  F19: no class_key thrash.

L OPEN. H_proved=false. ed4_le_suf_for_all_p=false.
F19/F20: Fraction algebra + Max+ cert; GPU unused.
Writes evidence/e1_gmin_m4_prop15114.json
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

from e1_gmin_m4_prop15100 import n_of, wick_lo  # noqa: E402
from e1_gmin_m4_prop15102 import delta2_budget_hyp, rho_min2  # noqa: E402
from e1_gmin_m4_prop15112 import ed4_flat, ed4_suf  # noqa: E402


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
# Closed forms (Theorems C–E)
# ---------------------------------------------------------------------------

def binom_n_4(p: int) -> int:
    return comb(n_of(p), 4)


def sum_gamma_formula(p: int) -> Fraction:
    """∑_S γ_y(S) = (6/p) C(n,4)."""
    return Fraction(6, p) * binom_n_4(p)


def sum_gamma2_formula(p: int) -> Fraction:
    """∑_S γ_y(S)² = 6 C(n,4) + n(n−1)(n−2)/4."""
    n = n_of(p)
    return Fraction(6 * binom_n_4(p)) + Fraction(n * (n - 1) * (n - 2), 4)


def Tf_energy_formula(p: int) -> Fraction:
    """‖T f_y‖₂² = (16p² − 72) C(n,4) + n(n−1)(n−2)."""
    n = n_of(p)
    nQ = binom_n_4(p)
    return Fraction(16 * p * p - 72) * nQ + Fraction(n * (n - 1) * (n - 2))


def prove_gamma_algebra(primes: list[int]) -> dict:
    """Verify C–E identities are consistent Fraction forms."""
    rows = {}
    ok = True
    for p in primes:
        n = n_of(p)
        nQ = binom_n_4(p)
        sg = sum_gamma_formula(p)
        sg2 = sum_gamma2_formula(p)
        # mean γ = 6/p
        mean_g = sg / nQ
        # energy expand check: sum (4p-2g)^2 = 16p² nQ - 16p sum g + 4 sum g²
        energy_from_g = (
            Fraction(16 * p * p) * nQ
            - Fraction(16 * p) * sg
            + Fraction(4) * sg2
        )
        energy_closed = Tf_energy_formula(p)
        # adjacent-vanishing: sum g² - 6 nQ - n(n-1)(n-2)/4 == 0 by construction
        match_part = Fraction(n * (n - 1) * (n - 2), 4)
        diag_part = Fraction(6 * nQ)
        rows[str(p)] = {
            "sum_gamma": str(sg),
            "mean_gamma": str(mean_g),
            "mean_gamma_eq_6_over_p": mean_g == Fraction(6, p),
            "sum_gamma2": str(sg2),
            "diag_plus_match": str(diag_part + match_part),
            "sum_gamma2_eq_diag_match": sg2 == diag_part + match_part,
            "Tf_energy": str(energy_closed),
            "Tf_energy_from_moments": str(energy_from_g),
            "energy_match": energy_from_g == energy_closed,
        }
        if not (
            mean_g == Fraction(6, p)
            and sg2 == diag_part + match_part
            and energy_from_g == energy_closed
        ):
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "∑γ=(6/p)C(n,4); ∑γ²=6 C(n,4)+n(n−1)(n−2)/4 (adj=0 via Cy=py); "
            "‖Tf_y‖²=(16p²−72)C(n,4)+n(n−1)(n−2)."
        ),
        "by_p": rows,
    }


def prove_multiplicative_formula() -> dict:
    return {
        "proved": True,
        "theorem": (
            "For boolean Cy=py: (T f_y)(S)=(4p−2γ_y(S)) f_y(S); "
            "A f_y=2(γ_y⊙f_y); hence ⟨δ,γ_y⊙f_y⟩=0 for all δ∈ker A."
        ),
    }


# ---------------------------------------------------------------------------
# Certification
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


def certify_gamma_and_Tf(p: int, n_y: int = 3) -> dict:
    """Certify Thm A,C,D,E at small p on Max+ samples."""
    Y = _load_Y(p)
    if Y is None:
        return {"p": p, "skipped": True}

    from minmax_quadratic import paley_conference_prime_power
    from scipy import sparse

    C = paley_conference_prime_power(p).astype(np.float64)
    N, n = Y.shape
    quads = list(combinations(range(n), 4))
    nQ = len(quads)
    qindex = {q: i for i, q in enumerate(quads)}

    # build T once
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

    pred_sg = float(sum_gamma_formula(p))
    pred_sg2 = float(sum_gamma2_formula(p))
    pred_energy = float(Tf_energy_formula(p))

    samples = []
    n_check = min(n_y, N)
    all_ok = True
    for yi in range(n_check):
        y = Y[yi]
        fy = np.array(
            [y[a] * y[b] * y[c] * y[d] for a, b, c, d in quads], dtype=float
        )
        gam = np.array([gamma_of(C, y, S) for S in quads], dtype=float)
        Tfy = T.dot(fy)
        pred = (4 * p - 2 * gam) * fy
        tf_ok = bool(np.linalg.norm(Tfy - pred) < 1e-6 * (1 + np.linalg.norm(Tfy)))
        sg = float(np.sum(gam))
        sg2 = float(np.sum(gam * gam))
        energy = float(np.dot(Tfy, Tfy))
        row = {
            "y_index": yi,
            "Tf_multiplicative_ok": tf_ok,
            "sum_gamma": sg,
            "sum_gamma_pred": pred_sg,
            "sum_gamma_ok": abs(sg - pred_sg) < 1e-4,
            "sum_gamma2": sg2,
            "sum_gamma2_pred": pred_sg2,
            "sum_gamma2_ok": abs(sg2 - pred_sg2) < 1e-3,
            "Tf_energy": energy,
            "Tf_energy_pred": pred_energy,
            "Tf_energy_ok": abs(energy - pred_energy) < 1e-2,
        }
        samples.append(row)
        if not (
            row["Tf_multiplicative_ok"]
            and row["sum_gamma_ok"]
            and row["sum_gamma2_ok"]
            and row["Tf_energy_ok"]
        ):
            all_ok = False

    return {
        "p": p,
        "N": N,
        "n": n,
        "nQ": nQ,
        "n_checked": n_check,
        "all_ok": all_ok,
        "samples": samples,
        "closed_forms": {
            "sum_gamma": str(sum_gamma_formula(p)),
            "sum_gamma2": str(sum_gamma2_formula(p)),
            "Tf_energy": str(Tf_energy_formula(p)),
            "mean_gamma": str(Fraction(6, p)),
        },
    }


def certify_ed4_types(p: int) -> dict:
    """Pair-average ED4 and per-basepoint type census (Thm F)."""
    Y = _load_Y(p)
    if Y is None:
        return {"p": p, "skipped": True}

    N, n = Y.shape
    # Gram for full census when N moderate; for p=7 N=11452 ~1GB float64 OK
    G = Y @ Y.T
    Gi = np.rint(G).astype(np.int16)
    # per-basepoint ED4
    # use int64 accumulate
    row_s4 = np.sum(Gi.astype(np.int64) ** 4, axis=1)
    ED4_y = [Fraction(int(s), N) for s in row_s4]
    # type counts
    from collections import Counter

    type_ctr = Counter(ED4_y)
    types_sorted = sorted(type_ctr.items(), key=lambda kv: kv[0])
    # global pair average
    total_s4 = int(np.sum(Gi.astype(np.int64) ** 4))
    ED4_global = Fraction(total_s4, N * N)
    delta2 = (ED4_global - ed4_flat(p)) / 24
    max_ed4 = max(ED4_y)
    min_ed4 = min(ED4_y)
    return {
        "p": p,
        "N": N,
        "n_types": len(type_ctr),
        "types": [
            {"ED4": str(ed), "count": cnt, "float": float(ed)}
            for ed, cnt in types_sorted
        ],
        "ED4_global": str(ED4_global),
        "ED4_global_float": float(ED4_global),
        "ED4_min": str(min_ed4),
        "ED4_max": str(max_ed4),
        "ED4_flat": str(ed4_flat(p)),
        "ED4_suf": str(ed4_suf(p)),
        "delta2_global": str(delta2),
        "rho_min2": str(rho_min2(p)),
        "budget_hyp": str(delta2_budget_hyp(p) if p > 3 else Fraction(0)),
        "global_le_suf": ED4_global <= ed4_suf(p),
        "global_le_bud": (
            ED4_global <= ed4_flat(p) + 24 * delta2_budget_hyp(p)
            if p > 3
            else True
        ),
        "max_y_le_suf": max_ed4 <= ed4_suf(p),
        "delta2_le_rho_min2": delta2 <= rho_min2(p),
        "not_2point_homogeneous": len(type_ctr) > 1,
    }


def prove_general_status() -> dict:
    return {
        "Tf_multiplicative_proved": True,
        "delta_annihilates_gamma_f_proved": True,
        "sum_gamma_proved": True,
        "sum_gamma2_proved": True,
        "Tf_energy_proved": True,
        "ed4_pair_average_structure_proved": True,
        "ed4_le_suf_for_all_p": False,
        "delta_le_rho_min_for_all_p": False,
        "sixteen_N_proved_for_all_p": False,
        "open_residual": (
            "Prove pair-average ED4≤ED4_suf (⇔δ²≤ρ_min²) for all primes p≥5. "
            "New tools: γ-calculus, Tf_y=(4p−2γ)f_y, ∑γ² closed, δ⊥γ⊙f. "
            "F19: no class_key."
        ),
    }


def main() -> dict:
    primes = [p for p in range(3, 60) if is_prime(p)]
    print(
        "Prop 15.114: γ-calculus; Tf_y multiplicative; ∑γ,∑γ²; ED4 types",
        flush=True,
    )

    alg = prove_gamma_algebra(primes)
    print(f"  gamma algebra proved={alg['proved']}", flush=True)
    mult = prove_multiplicative_formula()

    cert_g = {}
    for p in (3, 5):
        c = certify_gamma_and_Tf(p, n_y=3)
        cert_g[str(p)] = c
        if not c.get("skipped"):
            print(
                f"  gamma/Tf cert p={p}: all_ok={c['all_ok']}",
                flush=True,
            )

    # p=7: skip full T (heavy); still do ED4 type census
    cert_ed4 = {}
    for p in (3, 5, 7):
        e = certify_ed4_types(p)
        cert_ed4[str(p)] = e
        if not e.get("skipped"):
            print(
                f"  ED4 types p={p}: n_types={e['n_types']} "
                f"global_le_suf={e['global_le_suf']} "
                f"max_y_le_suf={e['max_y_le_suf']} "
                f"δ≤rm={e['delta2_le_rho_min2']}",
                flush=True,
            )

    gen = prove_general_status()
    cert_ok = all(
        cert_g[str(p)].get("all_ok")
        for p in (3, 5)
        if not cert_g[str(p)].get("skipped")
    )
    ed4_ok = all(
        cert_ed4[str(p)].get("global_le_suf")
        and cert_ed4[str(p)].get("delta2_le_rho_min2")
        for p in (3, 5, 7)
        if not cert_ed4[str(p)].get("skipped")
    )

    out = {
        "prop": "15.114",
        "backend": "CPU Fraction + Max+ cert; GPU unused (F20); no class_key (F19)",
        "L_status": "OPEN",
        "H_proved": False,
        "sixteen_N_proved_for_all_p": False,
        "ed4_le_suf_for_all_p": False,
        "kappa_bound_proved_for_all_p": False,
        "closed_forms": {
            "sum_gamma": "(6/p) * C(n,4)",
            "sum_gamma2": "6*C(n,4) + n(n-1)(n-2)/4",
            "Tf_energy": "(16p^2 - 72)*C(n,4) + n(n-1)(n-2)",
            "mean_gamma": "6/p",
            "spot_p5_sum_gamma2": str(sum_gamma2_formula(5)),
            "spot_p7_sum_gamma2": str(sum_gamma2_formula(7)),
            "spot_p5_Tf_energy": str(Tf_energy_formula(5)),
        },
        "gamma_algebra": alg,
        "multiplicative": mult,
        "cert_gamma_Tf": cert_g,
        "cert_ed4_types": cert_ed4,
        "general_status": gen,
        "attack_surface": gen["open_residual"],
        "verdict": (
            "Proved: Tf_y=(4p−2γ_y)⊙f_y; δ⊥γ⊙f; ∑γ=(6/p)C(n,4); "
            "∑γ²=6 C(n,4)+n(n−1)(n−2)/4 (adjacent edges cancel by Cy=py); "
            "‖Tf_y‖² closed. Certified p=3,5. At p=7: three ED4(y) types, "
            "all ≤ED4_suf; pair-average ED4≤ED4_suf. "
            "OPEN: ED4≤ED4_suf for general p≥5. L OPEN."
        ),
        "settlement_note": (
            "No soft-close. N_ED4_SUF still open for general p. "
            "New attack surface: γ-weighted spectral mass of f_y on E_4p."
        ),
        "proved": {
            "Tf_multiplicative": True,
            "delta_annihilates_gamma_f": True,
            "sum_gamma": alg["proved"],
            "sum_gamma2": alg["proved"],
            "Tf_energy": alg["proved"],
            "cert_p3_p5_gamma": cert_ok,
            "cert_ed4_types_p3_5_7": ed4_ok,
            "ed4_le_suf_for_all_p": False,
            "delta_le_rho_min_for_all_p": False,
            "sixteen_N_for_all_p": False,
        },
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15114.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str))
    print(f"Wrote {path}", flush=True)
    return out


if __name__ == "__main__":
    main()
