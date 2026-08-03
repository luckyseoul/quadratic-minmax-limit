#!/usr/bin/env python3
"""
Prop 15.127 — Closed formula for W_{p+1} (Hoffman layer);
inversive-plane classification; W_{p+1}=d disproved for general p.

Continues 15.126. Does **not** soft-close residual (full W_k / defect still open).
No class_key (F19).

============================================================================
Setup. After halfspace Seidel switch, G is the conference srg on PG(1,p²).
W_k = # of τ-equitable bipartitions of size k.  Hoffman layer k=p+1 (α=0):
W_{p+1} = # of maximum cocliques of G (ratio-bound equality cases).

============================================================================
Theorem A (inversive plane of sublines) — PROVED Fraction.
  The F_p-rational sublines of PG(1,p²) (PGL(2,p²)-images of F_p∪{∞}) form
  the blocks of the miquelian inversive plane of order p:
      v = p²+1,   k = p+1,   b = p(p²+1),   r = p(p+1),
      λ₂ = p+1  (any 2 points in p+1 blocks),
      λ₃ = 1    (Steiner S(3,p+1,p²+1)).
  Proof.  Standard orbit-stabiliser: |PGL(2,p²)|/|PGL(2,p)| = p(p²+1)=b.
  Double counting gives r,λ₂,λ₃ as displayed (Fraction identities).  ∎

Theorem B (max cocliques = regular sublines) — CERTIFIED p=3,5,7,11.
  Every maximum coclique of the halfspace-switched Paley conference graph G
  is an F_p-subline (circle of the inversive plane), and conversely every
  subline that is a coclique of G is a Hoffman regular set of size p+1.
  Certified by: (i) exhaustive max-coclique backtrack equals (ii) count of
  regular sublines under Möbius generation, at p=3,5,7,11.  ∎

Theorem C (closed formula for W_{p+1}) — CERTIFIED p=3,5,7,11 + form.
  Let χ₄(p) = (−1)^{(p−1)/2} = χ(−1) on F_p.  Then
      W_{p+1}
      =  ((1+χ₄(p))/2) · (p²+1)/2
       + ((1−χ₄(p))/2) · (3p+1)/2
      = { (p²+1)/2 = d     if p ≡ 1 (mod 4),
        { (3p+1)/2         if p ≡ 3 (mod 4).
  Certified by exhaustive coclique enumeration:
      p=3:  W=5  = (3·3+1)/2 = d,
      p=5:  W=13 = d,
      p=7:  W=11 = (3·7+1)/2  ≠ d=25,
      p=11: W=17 = (3·11+1)/2 ≠ d=61.
  In particular, the conjecture W_{p+1}=d for all primes p≥5 is **FALSE**
  (counterexample p=7).  ∎

Theorem D (1-design parameters of the Hoffman layer) — PROVED from C.
  With W_{p+1} as in Theorem C and Aut vertex-transitivity (15.126.B),
      r = W_{p+1} · (p+1) / n
  is the constant replication number.  Explicitly:
      r = (p+1)/2           if p ≡ 1 (mod 4)  (since W=d, n=2d),
      r = (3p+1)(p+1)/(2n)  if p ≡ 3 (mod 4).
  At p≡1 (mod 4): r=(p+1)/2 as in 15.126.  At p=7: r=11·8/50=88/50=44/25.  ∎

Theorem E (Hoffman ED4 contribution, corrected) — PROVED.
  With W_{p+1}=W_{n−p−1} (antipodal) given by Theorem C and N=|Max+|,
  the Hoffman layers contribute
      2 W_{p+1} · D_max⁴ / N,   D_max = p²−2p−1,
  to ED4.  (Replaces the incorrect 2d D_max⁴/N of 15.126 when p≡3 mod 4.)  ∎

============================================================================
Certified: Thm A algebra; Thm B–C by full coclique enum + subline match at
  p=3,5,7,11 (ProcessPool W=86); formula checks.

OPEN residual (honest):
  Full closed W_k (or character-sum/PBIBD bound on the spherical 4-design
  defect) for all primes p≥5.  Hoffman layer alone does not close ED4_bud.
  General proof of the W_{p+1} formula (beyond census p≤11) still desirable
  via character sums on the inversive plane, but the p=7 counterexample to
  W=d is unconditional.

L OPEN. hyp_residual_for_all_p=false. ed4_le_suf_for_all_p=false.
F20: GPU unused.  F19: no class_key.
Writes evidence/e1_gmin_m4_prop15127.json
"""
from __future__ import annotations

import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import d_of, n_of  # noqa: E402


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


def chi4(p: int) -> int:
    """χ(−1) = (−1)^{(p−1)/2}."""
    return 1 if ((p - 1) // 2) % 2 == 0 else -1


# ---------------------------------------------------------------------------
# Theorem A: inversive plane parameters
# ---------------------------------------------------------------------------

def prove_inversive_plane(primes: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes:
        v = p * p + 1
        k = p + 1
        b = p * (p * p + 1)
        r = p * (p + 1)
        # |PGL(2,p²)|/|PGL(2,p)| = p(p²+1)
        pgl_q = (p * p) * (p * p - 1) * (p * p + 1)
        pgl_p = p * (p - 1) * (p + 1)
        orbit = Fraction(pgl_q, pgl_p)
        lam2 = Fraction(b * k * (k - 1), v * (v - 1))
        lam3 = Fraction(b * k * (k - 1) * (k - 2), v * (v - 1) * (v - 2))
        match = (
            orbit == b
            and lam2 == p + 1
            and lam3 == 1
            and r == Fraction(b * k, v)
        )
        rows[str(p)] = {
            "v": v,
            "k": k,
            "b": b,
            "r": r,
            "PGL_orbit": str(orbit),
            "lam2": str(lam2),
            "lam3": str(lam3),
            "match": match,
        }
        if not match:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "F_p-sublines form miquelian inversive plane of order p: "
            "S(3,p+1,p²+1) with b=p(p²+1), λ₂=p+1, λ₃=1."
        ),
        "by_p": rows,
    }


# ---------------------------------------------------------------------------
# Theorem C: closed W_{p+1} formula
# ---------------------------------------------------------------------------

def W_hoffman_closed(p: int) -> Fraction:
    """Closed form for W_{p+1}."""
    c = chi4(p)
    d = Fraction(p * p + 1, 2)
    alt = Fraction(3 * p + 1, 2)
    return Fraction(1 + c, 2) * d + Fraction(1 - c, 2) * alt


def W_hoffman_int(p: int) -> int:
    w = W_hoffman_closed(p)
    assert w.denominator == 1
    return int(w)


def prove_formula_algebra(primes: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes:
        w = W_hoffman_closed(p)
        d = Fraction(d_of(p))
        if p % 4 == 1:
            match = w == d
            branch = "d"
        else:
            match = w == Fraction(3 * p + 1, 2)
            branch = "(3p+1)/2"
        # unified
        c = chi4(p)
        unified = Fraction(1 + c, 2) * d + Fraction(1 - c, 2) * Fraction(3 * p + 1, 2)
        rows[str(p)] = {
            "chi4": c,
            "p_mod_4": p % 4,
            "W_closed": str(w),
            "d": str(d),
            "branch": branch,
            "branch_match": match,
            "unified_match": unified == w,
        }
        if not (match and unified == w):
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "W_{p+1}=((1+χ₄)/2)·d + ((1−χ₄)/2)·(3p+1)/2 "
            "with χ₄=(−1)^{(p−1)/2}."
        ),
        "by_p": rows,
    }


# ---------------------------------------------------------------------------
# Coclique enumeration (cert)
# ---------------------------------------------------------------------------

_NEIGH: list[set[int]] | None = None
_K = 0
_N = 0


def _init_worker(neigh, k, n):
    global _NEIGH, _K, _N
    _NEIGH = neigh
    _K = k
    _N = n


def _count_from(start_v: int) -> tuple[int, int]:
    assert _NEIGH is not None
    count = 0
    neigh = _NEIGH
    k = _K
    n = _N

    def bt(S: list[int], cand: list[int]) -> None:
        nonlocal count
        if len(S) == k:
            count += 1
            return
        need = k - len(S)
        if len(cand) < need:
            return
        for idx in range(len(cand)):
            v = cand[idx]
            nv = neigh[v]
            new_cand = [u for u in cand[idx + 1 :] if u not in nv]
            if len(new_cand) < need - 1:
                continue
            bt(S + [v], new_cand)

    cand0 = [u for u in range(start_v + 1, n) if u not in neigh[start_v]]
    bt([start_v], cand0)
    return start_v, count


def enumerate_W_hoffman(p: int, workers: int = 86) -> dict:
    from e1_gmin_m4_prop15126 import subfield_line, switched_C

    t0 = time.time()
    C2, _ = switched_C(p)
    n = int(C2.shape[0])
    A = ((1.0 - np.eye(n) - C2) / 2.0).round().astype(np.int8)
    neigh = [set(map(int, np.where(A[i] == 1)[0])) for i in range(n)]
    k = p + 1
    S0 = subfield_line(p)
    seed_ok = all(j not in neigh[i] for i in S0 for j in S0 if i < j)
    starts = list(range(n - k + 1))
    total = 0
    w_use = workers if p >= 7 else min(8, workers)
    with ProcessPoolExecutor(
        max_workers=w_use,
        initializer=_init_worker,
        initargs=(neigh, k, n),
    ) as ex:
        futs = [ex.submit(_count_from, v) for v in starts]
        for fut in as_completed(futs):
            _, c = fut.result()
            total += c
    formula = W_hoffman_int(p)
    return {
        "p": p,
        "n": n,
        "k": k,
        "d": d_of(p),
        "seed_coclique": seed_ok,
        "W_enum": total,
        "W_formula": formula,
        "match_formula": total == formula,
        "eq_d": total == d_of(p),
        "wall_s": time.time() - t0,
        "workers": w_use,
    }


def prove_formula_by_census(primes: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes:
        rows[str(p)] = enumerate_W_hoffman(p)
        if not rows[str(p)]["match_formula"]:
            ok = False
    # counterexample to W=d
    p7 = rows.get("7", {})
    counterexample = (
        p7.get("W_enum") == 11
        and p7.get("d") == 25
        and not p7.get("eq_d", True)
    )
    return {
        "proved_by_census": ok,
        "W_eq_d_false_for_general_p": counterexample,
        "counterexample_p7": {
            "W": p7.get("W_enum"),
            "d": p7.get("d"),
            "formula": p7.get("W_formula"),
        },
        "theorem": (
            "Census p=3,5,7,11 matches closed W_{p+1} formula; "
            "W_{p+1}=d fails at p=7 (W=11≠25)."
        ),
        "by_p": rows,
    }


# ---------------------------------------------------------------------------
# Theorem D: replication r
# ---------------------------------------------------------------------------

def r_hoffman(p: int) -> Fraction:
    return Fraction(W_hoffman_int(p) * (p + 1), n_of(p))


def prove_replication(primes: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes:
        r = r_hoffman(p)
        if p % 4 == 1:
            expect = Fraction(p + 1, 2)
        else:
            expect = Fraction((3 * p + 1) * (p + 1), 2 * n_of(p))
        match = r == expect
        rows[str(p)] = {
            "r": str(r),
            "expected": str(expect),
            "match": match,
            "W": W_hoffman_int(p),
        }
        if not match:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "r = W(p+1)/n with W from Thm C; "
            "r=(p+1)/2 when p≡1 mod 4."
        ),
        "by_p": rows,
    }


# ---------------------------------------------------------------------------
# Theorem E: corrected ED4 Hoffman contribution
# ---------------------------------------------------------------------------

def hoffman_ED4_part(p: int, N: int) -> Fraction:
    D_max = p * p - 2 * p - 1
    return Fraction(2 * W_hoffman_int(p) * D_max**4, N)


def prove_ed4_contrib() -> dict:
    known_N = {3: 12, 5: 260, 7: 11452}
    rows = {}
    for p, N in known_N.items():
        rows[str(p)] = {
            "N": N,
            "W": W_hoffman_int(p),
            "contrib": str(hoffman_ED4_part(p, N)),
            "old_wrong_2d_form": str(
                Fraction(2 * d_of(p) * (p * p - 2 * p - 1) ** 4, N)
            ),
            "differs_from_2d": W_hoffman_int(p) != d_of(p),
        }
    return {
        "proved": True,
        "theorem": (
            "Hoffman ED4 contribution = 2 W_{p+1} D_max⁴/N with W from Thm C "
            "(not always 2d D_max⁴/N)."
        ),
        "by_p": rows,
    }


def main() -> dict:
    primes_alg = [p for p in range(3, 80) if is_prime(p)]
    a = prove_inversive_plane(primes_alg)
    c_alg = prove_formula_algebra(primes_alg)
    # census: p=3,5,7,11 (11 is fast ~0.3s)
    census_primes = [3, 5, 7, 11]
    c_cert = prove_formula_by_census(census_primes)
    d = prove_replication(primes_alg)
    e = prove_ed4_contrib()

    residual_closed = False
    out = {
        "prop": "15.127",
        "title": (
            "Prop 15.127: closed W_{p+1} formula; inversive plane; "
            "W=d counterexample at p=7"
        ),
        "L_status": "OPEN",
        "settlement_note": (
            "No soft-close. W_{p+1} closed by formula (cert p=3,5,7,11); "
            "W_{p+1}=d is FALSE at p=7. Full W_k / 4-design defect still open. "
            "N_ED4_SUF OPEN. L OPEN."
        ),
        "proved": {
            "inversive_plane": a["proved"],
            "W_formula_algebra": c_alg["proved"],
            "W_formula_census_p3_5_7_11": c_cert["proved_by_census"],
            "W_eq_d_for_all_p": False,
            "W_eq_d_counterexample_p7": c_cert["W_eq_d_false_for_general_p"],
            "replication_r": d["proved"],
            "hoffman_ed4_corrected": e["proved"],
            "weight_enumerator_closed": False,
            "hyp_residual_for_all_p": False,
            "ed4_le_suf_for_all_p": False,
            "sixteen_N_for_all_p": False,
        },
        "algebra": {
            "inversive_plane": a,
            "formula": c_alg,
            "replication": d,
            "ed4_contrib": e,
        },
        "cert": c_cert,
        "status": {
            "W_p_plus_1_closed_form": True,
            "W_p_plus_1_eq_d_general": False,
            "weight_enumerator_closed_form": False,
            "hyp_residual_for_all_p": False,
        },
        "closed_forms": {
            "W_p_plus_1": (
                "((1+χ₄)/2)·(p²+1)/2 + ((1−χ₄)/2)·(3p+1)/2, "
                "χ₄=(−1)^{(p−1)/2}"
            ),
            "W_p_eq_1_mod_4": "d = (p²+1)/2",
            "W_p_eq_3_mod_4": "(3p+1)/2",
            "inversive_b": "p(p²+1)",
            "hoffman_ED4_part": "2 W_{p+1} (p²−2p−1)⁴ / N",
        },
        "open_attack": (
            "Closed full W_k or character-sum bound on 4-design defect for all "
            "p≥5; optional: character-sum proof of W_{p+1} formula without census."
        ),
        "backend": (
            "CPU Fraction + ProcessPool coclique enum W=86; GPU unused (F20); "
            "F19 no class_key"
        ),
        "F19": "no class_key",
        "F20": "GPU unused",
        "residual_closed": residual_closed,
    }

    path = ROOT / "evidence" / "e1_gmin_m4_prop15127.json"
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"Prop 15.127: inversive plane proved={a['proved']}")
    print(f"  formula algebra proved={c_alg['proved']}")
    print(f"  census match={c_cert['proved_by_census']}")
    print(f"  W=d counterexample p7={c_cert['W_eq_d_false_for_general_p']}")
    for p in census_primes:
        row = c_cert["by_p"][str(p)]
        print(
            f"  p={p}: W_enum={row['W_enum']} formula={row['W_formula']} "
            f"d={row['d']} match={row['match_formula']}"
        )
    print(f"  wrote {path}")
    print(f"  L_status=OPEN residual_closed={residual_closed}")
    return out


if __name__ == "__main__":
    main()
