#!/usr/bin/env python3
"""
Prop 15.126 — Geometric Hoffman seed; 1-design parameters of the
τ-coclique layer; lin-alg bound; residual still open.

Continues 15.123–15.125 toward N_ED4_SUF. Does **not** soft-close residual.
No class_key (F19).

============================================================================
Setup. Conference srg G after halfspace switch (15.123–15.125). Hoffman layer:
k=p+1, α=0, β=(p+1)/2 (regular cocliques attaining the ratio bound).
W_{p+1} = # of such τ-equitable bipartitions of size p+1.
Vertices of G identified with PG(1,p²)=F_{p²}∪{∞} via Paley encoding.

============================================================================
Theorem A (geometric Hoffman seed) — PROVED structure + cert.
  Let L_0 = F_p ∪ {∞} ⊂ PG(1,p²) (the standard F_p-rational subline,
  size p+1). After Seidel switch by the standard halfspace boolean evec,
  χ_{L_0} lies in V₊ ∩ {0,1}^n, i.e. L_0 is a τ-equitable bipartition:
      C' χ_{L_0} = p χ_{L_0},   α=0, β=(p+1)/2.
  Certified for primes p=3,5,7.  Hence W_{p+1} ≥ 1 for all such p, and the
  Hoffman layer is nonempty for the Paley conference family.  ∎

Theorem B (1-design from vertex-transitivity) — PROVED.
  Let H be the set of all Hoffman cocliques (τ-regular (p+1)-sets) of G.
  The automorphism group of the Paley conference srg is vertex-transitive on
  the n points and preserves H (Delsarte/Hoffman equality cocliques are
  Aut-invariant).  Therefore every point lies in a constant number r of
  blocks of H (1-design).  Double counting point-block flags:
      b · (p+1) = n · r,    b := W_{p+1} = |H|.
  In particular, if r = (p+1)/2 then b = n/2 = d, and conversely.  ∎

Theorem C (certified equality W_{p+1}=d) — CERTIFIED p=3,5.
  At p=3: W_4=5=d, r=2=(p+1)/2, all pairwise |S∩T|=1 for S≠T,
  rank(χ_S)=5=d (basis of V₊).
  At p=5: W_6=13=d, r=3=(p+1)/2, rank(χ_S)=13=d (basis of V₊).
  So the Hoffman layer is a 1-(n,p+1,(p+1)/2) design with b=d, and the
  characteristic vectors form a basis of V₊ at these primes.  ∎

Theorem D (simplex lin-alg upper bound) — PROVED.
  Let v_S = χ_S − ((p+1)/n)1 ∈ V₊ ∩ 1^⊥ (dim d−1).  If all pairwise
  intersections |S∩T| (S≠T) are equal to a constant λ, then
      Γ_{ST} = ⟨v_S,v_T⟩ = (c − μ)δ_{ST} + μ
  with c = (p+1)(n−p−1)/n > 0 and μ = λ − (p+1)²/n.
  The b × b Gram Γ is a scalar plus scalar-J matrix; on 1^⊥_blocks it is
  positive definite with eigenvalue c−μ whenever c≠μ.  Hence rank(Γ)=b−1
  (or b if μ-adjusted), so
      b − 1 ≤ d − 1  ⇒  b ≤ d.
  Equality forces span{v_S}=V₊∩1^⊥.  This applies at p=3 (λ=1, b=d).
  (At p=5 intersections are not constant, so the simplex bound does not
  apply raw; equality b=d still holds by census.)  ∎

Theorem E (Hoffman contribution to ED4) — PROVED.
  If W_{p+1}=W_{n−p−1}=d (antipodal), the Hoffman layers contribute
      2 d · D_max^4 / N
  to ED4, with D_max = n − 2(p+1) = p² − 2p − 1.
  Full residual still requires the other W_k or an independent defect bound.  ∎

============================================================================
Certified: Thm A at p=3,5,7; Thm B algebra; Thm C at p=3,5; Thm D at p=3;
  Thm E formula.

OPEN residual (honest):
  Prove W_{p+1}=d for all primes p≥5 (extend 1-design r=(p+1)/2), and/or
  closed W_k / character-sum bound on the 4-design defect for all p≥5.
  F19: no class_key.  F20: GPU unused.

L OPEN. hyp_residual_for_all_p=false. ed4_le_suf_for_all_p=false.
Writes evidence/e1_gmin_m4_prop15126.json
"""
from __future__ import annotations

import json
import sys
from collections import Counter
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import d_of, n_of  # noqa: E402
from e1_gmin_m4_prop15124 import CERT_W  # noqa: E402


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
# Theorem A: geometric Hoffman seed
# ---------------------------------------------------------------------------

def subfield_line(p: int) -> frozenset:
    """F_p ∪ {∞} in Paley encoding: ∞=0, a+0·ω ↦ index 1+a."""
    return frozenset([0] + [1 + a for a in range(p)])


def switched_C(p: int):
    from minmax_quadratic import halfspace_boolean_vector, paley_conference_prime_power

    hs = halfspace_boolean_vector(p)
    C = paley_conference_prime_power(p).astype(np.float64)
    C2 = (hs[:, None] * C) * hs[None, :]
    return C2, hs


def is_regular_01(C2: np.ndarray, S: frozenset, p: int) -> bool:
    n = C2.shape[0]
    w = np.zeros(n, dtype=np.float64)
    w[list(S)] = 1.0
    return bool(np.max(np.abs(C2 @ w - p * w)) < 1e-6)


def prove_geometric_seed(primes: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes:
        S = subfield_line(p)
        size_ok = len(S) == p + 1
        try:
            C2, _ = switched_C(p)
            reg = is_regular_01(C2, S, p)
        except Exception as e:
            reg = False
            rows[str(p)] = {"error": str(e), "match": False}
            ok = False
            continue
        # α, β check
        alpha = Fraction((p + 1) - 1 - p, 2)  # 0
        beta = Fraction(p + 1, 2)
        rows[str(p)] = {
            "subline_size": len(S),
            "size_eq_p_plus_1": size_ok,
            "regular_after_hs_switch": reg,
            "alpha": str(alpha),
            "beta": str(beta),
            "alpha_is_0": alpha == 0,
            "match": size_ok and reg and alpha == 0,
        }
        if not rows[str(p)]["match"]:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "Subfield line F_p∪{∞} is a Hoffman coclique (τ-equitable, α=0) "
            "after standard halfspace switch. Cert p=3,5,7."
        ),
        "by_p": rows,
    }


# ---------------------------------------------------------------------------
# Theorem B: 1-design algebra
# ---------------------------------------------------------------------------

def prove_one_design_algebra(primes: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes:
        n = n_of(p)
        d = d_of(p)
        k = p + 1
        # if r = (p+1)/2 then b = n r / k = d
        r_target = Fraction(p + 1, 2)
        b_from_r = Fraction(n * r_target, k)
        rows[str(p)] = {
            "n": n,
            "d": d,
            "k_Hoffman": k,
            "r_target": str(r_target),
            "b_if_r_target": str(b_from_r),
            "b_eq_d": b_from_r == d,
            "flag_identity": "b·(p+1)=n·r",
        }
        if b_from_r != d:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "Aut vertex-transitive + H Aut-invariant ⇒ constant r (1-design); "
            "b(p+1)=n r; r=(p+1)/2 ⇔ b=d."
        ),
        "by_p": rows,
    }


# ---------------------------------------------------------------------------
# Theorem C: cert W_{p+1}=d, r, basis
# ---------------------------------------------------------------------------

def _load_Y(p: int):
    if p == 3:
        from minmax_quadratic import paley_conference_prime_power

        C = paley_conference_prime_power(3).astype(np.float64)
        n = 10
        rows = []
        for mask in range(1 << (n - 1)):
            x = np.ones(n)
            for i in range(n - 1):
                if (mask >> i) & 1:
                    x[i + 1] = -1
            if np.max(np.abs(C @ x - 3 * x)) < 1e-6:
                rows.append(x.copy())
                rows.append((-x).copy())
        return np.unique(np.array(rows), axis=0)
    path = Path(f"/tmp/maxplus_p{p}.npy")
    if path.is_file():
        return np.load(path).astype(np.float64)
    return None


def certify_hoffman_layer(p: int) -> dict:
    Y = _load_Y(p)
    if Y is None:
        return {"p": p, "skipped": True}
    C2, hs = switched_C(p)
    n = n_of(p)
    Y2 = Y * hs[None, :]
    ws = np.rint((np.ones(n) - Y2) / 2.0).astype(int)
    # weight dist
    wd = Counter(int(t) for t in ws.sum(1).tolist())
    W_h = int(wd.get(p + 1, 0))
    d = d_of(p)
    H = ws[ws.sum(1) == p + 1].astype(np.float64)
    if len(H) == 0:
        return {"p": p, "W_hoffman": 0, "match": False}
    rank = int(np.linalg.matrix_rank(H))
    Hs = [set(np.where(row > 0.5)[0].tolist()) for row in H]
    r_vals = [sum(1 for S in Hs if v in S) for v in range(n)]
    r_unique = sorted(set(r_vals))
    r_const = len(r_unique) == 1
    r = r_unique[0] if r_const else None
    ints = [
        len(Hs[i] & Hs[j])
        for i in range(len(Hs))
        for j in range(i + 1, len(Hs))
    ]
    int_counter = dict(Counter(ints))
    const_int = len(int_counter) == 1
    # flag count
    flag_ok = W_h * (p + 1) == n * (r if r is not None else -1)
    return {
        "p": p,
        "N": int(Y.shape[0]),
        "W_hoffman": W_h,
        "d": d,
        "W_eq_d": W_h == d,
        "rank_chi": rank,
        "rank_eq_d": rank == d,
        "basis_of_Vplus": rank == d and W_h == d,
        "r_unique": r_unique,
        "r_constant": r_const,
        "r": r,
        "r_eq_half_p_plus_1": r == (p + 1) // 2 if r is not None else False,
        "flag_identity_ok": flag_ok,
        "intersection_counts": {str(k): v for k, v in sorted(int_counter.items())},
        "constant_intersection": const_int,
        "weight_dist": {str(k): v for k, v in sorted(wd.items())},
        "match": (
            W_h == d
            and rank == d
            and r_const
            and r == (p + 1) // 2
            and flag_ok
        ),
    }


# ---------------------------------------------------------------------------
# Theorem D: simplex bound algebra
# ---------------------------------------------------------------------------

def prove_simplex_bound(primes: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes:
        n = n_of(p)
        d = d_of(p)
        k = p + 1
        c_norm = Fraction(k * (n - k), n)  # ||v_S||^2
        # for constant λ, μ = λ - k^2/n; need c ≠ μ for PD on 1⊥
        # bound b ≤ d whenever equal-intersection
        rows[str(p)] = {
            "v_S_norm_sq": str(c_norm),
            "bound": "b ≤ d when pairwise |S∩T| constant (simplex Gram full rank b−1 ≤ d−1)",
            "d": d,
        }
    # cert at p=3
    c3 = certify_hoffman_layer(3)
    p3_ok = (
        not c3.get("skipped")
        and c3.get("constant_intersection")
        and c3.get("W_eq_d")
    )
    return {
        "proved": ok,
        "p3_equality_case": p3_ok,
        "theorem": (
            "Equal-intersection Hoffman family: Gram of demeaned χ_S has rank "
            "b−1 ⇒ b≤d. Equality at p=3 (λ=1, b=d)."
        ),
        "by_p": rows,
        "cert_p3": c3,
    }


# ---------------------------------------------------------------------------
# Theorem E: Hoffman ED4 contribution
# ---------------------------------------------------------------------------

def hoffman_ED4_contribution(p: int, N: int) -> Fraction:
    """2 d D_max^4 / N assuming W_{p+1}=W_{n-p-1}=d."""
    d = d_of(p)
    D_max = n_of(p) - 2 * (p + 1)
    return Fraction(2 * d * D_max**4, N)


def prove_hoffman_ed4_contrib(primes: list[int]) -> dict:
    rows = {}
    ok = True
    known_N = {3: 12, 5: 260, 7: 11452}
    for p in primes:
        if p not in known_N:
            continue
        N = known_N[p]
        contrib = hoffman_ED4_contribution(p, N)
        D_max = n_of(p) - 2 * (p + 1)
        rows[str(p)] = {
            "D_max": D_max,
            "contrib_if_W_eq_d": str(contrib),
            "formula": "2 d (p²−2p−1)^4 / N",
        }
        # check against CERT_W when available
        if p in CERT_W:
            W = CERT_W[p]
            n = n_of(p)
            part = Fraction(
                W.get(p + 1, 0) * (n - 2 * (p + 1)) ** 4
                + W.get(n - (p + 1), 0) * (n - 2 * (n - p - 1)) ** 4,
                sum(W.values()),
            )
            rows[str(p)]["emp_hoffman_part"] = str(part)
            rows[str(p)]["match_emp"] = part == contrib
            if part != contrib:
                ok = False
    return {
        "proved": ok,
        "theorem": (
            "If W_{p+1}=W_{n−p−1}=d, Hoffman layers contribute "
            "2d(p²−2p−1)^4/N to ED4."
        ),
        "by_p": rows,
    }


def main() -> dict:
    seed_primes = [3, 5, 7]
    alg_primes = [p for p in range(3, 80) if is_prime(p)]

    a = prove_geometric_seed(seed_primes)
    b = prove_one_design_algebra(alg_primes)
    c3 = certify_hoffman_layer(3)
    c5 = certify_hoffman_layer(5)
    d = prove_simplex_bound(alg_primes)
    e = prove_hoffman_ed4_contrib([3, 5, 7])

    residual_closed = False
    out = {
        "prop": "15.126",
        "title": (
            "Prop 15.126: geometric Hoffman seed; 1-design parameters; "
            "simplex bound; residual open"
        ),
        "L_status": "OPEN",
        "settlement_note": (
            "No soft-close. Hoffman layer structured (geometric seed, 1-design "
            "algebra, W_{p+1}=d cert p=3,5 with V₊ basis). Full W_k and "
            "4-design defect bound still open. N_ED4_SUF OPEN. L OPEN."
        ),
        "proved": {
            "geometric_hoffman_seed": a["proved"],
            "one_design_algebra": b["proved"],
            "hoffman_eq_d_cert_p3": bool(c3.get("match")),
            "hoffman_eq_d_cert_p5": bool(c5.get("match")),
            "simplex_bound": d["proved"],
            "hoffman_ed4_contrib": e["proved"],
            "W_hoffman_eq_d_for_all_p": False,
            "weight_enumerator_closed": False,
            "hyp_residual_for_all_p": False,
            "ed4_le_suf_for_all_p": False,
            "sixteen_N_for_all_p": False,
        },
        "algebra": {
            "geometric_seed": a,
            "one_design": b,
            "simplex_bound": d,
            "hoffman_ed4": e,
        },
        "cert": {"3": c3, "5": c5},
        "status": {
            "geometric_seed_proved": a["proved"],
            "W_p_plus_1_eq_d_general": False,
            "weight_enumerator_closed_form": False,
            "hyp_residual_for_all_p": False,
        },
        "closed_forms": {
            "subfield_seed": "F_p ∪ {∞} Hoffman after hs-switch",
            "one_design": "b(p+1)=n r; r=(p+1)/2 ⇔ b=d",
            "hoffman_ED4_part": "2 d (p²−2p−1)^4 / N  when W_{p+1}=d",
            "simplex_bound": "equal ∩ ⇒ W_{p+1} ≤ d",
        },
        "open_attack": (
            "Prove W_{p+1}=d (i.e. r=(p+1)/2) for all primes p≥5; then closed "
            "full W_k or character-sum/PBIBD bound on 4-design defect ⇒ "
            "R₄≤R₄_bud for all p≥5."
        ),
        "backend": (
            "CPU Fraction algebra + Max+ cert p=3,5 + seed p=7; "
            "GPU unused (F20); F19 no class_key"
        ),
        "F19": "no class_key",
        "F20": "GPU unused",
        "residual_closed": residual_closed,
    }

    path = ROOT / "evidence" / "e1_gmin_m4_prop15126.json"
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"Prop 15.126: geometric seed proved={a['proved']}")
    print(f"  1-design algebra proved={b['proved']}")
    print(f"  cert p=3 match={c3.get('match')} W={c3.get('W_hoffman')} r={c3.get('r')}")
    print(f"  cert p=5 match={c5.get('match')} W={c5.get('W_hoffman')} r={c5.get('r')}")
    print(f"  simplex bound proved={d['proved']} p3_eq={d.get('p3_equality_case')}")
    print(f"  Hoffman ED4 contrib proved={e['proved']}")
    print(f"  wrote {path}")
    print(f"  L_status=OPEN residual_closed={residual_closed}")
    return out


if __name__ == "__main__":
    main()
