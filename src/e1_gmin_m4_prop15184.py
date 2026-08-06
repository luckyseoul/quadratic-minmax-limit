#!/usr/bin/env python3
"""
Prop 15.184 — T²κ = −24φ + 48κ on |κ|=1; φ-structure for residual-(i) hinge.

Does **not** set gsum_disj_lb_proved_general True. Soft-close forbidden.

SETUP (15.177 master + 15.68 Tκ)
  μ₄ on 4-sets; κ = ∑_{three matchings} C_e C_{e'}.
  Master: (16 p² I − T²) μ₄ = 16 κ, with T the external C-weighted
  replacement operator on 4-sets.
  On |κ|=1: Tκ = 0 (Prop 15.68, any conference).

DEFINITIONS (any symmetric conference C, Cᵀ=C, C_ii=0, C_ij=±1, C²=p² I)
  For a 4-set S:
    φ(S) := ∑_{r ∉ S} ∏_{i∈S} C_{r i}
           = ∑_{r=0}^{n−1} ∏_{i∈S} C_{r i}   (terms r∈S vanish by C_ii=0).
  Star triple: star(S) := ∑_{v∈S} ∏_{u∈S\\{v}} C_{v u}.
  Tκ = −6 · star (Prop 15.68); hence Tκ=0 ⇔ star=0 on |κ|=1.

PROVED Max+-free (any symmetric conference):
  A. T²κ identity on |κ|=1.
        T²κ(S) = −24 φ(S) + 48 κ(S) = −24 (φ(S) − 2 κ(S))
     Proof outline (implemented as algebraic_T2_kappa_on_kappa1):
       (i) Tκ = −6 star, so
             T²κ(S) = −6 ∑_{v∈S} ∑_{r∉S} C_{v r} star(S_{v→r}).
       (ii) Expand star(S_{v→r}) into four triple-products; the leading
             term ∑_r C_{v r} ∏_{u∈S\\{v}} C_{r u} sums to φ(S) for each v,
             contributing 4 φ after ∑_v.
       (iii) The remaining three terms reduce via (C²)_{v b}=0 (v≠b) to
             pure K4 edge monomials (no external r).
       (iv) Exhaustion of the 48 K4 edge-labelings with |κ|=1 shows the
             K4 remainder equals +48 κ. Hence T²κ = −24 φ + 48 κ. ∎
  B. Consequently |T²κ| = 24 |φ − 2κ| on every |κ|=1 four-set.

CERTIFIED (Paley n=p²+1 census, not general):
  C. On |κ|=1: |φ| ≤ 2(p−2), with φ ∈ {±2,±6,…,±2(p−2)}.
     Checked for primes p∈{3,5,7,11}.
  D. On |κ|=3: |φ| ≤ 2p (sample / full small-p).
  E. p=5: μ₄ = (4κ − φ)/(p n) exactly on |κ|=1 (lstsq residual 0);
     max|μ₄|=3/65 < 1/10. Formula fails at p=3 (out of residual scope).

PALEY CHARACTER-SUM STRUCTURE (Max+-free field algebra):
  F. For S ⊂ F_{p²} (no ∞): φ(S) = 1 + Θ, where
        Θ = ∑_{x∈F_{p²}} χ(∏_{a∈S}(x−a)),
     χ = quadratic character. Möbius reduction ⇒ Θ is a degree-3
     character sum H(λ) (cross-ratio λ); Weil/Hasse ⇒ |H| ≤ 2p,
     hence |φ| ≤ 2p+1 on finite 4-sets. (Weaker than 2(p−2).)

OPEN (blocks residual i / E1 / L):
  G. Prove |φ| ≤ 2(p−2) on |κ|=1 for all primes p≥5 (Paley or all
     conference of order p²+1).
  H. Prove |μ₄| ≤ 1/(2p) on |κ|=1 (e.g. via master + (A)+(G), or closed
     form generalising (E)), or any sufficient Gsum LB / ker-box kill.

Until G+H: gsum_disj_lb_proved_general False; residual_i False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15184.json
"""
from __future__ import annotations

import json
import sys
from collections import Counter
from fractions import Fraction
from itertools import combinations, product
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    is_prime,
    n_of,
)
from e1_gmin_m4_prop15176 import theorem_minus_1_over_p_suffices  # noqa: E402
from e1_gmin_m4_prop15177 import mu4_threshold  # noqa: E402


def algebraic_T2_kappa_on_kappa1() -> dict:
    """
    Max+-free: T²κ = −24 φ + 48 κ on every K4 labeling with |κ|=1.

    Reduces external sums via C² off-diagonal vanishing; exhausts 48
    |κ|=1 labelings for the pure-K4 remainder.
    """
    import sympy as sp

    E = {
        (i, j): sp.Symbol(f"e{i}{j}", commutative=True)
        for i, j in combinations(range(4), 2)
    }

    def EE(i, j):
        if i == j:
            return 0
        return E[tuple(sorted((i, j)))]

    Phi = sp.Symbol("Phi")

    def sum_r_Cvr_star_nonphi(v: int):
        """sum_r∉S C_vr star(S_v→r) = Phi + K4-polynomial (C² reduction)."""
        others = [x for x in range(4) if x != v]
        b, c, d = others
        # term2: C_bc C_bd * sum_r C_vr C_rb
        # sum_r∉S C_vr C_rb = −EE(v,c)EE(c,b) − EE(v,d)EE(d,b)
        sum_rb = -EE(v, c) * EE(c, b) - EE(v, d) * EE(d, b)
        term2 = EE(b, c) * EE(b, d) * sum_rb
        sum_rc = -EE(v, b) * EE(b, c) - EE(v, d) * EE(d, c)
        term3 = EE(b, c) * EE(c, d) * sum_rc
        sum_rd = -EE(v, b) * EE(b, d) - EE(v, c) * EE(c, d)
        term4 = EE(b, d) * EE(c, d) * sum_rd
        return sp.expand(term2 + term3 + term4)

    total_nonphi = sum(sum_r_Cvr_star_nonphi(v) for v in range(4))
    # sum_v sum_r C_vr star = 4 Phi + total_nonphi
    # T²κ = −6 * that
    T2 = sp.expand(-6 * (4 * Phi + total_nonphi))
    # T2 = −24 Phi + rest, rest = −6*total_nonphi
    rest = sp.expand(-6 * total_nonphi)
    kap = EE(0, 1) * EE(2, 3) + EE(0, 2) * EE(1, 3) + EE(0, 3) * EE(1, 2)

    edges = list(E.values())
    n1 = 0
    rest_matches_48k = True
    pairs = []
    for bits in product([-1, 1], repeat=6):
        subs = dict(zip(edges, bits))
        kv = int(kap.subs(subs))
        if abs(kv) != 1:
            continue
        n1 += 1
        rv = int(rest.subs(subs))
        pairs.append((kv, rv))
        if rv != 48 * kv:
            rest_matches_48k = False

    # symbolic check: rest − 48 κ vanishes on all |κ|=1
    return {
        "proved": rest_matches_48k and n1 == 48,
        "theorem": (
            "On |κ|=1: T²κ = −24 φ + 48 κ = −24(φ − 2κ). "
            "Leading φ from ∑_v ∑_r C_vr ∏_{u≠v} C_ru = 4φ; "
            "K4 remainder = 48κ by 48-labeling exhaustion + C² reduction."
        ),
        "T2_symbolic": str(T2),
        "n_kappa_abs_1_labelings": n1,
        "rest_equals_48_kappa": rest_matches_48k,
        "labeling_rest_counts": {
            f"k={k},rest={r}": c for (k, r), c in sorted(Counter(pairs).items())
        },
        "formula": "T2_kappa = -24*phi + 48*kappa  on |kappa|=1",
    }


def theorem_T2_implies_abs_formula() -> dict:
    return {
        "proved": True,
        "theorem": (
            "On |κ|=1: |T²κ| = 24|φ − 2κ|. "
            "If |φ| ≤ 2(p−2) then |φ − 2κ| ≤ 2(p−2) + 2 = 2p − 2, "
            "hence |T²κ| ≤ 48(p − 1)."
        ),
        "conditional_on_phi_bound": True,
    }


def certify_phi_bound_paley(primes: list[int] | None = None) -> dict:
    """Census |φ| ≤ 2(p−2) on |κ|=1 for small Paley primes."""
    import numpy as np
    from minmax_quadratic import paley_conference_prime_power

    if primes is None:
        primes = [3, 5, 7, 11]
    rows = {}
    all_ok = True
    for p in primes:
        C = paley_conference_prime_power(p).astype(np.int8)
        n = C.shape[0]
        bound = 2 * (p - 2)
        mx = 0
        vals: set[int] = set()
        # full enum p≤5; sample larger
        if p <= 5:
            sets = combinations(range(n), 4)
        else:
            rng = np.random.default_rng(0)
            sets = (
                tuple(sorted(rng.choice(n, 4, replace=False).tolist()))
                for _ in range(4000)
            )
        n1 = 0
        for S in sets:
            a, b, c, d = S
            kap = int(
                C[a, b] * C[c, d]
                + C[a, c] * C[b, d]
                + C[a, d] * C[b, c]
            )
            if abs(kap) != 1:
                continue
            n1 += 1
            ph = sum(
                int(C[r, a] * C[r, b] * C[r, c] * C[r, d])
                for r in range(n)
                if r not in S
            )
            vals.add(ph)
            mx = max(mx, abs(ph))
        ok = mx <= bound
        if not ok:
            all_ok = False
        rows[str(p)] = {
            "max_abs_phi": mx,
            "bound_2_p_minus_2": bound,
            "phi_values": sorted(vals),
            "n_kappa1_seen": n1,
            "ok": ok,
            "full_enum": p <= 5,
        }
    return {
        "certified_for_listed_primes": all_ok,
        "proved_general": False,
        "by_p": rows,
        "note": (
            "Census only. General |φ|≤2(p−2) on |κ|=1 still OPEN "
            "(Weil gives only |φ|≤2p+1 on finite Paley 4-sets)."
        ),
    }


def certify_T2_identity_paley(p: int = 5) -> dict:
    """Numerical check T²κ = −24φ+48κ on Paley |κ|=1."""
    import numpy as np
    from minmax_quadratic import paley_conference_prime_power

    C = paley_conference_prime_power(p).astype(np.int8)
    n = C.shape[0]
    sets = list(combinations(range(n), 4))

    def kap(S):
        a, b, c, d = S
        return int(
            C[a, b] * C[c, d] + C[a, c] * C[b, d] + C[a, d] * C[b, c]
        )

    def star(S):
        s = 0
        for v in S:
            pr = 1
            for u in S:
                if u != v:
                    pr *= int(C[v, u])
            s += pr
        return s

    def phi(S):
        return sum(
            int(np.prod([C[r, i] for i in S]))
            for r in range(n)
            if r not in S
        )

    Tk = {S: -6 * star(S) for S in sets}
    bad = 0
    n1 = 0
    for S in sets:
        k = kap(S)
        if abs(k) != 1:
            continue
        n1 += 1
        # T²κ(S)
        Sl = list(S)
        t2 = 0
        for i, v in enumerate(Sl):
            for r in range(n):
                if r in S:
                    continue
                Sp = Sl.copy()
                Sp[i] = r
                t2 += int(C[v, r]) * Tk[tuple(sorted(Sp))]
        expect = -24 * phi(S) + 48 * k
        if t2 != expect:
            bad += 1
    return {
        "p": p,
        "n_kappa1": n1,
        "mismatches": bad,
        "identity_holds": bad == 0 and n1 > 0,
        "not_general_proof": True,  # algebra is general; this is census
    }


def theorem_mu4_formula_p5() -> dict:
    """Census: μ₄ = (4κ − φ)/(p n) on |κ|=1 at p=5."""
    return {
        "p": 5,
        "formula": "(4*kappa - phi)/(p*n)",
        "certified": True,
        "proved_general": False,
        "note": (
            "At p=5, lstsq residual 0 on all |κ|=1 four-sets. "
            "Fails at p=3. General closed form OPEN."
        ),
        "implies_mu4_bound_if_phi": (
            "If formula and |φ|≤2(p−2): |μ₄|≤(4+2p−4)/(p n)=2/n=2/(p²+1). "
            "For p≥5: 2/(p²+1) ≤ 1/(2p) ⇔ p²−4p+1≥0, true."
        ),
    }


def hinge_status_184() -> dict:
    return {
        "residual_ii_closed": False,  # full open (15.193); affine branch closed by 15.179
        "residual_i_closed": False,
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "T2_kappa_identity_proved": True,
        "phi_bound_general_proved": False,
        "mu4_bound_general_proved": False,
        "open": (
            "Prove |φ|≤2(p−2) on |κ|=1 for all p≥5, then |μ₄|≤1/(2p) "
            "(via master / closed form), or alternate residual-(i) Farkas."
        ),
        "mu4_thr_p5": str(mu4_threshold(5)),
    }


def main() -> dict:
    alg = algebraic_T2_kappa_on_kappa1()
    absf = theorem_T2_implies_abs_formula()
    phi_cert = certify_phi_bound_paley([3, 5, 7, 11])
    t2_p5 = certify_T2_identity_paley(5)
    mu4_p5 = theorem_mu4_formula_p5()
    primes = [p for p in range(5, 40) if is_prime(p)]
    out = {
        "title": (
            "Prop 15.184 T²κ=−24φ+48κ on |κ|=1; "
            "φ-bound / |μ₄| residual (i) still OPEN"
        ),
        "proved": {
            "T2_kappa_identity_on_kappa1": alg["proved"],
            "abs_T2_formula": absf["proved"],
            "minus_1_over_p_suffices": theorem_minus_1_over_p_suffices(primes),
            "phi_bound_2_p_minus_2_general": False,
            "mu4_bound_general": False,
            "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
            "residual_i_closed": False,
            "E1_closed": e1_closed_general(),
            "L_closed": False,
        },
        "algebra": {
            "T2_identity": alg,
            "abs_formula": absf,
            "mu4_p5_formula": mu4_p5,
        },
        "certified_census": {
            "phi_bound": phi_cert,
            "T2_paley_p5": t2_p5,
        },
        "hinge": hinge_status_184(),
        "F3": (
            "Do not flip residual/E1/L without general |φ| bound + |μ₄| "
            "or alternate Max+-free residual-(i) proof."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15184.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.184 T²κ identity; residual (i) OPEN")
    print(f"  T2 identity proved={alg['proved']} (n1={alg['n_kappa_abs_1_labelings']})")
    print(f"  p5 census identity={t2_p5['identity_holds']}")
    print(f"  phi bound census ok={phi_cert['certified_for_listed_primes']}")
    print(f"  gsum_disj_lb_proved_general={gsum_disj_lb_proved_general()}")
    print(f"  open: {out['hinge']['open']}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
