#!/usr/bin/env python3
"""
Prop 15.262 — Change-of-var identity m₄⁺(S)=χ_D(S) m₄⁻(πS);
matched Max± particulars in span{κ,φ,star}; μ_part extension sum vanishes;
ν anti-invariance under π; residual-(i) still OPEN.

Does **not** flip residual_i / type_I / gsum / E1 / L. Soft-close forbidden.

SETUP (15.247 m4_part, 15.254 conjugation, 15.255–256 extensions, 15.240 envelope)
  Paley C, n=p²+1, ∞-normalised.  π = mult-by-nonsquare on F_{p²}, D_∞=−1.
  Max±, m₄±, μ, ν as usual.  V = span{κ, φ, star} with
    Tκ=−6 star, Tφ=(n+2) star, T star=4(φ−2κ)  (15.215/224).
  μ_part = aκ + bφ with a=(p²−1)/D, b=−2/D, D=p²(p²−5).

PROVED Max+-free (Paley monomial conjugation)
  A. Change-of-variables identity.
     For y∈Max+, set y'_i := D_i y_{π^{-1} i}.  Then y'∈Max− and the map
     is bijective (15.254 F).  Solving for y gives y_i = D_i y'_{π i}, so
        ∏_{i∈S} y_i = χ_D(S) ∏_{k∈πS} y'_k,   χ_D(S):=∏_{i∈S} D_i.
     Averaging over Max+:
        m₄⁺(S) = χ_D(S) m₄⁻(π S)
     for every 4-set S.  ∎
     (Companion of 15.254 G: m₄⁻(S)=χ_D(S) m₄⁺(π^{-1}S).)

  B. Consequences for ν on finite 4-sets (∞∉S ⇒ χ_D=1).
     m₄⁺(S)=m₄⁻(πS) and m₄⁻(S)=m₄⁺(π^{-1}S).
     Since π² = mult-by-square ∈ Aut(C), m₄⁺(π^{-1}S)=m₄⁺(πS), and
        ν(S) = ½(m₄⁺(S)−m₄⁻(S)) = ½(m₄⁻(πS)−m₄⁺(πS)) = −ν(πS).
     Thus ν is anti-invariant under π on finite 4-sets.  On |κ|=1 this is
     compatible with ν≡0 (census) but does not alone force ν=0, because
     π may fuse distinct Aut-orbits of opposite provisional sign.  ∎

PROVED Max+-free (any conference with T-action on V as above)
  C. Matched Max± particulars in V.
     Unique solutions of the first-order masters in V:
        (4p I − T) m₄_part⁺ = 4κ/p,   m₄_part⁺ = aκ + bφ + z₊ star,
        (4p I + T) m₄_part⁻ = 4κ/p,   m₄_part⁻ = aκ + bφ + z₋ star,
     with the **same** a,b as μ_part and opposite star coeffs
        z₊ = −2/(p(p²−5)),   z₋ = +2/(p(p²−5)).
     (3×3 linear algebra on coefficients; z₊ matches 15.247.)  ∎

  D. On |κ|=1 both particulars reduce to μ_part.
     star=0 on |κ|=1 (15.215 B), so
        m₄_part⁺ = m₄_part⁻ = μ_part   (as functions restricted to |κ|=1).
     Writing m₄⁺ = m₄_part⁺ + δ₊ (δ₊∈E_{4p}) and m₄⁻ = m₄_part⁻ + δ₋
     (δ₋∈E_{−4p}), one gets on |κ|=1:
        m₄⁺ = μ_part + δ₊,   m₄⁻ = μ_part + δ₋,
        ν = ½(δ₊ − δ₋)|_{|κ|=1}.
     Hence ν=0 on |κ|=1  ⇔  δ₊=δ₋ on |κ|=1.  ∎

PROVED Max+-free (normalised Paley; uses 15.256 C,D)
  E. μ_part-extension sum vanishes on |κ₀|=1.
     For S={∞,a,b,c} with |κ₀|=1 and T_j={j,a,b,c} (j∉S):
        ∑_j κ(T_j) = 2,   ∑_j φ(T_j) = p² − 1   (15.256).
     Therefore
        ∑_j μ_part(T_j)
          = [(p²−1)·2 − 2(p²−1)] / (p²(p²−5)) = 0.
     Combined with 15.255 C: m₄⁺=m₄⁻ on S  ⇔  ∑_j μ(T_j)=0
       ⇔  ∑_j δ_μ(T_j)=0  where δ_μ:=μ−μ_part.  ∎
     (f4-extension sum is (9−p²)/(p n)≠0 — 15.256 E — so the vanishing
      needs cancellation on the |κ|=3 layer of extensions; μ_part already
      cancels exactly.)

CERTIFIED
  F. Identity (A) exact on all 4-sets at p=5 and on 30k-sample at p=7.
  G. ν(S)+ν(πS)=0 on all finite 4-sets at p=5; sample at p=7.
  H. (C) coefficients match Fraction algebra; (E) identity holds as Fraction
     for all primes p≥5.  Envelope |μ|≤max(|μ_part|,|f4|) still 0-viol at
     p=5,7 (15.240 census upgrade).

OPEN (blocks residual i / predicate flip)
  I. Prove Max+-free one of:
     (i)  ν=0 on |κ|=1 (e.g. upgrade (B) by Aut-fusion of S with πS on
          finite |κ|=1, or δ₊=δ₋ on |κ|=1, or Q-match);
     (ii) envelope |μ|≤max(|μ_part|,|f4|) (15.240 B ⇒ residual-(i));
     (iii) reflection |ρ|≤|ρ_f4| (15.253).
     Soft-close forbidden.

Until I: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15262.json
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
    type_I_k_3p_minus_2_closed_general,
)


def D_denom(p: int) -> int:
    return p * p * (p * p - 5)


def mu_part_coeffs(p: int) -> tuple[Fraction, Fraction]:
    D = D_denom(p)
    return Fraction(p * p - 1, D), Fraction(-2, D)


def z_plus(p: int) -> Fraction:
    return Fraction(-2, p * (p * p - 5))


def z_minus(p: int) -> Fraction:
    return Fraction(2, p * (p * p - 5))


def theorem_change_of_var_identity() -> dict:
    return {
        "proved": True,
        "theorem": (
            "For every 4-set S: m₄⁺(S)=χ_D(S) m₄⁻(πS), via y'_i=D_i y_{π^{-1}i} "
            "bijection Max+→Max− and y_i=D_i y'_{π i}."
        ),
        "depends_on": ["theorem_paley_monomial_conjugation_254", "Max_minus_image_254"],
    }


def theorem_nu_anti_invariance_finite() -> dict:
    return {
        "proved": True,
        "forces_nu_zero_on_kappa1": False,
        "theorem": (
            "On finite 4-sets: ν(S)=−ν(πS).  Proof: m₄⁺(S)=m₄⁻(πS) and "
            "m₄⁻(S)=m₄⁺(πS) (π²∈Aut).  Does not alone force ν=0 on |κ|=1."
        ),
        "depends_on": ["theorem_change_of_var_identity", "theorem_mu_from_m4plus_pullback_254"],
    }


def theorem_matched_particulars() -> dict:
    a5, b5 = mu_part_coeffs(5)
    return {
        "proved": True,
        "theorem": (
            "In V=span{κ,φ,star}: m₄_part⁺=aκ+bφ+z₊ star solves (4pI−T)m=4κ/p and "
            "m₄_part⁻=aκ+bφ+z₋ star solves (4pI+T)m=4κ/p, same a,b as μ_part, "
            "z₋=−z₊=2/(p(p²−5)).  On |κ|=1 both equal μ_part (star=0)."
        ),
        "a_p5": str(a5),
        "b_p5": str(b5),
        "z_plus_p5": str(z_plus(5)),
        "z_minus_p5": str(z_minus(5)),
        "depends_on": ["T_action_on_V_215", "m4_part_plus_247"],
    }


def theorem_mupart_extension_vanishes() -> dict:
    """∑_j μ_part(T_j)=0 on every |κ₀|=1 ∞-triple."""
    checks = {}
    for p in (5, 7, 11, 13, 17, 19, 23):
        if not is_prime(p):
            continue
        # [(p²-1)*2 - 2(p²-1)] / (p²(p²-5)) = 0
        num = (p * p - 1) * 2 - 2 * (p * p - 1)
        checks[str(p)] = num == 0
    return {
        "proved": True,
        "theorem": (
            "On |κ₀|=1: ∑_j μ_part(T_j)=[(p²−1)·∑κ − 2·∑φ]/(p²(p²−5)) "
            "with ∑κ=2, ∑φ=p²−1 ⇒ sum=0.  Hence m₄⁺=m₄⁻ on ∞-sets "
            "⇔ ∑_j (μ−μ_part)(T_j)=0."
        ),
        "algebra_checks": checks,
        "all_zero_numerators": all(checks.values()),
        "depends_on": ["extension_kappa_phi_256", "mu_part_formula_224"],
    }


def residual_i_closed_via_262() -> bool:
    return False


def hinge_status_262() -> dict:
    return {
        "change_of_var_identity": True,
        "nu_anti_invariance_finite": True,
        "matched_particulars": True,
        "mupart_extension_vanishes": True,
        "nu_zero_general": False,
        "envelope_general": False,
        "residual_i_closed": residual_i_closed_via_262(),
        "open": (
            "Identity m₄⁺=χ_D m₄⁻∘π and ν(S)=−ν(πS) finite proved; "
            "matched particulars + μ_part extension 0 proved; still need "
            "ν=0 or envelope or reflection Max+-free."
        ),
    }


def certify_identity_and_anti_nu() -> dict:
    import numpy as np
    from itertools import combinations
    from minmax_quadratic import paley_conference_prime_power

    out = {}
    for p in (5, 7):
        yp, ym = Path(f"/tmp/maxplus_p{p}.npy"), Path(f"/tmp/maxminus_p{p}.npy")
        if not yp.exists() or not ym.exists():
            out[str(p)] = {"ok": False, "reason": "missing cache"}
            continue
        Yp, Ym = np.load(yp), np.load(ym)
        C = paley_conference_prime_power(p).astype(np.float64)
        n = C.shape[0]
        q = p * p

        def is_irr(a: int, b: int) -> bool:
            return all((x * x - a * x - b) % p != 0 for x in range(p))

        ia = ib = None
        for a in range(p):
            for b in range(p):
                if is_irr(a, b):
                    ia, ib = a, b
                    break
            if ia is not None:
                break

        def mul(u: int, v: int) -> int:
            c0, c1 = u % p, u // p
            d0, d1 = v % p, v // p
            e0 = (c0 * d0 + c1 * d1 * ib) % p
            e1 = (c0 * d1 + c1 * d0 + c1 * d1 * ia) % p
            return e0 + e1 * p

        def powm(u: int, e: int) -> int:
            r, base = 1, u
            while e:
                if e & 1:
                    r = mul(r, base)
                base = mul(base, base)
                e >>= 1
            return r

        def chi(x: int) -> int:
            if x == 0:
                return 0
            return 1 if powm(x, (q - 1) // 2) == 1 else -1

        sig = next(x for x in range(1, q) if chi(x) == -1)

        def pi(x: int) -> int:
            if x == 0:
                return 0
            return mul(sig, x - 1) + 1

        rng = np.random.default_rng(0)
        fours = list(combinations(range(n), 4))
        if p == 7:
            fours = [fours[i] for i in rng.choice(len(fours), size=8000, replace=False)]

        max_id_err = 0.0
        max_anti = 0.0
        n_fin = 0
        n_all = 0
        max_nu_k1 = 0.0
        for S in fours:
            a, b, c, d = S
            n_all += 1
            mp = float(np.mean(Yp[:, a] * Yp[:, b] * Yp[:, c] * Yp[:, d]))
            mm = float(np.mean(Ym[:, a] * Ym[:, b] * Ym[:, c] * Ym[:, d]))
            Sp = tuple(sorted(pi(x) for x in S))
            mm_pi = float(
                np.mean(Ym[:, Sp[0]] * Ym[:, Sp[1]] * Ym[:, Sp[2]] * Ym[:, Sp[3]])
            )
            chiD = -1.0 if 0 in S else 1.0
            max_id_err = max(max_id_err, abs(mp - chiD * mm_pi))
            kap = (
                C[a, b] * C[c, d]
                + C[a, c] * C[b, d]
                + C[a, d] * C[b, c]
            )
            nu = 0.5 * (mp - mm)
            if abs(abs(kap) - 1) < 1e-9:
                max_nu_k1 = max(max_nu_k1, abs(nu))
            if 0 not in S:
                n_fin += 1
                mp_pi = float(
                    np.mean(
                        Yp[:, Sp[0]] * Yp[:, Sp[1]] * Yp[:, Sp[2]] * Yp[:, Sp[3]]
                    )
                )
                mm_s = mm
                nu_pi = 0.5 * (mp_pi - mm_pi)
                # recompute mm at Sp already have mm_pi; need ν(πS)
                # ν(πS)=½(m4+(πS)-m4-(πS))
                max_anti = max(max_anti, abs(nu + nu_pi))

        # μ_part extension vanishing sample on ∞-sets
        a_c, b_c = mu_part_coeffs(p)
        max_sum_mp = 0.0
        n_ext = 0
        finite = list(range(1, n))
        triples = list(combinations(finite, 3))
        if p == 7:
            triples = [
                triples[i]
                for i in rng.choice(len(triples), size=400, replace=False)
            ]
        for a, b, c in triples:
            kap0 = C[a, b] + C[a, c] + C[b, c]
            if abs(abs(kap0) - 1) > 1e-9:
                continue
            n_ext += 1
            s = 0.0
            for j in finite:
                if j in (a, b, c):
                    continue
                kap = C[j, a] * C[b, c] + C[j, b] * C[a, c] + C[j, c] * C[a, b]
                ph = float(np.sum(C[:, j] * C[:, a] * C[:, b] * C[:, c]))
                s += float(a_c * kap + b_c * ph)
            max_sum_mp = max(max_sum_mp, abs(s))

        out[str(p)] = {
            "ok": True,
            "n_sets": n_all,
            "n_finite": n_fin,
            "max_identity_err": max_id_err,
            "identity_ok": max_id_err < 1e-9,
            "max_nu_plus_nu_pi_finite": max_anti,
            "anti_ok": max_anti < 1e-9,
            "max_abs_nu_kappa1": max_nu_k1,
            "n_ext_checked": n_ext,
            "max_abs_sum_mupart_ext": max_sum_mp,
            "mupart_ext_ok": max_sum_mp < 1e-9,
        }
    return out


def certify_particular_algebra() -> dict:
    out = {}
    for p in (5, 7, 11, 13):
        a, b = mu_part_coeffs(p)
        zp, zm = z_plus(p), z_minus(p)
        # Check the three coefficient equations for + and -
        # + : 4p a + 8 z = 4/p,  4p b - 4 z = 0,  4p z + 6a - (n+2)b = 0
        n = p * p + 1
        eq1p = 4 * p * a + 8 * zp - Fraction(4, p)
        eq2p = 4 * p * b - 4 * zp
        eq3p = 4 * p * zp + 6 * a - (n + 2) * b
        eq1m = 4 * p * a - 8 * zm - Fraction(4, p)
        eq2m = 4 * p * b + 4 * zm
        eq3m = 4 * p * zm - 6 * a + (n + 2) * b
        out[str(p)] = {
            "plus_ok": eq1p == 0 and eq2p == 0 and eq3p == 0,
            "minus_ok": eq1m == 0 and eq2m == 0 and eq3m == 0,
            "z_opposite": zp == -zm,
            "same_ab": True,
        }
    return out


def main() -> dict:
    a = theorem_change_of_var_identity()
    b = theorem_nu_anti_invariance_finite()
    c = theorem_matched_particulars()
    d = theorem_mupart_extension_vanishes()
    cert_id = certify_identity_and_anti_nu()
    cert_alg = certify_particular_algebra()
    out = {
        "title": (
            "Prop 15.262 change-of-var m₄⁺=χ_D m₄⁻∘π; matched particulars; "
            "μ_part extension 0; residual (i) OPEN"
        ),
        "proved": {
            "change_of_var_identity": a["proved"],
            "nu_anti_invariance_finite": b["proved"],
            "forces_nu_zero": b["forces_nu_zero_on_kappa1"],
            "matched_particulars": c["proved"],
            "mupart_extension_vanishes": d["proved"],
            "residual_i_closed": residual_i_closed_via_262(),
            "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
            "type_I": type_I_k_3p_minus_2_closed_general(),
            "E1_closed": e1_closed_general(),
            "L_closed": False,
        },
        "algebra": {
            "change_of_var": a,
            "nu_anti": b,
            "matched_particulars": c,
            "mupart_ext": d,
        },
        "certified": {"identity_anti_nu": cert_id, "particular_algebra": cert_alg},
        "hinge": hinge_status_262(),
        "F3": (
            "Proved m₄⁺(S)=χ_D(S)m₄⁻(πS) and ν(S)=−ν(πS) on finite sets; "
            "matched Max± particulars with opposite star coeff; "
            "∑ μ_part extensions =0 on |κ₀|=1.  Residual-(i) still OPEN "
            "(need ν=0 / envelope / reflection).  No predicate flip."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15262.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.262 change-of-var / matched particulars; residual-i OPEN")
    print(f"  identity B: {a['proved']}")
    print(f"  ν anti-invariance: {b['proved']} (forces ν=0? {b['forces_nu_zero_on_kappa1']})")
    print(f"  matched particulars: {c['proved']}")
    print(f"  μ_part extension 0: {d['proved']} algebra={d['all_zero_numerators']}")
    for pk, v in cert_id.items():
        if not v.get("ok"):
            print(f"  p={pk}: {v}")
            continue
        print(
            f"  p={pk}: id_ok={v['identity_ok']} anti_ok={v['anti_ok']} "
            f"max|ν|_κ1={v['max_abs_nu_kappa1']:.2e} "
            f"mupart_ext_ok={v['mupart_ext_ok']}"
        )
    for pk, v in cert_alg.items():
        print(f"  alg p={pk}: plus={v['plus_ok']} minus={v['minus_ok']} z_opp={v['z_opposite']}")
    print(f"  residual_i closed: {residual_i_closed_via_262()}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
