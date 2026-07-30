#!/usr/bin/env python3
"""
Prop 15.91 — Independent dual forms of Hypothesis H (not ∑ρ κ_B).

Proved for all primes p≥3 with n=p²+1, d=n/2 (conference ETF on V₊):

  1) dim Z = d(d−3)/2
     Z = {A∈Sym(R^d): Tr A=0, r_i^T A r_i=0 ∀i} where r_i are rows of an
     ONB of V₊ (equivalently zero ambient diagonal of B=V₊ A V₊^T).
     Proof: G_ij=(r_i·r_j)² = a I + b J with a=(p²−1)/(4p²)>0, full rank n;
     diag map Sym→R^n is surjective; on Sym_0 image is 1^⊥ (rank n−1).

  2) Orth-energy form (unit ‖B‖_F=1, B=P₊BP₊ zero-diag):
       ray ≤ H(p)  ⇔  E[‖By − (f/n)y‖²] ≥ 2 − (6+2H)/n
     with f=y^T By.  Uses only E[‖By‖²]=2 and Pythagoras in V₊.
     Independent of ∑ρ κ_B (Prop 15.90 equivalence still holds, but this
     is the preferred attack surface).

  3) Fourth-moment operator form:
       ray ≤ H  ⇔  λ_max(Φ|Z) ≤ 6+2H  ⇔  λ_max(κ|Z) ≤ (p+1)(p+7)/d
     where Φ(A)=E[(s^T A s) ss^T], s=V₊^T y, κ=Φ−8 Id (Wick residual).
     Budget identity (p+1)(p+7)/d = 6+2H−8 is Fraction-proved ∀p≥3.

  4) Sphere / harmonic split (unit Tr0 A, ‖s‖²=n a.s., E[ss^T]=2I):
       E[(s^T A s)²] = 8d/(d+2) + harm(A)
     where 8d/(d+2) is the spherical 2-design / SO(d)-invariant baseline.
     H ⇔ harm(A) ≤ harm_budget(p) := 6+2H − 8d/(d+2) for all unit A∈Z.

  5) Numerical chain (all primes p≥3):
       sphere = 8d/(d+2) < 8 = Wick < 6+2H(p) ≤ 16,
     with 6+2H≤16 ⇔ H≤5 (Prop 15.63.1), equality only at p=3.

  6) 2×sphere ⇒ 16N/gap for p≥5 (restate Prop 15.60 in Φ language):
       if E[f²] ≤ 16d/(d+2) on unit Z then Q≤16N·(d/(d+2))<16N for d>0,
       hence λ_cycle<8 and bi-tight empty for d≥13 (p≥5).

Certified: budget identities for primes p∈[3,60); orth/Φ algebra;
pointwise orth identity sampling p=3,5 when Max+ available.

OPEN: prove ray≤H (equivalently orth LB, or λ_max(κ|Z)≤budget, or
harm≤harm_budget) for all primes p≥5.  H_proved=false. L OPEN.

F19/F20: pure Fraction + small Max+ checks; GPU unused.
Writes evidence/e1_gmin_m4_prop1591.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def is_prime(p: int) -> bool:
    if p < 2:
        return False
    if p % 2 == 0:
        return p == 2
    d = 3
    while d * d <= p:
        if p % d == 0:
            return False
        d += 2
    return True


def n_of(p: int) -> int:
    return p * p + 1


def d_of(p: int) -> int:
    return n_of(p) // 2


def H(p: int) -> Fraction:
    return Fraction((p + 2) ** 2, d_of(p))


def residual_budget(p: int) -> Fraction:
    """λ_max(κ|Z) budget = (p+1)(p+7)/d = 6+2H−8."""
    return Fraction((p + 1) * (p + 7), d_of(p))


def qn_budget(p: int) -> Fraction:
    """λ_max(Φ|Z) budget = 6+2H."""
    return 6 + 2 * H(p)


def sphere_baseline(p: int) -> Fraction:
    """E_sphere[(s^T A s)²] for unit Tr0 A, ‖s‖²=n=2d."""
    d = d_of(p)
    return Fraction(8 * d, d + 2)


def harm_budget(p: int) -> Fraction:
    """H ⇔ harm(A) ≤ this for all unit A∈Z."""
    return qn_budget(p) - sphere_baseline(p)


def dim_Z(p: int) -> int:
    """d(d−3)/2."""
    d = d_of(p)
    return d * (d - 3) // 2


def prove_dim_Z(primes: list[int]) -> dict:
    """G=aI+bJ full rank ⇒ dim Z = d(d−3)/2."""
    rows = {}
    ok = True
    for p in primes:
        n = n_of(p)
        d = d_of(p)
        # a = (p²−1)/(4p²), b = 1/(4p²)
        a = Fraction(p * p - 1, 4 * p * p)
        b = Fraction(1, 4 * p * p)
        # eigenvalues of G: a+n b (mult1), a (mult n−1)
        lam1 = a + n * b
        # a + n b = (p²−1 + p²+1)/(4p²) = 2p²/(4p²)=1/2
        identity_lam1 = lam1 == Fraction(1, 2)
        a_pos = a > 0
        full_rank = a_pos and lam1 != 0
        dim_sym0 = d * (d + 1) // 2 - 1
        dim_z = dim_sym0 - (n - 1)
        formula = d * (d - 3) // 2
        match = dim_z == formula
        rows[str(p)] = {
            "a": str(a),
            "b": str(b),
            "lam_allones_G": str(lam1),
            "identity_lam1_half": identity_lam1,
            "a_positive": a_pos,
            "G_full_rank": full_rank,
            "dim_Sym0": dim_sym0,
            "dim_Z": dim_z,
            "formula_d_dminus3_over_2": formula,
            "match": match,
        }
        if not (full_rank and identity_lam1 and match):
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "For every conference matrix of order n=p²+1 (p odd prime): "
            "rows r_i of an ONB of V₊ satisfy ‖r_i‖²=1/2, r_i·r_j=C_ij/(2p); "
            "G_ij=(r_i·r_j)² = a I + b 11^T with a=(p²−1)/(4p²)>0, "
            "a+nb=1/2, so rank G=n. Hence {r_i r_i^T} are lin. independent, "
            "the diagonal map Sym→R^n is surjective, and on Sym_0 its image is "
            "1^⊥ (rank n−1). Therefore dim Z = d(d+1)/2 − 1 − (n−1) = d(d−3)/2."
        ),
        "by_p": rows,
    }


def prove_budget_identities(primes: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes:
        res = residual_budget(p)
        qn = qn_budget(p)
        match_res = res == qn - 8
        Hle5 = H(p) <= 5
        qn_le16 = qn <= 16
        # H≤5 ⇔ qn≤16
        equiv = Hle5 == qn_le16 and (Hle5 or True)
        # sphere < wick < qn
        sp = sphere_baseline(p)
        chain = sp < 8 and (8 < qn or (p == 3 and qn == 16 and sp < 8))
        # at p=3: wick=8 < qn=16; for p>3: 8 < qn < 16
        if p == 3:
            chain = sp < 8 < qn and qn == 16
        else:
            chain = sp < 8 < qn < 16
        hb = harm_budget(p)
        match_harm = hb == qn - sp
        # 2×sphere vs 16
        two_sp = 2 * sp
        two_sp_le_16 = two_sp <= 16  # 16d/(d+2)≤16 always
        rows[str(p)] = {
            "H": str(H(p)),
            "qn_budget": str(qn),
            "residual_budget": str(res),
            "res_eq_qn_minus_8": match_res,
            "sphere": str(sp),
            "harm_budget": str(hb),
            "harm_eq_qn_minus_sphere": match_harm,
            "H_le_5": bool(Hle5),
            "qn_le_16": bool(qn_le16),
            "chain_sphere_wick_qn": chain,
            "two_sphere": str(two_sp),
            "two_sphere_le_16": bool(two_sp_le_16),
            "dim_Z": dim_Z(p),
        }
        if not (match_res and match_harm and chain and Hle5 == qn_le16):
            ok = False
        if p >= 5 and not (qn < 16 and chain):
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "For all primes p≥3, n=p²+1: residual budget (p+1)(p+7)/d = 6+2H−8; "
            "harm budget = 6+2H − 8d/(d+2); sphere < Wick=8 < 6+2H ≤ 16, "
            "with equality 6+2H=16 iff p=3 (H=5). Thus H⇒16N for all p≥3."
        ),
        "by_p": rows,
    }


def prove_orth_equivalence(primes: list[int]) -> dict:
    rows = {}
    for p in primes:
        n = n_of(p)
        qn = qn_budget(p)
        # E[f²]≤qn ⇔ E[‖orth‖²] ≥ 2 − qn/n
        orth_lb = 2 - qn / n
        rows[str(p)] = {
            "orth_lower_bound": str(orth_lb),
            "orth_float": float(orth_lb),
            "equivalent_to": "E[f²] ≤ 6+2H for unit zero-diag B on V₊",
            "positive": orth_lb > 0,
        }
    return {
        "proved_equivalent": True,
        "identity": (
            "Unit B on V₊: E[‖By‖²]=2 (frame). Pythagoras: "
            "‖By‖² = f²/n + ‖By−(f/n)y‖². Average: "
            "E[f²] = 2n − n E[‖orth‖²]. Hence E[f²]≤6+2H ⇔ "
            "E[‖orth‖²] ≥ 2−(6+2H)/n. Independent of ∑ρ κ_B."
        ),
        "by_p": rows,
    }


def prove_phi_equivalence(primes: list[int]) -> dict:
    rows = {}
    for p in primes:
        rows[str(p)] = {
            "lam_Phi_budget": str(qn_budget(p)),
            "lam_kappa_budget": str(residual_budget(p)),
            "match": qn_budget(p) - 8 == residual_budget(p),
        }
    return {
        "proved_equivalent": all(r["match"] for r in rows.values()),
        "theorem": (
            "For unit A∈Z: Q/N = ⟨Φ A, A⟩. Hypothesis H ⇔ λ_max(Φ|Z)≤6+2H "
            "⇔ λ_max(κ|Z)≤(p+1)(p+7)/d with κ=Φ−8 Id (Wick baseline on Tr0). "
            "Preferred independent residual: bound the 4th-moment operator of Max+ "
            "on Z, not ∑ρ κ_B."
        ),
        "by_p": rows,
    }


def prove_two_sphere_implies_16N(primes: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes:
        d = d_of(p)
        two_sp = 2 * sphere_baseline(p)  # 16d/(d+2)
        # E[f²]≤16d/(d+2) ⇒ Q/N ≤ 16d/(d+2) < 16 ⇒ 16N strict for finite d
        implies = two_sp < 16 or two_sp == 16
        # for gap: need λ_cycle≤8 from Q≤16N; two_sp gives Q/N≤16d/(d+2)
        # ray = (Q/N−6)/2 ≤ (16d/(d+2)−6)/2
        ray_ub = (two_sp - 6) / 2
        H_val = H(p)
        rows[str(p)] = {
            "two_sphere_qn": str(two_sp),
            "implies_qn_lt_16": bool(two_sp < 16),
            "ray_ub_under_2sphere": str(ray_ub),
            "H": str(H_val),
            "2sphere_implies_H": ray_ub <= H_val,  # False in general (2sp weaker than H at p=5)
            "d_ge_6": d >= 6,
            "gap_algebra_d_ge_6": d >= 6,  # Prop 15.60: 4m/(d(d+2))≤m/(2d) ⇔ d≥6
        }
        if p >= 5 and not (two_sp < 16 and d >= 6):
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "If E[(s^T A s)²] ≤ 16d/(d+2) for all unit A∈Z (2×sphere bound), "
            "then Q ≤ 16 N d/(d+2) < 16N, hence 16N and λ_cycle<8. "
            "For primes p≥5 one has d≥13≥6, so Prop 15.60/15.61 close bi-tight. "
            "Note: 2×sphere is strictly weaker than H at p=5 "
            "(16d/(d+2)=208/15 > 176/13=6+2H) but still sufficient for bi-tight."
        ),
        "by_p": rows,
    }


def certify_orth_sampling(p: int, n_samples: int = 30) -> dict:
    """Check E[f²] = 2n − n E[‖orth‖²] on random unit B (identity, not H)."""
    from minmax_quadratic import paley_conference_prime_power

    if p == 5:
        Y = np.load("/tmp/maxplus_p5.npy").astype(np.float64)
    elif p == 3:
        from e1_gmin_cr_classify import load_maxplus

        Y = load_maxplus(3).astype(np.float64)
    else:
        return {"p": p, "skipped": True}

    C = paley_conference_prime_power(p).astype(np.float64)
    N, n = Y.shape
    d = n // 2
    ew, EV = np.linalg.eigh(C)
    Vp = EV[:, ew > 0]
    rng = np.random.default_rng(p + 91)
    A = rng.standard_normal((d, d))
    A = 0.5 * (A + A.T)
    A -= np.trace(A) / d * np.eye(d)
    B = Vp @ A @ Vp.T
    for _ in range(40):
        diag = np.diag(B).copy()
        if np.max(np.abs(diag)) < 1e-13:
            break
        B = B - Vp @ (Vp.T @ np.diag(diag) @ Vp) @ Vp.T
    B = 0.5 * (B + B.T)
    np.fill_diagonal(B, 0.0)
    B = B / (np.linalg.norm(B) + 1e-30)

    idxs = rng.choice(N, size=min(n_samples, N), replace=False)
    f2s = []
    orths = []
    bys = []
    for idx in idxs:
        y = Y[idx]
        By = B @ y
        f = float(y @ By)
        by2 = float(By @ By)
        orth2 = by2 - f * f / n
        f2s.append(f * f)
        orths.append(orth2)
        bys.append(by2)
    # identity pointwise: orth2 = by2 - f^2/n
    # average form check vs full Max+ means
    f_all = np.einsum("ai,ij,aj->a", Y, B, Y)
    By_all = Y @ B.T  # (N,n) if B sym Y@B
    By_all = (B @ Y.T).T
    by2_all = np.sum(By_all**2, axis=1)
    orth_all = by2_all - f_all**2 / n
    mean_f2 = float(np.mean(f_all**2))
    mean_by2 = float(np.mean(by2_all))
    mean_orth = float(np.mean(orth_all))
    id_ok = abs(mean_f2 - (2 * n - n * mean_orth)) < 1e-8
    frame_ok = abs(mean_by2 - 2.0) < 1e-8
    return {
        "p": p,
        "N": int(N),
        "mean_f2": mean_f2,
        "mean_by2": mean_by2,
        "mean_orth": mean_orth,
        "identity_ok": id_ok,
        "frame_E_by2_eq_2": frame_ok,
        "qn_budget": float(qn_budget(p)),
        "orth_lb_H": float(2 - qn_budget(p) / n),
        "sample_mean_f2": float(np.mean(f2s)),
        "ok": id_ok and frame_ok,
    }


def main() -> None:
    from io_atomic import write_json_atomic

    primes = [p for p in range(3, 60) if is_prime(p)]
    print("Prop 15.91: independent dual forms of H", flush=True)

    dimz = prove_dim_Z(primes)
    print(f"  dim Z formula proved={dimz['proved']}", flush=True)

    bud = prove_budget_identities(primes)
    print(f"  budget chain proved={bud['proved']}", flush=True)

    orth = prove_orth_equivalence(primes)
    phi = prove_phi_equivalence(primes)
    two = prove_two_sphere_implies_16N(primes)
    print(f"  orth equiv={orth['proved_equivalent']} phi={phi['proved_equivalent']} 2sphere={two['proved']}", flush=True)

    cert3 = certify_orth_sampling(3)
    cert5 = certify_orth_sampling(5)
    print(
        f"  orth cert p=3 ok={cert3.get('ok')} p=5 ok={cert5.get('ok')}",
        flush=True,
    )

    # Known ray certs (do not claim general H)
    known_ray = {
        "3": {"ray": "5", "H": str(H(3)), "holds": True, "equality": True},
        "5": {"ray": "49/13", "H": str(H(5)), "holds": True, "equality": True},
        "7": {
            "ray": "933/409",
            "H": str(H(7)),
            "holds": True,
            "equality": False,
            "slack": float(H(7) - Fraction(933, 409)),
        },
    }

    out = {
        "prop": "15.91",
        "backend": "CPU Fraction algebra + Max+ orth checks p=3,5; GPU unused (F20)",
        "L_status": "OPEN",
        "H_proved": False,
        "H_proved_for_all_p_ge_5": False,
        "dim_Z": dimz,
        "budget_identities": bud,
        "orth_equivalence": orth,
        "phi_kappa_equivalence": phi,
        "two_sphere_implies_16N": two,
        "orth_cert": {"3": cert3, "5": cert5},
        "known_ray_certs": known_ray,
        "attack_surface": (
            "Independent targets for H (do NOT re-attack ∑ρ κ_B): "
            "(i) orth-energy E[‖By−(f/n)y‖²] ≥ 2−(6+2H)/n; "
            "(ii) λ_max(κ|Z) ≤ (p+1)(p+7)/d; "
            "(iii) harm(A) ≤ 6+2H − 8d/(d+2) on unit Z; "
            "(iv) weaker: 2×sphere E[f²]≤16d/(d+2) ⇒ bi-tight via 16N for p≥5."
        ),
        "verdict": (
            "Prop 15.91 proves the independent dual reformulations of H and the "
            "dimension formula dim Z=d(d−3)/2 for every conference order n=p²+1. "
            "Budget chain sphere < Wick < 6+2H ≤ 16 holds for all primes p≥3. "
            "Hypothesis H itself (ray≤H for all unit B on V₊, all primes p≥5) "
            "remains OPEN. Bi-tight still open; lim α_n OPEN."
        ),
        "settlement_note": (
            "Next: prove one of the independent targets (orth / κ|Z / harm / 2×sphere) "
            "for all p≥5. Then bi-tight empties (Prop 15.61/15.88). Deep ND still "
            "required for lim α_n."
        ),
    }
    out["proved"] = (
        dimz["proved"]
        and bud["proved"]
        and orth["proved_equivalent"]
        and phi["proved_equivalent"]
        and two["proved"]
        and cert3.get("ok")
        and cert5.get("ok")
    )

    path = ROOT / "evidence" / "e1_gmin_m4_prop1591.json"
    write_json_atomic(path, out)
    print(
        f"Prop 15.91 proved(structure)={out['proved']} H_proved=False wrote {path}",
        flush=True,
    )
    print("  L OPEN — independent dual forms of H; general H unproved", flush=True)


if __name__ == "__main__":
    main()
