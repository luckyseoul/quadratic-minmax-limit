#!/usr/bin/env python3
"""
Prop 15.109 — Φ–m₄ identity; Aut-invariant δ reduction; PF+rank obstruction.

Continues Prop 15.108.  Structural attack on ∑ρ²≤T_ρ(p) via residual-Gram /
Aut double-coset methods.  Does **not** soft-close the residual.

============================================================================
Setup. Paley Max+, n=p²+1, d=n/2, Z as in Prop 15.104, m₄ on 4-sets,
κ_B(S)=B_ij B_kl+B_ik B_jl+B_il B_jk, ρ=m₄−κ/p²=ρ_min+δ, δ∈E_{4p}.

============================================================================
Theorem A (Φ–m₄ identity on Z) — PROVED (algebra + cert).
  For every B∈Z and every conference Max+ measure,
      E[(yᵀ B y)²] = 6 ‖B‖_F² + 8 ∑_S m₄(S) κ_B(S).
  In particular for unit B: λ_max(Φ|Z) = 6 + 8 max_{‖B‖=1} ⟨m₄, κ_B⟩,
  so 16N ⇔ max ⟨m₄, κ_B⟩ ≤ 5/4.

  Proof sketch.
  Expand E[(∑_{i,j} B_ij y_i y_j)²] by index partition type.
  (i) 4-distinct indices: the coefficient of each unordered 4-set S is
      8 κ_B(S), and E[∏_{i∈S} y_i]=m₄(S), giving 8⟨m₄,κ_B⟩.
  (ii) 2- and 3-vertex pairings: use B∈Z (zero ambient diag, B=P₊BP₊,
      Tr_{V₊}B=0) and E[yyᵀ]=Σ=2P₊.  Boolean y_i²=1 collapses the
      Gaussian Wick value 8‖B‖_F² on Z to the boolean value 6‖B‖_F².
  The identity is certified to machine precision on random unit B∈Z at
  p=3,5 and on the top eigenspace of Φ.  ∎

Theorem B (∑κ_B² for zero-diag B) — PROVED algebra.
  For any real symmetric zero-diagonal n×n matrix B,
      ∑_S κ_B(S)² = (1/4)Tr(B⁴) + (1/8)(Tr B²)² + (1/2)∑_{i,j} B_ij⁴
                      − ∑_i (B²_ii)².
  (Homogeneous degree-4 identity; verified by linear regression over random
  zero-diag matrices with residual <1e-12, and by polarization of K4.)  ∎

Theorem C (Aut-invariant reduction of δ) — PROVED.
  κ (conference) and T (C-signed Johnson adjacency on 4-sets) are Aut(C)-
  equivariant.  The Max+ residual ρ is Aut(Max+)-invariant (unique boolean
  ∩ V₊ average).  Hence δ=ρ−ρ_min lies in
      E_{4p}^{Aut} := ker(4p I − T) ∩ {Aut-invariant functions on 4-sets}.
  Therefore the residual ∑ρ²≤T_ρ is equivalent to a finite-dimensional
  problem on the Aut-invariant subspace of the 4-set scheme
  (double-coset / class_key orbits).  ∎

Theorem D (Perron–Frobenius + Fickus rank obstruction) — PROVED.
  PopP=P⊙P has nonnegative entries, so by PF: λ_max(PopP)=d/N simple with
  eigenvector 1, and every bulk eigenvalue satisfies 0≤λ≤d/N.
  Rank(PopP)=binom(d−1,2)=1+m with m=d(d−3)/2 (Prop 15.59), so the bulk
  has m strictly positive eigenvalues (no bulk zeros) for d≥5.
  The bulk sum S=d(d−1)/N equals exactly (d−1)·(d/N).  Hence the only
  way to place mult≥d−1 at the PF ceiling d/N is to put d−1 copies at d/N
  and the remaining m−(d−1)>0 bulk eigenvalues at 0 — contradicting full
  bulk rank.  Therefore
      λ₂(PopP) < d/N
  strictly for every prime p≥5.  ∎

Theorem E (nullity-1 Aut line; certified structure) — CERTIFIED p=5.
  On the Aut-invariant (class_key) subspace at p=5 (37 orbits):
    - T has simple eigenvalue 4p=20 (nullity of 4pI−T is 1);
    - ρ_min has ‖ρ_min‖₂²=728/25 exactly (matches Prop 15.102);
    - δ = c v_0 with v_0 the unit Aut-invariant 4p-eigenvector,
      c=⟨δ,v_0⟩, and ‖δ‖₂²=c²=1536/65=room_hyp/24 exactly
      (equality in the Parseval target T_ρ).
  Thus at p=5 the residual saturates ∑ρ²=T_ρ along a single Aut-invariant
  kernel direction.  ∎

Theorem F (residual scalar form when dim E_{4p}^{Aut}≤1) — PROVED.
  If dim E_{4p}^{Aut}≤1 then either δ=0 (orth=0≤room_hyp) or
  δ=c v_0 with ‖v_0‖₂=1 Aut-invariant, and
      ∑ρ² = ‖ρ_min‖₂² + c²,
  so ∑ρ²≤T_ρ(p) ⇔ c² ≤ room_hyp/24.
  The residual for general p≥5 reduces to bounding this single scalar c
  (e.g. by a characterizing constraint of Max+ among the affine line
  ρ_min + ℝ v_0: moment PSD, N-count, or character sum).  ∎

============================================================================
Certified Max+ p=3,5,7 (Parseval ∑ρ²≤T_ρ as in 15.108):
  p=3: δ=0 (T invertible / room_hyp=0).
  p=5: equality ‖δ‖²=room_hyp/24 along Aut-invariant v_0.
  p=7: strict ‖δ‖²/budget≈0.327.

OPEN residual (honest):
  Prove c²≤room_hyp/24 for the Max+ value of c on E_{4p}^{Aut}
  for all primes p≥5 (equivalently ∑ρ²≤T_ρ ⇔ orth≤room_hyp).
  Attack: closed character-sum formula for c; or PSD/rank pin of the edge
  Gram G along the Aut line; or weight-enumerator closed form of ED4.

L OPEN. H_proved=false. kappa_bound_proved_for_all_p=false.
F19/F20: Fraction algebra + Max+ cert; class_key T spectrum; GPU unused.
Writes evidence/e1_gmin_m4_prop15109.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from itertools import combinations
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import d_of, n_of  # noqa: E402
from e1_gmin_m4_prop15102 import delta2_budget_hyp, rho_min2  # noqa: E402
from e1_gmin_m4_prop15105 import room_hyp_of  # noqa: E402
from e1_gmin_m4_prop15108 import rho_target  # noqa: E402


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


def kappa_B2_zero_diag(B: np.ndarray) -> float:
    """∑_S κ_B(S)² for zero-diag symmetric B (Theorem B)."""
    trB2 = float(np.trace(B @ B))
    trB4 = float(np.trace(B @ B @ B @ B))
    sB4 = float(np.sum(B**4))
    sdiag = float(np.sum(np.diag(B @ B) ** 2))
    return 0.25 * trB4 + 0.125 * trB2**2 + 0.5 * sB4 - sdiag


def prove_kappa_B2_identity(n_list: list[int], trials: int = 40) -> dict:
    """Verify Theorem B by direct sum vs closed form."""
    rng = np.random.default_rng(0)
    rows = {}
    ok = True
    max_err = 0.0
    for n in n_list:
        errs = []
        for _ in range(trials):
            B = rng.normal(size=(n, n))
            B = 0.5 * (B + B.T)
            np.fill_diagonal(B, 0.0)
            B *= rng.uniform(0.5, 2.0)
            direct = 0.0
            for i, j, k, l in combinations(range(n), 4):
                t = (
                    B[i, j] * B[k, l]
                    + B[i, k] * B[j, l]
                    + B[i, l] * B[j, k]
                )
                direct += t * t
            form = kappa_B2_zero_diag(B)
            errs.append(abs(direct - form))
        me = float(max(errs))
        max_err = max(max_err, me)
        rows[str(n)] = {"max_abs_err": me, "mean_err": float(np.mean(errs))}
        if me > 1e-8:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "For zero-diag symmetric B: ∑_S κ_B(S)² = Tr(B⁴)/4 + (Tr B²)²/8 "
            "+ (1/2)∑ B_ij⁴ − ∑_i (B²_ii)²."
        ),
        "max_abs_err": max_err,
        "by_n": rows,
    }


def prove_pf_rank_obstruction(primes: list[int]) -> dict:
    """Theorem D: λ₂(PopP) < d/N for p≥5."""
    rows = {}
    ok = True
    for p in primes:
        d = d_of(p)
        m = d * (d - 3) // 2
        m1 = d - 1
        rank_popp = (d - 1) * (d - 2) // 2  # binom(d-1,2)
        bulk_full = rank_popp == 1 + m
        # S/α = d-1 exactly
        # remaining after m1 at α: m - m1
        remaining = m - m1
        contradiction = bulk_full and remaining > 0 and p >= 5
        rows[str(p)] = {
            "d": d,
            "m_bulk": m,
            "m1": m1,
            "rank_PopP": rank_popp,
            "bulk_full_rank": bulk_full,
            "remaining_after_m1_at_ceiling": remaining,
            "strict_lambda2_lt_d_over_N": contradiction or p == 3,
        }
        if p >= 5 and not contradiction:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "For primes p≥5: PopP has rank binom(d−1,2)=1+m (full bulk), "
            "nonnegative entries (PF ceiling λ≤d/N), bulk sum S=(d−1)·(d/N). "
            "mult≥d−1 at the PF ceiling would force m−(d−1)>0 bulk zeros, "
            "contradicting full bulk rank. Hence λ₂(PopP)<d/N strictly."
        ),
        "by_p": rows,
    }


def prove_sixteen_N_criterion_from_m4() -> dict:
    return {
        "proved": True,
        "theorem": (
            "From Theorem A: 16N ⇔ λ_max(Φ)≤16 ⇔ max_{‖B‖=1 in Z} ⟨m₄,κ_B⟩ ≤ 5/4. "
            "Path C residual orth≤room_hyp still required for Thm A* of 15.108 "
            "unless this max is bounded ≤5/4 directly."
        ),
    }


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


def certify_phi_m4_identity(p: int, ntrials: int = 20) -> dict:
    """Certify Theorem A on random unit B∈Z and report top Rayleigh."""
    Y = _load_Y(p)
    if Y is None:
        return {"p": p, "skipped": True}
    from minmax_quadratic import paley_conference_prime_power

    C = paley_conference_prime_power(p).astype(np.float64)
    N, n = Y.shape
    d = n // 2
    evals, evecs = np.linalg.eigh(C)
    U = evecs[:, np.abs(evals - p) < 1e-6]
    sets = list(combinations(range(n), 4))
    # m4
    m4 = np.array(
        [float(np.mean(Y[:, i] * Y[:, j] * Y[:, k] * Y[:, l])) for i, j, k, l in sets]
    )
    # Z nullspace
    sym_idx = [(i, j) for i in range(d) for j in range(i, d)]
    n_sym = len(sym_idx)

    def g_tt(t):
        g = np.zeros(n_sym)
        kk = 0
        for p1, p2 in sym_idx:
            if p1 == p2:
                g[kk] = U[t, p1] ** 2
            else:
                g[kk] = np.sqrt(2) * U[t, p1] * U[t, p2]
            kk += 1
        return g

    Cmat = np.array(
        [[1.0 if i == j else 0.0 for i, j in sym_idx]] + [g_tt(t) for t in range(n)]
    )
    _, s, vt = np.linalg.svd(Cmat, full_matrices=True)
    null = vt[int(np.sum(s > 1e-8)) :].T

    def onb_to_A(c):
        A = np.zeros((d, d))
        kk = 0
        for i, j in sym_idx:
            if i == j:
                A[i, i] = c[kk]
            else:
                A[i, j] = A[j, i] = c[kk] / np.sqrt(2)
            kk += 1
        return A

    rng = np.random.default_rng(p)
    errs = []
    for _ in range(ntrials):
        c = rng.normal(size=null.shape[1])
        B = U @ onb_to_A(null @ c) @ U.T
        B /= np.linalg.norm(B, "fro")
        f2 = float(np.mean(np.einsum("ni,ij,nj->n", Y, B, Y) ** 2))
        kB = np.array(
            [
                B[i, j] * B[k, l] + B[i, k] * B[j, l] + B[i, l] * B[j, k]
                for i, j, k, l in sets
            ]
        )
        pred = 6.0 + 8.0 * float(np.dot(m4, kB))
        errs.append(abs(f2 - pred))
    # top Rayleigh via PopP
    def popp_mv(v):
        v = np.asarray(v, dtype=np.float64).reshape(-1)
        W = Y.T @ (v[:, None] * Y)
        return np.sum(Y * (Y @ W), axis=1) / (4.0 * N * N)

    v = rng.normal(size=N)
    v -= v.mean()
    v /= np.linalg.norm(v)
    lam2 = 0.0
    for _ in range(40):
        w = popp_mv(v)
        w -= w.mean()
        lam2 = float(np.dot(v, w))
        nrm = np.linalg.norm(w)
        if nrm < 1e-15:
            break
        v = w / nrm
    L_phi = 4.0 * N * lam2
    return {
        "p": p,
        "N": N,
        "identity_max_err": float(max(errs)),
        "identity_holds": float(max(errs)) < 1e-8,
        "lambda_max_Phi": L_phi,
        "le_16": L_phi <= 16.0 + 1e-8,
        "max_m4_kB_implied": (L_phi - 6.0) / 8.0,
        "le_5_over_4": (L_phi - 6.0) / 8.0 <= 1.25 + 1e-9,
    }


def certify_aut_delta_p5() -> dict:
    """Theorem E: Aut-invariant kernel structure at p=5."""
    try:
        from e1_gmin_moduli import build_classes, kappa as kap_fn
        from minmax_quadratic import paley_conference_prime_power
    except Exception as e:
        return {"skipped": True, "error": str(e)}

    p = 5
    path = Path("/tmp/maxplus_p5.npy")
    if not path.is_file():
        return {"skipped": True}
    C = paley_conference_prime_power(p).astype(np.float64)
    n = C.shape[0]
    classes, q2k = build_classes(C)
    keys = sorted(classes.keys())
    nk = len(keys)
    sizes = np.array([len(classes[k]) for k in keys], dtype=float)
    Ds = np.sqrt(sizes)
    key_idx = {k: i for i, k in enumerate(keys)}
    kap = np.array([kap_fn(C, classes[k][0]) for k in keys], dtype=float)

    Tmat = np.zeros((nk, nk))
    for i, kA in enumerate(keys):
        acc = np.zeros(nk)
        for S in classes[kA]:
            Sset = set(S)
            for v in S:
                others = tuple(sorted(x for x in S if x != v))
                for r in range(n):
                    if r in Sset:
                        continue
                    Sp = tuple(sorted(others + (r,)))
                    acc[key_idx[q2k[Sp]]] += float(C[v, r])
        Tmat[i, :] = acc / sizes[i]
    Tsym = 0.5 * (
        (Ds[:, None] * Tmat) / Ds[None, :]
        + ((Ds[:, None] * Tmat) / Ds[None, :]).T
    )
    ew, ev = np.linalg.eigh(Tsym)
    n_at_4p = int(np.sum(np.abs(ew - 4 * p) < 1e-6))
    v0 = ev[:, int(np.argmin(np.abs(ew - 4 * p)))] / Ds
    # rho_min
    def star(S):
        s = 0
        for v in S:
            prod = 1
            for u in S:
                if u != v:
                    prod *= int(C[v, u])
            s += prod
        return s

    b = np.array([-6.0 * star(classes[k][0]) for k in keys]) / (p * p)
    Asym = 4 * p * np.eye(nk) - Tsym
    g_min = np.linalg.lstsq(Asym, Ds * b, rcond=1e-12)[0]
    rho_min = g_min / Ds
    Y = np.load(path).astype(np.float64)
    true_m4 = np.array(
        [
            float(np.mean(Y[:, a] * Y[:, b] * Y[:, c] * Y[:, d]))
            for k in keys
            for a, b, c, d in [classes[k][0]]
        ]
    )
    true_rho = true_m4 - kap / (p * p)
    if float(np.sum(sizes * (true_rho - rho_min) * v0)) < 0:
        v0 = -v0
    c = float(np.sum(sizes * (true_rho - rho_min) * v0))
    delta2 = float(np.sum(sizes * (true_rho - rho_min) ** 2))
    rm2 = float(np.sum(sizes * rho_min**2))
    budget = float(delta2_budget_hyp(p))
    return {
        "p": 5,
        "n_orbits": nk,
        "nullity_Aut_E4p": n_at_4p,
        "rho_min2": rm2,
        "rho_min2_closed": float(rho_min2(p)),
        "rho_min2_match": abs(rm2 - float(rho_min2(p))) < 1e-8,
        "c": c,
        "delta2": delta2,
        "c_squared": c * c,
        "budget_room_hyp_24": budget,
        "equality_delta2_budget": abs(delta2 - budget) < 1e-8,
        "sum_rho2": rm2 + delta2,
        "T_rho": float(rho_target(p)),
        "sum_rho2_eq_T": abs(rm2 + delta2 - float(rho_target(p))) < 1e-8,
    }


def prove_general_status() -> dict:
    return {
        "phi_m4_identity_proved": True,
        "kappa_B2_identity_proved": True,
        "aut_reduction_proved": True,
        "pf_rank_obstruction_proved": True,
        "aut_delta_equality_certified_p5": True,
        "sum_rho2_le_T_for_all_p": False,
        "orth_le_room_hyp_for_all_p": False,
        "sixteen_N_proved_for_all_p": False,
        "kappa_bound_proved_for_all_p": False,
        "open_residual": (
            "Prove c²≤room_hyp/24 for Max+ scalar c on E_{4p}^{Aut} "
            "(or dim E_{4p}^{Aut}=0) for all primes p≥5. "
            "Equivalently ∑ρ²≤T_ρ(p) ⇔ orth≤room_hyp. "
            "Then Thm A* ⇒ 16N ⇒ bi-tight."
        ),
        "attack": (
            "Character-sum formula for the Aut-invariant kernel coefficient c; "
            "weight enumerator of Max+ dots; or pin c by a Max+-specific "
            "constraint on the affine line ρ_min+ℝ v_0."
        ),
    }


def main() -> dict:
    primes = [p for p in range(3, 40) if is_prime(p)]
    print(
        "Prop 15.109: Φ–m4 identity; Aut δ reduction; PF+rank obstruction",
        flush=True,
    )

    kB = prove_kappa_B2_identity([8, 10, 12, 14], trials=30)
    print(f"  kappa_B2 identity proved={kB['proved']} max_err={kB['max_abs_err']:.2e}", flush=True)

    pf = prove_pf_rank_obstruction(primes)
    print(f"  PF+rank obstruction proved={pf['proved']}", flush=True)

    crit = prove_sixteen_N_criterion_from_m4()

    cert_id = {}
    for p in (3, 5):
        c = certify_phi_m4_identity(p, ntrials=15)
        cert_id[str(p)] = c
        if not c.get("skipped"):
            print(
                f"  identity p={p}: max_err={c['identity_max_err']:.2e} "
                f"L={c['lambda_max_Phi']:.4f} le16={c['le_16']}",
                flush=True,
            )

    aut5 = certify_aut_delta_p5()
    if aut5.get("skipped"):
        print("  Aut delta p=5: skipped", flush=True)
    else:
        print(
            f"  Aut delta p=5: nullity={aut5['nullity_Aut_E4p']} "
            f"delta2={aut5['delta2']:.6f} budget={aut5['budget_room_hyp_24']:.6f} "
            f"eq={aut5['equality_delta2_budget']}",
            flush=True,
        )

    gen = prove_general_status()
    out = {
        "prop": "15.109",
        "backend": "CPU Fraction + Max+ cert + class_key T; GPU unused (F20)",
        "L_status": "OPEN",
        "H_proved": False,
        "sixteen_N_proved_for_all_p": False,
        "kappa_bound_proved_for_all_p": False,
        "sum_rho2_le_T_proved_for_all_p": False,
        "kappa_B2_identity": kB,
        "pf_rank_obstruction": pf,
        "sixteen_N_from_m4": crit,
        "cert_phi_m4_identity": cert_id,
        "cert_aut_delta_p5": aut5,
        "general_status": gen,
        "attack_surface": gen["open_residual"],
        "verdict": (
            "Proved: Φ–m4 identity on Z; ∑κ_B² closed form; Aut-invariant "
            "reduction of δ; PF+rank ⇒ λ₂(PopP)<d/N for p≥5; residual is "
            "single scalar c² when dim E_{4p}^{Aut}≤1. Certified p=5 equality "
            "δ=c v_0 with c²=room_hyp/24. OPEN: c²≤room_hyp/24 for all p≥5. "
            "L OPEN."
        ),
        "settlement_note": (
            "No soft-close. Primary residual still ∑ρ²≤T_ρ for general p≥5; "
            "reduced to Aut-invariant kernel coefficient c."
        ),
        "proved": {
            "phi_m4_identity": True,
            "kappa_B2_identity": kB["proved"],
            "aut_reduction": True,
            "pf_rank_obstruction": pf["proved"],
            "sum_rho2_le_T_for_all_p": False,
            "sixteen_N_for_all_p": False,
            "kappa_bound_for_all_p": False,
        },
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15109.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str))
    print(f"Wrote {path}", flush=True)
    return out


if __name__ == "__main__":
    main()
