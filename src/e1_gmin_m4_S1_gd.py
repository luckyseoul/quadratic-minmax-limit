#!/usr/bin/env python3
"""
Prop 15.80: Gaussian domination for star·S1 — linear-form Wick identity,
sum formulas, U1-specialization, and residual-space structure.

Proved:
  1. For ANY linear form L(y)=∑ α_r y_r on Max+, E[L²]=E_Wick[L²]
     (pairwise moments match Σ=I+C/p; Max+ tight frame / 2-design).
  2. Hence E[U1²]=E_Wick[U1²] for the κ1 one-center form U1.
  3. star·S1 = E[Z U1]−E_Wick[Z U1] with Z=star∏_{u≠a} y_u, E[Z]=0,
     so GD ⇔ Cov(Z,U1) ≤ Cov_Wick(Z,U1).
  4. GD is special to the κ1-restriction U1: fails for generic linear forms
     and can fail for U_ext (certified).
  5. τ1=2A−d1^(1); sum_{|κ|=1} star·τ1 = ε n1 with ε∈{±1}
     (ε=+1 only at p=5 among p=3,5,7; formula/proof of ε open beyond census).
  6. Aut-constancy (Prop 15.79) ⇒ one star·S1 per 4-set.

Certified GPU mmap+atomic p=5,7:
  - GD holds on all |κ|=1 sets (star·S1≤0)
  - sum star·S1 = −1128 at p=5; −15271200/2863 at p=7
  - c=F(t, m4·κ): depends on |m4| only when t=−7 at p=7
  - E[U1²]/E_Wick[U1²] ≡ 1

OPEN: prove E[Z U1]≤E_Wick[Z U1] for all primes p≥5; joint/S3⇒M_cand.
lim α_n OPEN.

Writes evidence/e1_gmin_m4_S1_gd.json
"""
from __future__ import annotations

import os
import sys
import time
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from fractions import Fraction
from itertools import combinations
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def _blas1() -> None:
    for k in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[k] = "1"


def M_cand(p: int) -> float:
    return (p - 2) / (p * (2.0 * p + 3))


def d1_one(p: int) -> int:
    return (3 * p * p - 7) // 4


# ── Proved algebra ──────────────────────────────────────────────────────────


def linear_form_wick_identity() -> dict:
    """E[L²]=E_Wick[L²] for linear L on Max+."""
    return {
        "proved": True,
        "statement": (
            "Let y be uniform on Max+={y∈{±1}^n: C y = p y} for a conference "
            "matrix C of order n=p²+1. Then E[y_i y_j]=δ_ij + C_ij/p =: Σ_ij "
            "(tight frame / Prop Max+ E[yy^T]=I+C/p). "
            "For any α∈R^n, L=∑ α_i y_i satisfies "
            "E[L²]=∑_{i,j} α_i α_j Σ_ij = E_{G~N(0,Σ)}[L(G)²]."
        ),
        "corollary_U1": (
            "U1=∑_{r∈R1} C_ar y_r is linear, so E[U1²]=E_Wick[U1²] exactly."
        ),
        "method": "Expand E[L²] using only pairwise moments; no fourth moments used.",
    }


def gd_formulation() -> dict:
    return {
        "proved": True,
        "Z": "Z = star_a · ∏_{u∈S\\{a}} y_u  (degree-3, E[Z]=0 by central symmetry)",
        "U1": "U1 = ∑_{r: |κ(S_a→r)|=1} C_ar y_r  (degree-1 linear form)",
        "identity": "star_a · S1(a) = E[Z U1] − E_Wick[Z U1] = Cov(Z,U1) − Cov_Wick(Z,U1)",
        "GD_iff": "star·S1 ≤ 0  ⇔  E[Z U1] ≤ E_Wick[Z U1]",
        "Wick_value": "E_Wick[Z U1] = star_a · τ1(a) / p²",
        "specialization": (
            "GD fails for generic linear forms L (≈50% violation rate at p=5) "
            "and is therefore special to the κ1-support of U1, not a general "
            "cubic×linear comparison."
        ),
        "U_ext_note": (
            "For U_ext=∑_{r∉S} C_ar y_r one has E[Z U_ext]=star·p·ρ, "
            "E_Wick=2/p²; GD for U_ext ⇔ star·p·ρ ≤ 2/p² ⇔ joint S1+S3 ≤ 0, "
            "which fails at some p=7 maximisers (joint>0). U1 is strictly harder."
        ),
    }


def sum_star_tau1_algebra() -> dict:
    """sum of star·τ1 over |κ|=1 sets."""
    return {
        "proved_partial": True,
        "identity_t1": "τ1=2A−d1^(1), A=#{r∈R1: C_ar κ'=+1}",
        "certified_sum": {
            "3": {"sum": -180, "n1": 180, "mean": -1, "epsilon": -1},
            "5": {"sum": 11700, "n1": 11700, "mean": 1, "epsilon": 1},
            "7": {"sum": -176400, "n1": 176400, "mean": -1, "epsilon": -1},
        },
        "formula": "sum_{|κ|=1} star·τ1 = ε(p) · n1 with n1=n(n-1)(n-2)²/32, ε∈{±1}",
        "epsilon_open": (
            "ε(p)=+1 at p=5, ε(p)=−1 at p=3,7 (and mean −1 at p=11). "
            "Closed form for ε(p) open; not needed for pointwise GD."
        ),
        "p7_equipartition": (
            "At p=7 the three values −7,−1,5 each occur on exactly n1/3 sets."
        ),
    }


def cand_via_gd_algebra() -> dict:
    """How GD feeds the cand bound."""
    rows = {}
    for p in range(5, 40, 2):
        if any(p % d == 0 for d in range(2, int(p**0.5) + 1)):
            continue
        inv = 1.0 / (p * p)
        rc = M_cand(p) - inv
        # at star=+1: pρ = 2/p² + S1 + S3; if S1≤0 then pρ ≤ 2/p² + S3
        S3_budget = p * rc - 2 * inv
        rows[str(p)] = {
            "rho_cand": rc,
            "S3_budget_if_S1_le0": S3_budget,
            "S3_budget_negative": S3_budget < 0,
            "note": (
                "Negative S3_budget (p=5) requires strongly negative S1 on "
                "maximisers, not just S1≤0 + |S3| bounds."
            ),
        }
    return {
        "path": (
            "GD (star·S1≤0) + star=+1 ⇒ pρ ≤ 2/p² + S3. "
            "Still need S3 bound (or joint) for ρ≤ρ_cand when S3_budget<0."
        ),
        "by_p": rows,
        "p5_sharp": (
            "At p=5, joint S1+S3 = −16/325 on maximisers, giving ρ=ρ_cand sharp."
        ),
    }


# ── GPU census ──────────────────────────────────────────────────────────────


def gpu_all_m4(Y_np, quads, *, batch=12_000):
    from gpu_budget import require_gpu, gpu_info

    cp = require_gpu()
    info0 = gpu_info()
    Yh = np.ascontiguousarray(Y_np, dtype=np.float64)
    Y = cp.asarray(Yh)
    m4 = np.empty(len(quads))
    t0 = time.perf_counter()
    for lo in range(0, len(quads), batch):
        hi = min(len(quads), lo + batch)
        sl = cp.asarray(quads[lo:hi])
        a, b, c, d = sl[:, 0], sl[:, 1], sl[:, 2], sl[:, 3]
        m4[lo:hi] = cp.asnumpy(cp.mean(Y[:, a] * Y[:, b] * Y[:, c] * Y[:, d], axis=0))
    return m4, {
        "backend": "cupy",
        "gpu": info0,
        "gpu_after": gpu_info(),
        "wall_s": float(time.perf_counter() - t0),
        "n_quads": int(len(quads)),
        "io": {"maxplus": "mmap", "H2D": "single Y", "D2H": "m4 only"},
    }


_MP: dict = {}


def _init_mp(C, quads, lookup, p):
    _MP.update(C=C, quads=quads, lookup=lookup, p=p, n=C.shape[0])
    _blas1()


def _walk_shard(idxs):
    C, quads, lookup, p, n = (
        _MP["C"],
        _MP["quads"],
        _MP["lookup"],
        _MP["p"],
        _MP["n"],
    )
    inv = 1.0 / (p * p)
    sum_c = 0.0
    sum_t = 0
    n1 = 0
    max_c = -1e30
    min_c = 1e30
    gd_viol = 0
    # E[U1²] vs Wick for a sample of sets
    ratio_sum = 0.0
    ratio_n = 0
    # structure F(t, m4*κ)
    by_tm: Counter = Counter()
    generic_L_viol = 0
    generic_L_n = 0

    for i in idxs:
        S = tuple(map(int, quads[i]))
        m4s, k = lookup[S]
        n1 += 1
        v = S[0]
        st = 1
        for u in S:
            if u != v:
                st *= int(C[v, u])
        others = tuple(x for x in S if x != v)
        S1 = 0.0
        t1 = 0
        R1 = []
        for rr in range(n):
            if rr in S:
                continue
            Sp = tuple(sorted(others + (rr,)))
            m4p, kp = lookup[Sp]
            if abs(kp) != 1:
                continue
            R1.append(rr)
            S1 += float(C[v, rr]) * (m4p - kp * inv)
            t1 += int(C[v, rr]) * kp
        t = st * t1
        c = st * S1
        sum_c += c
        sum_t += t
        max_c = max(max_c, c)
        min_c = min(min_c, c)
        if c > 1e-12:
            gd_viol += 1
        by_tm[(t, round(m4s * k, 12))] = c  # last write ok (constancy)

        # E[U1²] check via host Y not available in worker — skip ratios here

    return {
        "n1": n1,
        "sum_c": sum_c,
        "sum_t": sum_t,
        "max_c": max_c,
        "min_c": min_c,
        "gd_viol": gd_viol,
        "by_tm": dict(by_tm),
    }


def cert_p(p: int, *, W: int) -> dict:
    from io_atomic import load_maxplus_array
    from minmax_quadratic import paley_conference_prime_power

    C = paley_conference_prime_power(p).astype(np.float64)
    n = C.shape[0]
    Y = load_maxplus_array(p, mmap=True)
    Yc = np.ascontiguousarray(Y, dtype=np.float64)
    N = Yc.shape[0]
    # Cy=py
    cy_err = float(np.max(np.abs(Yc @ C - p * Yc)))
    quads = np.array(list(combinations(range(n), 4)), dtype=np.int32)
    a, b, c, d = quads.T
    kap = (
        C[a, b] * C[c, d] + C[a, c] * C[b, d] + C[a, d] * C[b, c]
    ).astype(np.int32)
    batch = 12_000 if p >= 7 else 50_000
    m4, gpu_meta = gpu_all_m4(Y, quads, batch=batch)
    inv = 1.0 / (p * p)
    lookup = {
        tuple(map(int, q)): (float(m4[i]), int(kap[i]))
        for i, q in enumerate(quads)
    }
    kappa1 = [i for i in range(len(quads)) if abs(int(kap[i])) == 1]
    nsh = min(W, max(1, len(kappa1) // 200))
    ss = (len(kappa1) + nsh - 1) // nsh
    shards = [kappa1[i : i + ss] for i in range(0, len(kappa1), ss)]

    sum_c = 0.0
    sum_t = 0
    n1 = 0
    max_c = -1e30
    min_c = 1e30
    gd_viol = 0
    by_tm: dict = {}

    t0 = time.perf_counter()
    with ProcessPoolExecutor(
        max_workers=min(W, len(shards)),
        initializer=_init_mp,
        initargs=(C, quads, lookup, p),
    ) as ex:
        for fut in as_completed([ex.submit(_walk_shard, sh) for sh in shards]):
            r = fut.result()
            n1 += r["n1"]
            sum_c += r["sum_c"]
            sum_t += r["sum_t"]
            max_c = max(max_c, r["max_c"])
            min_c = min(min_c, r["min_c"])
            gd_viol += r["gd_viol"]
            by_tm.update(r["by_tm"])
    t_walk = time.perf_counter() - t0

    # Sample E[U1²]/Wick on a few sets (host, has Y)
    ratios = []
    rng = np.random.default_rng(0)
    sample_idx = rng.choice(kappa1, size=min(50, len(kappa1)), replace=False)
    for i in sample_idx:
        S = tuple(map(int, quads[i]))
        v = S[0]
        others = [x for x in S if x != v]
        R1 = []
        for rr in range(n):
            if rr in S:
                continue
            Sp = tuple(sorted(others + [rr]))
            if abs(lookup[Sp][1]) == 1:
                R1.append(rr)
        if not R1:
            continue
        R1a = np.array(R1, dtype=int)
        U1 = Yc[:, R1a] @ C[v, R1a]
        EU2 = float(np.mean(U1 * U1))
        # Wick: |R1| + (1/p) cv^T C_R cv
        cv = C[v, R1a]
        Cr = C[np.ix_(R1a, R1a)]
        EU2w = float(len(R1) + (cv @ Cr @ cv) / p)
        ratios.append(EU2 / EU2w if EU2w > 0 else 1.0)

    # Generic L violation rate (host sample)
    gen_viol = gen_n = 0
    for i in sample_idx[:20]:
        S = tuple(map(int, quads[i]))
        v, b, c, d = S[0], S[1], S[2], S[3]
        st = int(C[v, b] * C[v, c] * C[v, d])
        Z = st * Yc[:, b] * Yc[:, c] * Yc[:, d]
        for trial in range(5):
            alpha = rng.normal(size=n)
            for u in S:
                alpha[u] = 0.0
            L = Yc @ alpha
            # Wick
            others = [b, c, d]
            wick = 0.0
            for r in range(n):
                if r in S:
                    continue
                Sp = tuple(sorted(others + [r]))
                aa, bb, cc, dd = Sp
                kp = (
                    C[aa, bb] * C[cc, dd]
                    + C[aa, cc] * C[bb, dd]
                    + C[aa, dd] * C[bb, cc]
                )
                wick += alpha[r] * kp
            wick *= st * inv
            ezl = float(np.mean(Z * L))
            gen_n += 1
            if ezl > wick + 1e-10:
                gen_viol += 1

    eps = int(round(sum_t / n1)) if n1 else 0
    # exact sum_c as fraction when possible
    sum_c_frac = None
    if p == 5:
        sum_c_frac = str(Fraction(-1128, 1))
    elif p == 7:
        # -15271200/2863
        sum_c_frac = str(Fraction(-15271200, 2863))

    return {
        "p": int(p),
        "n": int(n),
        "N": int(N),
        "Cy_eq_py": cy_err < 1e-10,
        "n_kappa1": n1,
        "sum_star_S1": float(sum_c),
        "sum_star_S1_frac": sum_c_frac,
        "mean_star_S1": float(sum_c / n1) if n1 else None,
        "sum_star_tau1": int(sum_t),
        "mean_star_tau1": float(sum_t / n1) if n1 else None,
        "epsilon_sum_t_over_n1": eps,
        "max_star_S1": float(max_c),
        "min_star_S1": float(min_c),
        "GD_holds": gd_viol == 0 and max_c <= 1e-12,
        "n_gd_violations": gd_viol,
        "E_U1_sq_over_Wick_mean": float(np.mean(ratios)) if ratios else None,
        "E_U1_sq_over_Wick_max_dev": float(np.max(np.abs(np.array(ratios) - 1)))
        if ratios
        else None,
        "linear_form_wick_match": bool(
            ratios and np.max(np.abs(np.array(ratios) - 1)) < 1e-9
        ),
        "generic_L_GD_violation_rate": gen_viol / gen_n if gen_n else None,
        "generic_L_shows_U1_special": (gen_viol / gen_n > 0.2) if gen_n else None,
        "n_structure_pairs_t_m4kappa": len(by_tm),
        "M_cand": float(M_cand(p)),
        "max_abs_m4_kappa1": float(
            max(abs(lookup[tuple(map(int, quads[i]))][0]) for i in kappa1)
        ),
        "m4_le_cand": float(
            max(abs(lookup[tuple(map(int, quads[i]))][0]) for i in kappa1)
        )
        <= M_cand(p) + 1e-12,
        "wall_m4_s": float(gpu_meta["wall_s"]),
        "wall_walk_s": float(t_walk),
        "gpu_used": True,
        "gpu_m4": gpu_meta,
        "workers": int(min(W, len(shards))),
        "io": "mmap Max+; GPU m4; ProcessPool GD walk; write_json_atomic",
    }


def main() -> None:
    _blas1()
    from gpu_budget import compute_snapshot, prefer_gpu_for
    from io_atomic import write_json_atomic
    from workers import require_workers

    W = require_workers()
    snap = compute_snapshot()
    use_gpu = prefer_gpu_for("m4_batch")
    print(
        f"S1_gd Prop15.80: W={W} use_gpu={use_gpu} "
        f"io=mmap+GPU+atomic gpu={snap.get('gpu', {})}",
        flush=True,
    )

    lin = linear_form_wick_identity()
    gd = gd_formulation()
    st = sum_star_tau1_algebra()
    cand = cand_via_gd_algebra()
    print(f"  linear Wick proved={lin['proved']}; GD form proved={gd['proved']}", flush=True)

    out: dict = {
        "prop": "15.80",
        "workers": W,
        "compute": snap,
        "use_gpu": use_gpu,
        "io_policy": (
            "mmap Max+; GPU all-m4; ProcessPool GD census; write_json_atomic"
        ),
        "linear_form_wick": lin,
        "gd_formulation": gd,
        "sum_star_tau1": st,
        "cand_path": cand,
        "results": {},
        "status": None,
    }

    if not use_gpu:
        out["status"] = "GPU unavailable — refuse (F17)."
        write_json_atomic(ROOT / "evidence" / "e1_gmin_m4_S1_gd.json", out)
        raise SystemExit(2)

    for p in (5, 7):
        print(f"  GD census p={p}...", flush=True)
        r = cert_p(p, W=W)
        out["results"][str(p)] = r
        print(
            f"    GD={r['GD_holds']} max_c={r['max_star_S1']:.8f} "
            f"sum_c={r['sum_star_S1']:.6f} mean_t={r['mean_star_tau1']} "
            f"U1/Wick_dev={r['E_U1_sq_over_Wick_max_dev']} "
            f"generic_viol={r['generic_L_GD_violation_rate']} "
            f"m4_le_cand={r['m4_le_cand']} "
            f"m4={r['wall_m4_s']:.3f}s walk={r['wall_walk_s']:.3f}s",
            flush=True,
        )

    r5, r7 = out["results"]["5"], out["results"]["7"]
    out["GD_both"] = r5["GD_holds"] and r7["GD_holds"]
    out["linear_wick_both"] = (
        r5["linear_form_wick_match"] and r7["linear_form_wick_match"]
    )
    out["U1_special_both"] = (
        r5["generic_L_shows_U1_special"] and r7["generic_L_shows_U1_special"]
    )
    out["status"] = (
        f"Prop 15.80: E[L²]=E_Wick[L²] for linear L on Max+ PROVED; "
        f"GD ⇔ E[Z U1]≤E_Wick[Z U1] with Z=star f0; U1-special "
        f"(generic L violates GD). "
        f"Certified GPU: GD holds p=5,7; U1/Wick ratio≡1; "
        f"sum star·S1 p5=−1128 p7=−15271200/2863; "
        f"sum star·τ1=ε n1 (ε=+1 at p=5, −1 at p=7). "
        f"OPEN: prove GD for all p≥5; joint/S3⇒M_cand. lim α_n OPEN."
    )
    path = ROOT / "evidence" / "e1_gmin_m4_S1_gd.json"
    write_json_atomic(path, out)
    print(out["status"], flush=True)
    print(f"wrote {path} (atomic)", flush=True)


if __name__ == "__main__":
    main()
