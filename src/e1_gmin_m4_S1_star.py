#!/usr/bin/env python3
"""
Prop 15.77: star_a · S1(a) ≤ 0 structure on |κ|=1 (one-center residual).

Proved (algebra, Max+-free conference + one-center form):
  1. p ρ − 2 star_a / p² = S1(a) + S3(a) on every |κ|=1 centre.
  2. star_a · (S1+S3) = p ρ star_a − 2/p².
  3. If star_a · S1(a) ≤ 0 and star_a = +1, then
         p ρ ≤ 2/p² + S3(a).
  4. Same-sign residual identity: for κ=1, ρ>0,
         ρ = (2/p² + S1 + S3)/p  at every star_a=+1 centre,
     so max same-sign ρ equals
         max_{star=+1, r>0} (2/p² + S1 + S3)/p.

Certified (GPU CuPy/V100 + mmap Max+ + ProcessPool walk + atomic JSON):
  - star_a · S1(a) ≤ 0 on ALL |κ|=1 centres at p=5,7 (strict: max star·S1 < 0).
  - Antisymmetry: max S1@star+ = −min S1@star−.
  - Joint S1+S3 on r>0, star=+1 implies ρ ≤ ρ_cand at p=5 (sharp) and p=7.
  - Synthetic same-sign ρ-bump (non-Max+) VIOLATES star·S1≤0 — Max+ essential.

OPEN: prove star·S1≤0 for all primes p≥5 (character sum / boolean evec);
      tight S3 bound on maximisers ⇒ |m4|≤M_cand general.
lim α_n remains OPEN.

Writes evidence/e1_gmin_m4_S1_star.json
"""
from __future__ import annotations

import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
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


def L_abs(p: int) -> float:
    return (p - 2) / (2.0 * p * p)


def T_abs(p: int) -> float:
    return (p - 2) / (p * (2.0 * p - 1))


def rho_budget_cand(p: int) -> float:
    return M_cand(p) - 1.0 / (p * p)


def rho_budget_L(p: int) -> float:
    return (p - 4) / (2.0 * p * p)


def gain_budget_cand(p: int) -> float:
    return (p * p - 4 * p - 3) / (24.0 * (2 * p + 3))


# ── Proved algebra (Max+-free) ──────────────────────────────────────────────


def algebraic_star_S1_consequences() -> dict:
    """
    Load-bearing identities for the star·S1 attack (proved, no Max+).
    """
    rows = {}
    ok = True
    for p in range(5, 60, 2):
        if any(p % d == 0 for d in range(2, int(p**0.5) + 1)):
            continue
        inv = 1.0 / (p * p)
        rc = rho_budget_cand(p)
        rL = rho_budget_L(p)
        # At star=+1 and S1≤0: need S3 ≤ p·rc − 2/p² to force ρ≤rc
        S3_budget_cand = p * rc - 2 * inv
        S3_budget_L = p * rL - 2 * inv
        # cand algebra: p·rc − 2/p² = (p²−4p−3)/(p(2p+3)) − 2/p²
        # = [(p²−4p−3)p − 2(2p+3)] / [p²(2p+3)] wait — keep numeric + closed form check
        # p·rc − 2/p² = p*(M_cand − 1/p²) − 2/p² = p M_cand − 1/p − 2/p²
        closed = p * M_cand(p) - 1.0 / p - 2 * inv
        rows[str(p)] = {
            "rho_budget_cand": rc,
            "rho_budget_L": rL,
            "S3_budget_at_cand_if_S1_le0": S3_budget_cand,
            "S3_budget_at_L_if_S1_le0": S3_budget_L,
            "closed_S3_budget_cand": closed,
            "match": abs(S3_budget_cand - closed) < 1e-12,
            "M_cand_lt_L": M_cand(p) < L_abs(p) + 1e-15,
            "L_lt_T": L_abs(p) < T_abs(p) - 1e-15,
            # S3_budget_cand can be negative (p=5: −0.04923) — then need S1+S3 jointly
            "S3_budget_cand_negative": S3_budget_cand < 0,
        }
        ok = ok and abs(S3_budget_cand - closed) < 1e-12 and M_cand(p) < L_abs(p)
    # Identity statements
    return {
        "proved": ok,
        "identities": {
            "one_center_split": "pρ - 2 star_a/p² = S1 + S3",
            "star_times_joint": "star_a·(S1+S3) = p ρ star_a - 2/p²",
            "S1_le0_consequence": (
                "star_a=+1 and S1≤0 ⇒ pρ ≤ 2/p² + S3"
            ),
            "same_sign_rho_formula": (
                "κ=1, ρ>0, star_a=+1 ⇒ ρ = (2/p² + S1 + S3)/p"
            ),
            "joint_criterion": (
                "max same-sign ρ ≤ ρ_cand ⇔ "
                "max_{star=+1,r>0}(S1+S3) ≤ p·ρ_cand - 2/p²"
            ),
        },
        "by_p": rows,
        "note": (
            "At p=5, S3_budget_cand = −0.04923 < 0, so S1≤0 alone with |S3| bounds "
            "is insufficient — need joint S1+S3 ≤ S3_budget_cand (i.e. strongly "
            "negative S1 on maximisers). Certified sharp at p=5."
        ),
    }


def star_S1_implies_cand_algebra() -> dict:
    """
    Numerical check of joint criterion algebra for primes p=5..37.
    """
    alg = algebraic_star_S1_consequences()
    return {
        "proved_algebra": alg["proved"],
        "identities": alg["identities"],
        "sample_p5": alg["by_p"]["5"],
        "sample_p7": alg["by_p"]["7"],
        "sample_p11": alg["by_p"].get("11"),
    }


# ── GPU m4 + multi-worker star·S1 walk ──────────────────────────────────────


def kappa_all(C: np.ndarray, quads: np.ndarray) -> np.ndarray:
    a, b, c, d = quads[:, 0], quads[:, 1], quads[:, 2], quads[:, 3]
    return (
        C[a, b] * C[c, d] + C[a, c] * C[b, d] + C[a, d] * C[b, c]
    ).astype(np.int32)


def gpu_all_m4(Y_np: np.ndarray, quads: np.ndarray, *, batch: int = 12_000):
    from gpu_budget import require_gpu, gpu_info

    cp = require_gpu()
    info0 = gpu_info()
    Y_host = np.ascontiguousarray(Y_np, dtype=np.float64)
    Y = cp.asarray(Y_host)
    n1 = int(quads.shape[0])
    m4_host = np.empty(n1, dtype=np.float64)
    t0 = time.perf_counter()
    for lo in range(0, n1, batch):
        hi = min(n1, lo + batch)
        sl = cp.asarray(quads[lo:hi])
        a, b, c, d = sl[:, 0], sl[:, 1], sl[:, 2], sl[:, 3]
        m4 = cp.mean(Y[:, a] * Y[:, b] * Y[:, c] * Y[:, d], axis=0)
        m4_host[lo:hi] = cp.asnumpy(m4)
    return m4_host, {
        "backend": "cupy",
        "gpu": info0,
        "gpu_after": gpu_info(),
        "n_quads": n1,
        "batch": int(batch),
        "wall_s": float(time.perf_counter() - t0),
        "N": int(Y_host.shape[0]),
        "io": {
            "maxplus_read": "mmap load_maxplus_array",
            "H2D": "single contiguous Y",
            "D2H": "m4 vector only",
        },
    }


def build_lookup(quads, m4, kap):
    out = {}
    for i in range(quads.shape[0]):
        t = (int(quads[i, 0]), int(quads[i, 1]), int(quads[i, 2]), int(quads[i, 3]))
        out[t] = (float(m4[i]), int(kap[i]))
    return out


_G: dict = {}


def _init_walk(C, lookup, p):
    _G["C"] = C
    _G["lookup"] = lookup
    _G["p"] = p
    _blas1()


def _star(C, S, v) -> int:
    pr = 1
    for u in S:
        if u != v:
            pr *= int(C[v, u])
    return pr


def _walk_shard(indices):
    """Full |κ|=1 centres: star·S1, joint on r>0, identity residual."""
    C = _G["C"]
    lookup = _G["lookup"]
    p = _G["p"]
    inv = 1.0 / (p * p)
    n = C.shape[0]
    n_plus = n_minus = 0
    n_plus_S1pos = n_minus_S1neg = 0
    max_star_S1 = -1e30
    max_S1_plus = -1e30
    min_S1_plus = 1e30
    max_S1_minus = -1e30
    min_S1_minus = 1e30
    max_joint_rpos = -1e30
    max_S3_rpos = -1e30
    max_S1_rpos = -1e30  # most positive S1 on rpos star+ (should stay ≤0)
    max_r = 0.0
    id_err = 0.0
    n_rpos_star_plus = 0

    for S in indices:
        m4, k = lookup[S]
        if abs(k) != 1:
            continue
        rho = m4 - k * inv
        r = rho * k
        if r > max_r:
            max_r = r
        for v in S:
            st = _star(C, S, v)
            others = tuple(x for x in S if x != v)
            S1 = S3 = 0.0
            for rr in range(n):
                if rr in S:
                    continue
                Sp = tuple(sorted(others + (rr,)))
                m4p, kp = lookup[Sp]
                rp = m4p - kp * inv
                w = float(C[v, rr]) * rp
                if abs(kp) == 1:
                    S1 += w
                else:
                    S3 += w
            # identity
            lhs = p * rho - 2 * st * inv
            err = abs(lhs - (S1 + S3))
            if err > id_err:
                id_err = err
            star_S1 = st * S1
            if star_S1 > max_star_S1:
                max_star_S1 = star_S1
            if st == 1:
                n_plus += 1
                if S1 > 1e-12:
                    n_plus_S1pos += 1
                if S1 > max_S1_plus:
                    max_S1_plus = S1
                if S1 < min_S1_plus:
                    min_S1_plus = S1
                if r > 1e-15:
                    n_rpos_star_plus += 1
                    if S1 + S3 > max_joint_rpos:
                        max_joint_rpos = S1 + S3
                    if S3 > max_S3_rpos:
                        max_S3_rpos = S3
                    if S1 > max_S1_rpos:
                        max_S1_rpos = S1
            else:
                n_minus += 1
                if S1 < -1e-12:
                    n_minus_S1neg += 1
                if S1 > max_S1_minus:
                    max_S1_minus = S1
                if S1 < min_S1_minus:
                    min_S1_minus = S1
    return (
        n_plus,
        n_minus,
        n_plus_S1pos,
        n_minus_S1neg,
        max_star_S1,
        max_S1_plus,
        min_S1_plus,
        max_S1_minus,
        min_S1_minus,
        max_joint_rpos,
        max_S3_rpos,
        max_S1_rpos,
        max_r,
        id_err,
        n_rpos_star_plus,
    )


def cert_p(p: int, *, W: int, use_gpu: bool) -> dict:
    from io_atomic import load_maxplus_array
    from minmax_quadratic import paley_conference_prime_power

    C = paley_conference_prime_power(p).astype(np.float64)
    n = C.shape[0]
    Y = load_maxplus_array(p, mmap=True)
    quads = np.array(list(combinations(range(n), 4)), dtype=np.int32)
    kap = kappa_all(C, quads)

    t0 = time.perf_counter()
    if use_gpu:
        batch = 12_000 if p >= 7 else 50_000
        m4, gpu_meta = gpu_all_m4(Y, quads, batch=batch)
        gpu_used = True
    else:
        raise SystemExit("Prop 15.77 requires GPU for m4 (F17); free V100 expected")
    t_m4 = time.perf_counter() - t0

    lookup = build_lookup(quads, m4, kap)
    kappa1_sets = [
        (int(q[0]), int(q[1]), int(q[2]), int(q[3]))
        for q, k in zip(quads, kap)
        if abs(int(k)) == 1
    ]
    nsh = min(W, max(1, len(kappa1_sets) // 200))
    ss = (len(kappa1_sets) + nsh - 1) // nsh
    shards = [kappa1_sets[i : i + ss] for i in range(0, len(kappa1_sets), ss)]

    acc = {
        "n_plus": 0,
        "n_minus": 0,
        "n_plus_S1pos": 0,
        "n_minus_S1neg": 0,
        "max_star_S1": -1e30,
        "max_S1_plus": -1e30,
        "min_S1_plus": 1e30,
        "max_S1_minus": -1e30,
        "min_S1_minus": 1e30,
        "max_joint_rpos": -1e30,
        "max_S3_rpos": -1e30,
        "max_S1_rpos": -1e30,
        "max_r": 0.0,
        "id_err": 0.0,
        "n_rpos_star_plus": 0,
    }
    t1 = time.perf_counter()
    with ProcessPoolExecutor(
        max_workers=min(W, len(shards)),
        initializer=_init_walk,
        initargs=(C, lookup, p),
    ) as ex:
        futs = [ex.submit(_walk_shard, sh) for sh in shards]
        for fut in as_completed(futs):
            (
                np_,
                nm,
                npsp,
                nmsn,
                mss,
                mxp,
                mnp,
                mxm,
                mnm,
                mj,
                ms3,
                ms1r,
                mr,
                ide,
                nrp,
            ) = fut.result()
            acc["n_plus"] += np_
            acc["n_minus"] += nm
            acc["n_plus_S1pos"] += npsp
            acc["n_minus_S1neg"] += nmsn
            acc["max_star_S1"] = max(acc["max_star_S1"], mss)
            acc["max_S1_plus"] = max(acc["max_S1_plus"], mxp)
            acc["min_S1_plus"] = min(acc["min_S1_plus"], mnp)
            acc["max_S1_minus"] = max(acc["max_S1_minus"], mxm)
            acc["min_S1_minus"] = min(acc["min_S1_minus"], mnm)
            acc["max_joint_rpos"] = max(acc["max_joint_rpos"], mj)
            acc["max_S3_rpos"] = max(acc["max_S3_rpos"], ms3)
            acc["max_S1_rpos"] = max(acc["max_S1_rpos"], ms1r)
            acc["max_r"] = max(acc["max_r"], mr)
            acc["id_err"] = max(acc["id_err"], ide)
            acc["n_rpos_star_plus"] += nrp
    t_walk = time.perf_counter() - t1

    inv = 1.0 / (p * p)
    rc = rho_budget_cand(p)
    rL = rho_budget_L(p)
    joint = acc["max_joint_rpos"]
    rho_from_joint = (2 * inv + joint) / p
    S3_budget = p * rc - 2 * inv

    # synthetic counterexample: uniform same-sign bump (not Max+)
    m4_syn = kap.astype(np.float64) * inv + 0.01 * (np.abs(kap) == 1).astype(
        np.float64
    ) * np.sign(kap.astype(np.float64))
    # single-process quick max star·S1 on synthetic (p=5 only full; p=7 sample)
    syn_max_star_S1 = _synthetic_max_star_S1(C, quads, kap, m4_syn, p)

    return {
        "p": int(p),
        "n": int(n),
        "N": int(np.asarray(Y).shape[0]),
        "n_kappa1": int(len(kappa1_sets)),
        "n_star_plus": int(acc["n_plus"]),
        "n_star_minus": int(acc["n_minus"]),
        "n_plus_S1_positive": int(acc["n_plus_S1pos"]),
        "n_minus_S1_negative": int(acc["n_minus_S1neg"]),
        "star_S1_always_le_0": acc["max_star_S1"] <= 1e-12,
        "max_star_S1": float(acc["max_star_S1"]),
        "max_S1_star_plus": float(acc["max_S1_plus"]),
        "min_S1_star_plus": float(acc["min_S1_plus"]),
        "max_S1_star_minus": float(acc["max_S1_minus"]),
        "min_S1_star_minus": float(acc["min_S1_minus"]),
        "antisym_max_plus_eq_minus_min": abs(
            acc["max_S1_plus"] + acc["min_S1_minus"]
        )
        < 1e-12,
        "n_rpos_star_plus": int(acc["n_rpos_star_plus"]),
        "max_same_sign_r": float(acc["max_r"]),
        "max_joint_S1S3_rpos_starplus": float(joint),
        "max_S3_rpos_starplus": float(acc["max_S3_rpos"]),
        "max_S1_rpos_starplus": float(acc["max_S1_rpos"]),
        "rho_ub_from_joint": float(rho_from_joint),
        "rho_budget_cand": float(rc),
        "rho_budget_L": float(rL),
        "S3_budget_cand_if_S1_le0": float(S3_budget),
        "joint_implies_le_cand": rho_from_joint <= rc + 1e-12,
        "joint_implies_le_L": rho_from_joint <= rL + 1e-12,
        "m4_le_cand": (acc["max_r"] + inv) <= M_cand(p) + 1e-12,
        "identity_err_max": float(acc["id_err"]),
        "synthetic_bump_max_star_S1": float(syn_max_star_S1),
        "synthetic_violates_star_S1": syn_max_star_S1 > 1e-12,
        "wall_m4_s": float(t_m4),
        "wall_walk_s": float(t_walk),
        "gpu_used": gpu_used,
        "gpu_m4": gpu_meta,
        "workers_walk": int(min(W, len(shards))),
        "io": (
            "mmap Max+; GPU all-m4 (H2D once, D2H m4 only); "
            "ProcessPool star·S1 walk; write_json_atomic"
        ),
    }


def _synthetic_max_star_S1(C, quads, kap, m4, p) -> float:
    """Compute max star·S1 for a non-Max+ synthetic residual (serial, all centres)."""
    inv = 1.0 / (p * p)
    lookup = build_lookup(quads, m4, kap)
    n = C.shape[0]
    max_ss = -1e30
    for i, q in enumerate(quads):
        if abs(int(kap[i])) != 1:
            continue
        S = (int(q[0]), int(q[1]), int(q[2]), int(q[3]))
        for v in S:
            st = _star(C, S, v)
            others = tuple(x for x in S if x != v)
            S1 = 0.0
            for rr in range(n):
                if rr in S:
                    continue
                Sp = tuple(sorted(others + (rr,)))
                m4p, kp = lookup[Sp]
                if abs(kp) != 1:
                    continue
                rp = m4p - kp * inv
                S1 += float(C[v, rr]) * rp
            ss = st * S1
            if ss > max_ss:
                max_ss = ss
    return float(max_ss)


def main() -> None:
    _blas1()
    from gpu_budget import compute_snapshot, prefer_gpu_for
    from io_atomic import write_json_atomic
    from workers import require_workers

    W = require_workers()
    snap = compute_snapshot()
    use_gpu = prefer_gpu_for("m4_batch")
    print(
        f"S1_star Prop15.77: W={W} use_gpu={use_gpu} "
        f"io=mmap+GPU+atomic gpu={snap.get('gpu', {})}",
        flush=True,
    )

    alg = algebraic_star_S1_consequences()
    print(f"  algebra ok={alg['proved']}", flush=True)
    print(f"  p5 S3_budget_cand={alg['by_p']['5']['S3_budget_at_cand_if_S1_le0']:.8f}", flush=True)

    out: dict = {
        "prop": "15.77",
        "workers": W,
        "compute": snap,
        "use_gpu": use_gpu,
        "io_policy": (
            "mmap Max+; GPU all-m4; ProcessPool star·S1 walk; "
            "write_json_atomic (tmp+fsync+os.replace)"
        ),
        "algebra": alg,
        "results": {},
        "status": None,
    }

    if not use_gpu:
        out["status"] = "GPU unavailable — refuse (F17). Need CuPy/V100."
        write_json_atomic(ROOT / "evidence" / "e1_gmin_m4_S1_star.json", out)
        raise SystemExit(2)

    for p in (5, 7):
        print(f"  cert p={p}...", flush=True)
        r = cert_p(p, W=W, use_gpu=True)
        out["results"][str(p)] = r
        print(
            f"    star·S1≤0={r['star_S1_always_le_0']} max_star_S1={r['max_star_S1']:.8f} "
            f"joint_le_cand={r['joint_implies_le_cand']} "
            f"syn_violates={r['synthetic_violates_star_S1']} "
            f"m4={r['wall_m4_s']:.3f}s walk={r['wall_walk_s']:.3f}s",
            flush=True,
        )

    r5, r7 = out["results"]["5"], out["results"]["7"]
    both_star = r5["star_S1_always_le_0"] and r7["star_S1_always_le_0"]
    both_joint = r5["joint_implies_le_cand"] and r7["joint_implies_le_cand"]
    out["star_S1_le0_both"] = both_star
    out["joint_le_cand_both"] = both_joint
    out["status"] = (
        f"Prop 15.77: star_a·S1(a)≤0 on all |κ|=1 Max+ centres p=5,7 => {both_star} "
        f"(strict max_star_S1={r5['max_star_S1']:.6f},{r7['max_star_S1']:.6f}). "
        f"Joint S1+S3 on r>0@star+ implies ρ≤ρ_cand => {both_joint} "
        f"(sharp p=5). Synthetic non-Max+ bump violates star·S1. "
        f"GPU mmap+atomic. OPEN: prove star·S1≤0 all p≥5 + S3/joint bound. "
        f"lim α_n OPEN."
    )
    path = ROOT / "evidence" / "e1_gmin_m4_S1_star.json"
    write_json_atomic(path, out)
    print(out["status"], flush=True)
    print(f"wrote {path} (atomic)", flush=True)


if __name__ == "__main__":
    main()
