#!/usr/bin/env python3
"""
Prop 15.76 support: one-center extension degrees + residual bound attack.

Proved/certified:
  On every |κ|=1 4-set S and every a∈S (Paley), the one-vertex extension counts
    d1^(1) = (3p²-7)/4,   d3^(1) = (p²-5)/4
  are constant (equal to full d1/4, d3/4). Combined with σ_a=2·star_a (Prop 15.75)
  and one-center residual form.

Bound attempt:
  pρ - 2 star_a/p² = S1 + S3  (signed sums of ρ to κ1/κ3)
  |S3| ≤ d3^(1) · R3,  |S1| ≤ d1^(1) · R1
  with R1=max|ρ| on κ1, R3 on κ3.
  Absolute form does not contract (4p < d1), but one-center + star choice
  yields identities used in cand analysis.

GPU: optional Max+ residual moments via mmap+CuPy when certifying.
Atomic evidence write.

Writes evidence/e1_gmin_m4_onecenter_deg.json
"""
from __future__ import annotations

import os
import sys
import time
from collections import Counter
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


def L_abs(p: int) -> float:
    return (p - 2) / (2.0 * p * p)


def M_mid(p: int) -> float:
    return (p - 2) / (2.0 * p * (p + 1))


def M_cand(p: int) -> float:
    return (p - 2) / (p * (2.0 * p + 3))


def d1_one(p: int) -> int:
    return (3 * p * p - 7) // 4


def d3_one(p: int) -> int:
    return (p * p - 5) // 4


def d1_full(p: int) -> int:
    return 3 * p * p - 7


def d3_full(p: int) -> int:
    return p * p - 5


def algebra_degrees() -> dict:
    rows = {}
    ok = True
    for p in range(3, 80, 2):
        if any(p % d == 0 for d in range(2, int(p**0.5) + 1)):
            continue
        d1o, d3o = d1_one(p), d3_one(p)
        d1f, d3f = d1_full(p), d3_full(p)
        # divisibility
        div_ok = (3 * p * p - 7) % 4 == 0 and (p * p - 5) % 4 == 0
        sum_ok = d1o + d3o == p * p - 3  # n-4
        four_ok = 4 * d1o == d1f and 4 * d3o == d3f
        cand, mid, La = M_cand(p), M_mid(p), L_abs(p)
        T_abs = (p - 2) / (p * (2 * p - 1))
        # bound attempt: if R3 ≤ B3 and we had contraction...
        # document absolute fail: 4p - d1_full < 0
        gap = 4 * p - d1f
        rows[str(p)] = {
            "d1_one": d1o,
            "d3_one": d3o,
            "d1_full": d1f,
            "d3_full": d3f,
            "divisible_by_4": div_ok,
            "d1_one_plus_d3_one_eq_n_minus_4": sum_ok,
            "four_times_one_eq_full": four_ok,
            "4p_minus_d1_full": gap,
            "abs_bootstrap_contracts": gap > 0,
            "M_cand": cand,
            "M_mid": mid,
            "L_abs": La,
            "T_abs": T_abs,
            "cand_lt_T": cand < T_abs,
        }
        ok = ok and div_ok and sum_ok and four_ok and cand < T_abs
    return {
        "d1_one_closed": "(3p^2-7)/4",
        "d3_one_closed": "(p^2-5)/4",
        "identity_4_d1_one_eq_d1": True,
        "identity_d1_one_plus_d3_one": "n-4 = p^2-3",
        "by_p": rows,
        "proved_divisibility_algebra": ok,
        "note": (
            "Divisibility of 3p^2-7 and p^2-5 by 4 for odd p is elementary. "
            "Constancy of one-center degrees on Paley |κ|=1 certified multi-worker."
        ),
    }


def kappa(C, S) -> int:
    a, b, c, d = S
    return int(
        C[a, b] * C[c, d] + C[a, c] * C[b, d] + C[a, d] * C[b, c]
    )


_G: dict = {}


def _init_C(C: np.ndarray) -> None:
    _G["C"] = C
    _blas1()


def _deg_shard(quads: list) -> dict:
    C = _G["C"]
    n = C.shape[0]
    d1s: Counter = Counter()
    d3s: Counter = Counter()
    n1 = 0
    for S in quads:
        if abs(kappa(C, S)) != 1:
            continue
        n1 += 1
        Sset = set(S)
        for a in S:
            d1 = d3 = 0
            others = tuple(sorted(x for x in S if x != a))
            for r in range(n):
                if r in Sset:
                    continue
                Sp = tuple(sorted(others + (r,)))
                ak = abs(kappa(C, Sp))
                if ak == 1:
                    d1 += 1
                elif ak == 3:
                    d3 += 1
            d1s[d1] += 1
            d3s[d3] += 1
    return {"n1": n1, "d1s": dict(d1s), "d3s": dict(d3s)}


def certify_one_center_degrees(p: int, W: int) -> dict:
    from minmax_quadratic import paley_conference_prime_power

    t0 = time.perf_counter()
    C = paley_conference_prime_power(p).astype(np.float64)
    n = C.shape[0]
    quads = list(combinations(range(n), 4))
    n_shards = min(W, max(1, len(quads) // 400))
    ss = (len(quads) + n_shards - 1) // n_shards
    shards = [quads[i : i + ss] for i in range(0, len(quads), ss)]
    d1s: Counter = Counter()
    d3s: Counter = Counter()
    n1 = 0
    with ProcessPoolExecutor(
        max_workers=min(W, len(shards)), initializer=_init_C, initargs=(C,)
    ) as ex:
        futs = [ex.submit(_deg_shard, sh) for sh in shards]
        for fut in as_completed(futs):
            r = fut.result()
            n1 += r["n1"]
            d1s.update(r["d1s"])
            d3s.update(r["d3s"])
    e1, e3 = d1_one(p), d3_one(p)
    return {
        "p": int(p),
        "n_kappa1": n1,
        "n_center_checks": n1 * 4,
        "d1_one_unique": sorted(d1s.keys()),
        "d3_one_unique": sorted(d3s.keys()),
        "d1_one_expected": e1,
        "d3_one_expected": e3,
        "d1_constant_match": set(d1s.keys()) == {e1},
        "d3_constant_match": set(d3s.keys()) == {e3},
        "wall_s": float(time.perf_counter() - t0),
    }


def bound_attempt_algebra() -> dict:
    """
    Document the residual identity and what would close cand.

    For κ=1, ρ≥0 same-sign danger:
      pρ - 2 star_a/p² = S1 + S3
    With star_a=±1 and |S3|≤d3_one·R3, |S1|≤d1_one·R1.

    Choosing star_a = +1: pρ - 2/p² = S1+S3
    If S1+S3 ≤ d3_one·R3 + d1_one·R1 (always), no contraction.

    Sufficient for cand if we prove R1 ≤ rho_budget_cand and use energy...
    Record the identity target:
      rho_cand = (p-2)/(p(2p+3)) - 1/p² = (p²-4p-3)/(p²(2p+3))
    """
    rows = {}
    for p in range(5, 40, 2):
        if any(p % d == 0 for d in range(2, int(p**0.5) + 1)):
            continue
        rc = M_cand(p) - 1.0 / (p * p)
        # if S1+S3 could be bounded by α·rc with α < p when star=+1:
        # p rc - 2/p² ≤ α rc  ⇒ need bound on left vs right
        lhs_at_cand = p * rc - 2 / (p * p)
        rows[str(p)] = {
            "rho_budget_cand": rc,
            "lhs_star_plus_at_cand": lhs_at_cand,  # pρ-2/p² at ρ=rc
            "lhs_star_minus_at_cand": p * rc + 2 / (p * p),
            "d3_one": d3_one(p),
            "d1_one": d1_one(p),
            # if all mass on κ3 with |ρ|≤R3: need d3_one*R3 ≥ |lhs|
            "R3_needed_if_only_kappa3_star_plus": abs(lhs_at_cand) / d3_one(p),
            "wick_rho_on_kappa3_scale_3/p2": 3 / (p * p),  # not a bound on ρ
        }
    return {
        "identity_kappa1": "pρ - 2 star_a/p² = S1+S3 (Prop 15.75)",
        "cand_closes_if": (
            "max same-sign ρ ≤ rho_budget_cand, equivalent to |m4|≤M_cand"
        ),
        "by_p": rows,
        "note": (
            "Absolute |S1|+|S3| bounds do not contract. Need signed cancellation "
            "on S1 (κ1→κ1 transfer) or a Max+/design bound on R3 and S1 jointly."
        ),
    }


def gpu_residual_probe(p: int) -> dict | None:
    """GPU Max+ probe: max|ρ| on κ1/κ3 and max |S3| one-center (sample)."""
    try:
        from gpu_budget import require_gpu, gpu_info
        from io_atomic import load_maxplus_array
        from minmax_quadratic import paley_conference_prime_power
    except Exception:
        return None
    try:
        cp = require_gpu()
    except SystemExit:
        return None
    C = paley_conference_prime_power(p).astype(np.float64)
    Y = np.ascontiguousarray(load_maxplus_array(p, mmap=True), dtype=np.float64)
    n = C.shape[0]
    Yg = cp.asarray(Y)
    quads = np.array(list(combinations(range(n), 4)), dtype=np.int32)

    def kap_rows(idx):
        a, b, c, d = idx[:, 0], idx[:, 1], idx[:, 2], idx[:, 3]
        return C[a, b] * C[c, d] + C[a, c] * C[b, d] + C[a, d] * C[b, c]

    kap = kap_rows(quads)
    m4s = np.zeros(len(quads))
    batch = 12_000 if p >= 7 else 50_000
    t0 = time.perf_counter()
    for lo in range(0, len(quads), batch):
        hi = min(len(quads), lo + batch)
        sl = cp.asarray(quads[lo:hi])
        a, b, c, d = sl[:, 0], sl[:, 1], sl[:, 2], sl[:, 3]
        m4 = cp.mean(Yg[:, a] * Yg[:, b] * Yg[:, c] * Yg[:, d], axis=0)
        m4s[lo:hi] = cp.asnumpy(m4)
    wall = time.perf_counter() - t0
    rho = m4s - kap / (p * p)
    mask1 = np.abs(kap) == 1
    mask3 = np.abs(kap) == 3
    # same-sign r on κ1
    r1 = rho[mask1] * kap[mask1]
    return {
        "p": int(p),
        "gpu": gpu_info(),
        "wall_s": float(wall),
        "max_abs_m4_kappa1": float(np.max(np.abs(m4s[mask1]))),
        "max_abs_m4_kappa3": float(np.max(np.abs(m4s[mask3]))),
        "max_abs_rho_kappa1": float(np.max(np.abs(rho[mask1]))),
        "max_abs_rho_kappa3": float(np.max(np.abs(rho[mask3]))),
        "max_same_sign_r_kappa1": float(np.max(r1)),
        "M_cand": M_cand(p),
        "m4_le_cand": float(np.max(np.abs(m4s[mask1]))) <= M_cand(p) + 1e-12,
        "io": "mmap Max+; single H2D; device m4 batches",
        "gpu_used": True,
    }


def main() -> None:
    _blas1()
    sys.path.insert(0, str(ROOT / "src"))
    from gpu_budget import compute_snapshot, prefer_gpu_for
    from io_atomic import write_json_atomic
    from workers import require_workers

    W = require_workers()
    snap = compute_snapshot()
    print(f"onecenter_deg Prop15.76: W={W} gpu={snap['gpu'].get('available')}", flush=True)

    alg = algebra_degrees()
    print(f"  degree algebra ok={alg['proved_divisibility_algebra']}", flush=True)
    bound = bound_attempt_algebra()

    certs = {}
    for p in (3, 5, 7, 11):
        print(f"  one-center degrees p={p}...", flush=True)
        certs[str(p)] = certify_one_center_degrees(p, W)
        r = certs[str(p)]
        print(
            f"    d1={r['d1_one_unique']} exp={r['d1_one_expected']} ok={r['d1_constant_match']} "
            f"d3={r['d3_one_unique']} ok={r['d3_constant_match']} wall={r['wall_s']:.2f}s",
            flush=True,
        )

    gpu_res = {}
    if prefer_gpu_for("m4_batch"):
        for p in (5, 7):
            print(f"  GPU residual probe p={p}...", flush=True)
            gpu_res[str(p)] = gpu_residual_probe(p)
            r = gpu_res[str(p)]
            if r:
                print(
                    f"    max|m4|_κ1={r['max_abs_m4_kappa1']:.8f} κ3={r['max_abs_m4_kappa3']:.8f} "
                    f"R1={r['max_abs_rho_kappa1']:.6f} R3={r['max_abs_rho_kappa3']:.6f} "
                    f"le_cand={r['m4_le_cand']} wall={r['wall_s']:.2f}s",
                    flush=True,
                )
    else:
        print("  GPU unavailable — skip residual probe", flush=True)

    degs_ok = all(
        certs[k]["d1_constant_match"] and certs[k]["d3_constant_match"] for k in certs
    )
    gpu_ok = all(gpu_res[k]["m4_le_cand"] for k in gpu_res) if gpu_res else None

    out = {
        "workers": W,
        "compute": snap,
        "prop": "15.76",
        "algebra_degrees": alg,
        "bound_attempt": bound,
        "degree_certs": certs,
        "gpu_residual_probe": gpu_res,
        "one_center_degrees_constant": degs_ok,
        "gpu_cand_ok": gpu_ok,
        "io_policy": "mmap Max+; GPU device m4; write_json_atomic",
        "status": (
            f"Prop 15.76: one-center degrees d1^(1)=(3p^2-7)/4, d3^(1)=(p^2-5)/4 "
            f"constant on Paley |κ|=1 (certified p=3,5,7,11: {degs_ok}); "
            f"divisibility algebra for odd p proved. "
            f"GPU residual probe p=5,7 le_cand={gpu_ok}. "
            f"Absolute bootstrap still non-contracting (4p<d1). "
            f"OPEN: signed bound on S1+S3 ⇒ |m4|≤M_cand all p≥5. lim α_n OPEN."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_onecenter_deg.json"
    write_json_atomic(path, out)
    print(out["status"], flush=True)
    print(f"wrote {path} (atomic)", flush=True)


if __name__ == "__main__":
    main()
