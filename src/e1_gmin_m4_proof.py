#!/usr/bin/env python3
"""
P0 proof attack: |m4| ≤ (p-2)/(2p²) on |κ|=1 (⇔ g_min ≥ L(p)).

Uses master identity (Prop 15.66 lineage):
  m4(S) = κ(S)/p² + Ext(S)/(4p),  Ext = (T m4)(S)
  σ_sum = 4κ for every ±1 K4 labeling (combinatorial exhaustion).

Multi-worker only (F17). Shards 4-set work across W = nproc-2 processes.
Writes evidence/e1_gmin_m4_proof.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import Counter
from concurrent.futures import as_completed
from fractions import Fraction
from itertools import combinations, product
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from workers import cpu_count_reliable, pool, require_workers  # noqa: E402


def L_abs(p: int) -> Fraction:
    return Fraction(p - 2, 2 * p * p)


def T_abs(p: int) -> Fraction:
    return Fraction(p - 2, p * (2 * p - 1))


def kappa_pts(C: np.ndarray, S) -> int:
    a, b, c, d = S
    return int(C[a, b] * C[c, d] + C[a, c] * C[b, d] + C[a, d] * C[b, c])


def sigma_sum_identity() -> dict:
    """Combinatorial: σ_sum = 4κ on all 64 ±1 K4 edge labelings."""

    def check(bits):
        edges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
        Cv = {edges[i]: bits[i] for i in range(6)}

        def C(i, j):
            if i == j:
                return 0
            return Cv[tuple(sorted((i, j)))]

        def sigma(v):
            others = [x for x in range(4) if x != v]
            a, b, c = others
            return C(v, a) * C(b, c) + C(v, b) * C(a, c) + C(v, c) * C(a, b)

        kap = C(0, 1) * C(2, 3) + C(0, 2) * C(1, 3) + C(0, 3) * C(1, 2)
        ss = sum(sigma(v) for v in range(4))
        return ss == 4 * kap

    ok = sum(1 for bits in product([-1, 1], repeat=6) if check(bits))
    return {"n_labelings": 64, "n_ok": ok, "proved": ok == 64}


def _load_Y(p: int) -> np.ndarray:
    if p == 3:
        from e1_gmin_cr_classify import load_maxplus

        return load_maxplus(3).astype(np.float64)
    if p == 5:
        return np.load("/tmp/maxplus_p5.npy").astype(np.float64)
    if p == 7:
        return np.load("/tmp/e1_p7/maxplus.npy").astype(np.float64)
    raise ValueError(p)


def _shard_quads(n: int, n_shards: int):
    quads = list(combinations(range(n), 4))
    m = len(quads)
    sizes = [m // n_shards] * n_shards
    for i in range(m % n_shards):
        sizes[i] += 1
    out = []
    k = 0
    for s in sizes:
        out.append(quads[k : k + s])
        k += s
    return out


def _worker_m4_shard(args) -> dict:
    """Compute m4, kappa, Ext via identity Ext=4p m4 - 4κ/p on a shard of |κ|=1 quads.

    Full Ext expansion is O(n) per quad; identity evaluation is O(1) once m4 known.
    Here we only need m4 on the shard's quads for the direct |m4|≤L test.
    """
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[k] = "1"
    p, shard, Y_path, C = args
    Y = np.asarray(np.load(Y_path), dtype=np.float64)
    C = np.asarray(C, dtype=np.float64)

    max_abs_m4 = 0.0
    max_abs_ext_same = 0.0
    max_abs_ext_all = 0.0
    n_k1 = 0
    n_same = 0
    n_opp = 0
    n_le_L = 0
    worst = None
    thr_L = float(L_abs(p))
    thr_ext = 2.0 * (p - 4) / p
    unique_m4 = Counter()

    # Vectorized m4: for each quad, mean product
    for S in shard:
        a, b, c, d = S
        m = float(np.mean(Y[:, a] * Y[:, b] * Y[:, c] * Y[:, d]))
        kap = kappa_pts(C, S)
        if abs(kap) != 1:
            continue
        n_k1 += 1
        am = abs(m)
        if am > max_abs_m4:
            max_abs_m4 = am
            worst = {"S": list(S), "m4": m, "kappa": kap}
        if am <= thr_L + 1e-12:
            n_le_L += 1
        # Ext from identity (exact up to float)
        ext = 4 * p * m - 4 * kap / p
        ae = abs(ext)
        if ae > max_abs_ext_all:
            max_abs_ext_all = ae
        same = np.sign(ext) * np.sign(kap) > 0
        if same:
            n_same += 1
            if ae > max_abs_ext_same:
                max_abs_ext_same = ae
        elif np.sign(ext) * np.sign(kap) < 0:
            n_opp += 1
        # rational bucket
        try:
            fr = Fraction(m).limit_denominator(10**6)
            unique_m4[str(fr)] += 1
        except Exception:
            unique_m4[f"{m:.12f}"] += 1

    return {
        "n_k1": n_k1,
        "n_le_L": n_le_L,
        "max_abs_m4": max_abs_m4,
        "max_abs_ext_same": max_abs_ext_same,
        "max_abs_ext_all": max_abs_ext_all,
        "n_same": n_same,
        "n_opp": n_opp,
        "worst": worst,
        "unique_m4": dict(unique_m4),
        "thr_L": thr_L,
        "thr_ext": thr_ext,
        "ext_same_le_thr": max_abs_ext_same <= thr_ext + 1e-9,
        "m4_le_L": max_abs_m4 <= thr_L + 1e-12,
    }


def _merge_shard_results(parts: list[dict], p: int) -> dict:
    if not parts:
        return {"p": p, "empty": True}
    max_abs_m4 = max(r["max_abs_m4"] for r in parts)
    max_abs_ext_same = max(r["max_abs_ext_same"] for r in parts)
    max_abs_ext_all = max(r["max_abs_ext_all"] for r in parts)
    n_k1 = sum(r["n_k1"] for r in parts)
    n_le_L = sum(r["n_le_L"] for r in parts)
    n_same = sum(r["n_same"] for r in parts)
    n_opp = sum(r["n_opp"] for r in parts)
    thr_L = parts[0]["thr_L"]
    thr_ext = parts[0]["thr_ext"]
    worst = max(parts, key=lambda r: r["max_abs_m4"])["worst"]
    uniq: Counter = Counter()
    for r in parts:
        uniq.update(r["unique_m4"])
    # algebra
    lab = L_abs(p)
    tab = T_abs(p)
    cand_den = Fraction(p - 2, p * (2 * p + 3))  # exact at p=5
    return {
        "p": p,
        "n_kappa1": n_k1,
        "n_le_L": n_le_L,
        "all_m4_le_L": n_k1 == n_le_L and max_abs_m4 <= thr_L + 1e-12,
        "max_abs_m4": max_abs_m4,
        "max_abs_m4_frac": str(Fraction(max_abs_m4).limit_denominator(10**6)),
        "L_abs": float(lab),
        "L_abs_frac": str(lab),
        "T_abs": float(tab),
        "gap_L_minus_maxm4": float(lab) - max_abs_m4,
        "max_abs_Ext_same_sign": max_abs_ext_same,
        "max_abs_Ext_all": max_abs_ext_all,
        "thr_Ext_same_sign": thr_ext,
        "Ext_same_le_thr": max_abs_ext_same <= thr_ext + 1e-9,
        "n_same_sign": n_same,
        "n_opp_sign": n_opp,
        "worst": worst,
        "unique_m4_fracs": dict(uniq),
        "n_unique_m4": len(uniq),
        "candidate_ub_(p-2)/(p(2p+3))": float(cand_den),
        "max_m4_le_candidate": max_abs_m4 <= float(cand_den) + 1e-12,
        "L_gt_T_signed": lab < tab,
    }


def run_p(p: int, W: int) -> dict:
    from minmax_quadratic import paley_conference_prime_power

    C = paley_conference_prime_power(p).astype(np.float64)
    n = C.shape[0]
    if p == 7:
        Y_path = "/tmp/e1_p7/maxplus.npy"
        # memory: workers load from path
        y_arg = Y_path
    elif p == 5:
        y_arg = "/tmp/maxplus_p5.npy"
    else:
        Y = _load_Y(p)
        # save temp for workers
        path = f"/tmp/maxplus_p{p}_work.npy"
        np.save(path, Y)
        y_arg = path

    shards = _shard_quads(n, W)
    print(f"  p={p}: n={n} shards={len(shards)} total_quads={sum(len(s) for s in shards)}", flush=True)

    parts = []
    with pool(W) as ex:
        futs = [
            ex.submit(_worker_m4_shard, (p, shard, y_arg, C))
            for shard in shards
            if shard
        ]
        for i, fut in enumerate(as_completed(futs)):
            parts.append(fut.result())
            if (i + 1) % max(1, len(futs) // 10) == 0:
                print(f"    p={p} shards done {i+1}/{len(futs)}", flush=True)

    return _merge_shard_results(parts, p)


def algebra_note() -> dict:
    """Max+-free algebra that reduces P0 to an Ext / direct bound."""
    rows = {}
    for p in (5, 7, 11, 13, 17, 19, 23):
        lab = L_abs(p)
        rows[str(p)] = {
            "L_abs": str(lab),
            "wick_abs": str(Fraction(1, p * p)),
            "resid_budget_(p-4)/(2p^2)": str(Fraction(p - 4, 2 * p * p)) if p > 4 else None,
            "Ext_same_sign_thr_2(p-4)/p": str(Fraction(2 * (p - 4), p)) if p > 4 else None,
            # signed: L=-(p-2)/(2p²) > T=-(p-2)/(p(2p-1)) iff L_abs < T_abs
            "L_gt_T_signed": bool(L_abs(p) < T_abs(p)),
            "note": (
                "On |κ|=1: if sign(Ext)=sign(κ) and |Ext|≤2(p-4)/p then |m4|≤L_abs; "
                "or prove |m4|≤L_abs directly. Triangle |m4-κ/p²|≤(p-4)/(2p²) fails at p=5."
            ),
        }
    return {
        "master_identity": "m4=κ/p²+Ext/(4p), Ext=(T m4), σ_sum=4κ (K4 exhaustion)",
        "criterion_same_sign_Ext": "|Ext|≤2(p-4)/p on same-sign |κ|=1 ⇒ |m4|≤L_abs",
        "criterion_direct": "|m4|≤(p-2)/(2p²) on all |κ|=1 ⇒ g_min≥L(p)>T(p) ⇒ bi-tight empty (p≥5)",
        "by_p": rows,
    }


def main() -> None:
    W = require_workers(min_workers=4)
    print(f"m4_proof: W={W} (F17 multi-worker)", flush=True)

    out: dict = {
        "workers": W,
        "cpu_count": cpu_count_reliable(),
        "sigma_identity": sigma_sum_identity(),
        "algebra": algebra_note(),
        "certs": {},
        "status": None,
    }
    print(f"  sigma identity: {out['sigma_identity']}", flush=True)

    for p in (5, 7):
        print(f"running full |κ|=1 census p={p} ...", flush=True)
        out["certs"][str(p)] = run_p(p, W)
        r = out["certs"][str(p)]
        print(
            f"  p={p}: max|m4|={r['max_abs_m4']:.8f} L={r['L_abs']:.8f} "
            f"all_le_L={r['all_m4_le_L']} Ext_same_le={r['Ext_same_le_thr']} "
            f"gap={r['gap_L_minus_maxm4']:.8f}",
            flush=True,
        )

    c5 = out["certs"].get("5", {})
    c7 = out["certs"].get("7", {})
    ok = (
        out["sigma_identity"]["proved"]
        and c5.get("all_m4_le_L")
        and c7.get("all_m4_le_L")
        and c5.get("Ext_same_le_thr")
        and c7.get("Ext_same_le_thr")
    )
    out["status"] = (
        "CERTIFIED p=5,7: |m4|≤L_abs and same-sign Ext≤thr on all |κ|=1; "
        "σ_sum=4κ proved by exhaustion; master identity algebra ready. "
        "OPEN for general primes p≥5: prove |m4|≤L_abs (or Ext criterion) "
        "from Max+/boolean/C structure without per-p Max+ census."
        if ok
        else "FAILED cert or incomplete — inspect certs."
    )
    path = ROOT / "evidence" / "e1_gmin_m4_proof.json"
    path.write_text(json.dumps(out, indent=2))
    print(out["status"], flush=True)
    print(f"wrote {path}", flush=True)


if __name__ == "__main__":
    main()
