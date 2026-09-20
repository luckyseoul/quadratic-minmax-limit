#!/usr/bin/env python3
"""Triangle-free minus-graph census for plus-I Case B at p=5, m=13.

D10 uses a maximizer of Q (not |Q|). After switching that vector to 1,
M+ = Q(1) = 78 - 2 e_minus. Hamming-1 of a Q-maximizer is plus min-degree
≥6, i.e. minus Δ≤6. The Paley-gap Extra bound at n=26 is
Γ ≥ 48 - 2 M+, so Case B is M+ ≤ 22, i.e. e_minus ≥ 28.

A triangle-free minus graph has Γ=0 (Mantel), so it is a Case B
counterexample iff 1 is a Q-maximizer. That is the single cut condition
2 e_+^{cut} ≥ |U||U^c| for every U (equivalently e_minus^{cut} ≤ |U||U^c|/2).
Same-sign / |Q| constraints are not part of D10.

geng -t -D6 13 28:39 enumerates the unlabeled class completely
(163477 graphs). Isomorphism types suffice for existence.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

M = 13
BINOM = M * (M - 1) // 2
GENG = os.environ.get("GENG", "/tmp/nauty2_9_1/geng")
G6_PAIRS = [(i, j) for j in range(1, M) for i in range(j)]
assert len(G6_PAIRS) == BINOM


def cut_tables() -> tuple[np.ndarray, np.ndarray]:
    nmask = 1 << (M - 1)
    cut0 = np.zeros(nmask, dtype=np.uint64)
    cut1 = np.zeros(nmask, dtype=np.uint64)
    cut_sz = np.zeros(nmask, dtype=np.int32)
    for xb in range(nmask):
        sign = np.ones(M, dtype=np.int8)
        for i in range(M - 1):
            if not ((xb >> i) & 1):
                sign[i + 1] = -1
        nU = int(np.sum(sign < 0))
        cut_sz[xb] = nU * (M - nU)
        c0 = np.uint64(0)
        c1 = np.uint64(0)
        for b, (i, j) in enumerate(G6_PAIRS):
            if sign[i] == sign[j]:
                continue
            if b < 64:
                c0 |= np.uint64(1) << np.uint64(b)
            else:
                c1 |= np.uint64(1) << np.uint64(b - 64)
        cut0[xb], cut1[xb] = c0, c1
    return np.stack([cut0, cut1], axis=1), cut_sz


def parse_g6_bits(line: str) -> tuple[int, int] | None:
    s = line.strip()
    if not s or s[0] == ">":
        return None
    if ord(s[0]) - 63 != M:
        return None
    bits = 0
    k = 0
    for ch in s[1:]:
        v = ord(ch) - 63
        for sh in range(5, -1, -1):
            if k >= BINOM:
                break
            if (v >> sh) & 1:
                bits |= 1 << k
            k += 1
    return bits & ((1 << 64) - 1), bits >> 64


def is_q_maximizer(g0: int, g1: int, cut: np.ndarray, cut_sz: np.ndarray) -> bool:
    e_cut = np.bitwise_count(cut[:, 0] & np.uint64(g0)) + np.bitwise_count(
        cut[:, 1] & np.uint64(g1)
    )
    return bool(np.all(cut_sz >= 2 * e_cut))


def seidel_from_minus_bits(g0: int, g1: int) -> np.ndarray:
    B = np.ones((M, M), dtype=np.float64)
    np.fill_diagonal(B, 0.0)
    for b, (i, j) in enumerate(G6_PAIRS):
        word = g0 if b < 64 else g1
        bit = b if b < 64 else b - 64
        if (word >> bit) & 1:
            B[i, j] = B[j, i] = -1.0
    return B


def brute_q_range(B: np.ndarray) -> tuple[float, float]:
    n = B.shape[0]
    qs = []
    for xb in range(1 << (n - 1)):
        x = np.ones(n)
        for i in range(n - 1):
            x[i + 1] = 1.0 if (xb >> i) & 1 else -1.0
        qs.append(0.5 * float(x @ B @ x))
    return max(qs), min(qs)


def gamma_of(B: np.ndarray) -> tuple[int, int, int]:
    """Γ at the all-ones switching. Returns (Gamma, t_star, eplus_star)."""
    best = -10**9
    t_star = 0
    ep_star = 0
    plus = B > 0
    for mask in range(1 << M):
        T = [i for i in range(M) if mask & (1 << i)]
        t = len(T)
        if t <= 1:
            extra = 0 if t == 0 else -2
        else:
            eplus = 0
            for a in range(t):
                for b in range(a + 1, t):
                    if plus[T[a], T[b]]:
                        eplus += 1
            extra = 2 * t * (t - 2) - 8 * eplus
        if extra > best:
            best, t_star, ep_star = extra, t, eplus if t >= 2 else 0
    return int(best), t_star, ep_star


def self_test(cut, cut_sz) -> dict:
    rec = {}
    # All-plus Seidel (empty minus graph): all-ones is the unique Q-max, Q=78.
    rec["all_plus_is_qmax"] = is_q_maximizer(0, 0, cut, cut_sz)
    rec["all_plus_q"] = float(brute_q_range(seidel_from_minus_bits(0, 0))[0])
    # Complete minus: all-ones is a Q-min, not a Q-max.
    all_m0 = (1 << 64) - 1
    all_m1 = (1 << (BINOM - 64)) - 1
    rec["all_minus_is_qmax"] = is_q_maximizer(all_m0, all_m1, cut, cut_sz)
    # Known m_13=20 block: switch a Q-max to 1 and confirm the cut checker.
    root = Path(__file__).resolve().parents[1]
    B = np.loadtxt(root / "evidence" / "coherent_BI_lift_20260919" / "B13.txt")
    n = 13
    best_q, best_x = -1e9, None
    for xb in range(1 << (n - 1)):
        x = np.ones(n)
        for i in range(n - 1):
            x[i + 1] = 1.0 if (xb >> i) & 1 else -1.0
        q = 0.5 * float(x @ B @ x)
        if q > best_q:
            best_q, best_x = q, x.copy()
    Bp = (best_x[:, None] * B) * best_x[None, :]
    np.fill_diagonal(Bp, 0.0)
    bits = 0
    e_minus = 0
    for b, (i, j) in enumerate(G6_PAIRS):
        if Bp[i, j] < 0:
            bits |= 1 << b
            e_minus += 1
    g0, g1 = bits & ((1 << 64) - 1), bits >> 64
    rec["B13_Mplus"] = float(best_q)
    rec["B13_Q1_after_switch"] = 0.5 * float(np.ones(n) @ Bp @ np.ones(n))
    rec["B13_e_minus"] = e_minus
    rec["B13_checker_qmax"] = is_q_maximizer(g0, g1, cut, cut_sz)
    rec["B13_minus_triangle"] = bool(
        any(
            Bp[i, j] < 0 and Bp[i, k] < 0 and Bp[j, k] < 0
            for i in range(n)
            for j in range(i + 1, n)
            for k in range(j + 1, n)
        )
    )
    rec["B13_gamma"] = gamma_of(Bp)[0]
    rec["B13_needed"] = 48 - 2 * int(round(best_q))
    assert rec["all_plus_is_qmax"]
    assert rec["all_plus_q"] == 78.0
    assert rec["all_minus_is_qmax"] is False
    assert rec["B13_checker_qmax"] is True
    assert abs(rec["B13_Mplus"] - rec["B13_Q1_after_switch"]) < 1e-9
    return rec


def run_shard(emin: int, emax: int, res: int, mod: int, cut, cut_sz) -> dict:
    cmd = [GENG, "-t", "-q", "-D6", str(M), f"{emin}:{emax}", f"{res}/{mod}"]
    t0 = time.time()
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    n_seen = 0
    by_e = {}
    hits = []
    assert proc.stdout is not None
    for line in proc.stdout:
        parsed = parse_g6_bits(line)
        if parsed is None:
            continue
        g0, g1 = parsed
        n_seen += 1
        e = g0.bit_count() + g1.bit_count()
        by_e[e] = by_e.get(e, 0) + 1
        if is_q_maximizer(g0, g1, cut, cut_sz):
            B = seidel_from_minus_bits(g0, g1)
            mx, mn = brute_q_range(B)
            gam, t_star, ep_star = gamma_of(B)
            mplus = BINOM - 2 * e
            hits.append(
                {
                    "g6": line.strip(),
                    "e_minus": e,
                    "Mplus": mplus,
                    "brute_maxQ": mx,
                    "brute_minQ": mn,
                    "Gamma": gam,
                    "T_star": t_star,
                    "eplus_star": ep_star,
                    "needed": 48 - 2 * mplus,
                }
            )
    stderr = proc.stderr.read() if proc.stderr else ""
    rc = proc.wait()
    return {
        "emin": emin,
        "emax": emax,
        "res": res,
        "mod": mod,
        "n_seen": n_seen,
        "by_e": {str(k): by_e[k] for k in sorted(by_e)},
        "n_hits": len(hits),
        "hits": hits,
        "wall_s": round(time.time() - t0, 3),
        "geng_rc": rc,
        "geng_stderr": stderr.strip()[-500:],
        "cmd": cmd,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--emin", type=int, default=28)
    ap.add_argument("--emax", type=int, default=39)
    ap.add_argument("--res", type=int, default=0)
    ap.add_argument("--mod", type=int, default=1)
    ap.add_argument("--skip-self-test", action="store_true")
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()
    cut, cut_sz = cut_tables()
    rec = {
        "m": 13,
        "p": 5,
        "note": "Q-maximizer census of triangle-free minus graphs, Δ≤6, e≥28",
    }
    if not args.skip_self_test:
        rec["self_test"] = self_test(cut, cut_sz)
        print("SELF_TEST", json.dumps(rec["self_test"]), flush=True)
    shard = run_shard(args.emin, args.emax, args.res, args.mod, cut, cut_sz)
    rec.update(shard)
    rec["tf_case_b_counterexamples"] = rec["n_hits"]
    print(json.dumps({k: rec[k] for k in rec if k != "hits"}, indent=2), flush=True)
    if rec["hits"]:
        print("HITS", json.dumps(rec["hits"][:20]), flush=True)
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(rec, indent=2))


if __name__ == "__main__":
    main()
