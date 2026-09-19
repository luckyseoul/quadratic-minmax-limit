#!/usr/bin/env python3
"""Extend the two-half orientation minimax from n=5..8 to exact m_9=12, m_10=13.

Implication: for a recorded exact minimizer A, does there exist a skew R with
    B(A,R) := max_{x,y}(|Q_A(x)+Q_A(y)| + |x^T R y|)  <=  2 sqrt(2) Phi(A) ?
That is the zero-error diamond (6.13). Paley-R failed it at these orders in
the 2026-09-19 census; the stored n=7,8 orientations already pass. A passing
R is a finite feasibility certificate. A CP-SAT OPTIMAL B above the target
is a finite failure. Neither is an all-orders doubling proof.

Backend: ProcessPool over ILS seeds (W=86). Inner GEMM is 256^2 / 512^2.
GPU unused (driver down). Serial CP-SAT only if ILS does not find a passing R.
"""
from __future__ import annotations

import json
import math
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from original_mo_two_half_geometry import (  # noqa: E402
    analyze_orientation,
    bits_from_skew,
    boolean_states,
    edge_pairs,
    quadratic_values,
    skew_from_bits,
)

OUT = ROOT / "evidence" / "two_half_n9_n10_20260919"
TARGETS = {9: 12, 10: 13}


def random_signing(n: int, rng: np.random.Generator) -> np.ndarray:
    u = rng.integers(0, 2, size=(n, n)) * 2 - 1
    a = np.triu(u, 1)
    a = a + a.T
    np.fill_diagonal(a, 0)
    return a.astype(np.int64)


def energies_and_bits(a: np.ndarray) -> tuple[np.ndarray, list[np.ndarray]]:
    n = a.shape[0]
    m = 1 << (n - 1)
    u = np.arange(m, dtype=np.int64)
    bits = [np.ones(m, dtype=np.int8)]
    for k in range(1, n):
        bits.append((1 - 2 * ((u >> (k - 1)) & 1)).astype(np.int8))
    e = np.zeros(m, dtype=np.int32)
    for i in range(n):
        for j in range(i + 1, n):
            e += np.int32(a[i, j]) * (bits[i] * bits[j]).astype(np.int32)
    return e, bits


def ils_signing(n: int, target: int, seed: int, steps: int = 200) -> tuple[int, np.ndarray]:
    rng = np.random.default_rng(seed)
    a = random_signing(n, rng)
    e, bits = energies_and_bits(a)
    best = int(np.abs(e).max())
    for _ in range(steps):
        bi = bj = None
        bv = best
        for i in range(n):
            for j in range(i + 1, n):
                s = (bits[i] * bits[j]).astype(np.int32)
                v = int(np.abs(e - 2 * np.int32(a[i, j]) * s).max())
                if v < bv:
                    bv = v
                    bi, bj = i, j
        if bi is None:
            break
        e = e - 2 * np.int32(a[bi, bj]) * (bits[bi] * bits[bj]).astype(np.int32)
        a[bi, bj] *= -1
        a[bj, bi] *= -1
        best = bv
        if best <= target:
            break
    return best, a


def ils_orientation(a: np.ndarray, seed: int, steps: int = 400) -> tuple[int, np.ndarray]:
    """Greedy one-edge flips of R minimizing B = max(U+|xRy|)."""
    n = a.shape[0]
    states = boolean_states(n, fix_first=True).astype(np.int64)
    q = quadratic_values(a, states)
    u_mat = np.abs(q[:, None] + q[None, :])
    rng = np.random.default_rng(seed)
    r = np.zeros((n, n), dtype=np.int64)
    for i, j in edge_pairs(n):
        s = 1 if rng.integers(0, 2) else -1
        r[i, j] = s
        r[j, i] = -s
    # incremental xRy: states @ R @ states.T
    w = states @ r @ states.T
    b = int(np.max(u_mat + np.abs(w)))
    edges = edge_pairs(n)
    for _ in range(steps):
        improved = False
        rng.shuffle(edges)
        for i, j in edges:
            # flipping R_ij by -2 R_ij changes xRy by -2 R_ij (x_i y_j - x_j y_i)
            delta = -2 * int(r[i, j]) * (
                states[:, i : i + 1] * states[:, j] - states[:, j : j + 1] * states[:, i]
            )
            # delta[p,q] = -2 R_ij (x_p_i y_q_j - x_p_j y_q_i)
            d = -2 * int(r[i, j]) * (
                np.outer(states[:, i], states[:, j]) - np.outer(states[:, j], states[:, i])
            )
            w2 = w + d
            b2 = int(np.max(u_mat + np.abs(w2)))
            if b2 < b:
                w = w2
                b = b2
                r[i, j] *= -1
                r[j, i] *= -1
                improved = True
                break
        if not improved:
            # random kick
            i, j = edges[int(rng.integers(0, len(edges)))]
            d = -2 * int(r[i, j]) * (
                np.outer(states[:, i], states[:, j]) - np.outer(states[:, j], states[:, i])
            )
            w = w + d
            r[i, j] *= -1
            r[j, i] *= -1
            b = int(np.max(u_mat + np.abs(w)))
    return b, r


def hunt_signing(n: int, target: int, seeds: list[int]) -> np.ndarray:
    with ProcessPoolExecutor(max_workers=min(86, len(seeds))) as ex:
        futs = [ex.submit(ils_signing, n, target, s) for s in seeds]
        for fut in as_completed(futs):
            phi, a = fut.result()
            if phi <= target:
                for f in futs:
                    f.cancel()
                return a
    raise RuntimeError(f"no signing with Phi<={target} at n={n}")


def hunt_orientation(a: np.ndarray, seeds: list[int]) -> tuple[int, np.ndarray]:
    best_b = 10**9
    best_r = None
    with ProcessPoolExecutor(max_workers=min(86, len(seeds))) as ex:
        futs = [ex.submit(ils_orientation, a, s) for s in seeds]
        for fut in as_completed(futs):
            b, r = fut.result()
            if b < best_b:
                best_b = b
                best_r = r
                print(f"  ILS R B={b}", flush=True)
    assert best_r is not None
    return best_b, best_r


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    payload = {
        "status": "finite orientation extension; original limit OPEN",
        "zero_error_target": "B <= 2 sqrt(2) M",
        "backend": "ProcessPool ILS W<=86, GPU unused",
        "instances": {},
    }
    t_all = time.time()
    for n, m in TARGETS.items():
        t0 = time.time()
        print(f"=== n={n} target m={m} ===", flush=True)
        a = hunt_signing(n, m, list(range(200)))
        e, _ = energies_and_bits(a)
        phi = int(np.abs(e).max())
        print(f"  found A Phi={phi} in {time.time()-t0:.2f}s", flush=True)
        assert phi == m, (phi, m)
        b, r = hunt_orientation(a, list(range(86)))
        rec = analyze_orientation(a, m, r)
        target = 2 * math.sqrt(2) * m
        rec_out = {
            "n": n,
            "M": m,
            "Phi": rec["phi_A"],
            "B_ils": b,
            "B_checked": rec["B"],
            "max_W_minus_D": rec["max_W_minus_D"],
            "target_2sqrt2_M": target,
            "zero_error_passes": rec["B"] <= target + 1e-9,
            "excess_over_2sqrt2_M": rec["B"] - target,
            "normalized_excess": (rec["B"] - target) / (n ** 1.5),
            "r_bits": bits_from_skew(r),
            "A": a.tolist(),
            "seconds": round(time.time() - t0, 3),
            "stored_orientation_summary": {
                "B": rec["B"],
                "zero_error_doubling_passes": rec["zero_error_doubling_passes"],
                "target_margin": rec["target_margin"],
                "directed_halfcut_max_phi": rec["directed_halfcut"]["max_phi"],
            },
        }
        payload["instances"][str(n)] = rec_out
        np.save(OUT / f"A{n}_exact.npy", a)
        print(json.dumps({k: rec_out[k] for k in rec_out if k != "A"}, sort_keys=True), flush=True)
    payload["seconds"] = round(time.time() - t_all, 3)
    (OUT / "summary.json").write_text(json.dumps(payload, indent=2) + "\n")
    print("WROTE", OUT / "summary.json", flush=True)


if __name__ == "__main__":
    main()
