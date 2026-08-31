#!/usr/bin/env python3
"""Min-conflicts search for an integral affine residual-(ii) candidate.

This is a heuristic companion to ``residual_affine_johnson_milp.py``.  It
searches only the exact affine Johnson cardinality cuts and, if it finds a
candidate, reports its violations on the complete cached shells.  A found
candidate is evidence; failure to find one is not a certificate.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from minmax_quadratic import paley_conference_prime_power  # noqa: E402
from residual_affine_johnson_milp import (  # noqa: E402
    affine_shell,
    feature_rows,
    unique_rows,
)


def cardinality_model(p: int, shell_mode: str):
    if shell_mode not in ("affine", "full"):
        raise ValueError("shell_mode must be affine or full")
    C = np.rint(paley_conference_prime_power(p)).astype(np.int8)
    blocks: list[np.ndarray] = []
    caps: list[int] = []
    edges = None
    ei = None
    k = 4 * p
    for shell, eigen_sign in (("plus", 1), ("minus", -1)):
        if shell_mode == "affine":
            Y = affine_shell(p, eigen_sign, C)
        else:
            path = Path(f"/tmp/max{shell}_p{p}.npy")
            if not path.exists():
                raise FileNotFoundError(path)
            Y = np.sign(np.load(path, mmap_mode="r")).astype(np.int8)
        these_edges, F = feature_rows(Y, C)
        if edges is None:
            edges = these_edges
            ei = edges.index((0, 1))
        elif these_edges != edges:
            raise RuntimeError("edge orders differ")
        assert ei is not None
        fe = F[:, ei]
        for sign in (-1, 1):
            part = unique_rows(F[fe == sign])
            if shell == "plus":
                bound = 4 if sign == -1 else 2
                blocks.append((part < 0).astype(np.uint8))
                caps.extend([(k - bound) // 2] * len(part))
            else:
                bound = -2 if sign == -1 else -4
                blocks.append((part > 0).astype(np.uint8))
                caps.extend([(k + bound) // 2] * len(part))
    assert edges is not None and ei is not None
    return C, edges, ei, np.concatenate(blocks), np.asarray(caps, dtype=np.int16)


def lp_seed(path: Path | None, edges: list[tuple[int, int]]) -> np.ndarray:
    seed = np.zeros(len(edges), dtype=np.float64)
    if path is None:
        return seed
    data = json.loads(path.read_text())
    idx = {edge: j for j, edge in enumerate(edges)}
    for a, b, value in data.get("fractional_edges", []):
        seed[idx[(int(a), int(b))]] = float(value)
    return seed


def penalty(counts: np.ndarray, caps: np.ndarray) -> int:
    excess = np.maximum(counts - caps, 0).astype(np.int32)
    return int(np.sum(excess * excess * 100 + excess))


def search(
    p: int,
    shell_mode: str,
    seconds: float,
    seed_path: Path | None,
    random_seed: int,
    add_sample: int,
    restarts: int,
) -> dict:
    t0 = time.time()
    C, edges, ei, A, caps = cardinality_model(p, shell_mode)
    nvar = len(edges)
    k = 4 * p
    rng = np.random.default_rng(random_seed)
    frac = lp_seed(seed_path, edges)
    allowed = np.ones(nvar, dtype=bool)
    allowed[ei] = False
    best_penalty = 10**18
    best_x = None
    steps = 0

    for restart in range(restarts):
        if time.time() - t0 >= seconds:
            break
        noise = rng.gumbel(size=nvar)
        priority = frac + (0.08 + 0.01 * min(restart, 12)) * noise
        priority[~allowed] = -1e100
        selected = np.argpartition(priority, -k)[-k:]
        x = np.zeros(nvar, dtype=np.uint8)
        x[selected] = 1
        counts = A.astype(np.int16) @ x.astype(np.int16)
        cur = penalty(counts, caps)

        for local_step in range(20000):
            steps += 1
            if cur < best_penalty:
                best_penalty = cur
                best_x = x.copy()
            if cur == 0:
                elapsed = time.time() - t0
                out = finish(
                    p, C, edges, ei, A, caps, x, elapsed, steps, restart + 1
                )
                out["shell_mode"] = shell_mode
                return out
            if time.time() - t0 >= seconds:
                break

            excess = counts - caps
            violated = np.flatnonzero(excess > 0)
            weights = excess[violated].astype(np.float64)
            c = int(rng.choice(violated, p=weights / weights.sum()))
            rem = np.flatnonzero((x == 1) & (A[c] == 1))
            add = np.flatnonzero((x == 0) & allowed & (A[c] == 0))
            if len(rem) == 0 or len(add) == 0:
                break
            if len(rem) > 12:
                rem = rng.choice(rem, 12, replace=False)
            if len(add) > add_sample:
                add = rng.choice(add, add_sample, replace=False)

            chosen = None
            chosen_penalty = 10**18
            for i in rem:
                trial = (
                    counts[:, None].astype(np.int16)
                    - A[:, i, None].astype(np.int16)
                    + A[:, add].astype(np.int16)
                )
                ex = np.maximum(trial - caps[:, None], 0).astype(np.int32)
                scores = np.sum(ex * ex * 100 + ex, axis=0)
                jpos = int(np.argmin(scores))
                score = int(scores[jpos])
                if score < chosen_penalty:
                    chosen_penalty = score
                    chosen = (int(i), int(add[jpos]))
            if chosen is None:
                break
            i, j = chosen
            # Accept descent and plateaux.  At a strict local minimum, kick.
            if chosen_penalty <= cur or rng.random() < 0.003:
                x[i] = 0
                x[j] = 1
                counts += A[:, j].astype(np.int16) - A[:, i].astype(np.int16)
                cur = chosen_penalty
            else:
                # Two random swaps preserve cardinality and leave the fixed edge out.
                for _ in range(2):
                    si = int(rng.choice(np.flatnonzero(x == 1)))
                    uj = int(rng.choice(np.flatnonzero((x == 0) & allowed)))
                    x[si] = 0
                    x[uj] = 1
                    counts += A[:, uj].astype(np.int16) - A[:, si].astype(np.int16)
                cur = penalty(counts, caps)

    return {
        "p": p,
        "shell_mode": shell_mode,
        "found": False,
        "best_penalty": int(best_penalty),
        "steps": steps,
        "restarts": min(restarts, restart + 1 if restarts else 0),
        "seconds": round(time.time() - t0, 3),
        "best_edges": (
            [list(edges[j]) for j in np.flatnonzero(best_x)]
            if best_x is not None
            else []
        ),
        "not_a_certificate": True,
    }


def finish(p, C, edges, ei, A, caps, x, elapsed, steps, restarts):
    out = {
        "p": p,
        "found": True,
        "affine_constraints": int(len(A)),
        "affine_max_excess": int(np.max(A.astype(np.int16) @ x - caps)),
        "steps": steps,
        "restarts": restarts,
        "seconds": round(elapsed, 3),
        "chosen_edges": [list(edges[j]) for j in np.flatnonzero(x)],
    }
    paths = {
        5: (Path("/tmp/maxplus_p5.npy"), Path("/tmp/maxminus_p5.npy")),
        7: (Path("/tmp/maxplus_p7.npy"), Path("/tmp/maxminus_p7.npy")),
    }
    if p in paths and all(path.exists() for path in paths[p]):
        Yp = np.sign(np.load(paths[p][0], mmap_mode="r")).astype(np.int8)
        Ym = np.sign(np.load(paths[p][1], mmap_mode="r")).astype(np.int8)
        _, Fp = feature_rows(Yp, C)
        _, Fm = feature_rows(Ym, C)
        Sp = Fp.astype(np.int16) @ x.astype(np.int16)
        Sm = Fm.astype(np.int16) @ x.astype(np.int16)
        bp = np.where(Fp[:, ei] == -1, 4, 2)
        bm = np.where(Fm[:, ei] == -1, -2, -4)
        out.update(
            {
                "full_plus_violations": int(np.count_nonzero(Sp < bp)),
                "full_minus_violations": int(np.count_nonzero(Sm > bm)),
                "full_min_plus_margin": int(np.min(Sp - bp)),
                "full_min_minus_margin": int(np.min(bm - Sm)),
            }
        )
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--p", type=int, default=7)
    ap.add_argument("--shell-mode", choices=("affine", "full"), default="affine")
    ap.add_argument("--seconds", type=float, default=300.0)
    ap.add_argument("--lp-seed", type=Path)
    ap.add_argument("--seed", type=int, default=628)
    ap.add_argument("--add-sample", type=int, default=64)
    ap.add_argument("--restarts", type=int, default=1000)
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()
    out = search(
        args.p,
        args.shell_mode,
        args.seconds,
        args.lp_seed,
        args.seed,
        args.add_sample,
        args.restarts,
    )
    print(json.dumps(out, indent=2), flush=True)
    if args.output:
        args.output.write_text(json.dumps(out, indent=2) + "\n")


if __name__ == "__main__":
    main()
