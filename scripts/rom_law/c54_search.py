#!/usr/bin/env python3
"""Lower bounds for Phi(C54) by batched 1-opt ascent (ILS).

Phi(C54) <= 196 by the spectral ceiling 54*sqrt(53)/2 = 196.56 and integrality
(Q is a sum of 1431 terms of +-1, so Phi is odd). A witness at 188 or more
would already decide the source note's prediction 5 ("187 or 185") in the
negative. This script only produces witnesses: it never claims the maximum.
A witness x with |Q(x)| = v proves Phi >= v.

Q is the halved form Q(x) = (1/2) x^T C x, so a coordinate flip changes Q by
-2 x_j h_j with h = C x. Rows ascend in lockstep; every improving coordinate
flips at once and Q is recomputed by matmul, so nothing is inferred from the
per-coordinate deltas.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from paley_structural_lift_family import paley_conference  # noqa: E402


def descend(x, A, iters=400):
    for _ in range(iters):
        h = x @ A.T
        d = -2.0 * x * h
        m = d > 1e-9
        if not m.any():
            break
        x = np.where(m, -x, x)
    h = x @ A.T
    return x, 0.5 * (x * h).sum(axis=1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--q", type=int, default=53)
    ap.add_argument("--batch", type=int, default=131072)
    ap.add_argument("--rounds", type=int, default=24)
    ap.add_argument("--kick-rounds", type=int, default=8)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    C = paley_conference(args.q).astype(np.float32)
    n = C.shape[0]
    ceiling = n * float(np.sqrt(n - 1)) / 2.0
    print(f"[c54] n={n} spectral ceiling={ceiling:.4f} -> Phi<={int(ceiling)}",
          flush=True)

    rng = np.random.default_rng(11)
    best = {"Q": -1.0, "x": None, "sign": 1.0, "round": None}

    for sign in (1.0, -1.0):
        A = (sign * C).astype(np.float32)
        for rnd in range(args.rounds):
            x = (2.0 * rng.integers(0, 2, size=(args.batch, n)) - 1.0).astype(np.float32)
            x, Q = descend(x, A)
            k = int(Q.argmax())
            v = float(Q[k])
            if v > best["Q"]:
                best = {"Q": v, "x": x[k].copy(), "sign": sign, "round": rnd}
                print(f"[sign {sign:+.0f} rnd {rnd}] best Q={v}", flush=True)
        x0 = best["x"] * best["sign"]
        for kr in range(args.kick_rounds):
            x = np.tile(x0, (args.batch, 1)).astype(np.float32)
            idx = rng.integers(0, n, size=(args.batch, 3))
            rows = np.arange(args.batch)
            for t in range(3):
                x[rows, idx[:, t]] *= -1.0
            x, Q = descend(x, A)
            k = int(Q.argmax())
            v = float(Q[k])
            if v > best["Q"]:
                best = {"Q": v, "x": x[k].copy(), "sign": sign,
                        "round": f"kick{kr}"}
                print(f"[sign {sign:+.0f} kick {kr}] best Q={v}", flush=True)

    xw = (best["x"] * best["sign"]).astype(np.int64)
    Ci = C.astype(np.int64)
    Qint = int(0.5 * (xw @ Ci @ xw))
    assert abs(Qint) == abs(int(round(best["Q"]))), (Qint, best["Q"])
    out = {
        "q": args.q, "n": n, "spectral_ceiling": ceiling,
        "Phi_upper": int(ceiling),
        "witness_Q_exact_int": abs(Qint),
        "witness_Q_signed": Qint,
        "witness_x": xw.tolist(),
        "batch": args.batch, "rounds": args.rounds,
        "prediction_5": "187 or 185",
        "decides_prediction": bool(abs(Qint) >= 188),
    }
    print(json.dumps({k: v for k, v in out.items() if k != "witness_x"}))
    if args.out:
        Path(args.out).write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()
