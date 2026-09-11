"""Multi-restart tabu search for the structural Paley lift family B=C+I.

Any state found is a rigorous LOWER bound on Phi(K), which is what the
divergence question needs: if Phi(K) grows like N^2 then the family cannot
supply asymptotic upper bounds on m_N. Exactness is only claimed where
scripts/paley_structural_lift_family.py --exact can enumerate.
"""

import argparse
import json
import sys

import numpy as np

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from paley_structural_lift_family import lift, paley_conference  # noqa: E402


def search(K, restarts, sweeps, seed):
    rng = np.random.default_rng(seed)
    m = K.shape[0]
    best = 0
    best_z = None
    for _ in range(restarts):
        z = rng.choice([-1, 1], size=m).astype(np.int64)
        g = K @ z  # gradient: flipping i changes the form by -z_i * g_i
        val = int(z @ g) // 2
        tabu = np.zeros(m, dtype=np.int64)
        cur_best = abs(val)
        for t in range(sweeps):
            sign = 1 if val >= 0 else -1
            delta = -sign * 2 * z * g  # change in |form| if flipped (first order)
            delta = np.where(tabu > t, -10**9, delta)
            i = int(np.argmax(delta))
            if delta[i] <= 0 and rng.random() < 0.3:
                i = int(rng.integers(m))
            val -= 2 * int(z[i] * g[i])
            z[i] = -z[i]
            g += 2 * z[i] * K[:, i]
            tabu[i] = t + 3
            if abs(val) > cur_best:
                cur_best = abs(val)
                if cur_best > best:
                    best, best_z = cur_best, z.copy()
    return best, best_z


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("q", type=int, nargs="+")
    ap.add_argument("--restarts", type=int, default=200)
    ap.add_argument("--sweeps", type=int, default=4000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()

    rows = []
    for q in args.q:
        C = paley_conference(q)
        K = lift(C)
        N = K.shape[0]
        val, z = search(K, args.restarts, args.sweeps, args.seed + q)
        assert abs(int(z @ K @ z) // 2) == val, "witness replay failed"
        predicted = N * (N + 12) // 16
        row = {
            "q": q,
            "lifted_order": N,
            "best_found": val,
            "predicted": predicted,
            "matches": val == predicted,
            "exceeds": val > predicted,
            "normalized": round(val / N**1.5, 6),
            "seed": args.seed + q,
            "witness": [int(t) for t in z],
        }
        rows.append(row)
        print(json.dumps({k: v for k, v in row.items() if k != "witness"}), flush=True)

    if args.out:
        with open(args.out, "w") as fh:
            json.dump(rows, fh, indent=2)


if __name__ == "__main__":
    main()
