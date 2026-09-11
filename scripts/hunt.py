#!/usr/bin/env python3
"""Local search for low-Phi signings via the linear-response representation.

Q_K(s) = sum_{i<j} K_ij * s_i s_j is LINEAR in K. With G_a(s) = s_i s_j
precomputed per pair a, a move that flips a few entries evaluates as
Q' = Q - 2*sum_{a in move} K_a * G_a  -- pure vector arithmetic on the
Q vector, no per-move re-enumeration. This makes deep neighborhoods
(1-, 2-, 3-flip) exact and cheap.

Modes:
  general:    all n(n-1)/2 pairs free.
  --fixed-block PATH: pin A (top-left) and -A (bottom-right); search only
              the m x m cross block B (the completion family).
"""
import argparse
import json
import time

import numpy as np


def from_upper_string(s):
    L = len(s)
    n = 2
    while n * (n - 1) // 2 < L:
        n += 1
    A = np.zeros((n, n), dtype=np.float32)
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            A[i, j] = A[j, i] = 1.0 if s[k] == "1" else -1.0
            k += 1
    return A


def load_json_A(path):
    j = json.load(open(path))
    if "A" in j:
        return np.array(j["A"], dtype=np.float32), j.get("Phi")
    if "source_upper" in j:
        return from_upper_string(j["source_upper"]), j.get("Phi")
    raise ValueError(f"no A/source_upper in {path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", action="append", default=[])
    ap.add_argument("--fixed-block", default=None)
    ap.add_argument("--order", type=int, default=None)
    ap.add_argument("--minutes", type=float, default=10.0)
    ap.add_argument("--kick-seed", type=int, default=1)
    ap.add_argument("--out", default=None)
    ap.add_argument("--no-random", action="store_true")
    args = ap.parse_args()

    rng = np.random.default_rng(args.kick_seed)

    try:
        import cupy as xp
        dev = xp.cuda.runtime.getDeviceProperties(0)["name"].decode()
        backend = f"cupy {xp.__version__} on {dev}"
    except Exception:
        import numpy as xp
        backend = "numpy " + np.__version__
    deep = backend.startswith("cupy")

    mode = "fixed-block" if args.fixed_block else "general"

    if args.fixed_block:
        A = np.array(json.load(open(args.fixed_block))["A"], dtype=np.float32)
        m = A.shape[0]
        n = 2 * m
        pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        free_idx = np.array([a for a, (i, j) in enumerate(pairs)
                             if (i < m) != (j < m)], dtype=np.int64)
    else:
        if args.seed:
            K0, _ = load_json_A(args.seed[0])
            n = K0.shape[0]
        elif args.order:
            n = args.order
        else:
            raise SystemExit("need --seed, --order, or --fixed-block")
        pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        free_idx = np.arange(len(pairs), dtype=np.int64)

    P = len(free_idx)
    N = 1 << (n - 1)
    print(f"[hunt] order={n} mode={mode} pairs={P} states={N} backend={backend}",
          flush=True)

    u = xp.arange(N, dtype=xp.int64)
    S = xp.ones((N, n), dtype=xp.float32)
    for c in range(1, n):
        S[:, c] = 1 - 2 * ((u >> (c - 1)) & 1)
    G = xp.empty((len(pairs), N), dtype=xp.float32)
    for a, (i, j) in enumerate(pairs):
        G[a] = S[:, i] * S[:, j]
    Gf = G[free_idx]
    del u, S

    def K_to_Kv(K):
        return xp.asarray([K[i, j] for i, j in pairs], dtype=xp.float32)

    def Kv_to_K(Kv):
        K = np.zeros((n, n), dtype=np.float32)
        Kv_c = xp.asnumpy(Kv)
        for a, (i, j) in enumerate(pairs):
            K[i, j] = K[j, i] = Kv_c[a]
        return K

    def evalK(Kv):
        Q = Kv @ G
        return Q, float(xp.max(xp.abs(Q)))

    def descend1(Kv, Q, val):
        steps = 0
        while True:
            M = xp.abs(Q[None, :] - 2.0 * Kv[free_idx][:, None] * Gf)
            mv = M.max(axis=1)
            j = int(mv.argmin())
            if float(mv[j]) < val - 1e-7:
                a = int(free_idx[j])
                Q = Q - (2.0 * Kv[a]) * G[a]
                Kv[a] = -Kv[a]
                val = float(mv[j])
                steps += 1
                del M, mv
            else:
                del M, mv
                break
        return Kv, Q, val, steps

    chunk = max(64, (1 << 26) // N)

    def best2(Kv, Q, val):
        idxs = [(a, b) for a in range(P) for b in range(a + 1, P)]
        best_val, best = val, None
        for s in range(0, len(idxs), chunk):
            part = idxs[s:s + chunk]
            ia = xp.asarray([int(free_idx[a]) for a, b in part])
            ib = xp.asarray([int(free_idx[b]) for a, b in part])
            R = (Q[None, :] - 2.0 * Kv[ia][:, None] * G[ia]
                 - 2.0 * Kv[ib][:, None] * G[ib])
            mv = xp.max(xp.abs(R), axis=1)
            k = int(mv.argmin())
            v = float(mv[k])
            if v < best_val - 1e-7:
                best_val = v
                best = (int(free_idx[part[k][0]]), int(free_idx[part[k][1]]))
            del R, mv
        return best_val, best

    def best3(Kv, Q, val):
        best_val, best = val, None
        idxs = [(a, b, c) for a in range(P) for b in range(a + 1, P)
                for c in range(b + 1, P)]
        for s in range(0, len(idxs), chunk):
            part = idxs[s:s + chunk]
            ia = xp.asarray([int(free_idx[a]) for a, b, c in part])
            ib = xp.asarray([int(free_idx[b]) for a, b, c in part])
            ic = xp.asarray([int(free_idx[c]) for a, b, c in part])
            R = (Q[None, :] - 2.0 * Kv[ia][:, None] * G[ia]
                 - 2.0 * Kv[ib][:, None] * G[ib]
                 - 2.0 * Kv[ic][:, None] * G[ic])
            mv = xp.max(xp.abs(R), axis=1)
            k = int(mv.argmin())
            v = float(mv[k])
            if v < best_val - 1e-7:
                best_val = v
                best = (int(free_idx[part[k][0]]), int(free_idx[part[k][1]]),
                        int(free_idx[part[k][2]]))
            del R, mv
        return best_val, best

    def flip(Q, Kv, move):
        for a in move:
            Q = Q - (2.0 * Kv[a]) * G[a]
            Kv[a] = -Kv[a]
        return Q, Kv

    state = {"bestVal": None, "bestKv": None}

    def record(Kv, val, tag):
        if state["bestVal"] is None or val < state["bestVal"] - 1e-9:
            state["bestVal"] = val
            state["bestKv"] = Kv.copy()
            K = Kv_to_K(Kv)
            val2 = float(xp.max(xp.abs(K_to_Kv(K) @ G)))
            out = dict(A=K.tolist(), Phi=val, Phi_fresh=val2, order=n, mode=mode,
                       backend=backend, kick_seed=args.kick_seed, tag=tag,
                       elapsed_s=round(time.time() - t0, 1))
            path = args.out or f"hunt{n}_{mode}_{args.kick_seed}.json"
            json.dump(out, open(path, "w"))
            print(f"[NEW BEST] Phi={val:.6f} (fresh {val2:.6f}) tag={tag} "
                  f"t={out['elapsed_s']}s -> {path}", flush=True)

    t0 = time.time()
    budget = args.minutes * 60

    # build seeds
    Kv_seeds = []
    if args.fixed_block:
        m8 = A
        Kv_seeds.append(K_to_Kv(np.block([[m8, m8], [m8.T, -m8]])))
        Kv_seeds.append(K_to_Kv(np.block([[m8, -m8], [-m8.T, -m8]])))
        for _ in range(6):
            B = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), (m, m))
            Kb = np.block([[m8, B], [B.T, -m8]])
            Kv_seeds.append(K_to_Kv(Kb))
    else:
        for p in args.seed:
            K, phi = load_json_A(p)
            Kv = K_to_Kv(K)
            Q, val = evalK(Kv)
            print(f"[seed] {p} phi_file={phi} phi_eval={val:.6f}", flush=True)
            Kv_seeds.append(Kv)
    if not args.no_random:
        for _ in range(3 if not args.fixed_block else 0):
            K = np.zeros((n, n), dtype=np.float32)
            ku = np.triu_indices(n, 1)
            K[ku] = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), len(ku[0]))
            K = K + K.T
            Kv_seeds.append(K_to_Kv(K))

    rounds = 0
    for Kv0 in Kv_seeds:
        if time.time() - t0 > budget:
            break
        Q, val = evalK(Kv0)
        Kv, Q, val, st = descend1(Kv0, Q, val)
        record(Kv, val, f"descent{rounds}")
        for _ in range(2):
            v2, mv2 = best2(Kv, Q, val)
            if mv2 is None:
                break
            Q, Kv = flip(Q, Kv, mv2)
            val = float(xp.max(xp.abs(Q)))
            Kv, Q, val, _ = descend1(Kv, Q, val)
            record(Kv, val, "2flip")
        rounds += 1

    print(f"[hunt] seed phase done, rounds={rounds}, best={state['bestVal']} "
          f"t={time.time()-t0:.1f}s", flush=True)

    if state["bestVal"] is None:
        print("[hunt] no seeds evaluated", flush=True)
        return

    Kv = state["bestKv"].copy()
    Q, val = evalK(Kv)

    while time.time() - t0 < budget:
        rounds += 1
        if rounds % 40 == 0 and not args.no_random:
            if args.fixed_block:
                B = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), (m, m))
                K = np.block([[A, B], [B.T, -A]])
            else:
                K = np.zeros((n, n), dtype=np.float32)
                ku = np.triu_indices(n, 1)
                K[ku] = rng.choice(np.array([-1.0, 1.0], dtype=np.float32),
                                   len(ku[0]))
                K = K + K.T
            Kv = K_to_Kv(K)
        else:
            Kv = Kv.copy()
            k = int(rng.integers(3, 10))
            for _ in range(k):
                a = int(free_idx[int(rng.integers(0, P))])
                Q = Q - (2.0 * Kv[a]) * G[a]
                Kv[a] = -Kv[a]
        Q, val = evalK(Kv)
        Kv, Q, val, st = descend1(Kv, Q, val)
        record(Kv, val, f"ils{rounds}")
        for _ in range(3):
            v2, mv2 = best2(Kv, Q, val)
            if mv2 is None:
                break
            Q, Kv = flip(Q, Kv, mv2)
            val = float(xp.max(xp.abs(Q)))
            Kv, Q, val, _ = descend1(Kv, Q, val)
            record(Kv, val, "2flip")
        if deep:
            v3, mv3 = best3(Kv, Q, val)
            if mv3 is not None:
                Q, Kv = flip(Q, Kv, mv3)
                val = float(xp.max(xp.abs(Q)))
                Kv, Q, val, _ = descend1(Kv, Q, val)
                record(Kv, val, "3flip")
        if rounds % 5 == 0:
            print(f"[ils] round={rounds} val={val:.6f} "
                  f"best={state['bestVal']:.6f} t={time.time()-t0:.0f}s",
                  flush=True)

    print(f"[hunt] done rounds={rounds} best={state['bestVal']} "
          f"t={time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
