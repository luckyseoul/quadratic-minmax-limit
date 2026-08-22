#!/usr/bin/env python3
"""p=13 k=5 (CPU gauged) and k=6 (GPU gauged) drivers from p13p17enum.md.

Outputs go to /mnt/storage/e1work/maxplus_p13/ (not tmpfs).
Usage:
  python3 -u run_p13_stratum.py validate
  python3 -u run_p13_stratum.py k5
  python3 -u run_p13_stratum.py k6
  python3 -u run_p13_stratum.py probe
"""
from __future__ import annotations

import itertools
import os
import sys
import time
import multiprocessing as mp
from multiprocessing import Pool
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

TOOLKIT = Path(__file__).resolve().parent
sys.path.insert(0, str(TOOLKIT))
# modules still insert /tmp/e1work; keep a symlink farm there.

from kgen import square_coords  # noqa: E402
from kgen3 import prep_subset  # noqa: E402
from kgen4 import enum_chunk  # noqa: E402
from kgen6 import enum_gauged_task, translation_tables  # noqa: E402

OUT = Path("/mnt/storage/e1work/maxplus_p13")
OUT.mkdir(parents=True, exist_ok=True)

# GPU worker globals (initializer); Python 3.14 defaults away from fork.
_GPU_CTXS = None
_GPU_P = None
_GPU_K = None
_GPU_Q = None


def _gpu_init(ctxs, p, k):
    global _GPU_CTXS, _GPU_P, _GPU_K, _GPU_Q
    _GPU_CTXS = ctxs
    _GPU_P = p
    _GPU_K = k
    _GPU_Q = p * p


def _gpu_worker(args):
    import numpy as np
    from gpu_inner import GpuTester, process_outer_gpu
    from kgen5 import _activity_filter

    phase, sub, cf, tvidx = args
    ctx = _GPU_CTXS[sub]
    tester = GpuTester(_GPU_P, _GPU_K, ctx["Tm"], ctx["UU"])
    s_ar = np.arange(_GPU_P, dtype=np.int64)
    upper = np.zeros((_GPU_K, _GPU_P), dtype=np.int64)
    for d, vec in cf.items():
        upper = (upper + np.outer(vec, (s_ar ** d) % _GPU_P)) % _GPU_P
    sols = []
    t1 = time.time()
    nc = process_outer_gpu(
        _GPU_P, _GPU_K, _GPU_Q, upper, ctx["UU"], ctx["Tm"], ctx["c0"], 1, tester, sols
    )
    reps = _activity_filter(sols, ctx["Tm"], _GPU_P, _GPU_K, 1)
    R = np.stack(reps).astype(np.int8) if reps else np.zeros((0, _GPU_Q), np.int8)
    return phase, sub, tvidx, R, nc, time.time() - t1


def _stack_save(path: Path, rows: list) -> int:
    if not rows:
        print(f"  empty -> skip save {path}", flush=True)
        return 0
    A = __import__("numpy").stack(rows)
    npy = __import__("numpy")
    view = A.reshape(A.shape[0], -1).view(npy.dtype((npy.void, A.dtype.itemsize * A.shape[1])))
    _, idx = npy.unique(view, return_index=True)
    B = A[npy.sort(idx)]
    npy.save(path, B.astype(npy.int8))
    with open(path, "rb") as fh:
        os.fsync(fh.fileno())
    print(f"  wrote {path} rows={len(B)} distinct_from={len(A)}", flush=True)
    return int(len(B))


def run_cpu_chunks(p: int, k: int, workers: int, out_name: str) -> int:
    """Ungauged staged CPU enum (kgen4). Ground-truth path."""
    import numpy as np

    dirs, forms, coords = square_coords(p)
    tasks = []
    for sub in itertools.combinations(range(len(dirs)), k):
        ctx = prep_subset(p, list(sub), forms, coords)
        tot = ctx["outer_total"]
        step = max(1, tot // (workers * 4))
        for lo in range(0, tot, step):
            tasks.append((ctx, lo, min(lo + step, tot), 1))
    print(
        f"CPU ungauged p={p} k={k} subsets={len(list(itertools.combinations(range(len(dirs)), k)))} "
        f"tasks={len(tasks)} workers={workers}",
        flush=True,
    )
    t0 = time.time()
    allsols: list = []
    per_sub: dict = {}
    with Pool(workers) as pool:
        for sub, lo, hi, sols in pool.imap_unordered(enum_chunk, tasks):
            per_sub.setdefault(sub, []).extend(sols)
    for sub in sorted(per_sub):
        print(f"  subset {sub}: {len(per_sub[sub])}", flush=True)
        allsols.extend(per_sub[sub])
    print(f"p={p} k={k} TOTAL={len(allsols)} time={time.time() - t0:.1f}s", flush=True)
    return _stack_save(OUT / out_name, allsols)


def run_cpu_gauged(p: int, k: int, workers: int, out_name: str) -> int:
    """Translation-gauged CPU enum (kgen6). k>=5."""
    dirs, forms, coords = square_coords(p)
    TT = translation_tables(p)
    tasks = []
    nsub = 0
    for sub in itertools.combinations(range(len(dirs)), k):
        nsub += 1
        ctx = prep_subset(p, list(sub), forms, coords)
        deg = k - 2
        tsize = p - 1
        for d in range(deg - 2, 1, -1):
            tsize *= p ** len(ctx["kern"][d])
        step = max(1, tsize // max(workers, 1))
        for lo in range(0, tsize, step):
            tasks.append((ctx, "T", lo, min(lo + step, tsize), 1, TT))
        lsize = 1
        for d in range(deg - 1, 1, -1):
            lsize *= p ** len(ctx["kern"][d])
        step = max(1, lsize // max(workers, 1))
        for lo in range(0, lsize, step):
            tasks.append((ctx, "L", lo, min(lo + step, lsize), 1, TT))
    print(
        f"CPU gauged p={p} k={k} subsets={nsub} tasks={len(tasks)} workers={workers}",
        flush=True,
    )
    t0 = time.time()
    per: dict = {}
    with Pool(workers) as pool:
        for sub, phase, lo, hi, sols, nre in pool.imap_unordered(enum_gauged_task, tasks):
            per.setdefault(sub, []).extend(sols)
    allsols: list = []
    for sub in sorted(per):
        print(f"  subset {sub}: {len(per[sub])}", flush=True)
        allsols.extend(per[sub])
    print(f"p={p} k={k} gauged TOTAL={len(allsols)} time={time.time() - t0:.1f}s", flush=True)
    return _stack_save(OUT / out_name, allsols)


def run_gpu_gauged(p: int, k: int, nw: int) -> int:
    """Dilation+translation gauged GPU enum. Copied from run_kgauged.py with p/outdir."""
    import numpy as np
    from dilation import build_group, orbits
    from kgen6 import translation_tables as ttables

    q = p * p
    outdir = OUT / f"k{k}_gpu_out"
    outdir.mkdir(parents=True, exist_ok=True)

    def lattice(basis):
        if len(basis) == 0:
            return [np.zeros(0, dtype=np.int64)]
        basis = np.array(basis)
        dim = len(basis)
        out = []
        for combo in itertools.product(range(p), repeat=dim):
            v = np.zeros(basis.shape[1], dtype=np.int64)
            for c, b in zip(combo, basis):
                v = (v + c * b) % p
            out.append(v)
        return out

    dirs, forms, coords = square_coords(p)
    m = len(dirs)
    subsets = list(itertools.combinations(range(m), k))
    ctxs = {sub: prep_subset(p, list(sub), forms, coords) for sub in subsets}
    group, gcoords = build_group(p)
    TT = ttables(p)
    deg = k - 2
    T, L = [], []
    for sub in subsets:
        kern = ctxs[sub]["kern"]
        low_levels = [d for d in range(deg - 2, 1, -1)]
        lows = [lattice(kern[d]) for d in low_levels]
        for lam in range(1, p):
            for combo in itertools.product(*lows):
                cf = {deg: (lam * kern[deg][0]) % p}
                if deg - 1 >= 2:
                    cf[deg - 1] = np.zeros(k, dtype=np.int64)
                for d, vec in zip(low_levels, combo):
                    cf[d] = vec
                T.append((sub, cf))
        levs = [d for d in range(deg - 1, 1, -1)]
        lats = [lattice(kern[d]) for d in levs]
        for combo in itertools.product(*lats):
            cf = {deg: np.zeros(k, dtype=np.int64)}
            for d, vec in zip(levs, combo):
                cf[d] = vec
            L.append((sub, cf))
    print(f"states: T={len(T)} L={len(L)} p={p} k={k} subsets={len(subsets)}", flush=True)
    t0 = time.time()
    orbT = orbits(T, group, p)
    orbL = orbits(L, group, p)
    print(f"orbits: T={len(orbT)} L={len(orbL)}  ({time.time() - t0:.0f}s)", flush=True)

    tasks_all = []
    tvstore = []
    for phase, orb in (("T", orbT), ("L", orbL)):
        for (sub, cf), tv in orb:
            tvstore.append(tv)
            tasks_all.append((phase, sub, cf, len(tvstore) - 1))
    total = 0
    done_prior = 0
    tasks = []
    for t in tasks_all:
        tvidx = t[3]
        fp = outdir / f"orb{tvidx}.npy"
        if fp.exists():
            total += len(np.load(fp))
            done_prior += 1
        else:
            tasks.append(t)
    if done_prior:
        print(
            f"RESUMED: {done_prior}/{len(tasks_all)} outers on disk "
            f"({total} sols); {len(tasks)} remaining",
            flush=True,
        )
    done = done_prior
    # GPU_WORKERS env caps CuPy pool; keep nw == GPU_WORKERS
    with Pool(nw, initializer=_gpu_init, initargs=(ctxs, p, k)) as pool:
        for phase, sub, tvidx, R, nc, dt in pool.imap_unordered(_gpu_worker, tasks):
            tv = tvstore[tvidx]
            outs = []
            if len(R):
                for g in tv:
                    pi = np.argsort(g[0])
                    Rg = R[:, pi]
                    if phase == "T":
                        for y in Rg:
                            outs.append(y[TT])
                    else:
                        outs.append(Rg)
            fp = outdir / f"orb{tvidx}.npy"
            if outs:
                A = np.concatenate([o if o.ndim == 2 else o[None, :] for o in outs])
                cnt_here = len(A)
                np.save(fp, A.astype(np.int8))
            else:
                cnt_here = 0
                np.save(fp, np.zeros((0, q), np.int8))
            with open(fp, "rb") as fh:
                os.fsync(fh.fileno())
            total += cnt_here
            done += 1
            print(
                f"[{done}/{len(tasks_all)}] {phase} tvidx={tvidx} reps={len(R)} |tv|={len(tv)} "
                f"-> {cnt_here}  cand={nc} {dt:.0f}s  cum={total}",
                flush=True,
            )
    print(f"k={k} gauged GPU TOTAL raw={total}", flush=True)
    parts = [np.load(outdir / f"orb{i}.npy") for i in range(len(tasks_all))]
    nonempty = [a for a in parts if len(a)]
    if not nonempty:
        print("no solutions", flush=True)
        return 0
    A = np.concatenate(nonempty)
    return _stack_save(OUT / f"k{k}_p{p}_full.npy", list(A))


def probe(p: int, k: int) -> None:
    dirs, forms, coords = square_coords(p)
    m = len(dirs)
    nsub = 0
    tot_T = tot_L = 0
    for sub in itertools.combinations(range(m), k):
        nsub += 1
        ctx = prep_subset(p, list(sub), forms, coords)
        deg = k - 2
        tsize = p - 1
        for d in range(deg - 2, 1, -1):
            tsize *= p ** len(ctx["kern"][d])
        lsize = 1
        for d in range(deg - 1, 1, -1):
            lsize *= p ** len(ctx["kern"][d])
        tot_T += tsize
        tot_L += lsize
        if nsub <= 3:
            print(
                f"  sub{sub} outer_total={ctx['outer_total']} T={tsize} L={lsize} "
                f"kern_dims={{d:len(ctx['kern'][d]) for d in ctx['kern']}}".replace(
                    "{d:len(ctx['kern'][d]) for d in ctx['kern']}",
                    str({d: int(ctx["kern"][d].shape[0]) for d in ctx["kern"]}),
                ),
                flush=True,
            )
    print(
        f"p={p} k={k} m={m} subsets={nsub} sum_T={tot_T} sum_L={tot_L} "
        f"sum_outer={tot_T + tot_L}",
        flush=True,
    )


def validate() -> None:
    """p=5 -> 130 eps=+1; p=7 per-k {1:140, 3:1176, 4:4410}."""
    # k=1 is closed-form in the ledger; enum k=3,4 which exercise the pipeline.
    n5_3 = run_cpu_chunks(5, 3, workers=8, out_name="gt_k3_p5.npy")
    n7_3 = run_cpu_chunks(7, 3, workers=16, out_name="gt_k3_p7.npy")
    n7_4 = run_cpu_chunks(7, 4, workers=32, out_name="gt_k4_p7.npy")
    ok = n5_3 == 100 and n7_3 == 1176 and n7_4 == 4410
    # k=1 closed form
    def k1(p):
        m = (p + 1) // 2
        from math import comb

        return m * comb(p, m)

    print(f"k=1 closed form p=5 {k1(5)} expect 30; p=7 {k1(7)} expect 140", flush=True)
    print(
        f"p=5 k=3+k=1 = {n5_3 + k1(5)} expect 130; "
        f"p=7 k=1+3+4 = {k1(7)+n7_3+n7_4} expect 5726",
        flush=True,
    )
    if not ok:
        raise SystemExit(
            f"GROUND TRUTH FAIL k3_p5={n5_3} k3_p7={n7_3} k4_p7={n7_4}"
        )
    if n5_3 + k1(5) != 130 or k1(7) + n7_3 + n7_4 != 5726:
        raise SystemExit("GROUND TRUTH FAIL half-set totals")
    print("GROUND TRUTH OK", flush=True)


def main() -> None:
    try:
        mp.set_start_method("fork")
    except RuntimeError:
        pass
    if len(sys.argv) < 2:
        print(__doc__)
        raise SystemExit(2)
    cmd = sys.argv[1]
    if cmd == "validate":
        validate()
    elif cmd == "probe":
        for k in (4, 5, 6, 7):
            probe(13, k)
    elif cmd == "k5":
        # p=11 k=5 gauged was 1825s / 84w. p=13 has 21 vs 6 subsets.
        # Estimate ~3–6 h at 50 workers. Split CPU with k=6 GPU (3 workers).
        w = int(os.environ.get("K5_WORKERS", "50"))
        print(f"LAUNCH p=13 k=5 CPU gauged workers={w}  est. 3-6 hours", flush=True)
        run_cpu_gauged(13, 5, workers=w, out_name="k5_p13_full.npy")
    elif cmd == "k6":
        nw = int(os.environ.get("GPU_WORKERS", "3"))
        print(
            f"LAUNCH p=13 k=6 GPU gauged GPU_WORKERS={nw}  "
            f"est. scale from p=11 ~16h (1715 outers * ~100s / 3); "
            f"p=13 has 7 subsets vs 1, so several hours to ~1-2 days",
            flush=True,
        )
        run_gpu_gauged(13, 6, nw=nw)
    else:
        raise SystemExit(f"unknown cmd {cmd}")


if __name__ == "__main__":
    main()
