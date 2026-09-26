#!/usr/bin/env python3
"""
Exact discriminator for the live Banaszczyk one-vertex route.

For a supplied signing A, compute the minimax depth of the ZERO-DEFICIT,
SAME-PHASE compatibility recursion.

At a current principal block B and phase sigma, take all projective global
maximizers of sigma*Q_B.  For a chosen pivot coordinate, orient every
maximizer so its pivot spin is +1.  Every pair produces the disagreement
principal block T; the stable-pair factorization says its restricted state is
again an exact sigma-maximizer on B[T].  Recursing gives the number of genuine
compatibility generations that can consume one raw unit of slab width without
paying any stable deficit.

We minimize the worst child depth over the pivot choice.  Thus:
    depth < r
means the pure zero-deficit same-phase branch cannot exhaust an initial
extension buffer r under an optimal elimination order.
    depth >= r
kills that hoped-for closure mechanism on this matrix.

This is a finite discriminator, not an asymptotic proof.
"""

from __future__ import annotations
import argparse, json, math, os, sys
from functools import lru_cache
from pathlib import Path

import numpy as np

def load_matrix(path: str) -> np.ndarray:
    p = Path(path)
    txt = p.read_text()
    try:
        obj = json.loads(txt)
        if isinstance(obj, dict):
            for key in ("A", "matrix", "M"):
                if key in obj:
                    arr = np.asarray(obj[key], dtype=np.int16)
                    break
            else:
                raise ValueError("JSON has no A/matrix/M field")
        else:
            arr = np.asarray(obj, dtype=np.int16)
    except Exception:
        rows = []
        for line in txt.splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            vals = [int(x) for x in line.replace(",", " ").split()]
            rows.append(vals)
        arr = np.asarray(rows, dtype=np.int16)
    if arr.ndim != 2 or arr.shape[0] != arr.shape[1]:
        raise ValueError(f"matrix must be square, got {arr.shape}")
    if np.any(arr != arr.T):
        raise ValueError("matrix is not symmetric")
    if np.any(np.diag(arr) != 0):
        raise ValueError("diagonal must be zero")
    off = arr[~np.eye(arr.shape[0], dtype=bool)]
    if not np.all(np.isin(off, (-1, 1))):
        raise ValueError("off-diagonal entries must be +/-1")
    return arr

class Engine:
    def __init__(self, A: np.ndarray, chunk: int, gpu: bool):
        self.A = A
        self.n = A.shape[0]
        self.chunk = chunk
        self.gpu = False
        self.cp = None
        if gpu:
            try:
                import cupy as cp
                self.cp = cp
                self.gpu = True
            except Exception as e:
                print(f"[warn] CuPy unavailable, using NumPy: {e}", file=sys.stderr)

        self.max_cache = {}
        self.depth_cache = {}
        self.calls = 0

    def _enumerate_extrema(self, subset: int):
        """Return (P,N,pos_masks,neg_masks) on this principal block.
        Masks are local m-bit projective reps with local bit 0 fixed to 0 (+1).
        """
        if subset in self.max_cache:
            return self.max_cache[subset]
        verts = [i for i in range(self.n) if (subset >> i) & 1]
        m = len(verts)
        if m <= 1:
            out = (0, 0, [0], [0])
            self.max_cache[subset] = out
            return out

        B_np = self.A[np.ix_(verts, verts)].astype(np.int16, copy=False)
        total = 1 << (m - 1)
        P = -10**18
        Qmin = 10**18
        pos_masks = []
        neg_masks = []

        xp = self.cp if self.gpu else np
        B = xp.asarray(B_np) if self.gpu else B_np
        shifts = xp.arange(m - 1, dtype=xp.uint64)

        for start in range(0, total, self.chunk):
            stop = min(total, start + self.chunk)
            idx = xp.arange(start, stop, dtype=xp.uint64)
            bits = ((idx[:, None] >> shifts[None, :]) & 1).astype(xp.int16)
            X = xp.empty((stop - start, m), dtype=xp.int16)
            X[:, 0] = 1
            X[:, 1:] = 1 - 2 * bits
            H = X @ B
            q = (xp.sum(H * X, axis=1) // 2).astype(xp.int64)

            qmax = int(q.max().item())
            qmin = int(q.min().item())

            if qmax > P:
                P = qmax
                pos_masks = []
            if qmax == P:
                loc = xp.where(q == qmax)[0]
                loc = xp.asnumpy(loc) if self.gpu else loc
                for z in loc.tolist():
                    pos_masks.append(start + int(z))

            if qmin < Qmin:
                Qmin = qmin
                neg_masks = []
            if qmin == Qmin:
                loc = xp.where(q == qmin)[0]
                loc = xp.asnumpy(loc) if self.gpu else loc
                for z in loc.tolist():
                    neg_masks.append(start + int(z))

        # Enumeration index uses bits 0..m-2 for local coordinates 1..m-1.
        pos_masks = [int(v) << 1 for v in pos_masks]
        neg_masks = [int(v) << 1 for v in neg_masks]
        out = (int(P), int(-Qmin), pos_masks, neg_masks)
        self.max_cache[subset] = out
        return out

    def phi(self, subset: int):
        P,N,_,_ = self._enumerate_extrema(subset)
        return max(P,N),P,N

    @staticmethod
    def _local_to_global(local_mask: int, verts: list[int]) -> int:
        g = 0
        while local_mask:
            lsb = local_mask & -local_mask
            j = lsb.bit_length() - 1
            g |= 1 << verts[j]
            local_mask ^= lsb
        return g

    def depth(self, subset: int, sigma: int):
        key=(subset,sigma)
        if key in self.depth_cache:
            return self.depth_cache[key]
        self.calls += 1
        verts=[i for i in range(self.n) if (subset>>i)&1]
        m=len(verts)
        if m<=1:
            ans=(0,None,None)
            self.depth_cache[key]=ans
            return ans

        P,N,pos,neg=self._enumerate_extrema(subset)
        states = pos if sigma==1 else neg
        states = list(dict.fromkeys(states))
        if len(states)<=1:
            ans=(0,None,None)
            self.depth_cache[key]=ans
            return ans

        full=(1<<m)-1
        best_depth=10**9
        best_pivot=None
        best_child=None

        for p in range(m):
            oriented=[]
            for s in states:
                if (s>>p)&1:
                    s ^= full
                oriented.append(s)
            oriented=list(dict.fromkeys(oriented))
            child_masks=set()
            L=len(oriented)
            for a in range(L):
                sa=oriented[a]
                for b in range(a+1,L):
                    t=sa ^ oriented[b]
                    if t:
                        # By construction pivot p agrees, so t is proper.
                        child_masks.add(t)
            if not child_masks:
                worst=0
                worst_child=None
            else:
                worst=-1
                worst_child=None
                # Large children first: they are usually the dangerous ones.
                ordered=sorted(child_masks, key=int.bit_count, reverse=True)
                for tloc in ordered:
                    tg=self._local_to_global(tloc,verts)
                    d,_,_=self.depth(tg,sigma)
                    cand=1+d
                    if cand>worst:
                        worst=cand
                        worst_child=tg
                    if worst>=best_depth:
                        break
            if worst<best_depth:
                best_depth=worst
                best_pivot=verts[p]
                best_child=worst_child

        ans=(int(best_depth),best_pivot,best_child)
        self.depth_cache[key]=ans
        return ans

    def witness_chain(self, subset: int, sigma: int):
        chain=[]
        while subset and subset.bit_count()>1:
            d,p,ch=self.depth(subset,sigma)
            if d<=0 or ch is None:
                break
            chain.append({
                "depth_remaining": d,
                "order": subset.bit_count(),
                "pivot": p,
                "child_order": ch.bit_count(),
                "subset_hex": hex(subset),
                "child_hex": hex(ch),
            })
            subset=ch
        return chain

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("matrix", help="JSON containing A, or whitespace square matrix")
    ap.add_argument("--chunk", type=int, default=262144)
    ap.add_argument("--gpu", action="store_true", help="use CuPy if available")
    ap.add_argument("--r", type=float, default=None, help="override extension buffer")
    args=ap.parse_args()

    A=load_matrix(args.matrix)
    n=A.shape[0]
    E=Engine(A,args.chunk,args.gpu)
    root=(1<<n)-1
    F,P,N=E.phi(root)
    r = args.r if args.r is not None else F*((1.0+1.0/n)**1.5-1.0)

    phases=[]
    for sigma,name,edge in ((1,"+",P),(-1,"-",N)):
        if edge != F:
            phases.append({
                "phase": name,
                "root_edge": edge,
                "root_deficit": F-edge,
                "zero_deficit_present": False,
                "depth": 0,
                "margin_r_minus_depth": r,
                "chain": [],
            })
            continue
        d,p,ch=E.depth(root,sigma)
        phases.append({
            "phase": name,
            "root_edge": edge,
            "root_deficit": 0,
            "zero_deficit_present": True,
            "depth": d,
            "margin_r_minus_depth": r-d,
            "safe_for_pure_zero_deficit_branch": bool(d < r),
            "best_first_pivot": p,
            "chain": E.witness_chain(root,sigma),
        })

    out={
        "matrix": args.matrix,
        "n": n,
        "Phi": F,
        "P": P,
        "N": N,
        "neutral_r": r,
        "criterion": "pure zero-deficit same-phase branch is safe iff minimax depth < r",
        "phases": phases,
        "principal_blocks_evaluated": len(E.max_cache),
        "depth_states_evaluated": len(E.depth_cache),
        "gpu": E.gpu,
    }
    print(json.dumps(out,indent=2))

if __name__=="__main__":
    main()
