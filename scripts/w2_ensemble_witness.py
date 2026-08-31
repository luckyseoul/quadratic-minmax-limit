#!/usr/bin/env python3
"""Find U-differences coprime to g from the Max- ensemble (p=5,7)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15406 import gf2_rref, load_minus  # noqa: E402
from e1_gmin_m4_prop15598 import field_ctx  # noqa: E402
from e1_gmin_m4_prop15604 import _primitive  # noqa: E402
from e1_gmin_m4_prop15610 import _dil_fn  # noqa: E402
from e1_gmin_m4_prop15611 import _v2  # noqa: E402
from e1_gmin_m4_prop15612 import _f2_divmod, _f2_factors  # noqa: E402


def nrm(p):
    p = list(p)
    while p and p[-1] == 0:
        p.pop()
    return p or [0]


def poly_gcd(a, b):
    a, b = nrm(a), nrm(b)
    while b and b != [0]:
        _, r = _f2_divmod(a, b)
        a, b = b, nrm(r)
    return nrm(a)


def annihilator(wfn, mul, q):
    N = (q - 1) // 2
    omega = _primitive(mul, q)
    gen = mul(omega, omega)
    cols = []
    cur = wfn.copy()
    for _ in range(N):
        cols.append(cur.copy())
        cur = _dil_fn(cur, mul, gen, q)
        M = np.stack(cols + [cur], axis=1)
        r = gf2_rref(M.copy())[2]
        if r <= len(cols):
            A = np.stack(cols, axis=1)
            Aug = np.concatenate([A, cur.reshape(-1, 1)], axis=1)
            R, pivots, _ = gf2_rref(Aug.copy())
            cof = np.zeros(len(cols), dtype=np.uint8)
            for i, pv in enumerate(pivots):
                if pv < len(cols):
                    cof[pv] = R[i, len(cols)]
            return list(map(int, cof)) + [1]
    return None


def g_poly(p):
    N = (p * p - 1) // 2
    m = N >> _v2(N)
    xm = [0] * m + [1]
    xm[0] = 1
    g, _ = _f2_divmod(xm, [1, 1])
    return g


def run(p, max_test=80):
    Y, C = load_minus(p)
    Y = np.sign(Y.astype(np.float64)).astype(np.int8)
    B = ((1 - Y) // 2).astype(np.uint8)
    fe = np.rint(C[0, 1]).astype(np.int64) * Y[:, 0] * Y[:, 1]
    BU = B[fe < 0]
    q, mul, add, chi, frob, norm, ia, ib = field_ctx(p)
    g = g_poly(p)
    print(f"p={p} |U|={len(BU)} g_deg={len(g)-1}", flush=True)
    y0 = BU[0]
    n_cop = 0
    examples = []
    rng = np.random.default_rng(0)
    idx = rng.choice(len(BU), size=min(max_test, len(BU)), replace=False)
    for j in idx:
        d = (y0 ^ BU[j]) & 1
        wfn = d[1 : 1 + q].copy()
        if d[0]:
            wfn ^= 1
        if not wfn.max():
            continue
        ann = annihilator(wfn, mul, q)
        gg = poly_gcd(ann, g) if ann else None
        if gg == [1]:
            n_cop += 1
            wt = int(wfn.sum())
            # φ-profile
            sinv = None
            sig = next(e for e in range(1, q) if chi(e) == -1)
            from e1_gmin_m4_prop15613 import _finv

            sinv = _finv(mul, q, sig)
            prof = [0] * p
            for x in range(q):
                if wfn[x]:
                    prof[mul(sinv, x) // p] += 1
            examples.append(
                {
                    "j": int(j),
                    "wt": wt,
                    "w0": int(wfn[0]),
                    "phi_counts": prof,
                    "ann_deg": len(ann) - 1,
                }
            )
            print(f"  COPRIME j={j} wt={wt} w0={wfn[0]} phi={prof}", flush=True)
            if n_cop >= 5:
                break
    print(f"  n_coprime among {len(idx)} = {n_cop}", flush=True)
    return {"p": p, "nU": int(len(BU)), "tested": int(len(idx)), "n_cop": n_cop, "ex": examples}


def main():
    recs = {}
    for p, n in ((5, 80), (7, 40)):
        recs[str(p)] = run(p, max_test=n)
    dest = ROOT / "evidence" / "w2_ensemble_witness.json"
    import json

    dest.write_text(json.dumps(recs, indent=2, default=str) + "\n")
    print("wrote", dest, flush=True)


if __name__ == "__main__":
    main()
