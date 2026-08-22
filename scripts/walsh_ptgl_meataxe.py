#!/usr/bin/env python3
"""PΓL-submodules of H0=ker S, and whether 1_QR is constant on U.

If H0/⟨1⟩ is irreducible as a PΓL-module, Max- spans H0.
If 1_QR is nonconstant on U, it is not a Walsh character.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15406 import gf2_nullspace, gf2_rref, load_minus  # noqa: E402
from e1_gmin_m4_prop15598 import field_ctx  # noqa: E402
from walsh_linecode_rank import (  # noqa: E402
    _mobius_perm,
    aut_e_generators,
    square_line_matrix,
)
from walsh_rank_equality import psl_generators  # noqa: E402


def action_matrices(H0, gens):
    """Each gen permutes coordinates; induce F2-linear map on H0 columns."""
    n, k = H0.shape
    # H0 is n x k, columns a basis. g acts on vectors by permuting coords:
    # v |-> v[g]. In basis: H0 t |-> H0[g,:] t = H0 (M t), so H0 M = H0[g,:].
    mats = []
    for g in gens:
        Hg = H0[g, :]
        # solve H0 M = Hg, i.e. each column of Hg is H0 @ M[:,j]
        M = np.zeros((k, k), dtype=np.uint8)
        for j in range(k):
            # rref augment
            A = np.concatenate([H0, Hg[:, j : j + 1]], axis=1)
            R, pivots, rank = gf2_rref(A)
            # H0 has full column rank k
            x = np.zeros(k, dtype=np.uint8)
            prow = {pivots[i]: i for i in range(min(rank, k))}
            for col, i in prow.items():
                if col < k:
                    x[col] = R[i, k]
            M[:, j] = x
        mats.append(M)
    return mats


def cyclic_span(mats, v):
    """Dim of F2⟨G v⟩ by spinning: at most k·|gens| reductions, not |G|."""
    basis = []

    def reduce_add(x):
        w = x.copy()
        for b in basis:
            nz = np.flatnonzero(b)
            if nz.size and w[int(nz[0])]:
                w ^= b
        if np.any(w):
            basis.append(w)
            return True
        return False

    reduce_add(v)
    changed = True
    while changed:
        changed = False
        current = list(basis)
        for b in current:
            for M in mats:
                y = (M.astype(np.int32) @ b.astype(np.int32)) % 2
                if reduce_add(y.astype(np.uint8)):
                    changed = True
    return len(basis)


def qr_indicator(p: int) -> np.ndarray:
    q, mul, add, chi, frob, norm, ia, ib = field_ctx(p)
    n = q + 1
    v = np.zeros(n, dtype=np.uint8)
    for e in range(1, q):
        if chi(e) == 1:
            v[1 + e] = 1
    return v


def analyse(p: int) -> dict:
    print(f"\n======== p={p} ========", flush=True)
    S = square_line_matrix(p)
    H0, _ = gf2_nullspace(S)
    n, k = H0.shape[0], H0.shape[1]
    ones = np.ones(n, dtype=np.uint8)
    gens = psl_generators(p)
    mats = action_matrices(H0, gens)
    rng = np.random.default_rng(p)
    dims = []
    for _ in range(12):
        t = rng.integers(0, 2, size=k, dtype=np.uint8)
        if t.sum() == 0:
            continue
        d = cyclic_span(mats, t)
        dims.append(d)
    print(f"  dim H0={k}  PΓL cyclic dims={sorted(set(dims))}  all_full={all(d==k for d in dims)}", flush=True)

    # 1_QR vs H0 and vs Max- slices
    qr = qr_indicator(p)
    # is QR in rowspan(S)? 
    rS = gf2_rref(S)[2]
    rQR = gf2_rref(np.vstack([S, qr]))[2]
    ell = np.zeros(n, dtype=np.uint8)
    ell[0] = 1
    ell[1] = 1
    rEll = gf2_rref(np.vstack([S, ell]))[2]
    rBoth = gf2_rref(np.vstack([S, ell, qr]))[2]
    print(f"  QR in rowspan S? {rQR==rS}  ell extra? {rEll==rS+1}  S+ell+QR rank={rBoth}", flush=True)

    rec = {
        "p": p,
        "dim_H0": int(k),
        "ptgl_dims": dims,
        "ptgl_all_full": bool(dims) and all(d == k for d in dims),
        "qr_in_rowspan": rQR == rS,
        "ell_independent": rEll == rS + 1,
        "S_ell_QR_rank": int(rBoth),
    }
    if p in (5, 7):
        Y, C = load_minus(p)
        Y = np.sign(Y.astype(np.float64)).astype(np.int8)
        B = ((1 - Y) // 2).astype(np.uint8)
        fe = np.rint(C[0, 1]).astype(np.int64) * Y[:, 0] * Y[:, 1]
        U = fe < 0
        ip = (B.astype(np.int32) @ qr.astype(np.int32)) % 2
        rec["QR_on_U"] = sorted(set(ip[U].tolist()))
        rec["QR_on_Uc"] = sorted(set(ip[~U].tolist()))
        rec["QR_walsh_character"] = len(rec["QR_on_U"]) == 1
        rec["QR_const_on_Uc"] = len(rec["QR_on_Uc"]) == 1
        rec["QR_same_both"] = rec["QR_on_U"] == rec["QR_on_Uc"]
        print(
            f"  QR on U={rec['QR_on_U']} Uc={rec['QR_on_Uc']} "
            f"constU={rec['QR_walsh_character']} constUc={rec['QR_const_on_Uc']}",
            flush=True,
        )
    return rec


def main():
    out = {}
    for p in (5, 7):
        out[str(p)] = analyse(p)
    dest = ROOT / "evidence" / "walsh_ptgl_meataxe.json"
    dest.write_text(__import__("json").dumps(out, indent=2, default=str) + "\n")
    print(f"\nwrote {dest}", flush=True)


if __name__ == "__main__":
    main()
