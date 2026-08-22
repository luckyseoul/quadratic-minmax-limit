#!/usr/bin/env python3
"""GF(2) dual / affine geometry of Max- and of a pair-slice U.

Goal: a general-p reason for 15.406 Theorem C (Walsh: characters constant
on U are constant on U^c), including the p=11 case rank(B_U)=n/2-1.

No flag flip. Not an identity file. Serial GE is inherently sequential;
p=3,5,7 matrices are tiny (the inner numpy is vectorized).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15406 import (  # noqa: E402
    gf2_nullspace,
    gf2_rref,
    gf2_solve,
    load_minus,
)
from e1_gmin_m4_prop15590 import paley_conference  # noqa: E402


def bits_of(Y: np.ndarray) -> np.ndarray:
    return ((1 - Y.astype(np.int8)) // 2).astype(np.uint8)


def affine_rank(B: np.ndarray) -> int:
    """dim of affine span = rank of {rows} after homogenizing with a 1-column."""
    if B.size == 0:
        return 0
    A = np.concatenate([B, np.ones((B.shape[0], 1), dtype=np.uint8)], axis=1)
    return gf2_rref(A)[2]


def dual_right(B: np.ndarray) -> np.ndarray:
    """Columns a with B a = 0 over F2 (linear dual of the row-span)."""
    N, n = gf2_nullspace(B)
    if N.size == 0:
        return np.zeros((B.shape[1], 0), dtype=np.uint8)
    return N


def affine_chars(B: np.ndarray) -> dict:
    """Linear dual + whether 1 is in the column space (constant-1 solvable)."""
    N = dual_right(B)
    x1 = gf2_solve(B, np.ones(B.shape[0], dtype=np.uint8))
    return {
        "lin_rank": int(B.shape[1] - N.shape[1]),
        "ker_dim": int(N.shape[1]),
        "aff_rank": int(affine_rank(B)),
        "one_solvable": x1 is not None,
        "dual": N,
        "x1": x1,
    }


def support_wt(v: np.ndarray) -> list[int]:
    return [int(i) for i in np.flatnonzero(v)]


def paley_adj_mod2(C: np.ndarray) -> np.ndarray:
    """A = (J - I - C)/2, 0-1 Paley (plus infinity star). Then mod 2."""
    n = C.shape[0]
    J = np.ones((n, n), dtype=np.int64)
    np.fill_diagonal(J, 0)
    Ci = np.rint(C).astype(np.int64)
    A = (J - Ci) // 2
    return np.mod(A, 2).astype(np.uint8)


def describe_dual(N: np.ndarray, C: np.ndarray, tag: str) -> dict:
    n = C.shape[1] if C.ndim == 2 else C.shape[0]
    rows = []
    ones = np.ones(n, dtype=np.uint8)
    for c in range(N.shape[1]):
        v = N[:, c]
        wt = int(v.sum())
        supp = support_wt(v)
        is_ones = bool(np.array_equal(v, ones))
        # match against Paley-adjacency rows and C-rows mod 2
        A = paley_adj_mod2(C)
        match_A = [int(i) for i in range(n) if np.array_equal(v, A[i])]
        C2 = np.mod(np.rint(C).astype(np.int64), 2).astype(np.uint8)
        # C ≡ J-I, so C2[i] = 1 - e_i
        match_C2 = [int(i) for i in range(n) if np.array_equal(v, C2[i])]
        rows.append(
            {
                "wt": wt,
                "supp_head": supp[:24],
                "is_allones": is_ones,
                "match_A_rows": match_A[:8],
                "match_C2_rows": match_C2[:8],
            }
        )
    wts = sorted({r["wt"] for r in rows})
    n_ones = sum(r["is_allones"] for r in rows)
    print(
        f"  {tag}: ker_dim={N.shape[1]}  weights={wts}  #allones={n_ones}",
        flush=True,
    )
    for r in rows[:8]:
        print(
            f"    wt={r['wt']:3d} ones={r['is_allones']} "
            f"A={r['match_A_rows']} C2={r['match_C2_rows']} "
            f"supp={r['supp_head'][:12]}",
            flush=True,
        )
    return {"n_dual": int(N.shape[1]), "weights": wts, "rows_sample": rows[:12]}


def xor_cut(Y: np.ndarray, C: np.ndarray, i: int = 0, j: int = 1):
    """U = {C_ij y_i y_j < 0}. Over F2 this is x_i XOR x_j = 0 or 1."""
    fe = C[i, j] * Y[:, i] * Y[:, j]
    U = fe < 0
    B = bits_of(Y)
    xor = (B[:, i] ^ B[:, j]).astype(np.uint8)
    # y_i y_j = -1  <=>  x_i XOR x_j = 1, independently of C
    yprod_neg = (Y[:, i] * Y[:, j] < 0)
    return U, ~U, B, xor, yprod_neg


def difference_span(B: np.ndarray) -> int:
    """Rank of {row - row0} = affine direction dimension."""
    if len(B) == 0:
        return 0
    D = (B ^ B[0]) & 1
    return gf2_rref(D)[2]


def is_linear_code(B: np.ndarray) -> dict:
    """Is the set of rows a coset of a linear code? Sample-closed under XOR."""
    nR, n = B.shape
    # translate so a point is 0
    D = (B ^ B[0]) & 1
    # 0 in D?
    has0 = bool((D.sum(axis=1) == 0).any())
    # random triple closure: d1^d2 in D
    rng = np.random.default_rng(0)
    n_tests = min(200, nR * (nR - 1) // 2)
    hits = 0
    # hash rows as bytes for membership
    packed = np.packbits(D, axis=1)
    keys = {row.tobytes() for row in packed}
    for _ in range(n_tests):
        a, b = rng.integers(0, nR, size=2)
        s = (D[a] ^ D[b]).astype(np.uint8)
        if np.packbits(s).tobytes() in keys:
            hits += 1
    return {
        "translated_has0": has0,
        "xor_closure_hits": int(hits),
        "xor_closure_tests": int(n_tests),
        "xor_closed_sample": hits == n_tests,
        "dir_rank": int(gf2_rref(D)[2]),
    }


def analyse(p: int) -> dict:
    print(f"\n======== p={p} ========", flush=True)
    Y, C = load_minus(p)
    Y = np.sign(Y.astype(np.float64)).astype(np.int8)
    n = C.shape[0]
    B = bits_of(Y)
    print(f"  |Max-|={len(Y)} n={n} n/2={n // 2}", flush=True)

    full = affine_chars(B)
    print(
        f"  Max- lin_rank={full['lin_rank']} aff_rank={full['aff_rank']} "
        f"dir={difference_span(B)} ker={full['ker_dim']} "
        f"1-solvable={full['one_solvable']}",
        flush=True,
    )
    even = int((B.sum(axis=1) % 2).min()), int((B.sum(axis=1) % 2).max())
    print(f"  weight parity of x: min/max {even}", flush=True)

    dual_info = describe_dual(full["dual"], C, "Max- dual")
    lin = is_linear_code(B)
    print(f"  linear-coset sample: {lin}", flush=True)

    U, Uc, B, xor, yneg = xor_cut(Y, C, 0, 1)
    print(
        f"  |U|={int(U.sum())} |Uc|={int(Uc.sum())} "
        f"P(U)={U.mean():.4f}  (p+1)/(2p)={(p + 1) / (2 * p):.4f}",
        flush=True,
    )
    print(
        f"  y_i y_j=-1 iff xor=1? {(yneg == (xor == 1)).all()}  "
        f"U vs xor: U_xor_values={sorted(set(xor[U].tolist()))} "
        f"Uc_xor={sorted(set(xor[Uc].tolist()))}",
        flush=True,
    )

    BU, BC = B[U], B[Uc]
    uinfo, cinfo = affine_chars(BU), affine_chars(BC)
    print(
        f"  U  lin={uinfo['lin_rank']} aff={uinfo['aff_rank']} "
        f"dir={difference_span(BU)} ker={uinfo['ker_dim']} "
        f"1solv={uinfo['one_solvable']}",
        flush=True,
    )
    print(
        f"  Uc lin={cinfo['lin_rank']} aff={cinfo['aff_rank']} "
        f"dir={difference_span(BC)} ker={cinfo['ker_dim']} "
        f"1solv={cinfo['one_solvable']}",
        flush=True,
    )

    # containment: dual of U ⊆ dual of (U union something)? Walsh dual view:
    # every a with BU a = const should have BC a = const.
    NU = uinfo["dual"]
    ker_mixed = 0
    if NU.shape[1]:
        KN = (BC.astype(np.int32) @ NU.astype(np.int32)) % 2
        ker_mixed = int(sum(KN[:, c].min() != KN[:, c].max() for c in range(KN.shape[1])))
    aff_mixed = 0
    x1 = uinfo["x1"]
    if x1 is not None:
        w1 = (BC.astype(np.int32) @ x1.astype(np.int32)) % 2
        if w1.size and w1.min() != w1.max():
            aff_mixed = 1
    print(f"  Walsh ker_mixed={ker_mixed} aff_mixed={aff_mixed}", flush=True)

    # is e_i + e_j in the dual of U? (linear: xor constant 0 on U)
    eij = np.zeros(n, dtype=np.uint8)
    eij[0] = 1
    eij[1] = 1
    # check BU @ eij constant
    u_eij = (BU.astype(np.int32) @ eij) % 2
    c_eij = (BC.astype(np.int32) @ eij) % 2
    f_eij = (B.astype(np.int32) @ eij) % 2
    print(
        f"  e0+e1 on Max- values={sorted(set(f_eij.tolist()))} "
        f"on U={sorted(set(u_eij.tolist()))} on Uc={sorted(set(c_eij.tolist()))}",
        flush=True,
    )

    # Paley A rank over F2
    A = paley_adj_mod2(C)
    rA = gf2_rref(A)[2]
    print(f"  Paley-A rank_F2={rA}  (n={n})", flush=True)

    # Does the Max- dual contain the all-ones?
    ones = np.ones(n, dtype=np.uint8)
    has_ones = bool(
        any(np.array_equal(full["dual"][:, c], ones) for c in range(full["dual"].shape[1]))
    )
    print(f"  dual contains 1? {has_ones}", flush=True)

    # Extra dual of U not in dual of Max-: these cut U vs the span
    # Represent duals as sets of packed columns
    def pack_cols(N):
        if N.size == 0:
            return set()
        return {N[:, c].tobytes() for c in range(N.shape[1])}

    extra_U = pack_cols(NU) - pack_cols(full["dual"])
    print(f"  extra dual vectors of U beyond Max-: {len(extra_U)}", flush=True)
    extra_rows = []
    for col in range(NU.shape[1]):
        v = NU[:, col]
        if v.tobytes() in extra_U:
            extra_rows.append({"wt": int(v.sum()), "supp": support_wt(v)[:32]})
    for r in extra_rows[:6]:
        print(f"    extra U-dual wt={r['wt']} supp={r['supp']}", flush=True)

    # Hypothesis: extra dual is exactly e_i+e_j, possibly plus a Max- dual vector
    extra_is_eij = any(
        np.array_equal(NU[:, c], eij) for c in range(NU.shape[1])
    )
    print(f"  e0+e1 in dual(U)? {extra_is_eij}", flush=True)

    return {
        "p": p,
        "n": n,
        "nMax": int(len(Y)),
        "full": {k: v for k, v in full.items() if k not in ("dual", "x1")},
        "U": {k: v for k, v in uinfo.items() if k not in ("dual", "x1")},
        "Uc": {k: v for k, v in cinfo.items() if k not in ("dual", "x1")},
        "walsh": {"ker_mixed": ker_mixed, "aff_mixed": aff_mixed},
        "dual_info": dual_info,
        "linear_coset": lin,
        "paley_A_rank": int(rA),
        "dual_has_ones": has_ones,
        "n_extra_U_dual": len(extra_U),
        "extra_U_sample": extra_rows[:8],
        "eij_in_dual_U": extra_is_eij,
        "weight_parity": list(even),
        "P_U": float(U.mean()),
    }


def main():
    out = {"primes": {}}
    for p in (3, 5, 7):
        out["primes"][str(p)] = analyse(p)
    dest = ROOT / "evidence" / "walsh_gf2_dual.json"
    dest.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"\nwrote {dest}", flush=True)
    return out


if __name__ == "__main__":
    main()
