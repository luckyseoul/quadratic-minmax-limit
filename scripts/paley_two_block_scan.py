#!/usr/bin/env python3
"""Verified two-block decompositions of Paley conference witnesses.

- Complete switching+permutation canonical form (self-tested).
- Full partition scans for C6 (m=3), C10 (m=5), C14 (m=7).
- Brute-force switching extraction + explicit witness for hits; exact Phi checks.
- Squares/nonsquares character-split test for q=5,9,13,17,25.
- Phi(C18), Phi(C26) context values.
No all-orders claim; exact finite computations only.
"""
import itertools
import json
import sys
import time

import numpy as np

sys.path.insert(0, "/home/nick/quadratic-minmax-limit/scripts")
from paley_structural_lift_family import paley_conference, gf_elements  # noqa: E402

T0 = time.time()


def log(*a):
    print(f"[{time.time()-T0:6.1f}s]", *a, flush=True)


def phi_signing(A):
    n = A.shape[0]
    m = 1 << (n - 1)
    u = np.arange(m, dtype=np.int64)
    X = np.ones((m, n), dtype=np.int64)
    for c in range(1, n):
        X[:, c] = 1 - 2 * ((u >> (c - 1)) & 1)
    Q = np.einsum('mi,ij,mj->m', X, A, X) * 0.5
    return int(np.abs(Q).max())


def phi_big_pairs(A):
    """Phi for n up to 26 via pairwise accumulation (memory-light)."""
    n = A.shape[0]
    m = 1 << (n - 1)
    chunk = 1 << 21
    best = 0
    pairs = [(i, j, int(A[i, j])) for i in range(n) for j in range(i + 1, n)]
    for lo in range(0, m, chunk):
        hi = min(lo + chunk, m)
        u = np.arange(lo, hi, dtype=np.int64)
        X = np.ones((hi - lo, n), dtype=np.int8)
        for c in range(1, n):
            X[:, c] = 1 - 2 * ((u >> (c - 1)) & 1).astype(np.int8)
        acc = np.zeros(hi - lo, dtype=np.int32)
        for i, j, a in pairs:
            if a:
                acc += a * (X[:, i] * X[:, j]).astype(np.int32)
        v = int(np.abs(acc).max())
        if v > best:
            best = v
    return best


PERMS = {k: list(itertools.permutations(range(k - 1))) for k in range(2, 8)}


def canon_complete(A):
    """Complete invariant of switching+permutation for zero-diagonal signings (k<=7)."""
    A = np.asarray(A)
    k = A.shape[0]
    assert k <= 7, "canon_complete capped at k=7"
    best = None
    perms = PERMS[k]
    for v in range(k):
        d = A[v].copy()
        d[v] = 1
        N = (d[:, None] * d[None, :]) * A
        others = [x for x in range(k) if x != v]
        for pm in perms:
            order = [v] + [others[t] for t in pm]
            b = N[np.ix_(order, order)].astype(np.int8).tobytes()
            if best is None or b < best:
                best = b
    return best


def scan_partitions(C, m):
    n2 = C.shape[0]
    coords = list(range(n2))
    hits = []
    for T in itertools.combinations(coords, m):
        if T[0] != 0:
            continue
        Tc = tuple(c for c in coords if c not in T)
        A = C[np.ix_(T, T)]
        Cc = C[np.ix_(Tc, Tc)]
        if canon_complete(A) == canon_complete(-Cc):
            hits.append(T)
    return hits


def realize(C, T, m):
    """Brute switching (d,e) with dAd = -eCe; return witness info or None."""
    Tc = tuple(c for c in range(C.shape[0]) if c not in T)
    A = C[np.ix_(T, T)]
    Cc = C[np.ix_(Tc, Tc)]
    chosen = None
    for dmask in range(1 << m):
        d = np.array([1 if (dmask >> i) & 1 else -1 for i in range(m)], dtype=np.int64)
        Ad = (d[:, None] * d[None, :]) * A
        for emask in range(1 << m):
            e = np.array([1 if (emask >> i) & 1 else -1 for i in range(m)], dtype=np.int64)
            if (Ad == -((e[:, None] * e[None, :]) * Cc)).all():
                chosen = (d.copy(), e.copy())
                break
        if chosen:
            break
    if chosen is None:
        return None
    d, e = chosen
    g = np.ones(C.shape[0], dtype=np.int64)
    g[list(T)] = d
    g[list(Tc)] = e
    S2 = (g[:, None] * g[None, :]) * C
    assert (S2[np.ix_(T, T)] == -S2[np.ix_(Tc, Tc)]).all()
    return dict(T=[int(x) for x in T], d=d.tolist(), e=e.tolist(), S2=S2)


def split_sets(q):
    els, add, mul = gf_elements(q)
    zero = els[0] if not isinstance(els[0], tuple) else tuple([0] * len(els[0]))
    squares = {mul(e, e) for e in els if e != zero}
    idx = {e: i + 1 for i, e in enumerate(els)}
    sq = sorted(idx[e] for e in squares if e != zero)
    ns = sorted(i for i in range(1, q + 1) if i not in sq and i != idx[zero])
    return sq, ns


def main():
    out = {}

    # canon self-test
    rng = np.random.default_rng(1)
    ok = True
    for _ in range(20):
        A = (rng.random((7, 7)) < 0.5).astype(np.int64) * 2 - 1
        A = np.triu(A, 1)
        A = A + A.T
        P = rng.permutation(7)
        d = rng.choice([-1, 1], size=7)
        B = ((A * d[:, None] * d[None, :])[np.ix_(P, P)]).astype(np.int64)
        if canon_complete(A) != canon_complete(B):
            ok = False
    log("canon_complete self-test:", "PASS" if ok else "FAIL")
    out["canon_selftest"] = ok
    assert ok

    report = {}
    for q, m in ((5, 3), (9, 5), (13, 7)):
        C = paley_conference(q)
        phiC = phi_signing(C)
        hits = scan_partitions(C, m)
        block_phis = sorted({phi_signing(C[np.ix_(T, T)]) for T in hits})
        report[q] = dict(phi=phiC, n_hits=len(hits), total_partitions=len(list(itertools.combinations(range(q + 1), m))) // 2,
                         block_phis=block_phis)
        log(f"q={q}: Phi(C{q+1})={phiC}, hits={len(hits)}, block Phis={block_phis}")

    # explicit witness for C14, prefer a block-Phi 9 hit
    C14 = paley_conference(13)
    hits14 = scan_partitions(C14, 7)
    wit = None
    for T in hits14:
        w = realize(C14, T, 7)
        if w is None:
            continue
        A2 = w["S2"][np.ix_(w["T"], w["T"])]
        phA = phi_signing(A2)
        if phA == 9 and wit is None:
            wit = dict(T=w["T"], d=w["d"], e=w["e"],
                       phi_S=phi_signing(w["S2"]), phi_block=phA,
                       S2=w["S2"].tolist())
    log("C14 witness:", None if wit is None else (wit["T"], wit["phi_S"], wit["phi_block"]))
    out["C14_witness"] = wit

    # compare hit partitions with character splits
    sq13, ns13 = split_sets(13)
    char_part = sorted([0] + sq13)
    char_part2 = sorted([0] + ns13)
    hit_sets = {tuple(sorted(T)) for T in hits14}
    log("char split inf+sq hit:", tuple(char_part) in hit_sets, "| inf+ns hit:", tuple(char_part2) in hit_sets)
    log("C14 hit partitions:", [list(T) for T in hits14][:10])
    out["C14_hit_sets"] = [list(T) for T in hits14]
    out["C14_char_splits"] = dict(inf_sq=[0] + sq13, inf_ns=[0] + ns13,
                                  inf_sq_hit=tuple(char_part) in hit_sets,
                                  inf_ns_hit=tuple(char_part2) in hit_sets)

    # pattern test for q=17, q=25 (character splits; direct identity check)
    for q in (17, 25):
        C = paley_conference(q)
        m = (q + 1) // 2
        sq, ns = split_sets(q)
        res = {}
        for name, T in (("inf_sq", [0] + sq), ("inf_ns", [0] + ns)):
            Tc = [i for i in range(q + 1) if i not in set(T)]
            # order Tc as [zero_index] + [phi(x) for x in T-order] via nonsquare multiplier
            els, add, mul = gf_elements(q)
            zero = els[0] if not isinstance(els[0], tuple) else tuple([0] * len(els[0]))
            idx = {e: i + 1 for i, e in enumerate(els)}
            zidx = idx[zero]
            nsq_el = None
            squares = {mul(e, e) for e in els if e != zero}
            for e in els:
                if e != zero and e not in squares:
                    nsq_el = e
                    break
            if name == "inf_sq":
                pair = [zidx] + [idx[mul(nsq_el, els[i - 1])] for i in T if i != 0]
            else:
                # inverse multiplier still a nonsquare
                pair = [zidx] + [idx[mul(nsq_el, els[i - 1])] for i in T if i != 0]
            A = C[np.ix_(T, T)]
            Cc = C[np.ix_(pair, pair)]
            # find switching (vectorized per d, derive e from first row)
            k = len(T)
            found = False
            for dmask in range(1 << min(k, 14)):
                d = np.array([1 if (dmask >> i) & 1 else -1 for i in range(k)], dtype=np.int64)
                Ad = (d[:, None] * d[None, :]) * A
                e = np.ones(k, dtype=np.int64)
                # e_0 = 1; e_j from row 0: Ad_0j = -e_0 e_j Cc_0j
                e[1:] = -Ad[0, 1:] * Cc[0, 1:]
                Ee = (e[:, None] * e[None, :]) * Cc
                if (Ad == -Ee).all():
                    found = True
                    break
            res[name] = dict(found=found, k=k)
            log(f"q={q} split {name}: switch-reali zation found = {found} (k={k})")
        # Phi values
        if q == 17:
            phiC = phi_signing(C)
        else:
            phiC = phi_big_pairs(C)
        res["phi_C"] = phiC
        res["m"] = m
        out[f"q{q}"] = res
        log(f"q={q}: Phi(C{q+1}) = {phiC}")

    json.dump(out, open("/home/nick/scratch/paley_two_block_results.json", "w"), indent=2)
    log("-> /home/nick/scratch/paley_two_block_results.json")


if __name__ == "__main__":
    main()
