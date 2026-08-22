#!/usr/bin/env python3
"""Involution E[s²] index split, and Kloosterman/Bessel vs dilation Γ.

No flag flip.  Not an identity file unless involution is actually
Max+-free from 2-point (it may not be).
"""
from __future__ import annotations

import json
import math
import sys
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from e1_gmin_m4_prop15590 import MuLab, field_ops, paley_conference  # noqa: E402


def frac(x, cap=20_000_000):
    return Fraction(float(x)).limit_denominator(cap)


def primitive_root(q, fmul, one):
    def order_of(e):
        x, o = e, 1
        while x != one:
            x = fmul(x, e)
            o += 1
            if o > q:
                return 0
        return o

    return next(e for e in range(2, q) if order_of(e) == q - 1)


def dlog_table(gen, q, fmul, one):
    tab = [-1] * q
    x = one
    for k in range(q - 1):
        tab[x] = k
        x = fmul(x, gen)
    return tab


def inv_table(q, fmul, one):
    finv = np.zeros(q, dtype=np.int32)
    for e in range(1, q):
        r, b, ee = one, e, q - 2
        while ee:
            if ee & 1:
                r = fmul(r, b)
            b = fmul(b, b)
            ee >>= 1
        finv[e] = r
    return finv


def signed_lift(pi, C):
    n = C.shape[0]
    for s in (1, -1):
        d = np.zeros(n, dtype=np.int64)
        d[0] = 1
        d[1:] = s * C[pi[0], pi[1:]] * C[0, 1:]
        ok = s * d[:, None] * d[None, :] * C == C[np.ix_(pi, pi)]
        np.fill_diagonal(ok, True)
        if ok.all():
            return d, s
    return None, 0


def dilation_pi(t, q, n, fmul):
    pi = np.zeros(n, dtype=np.int64)
    for e in range(q):
        pi[1 + e] = 1 + fmul(e, t)
    return pi


def involution_split(p: int) -> dict:
    """Split E[s²] by |{i,πi,j,πj}|.  2-point prediction from E[yyᵀ]=I+C/p."""
    q = p * p
    n = q + 1
    fmul, fadd, fneg, one = field_ops(p)
    C = paley_conference(p)
    neg1 = fneg(one)
    pi = dilation_pi(neg1, q, n, fmul)
    d, s = signed_lift(pi, C)
    assert s == 1
    Y = MuLab(p, with_deg6=False).Yp.astype(np.int32)
    M = len(Y)
    # per-index pair values u_i = d_i y_i y_{πi}
    u = (d[None, :] * Y * Y[:, pi]).astype(np.float64)  # (M, n)
    svec = u.sum(axis=1)
    Es2 = float((svec * svec).mean())
    # pair (i,j) contribution classified by orbit size of {i,πi,j,πj} under π
    # Use one Max+ row? No, ensemble.  Vectorize over i<=j types via numpy.
    # Collision type from index sets:
    #  fix = {k: πk=k} = {0, 1}  (∞ and field 0)
    #  pairs {x,-x} for x in F_q^×/{±1}
    print(f"\n=== involution index split p={p} ===", flush=True)
    print(f"  E[s²]={Es2:.6f}  4(n-1)={4*(n-1)}  Γ={Es2-2*n:.6f}  2(n-2)={2*(n-2)}")
    print(f"  unique s (up to 12): {sorted(set(svec.astype(int)))[:12]}...")

    # 2-point only: E[u_i u_j] for i≠j using E[y_a y_b]=C_ab/p, Boolean y^2=1
    # u_i = d_i y_i y_πi
    # i=j: E[u_i²]=1
    # i≠j: E[y_i y_πi y_j y_πj] d_i d_j
    #   if 4 distinct: unknown m4
    #   if 3 distinct: 3-point
    #   if 2 distinct: either {i,πi}={j,πj} so j=πi, or {i,j} with π fixing both, etc.
    Cf = C.astype(np.float64)
    two_pt = 0.0  # sum_{i,j} d_i d_j * (Wick/boolean-2 from Σ)
    # Boolean 2-point contraction treating four labels with possible repeats:
    # use actual ensemble partition empirically
    n_eq = n  # diag
    # Classify all (i,j) by cardinality of {i, πi, j, πj}
    card_mass = defaultdict(float)  # empirical E[u_i u_j] summed
    card_count = defaultdict(int)
    card_wick = defaultdict(float)  # 2-point prediction
    Pi = pi
    D = d
    # Wick from pairwise: for a 4-tuple (possibly with repeats) of coordinates
    # Boolean: E[prod y] = 1 if all paired, C/p if one pair of distinct, 0 if odd, m4 if 4 distinct
    def two_point_moment(a, b, c, e):
        """E[y_a y_b y_c y_e] using only y^2=1 and E[y_i y_j]=C_ij/p.
        Returns None if 4 distinct (unknown)."""
        idx = [a, b, c, e]
        # fold squares
        from collections import Counter
        cnt = Counter(idx)
        # odd multiplicity leftover
        leftover = []
        for v, k in cnt.items():
            if k % 2:
                leftover.append(v)
            # k>=2 contributes 1 per full pair
        if len(leftover) == 0:
            return 1.0
        if len(leftover) == 2:
            x, y_ = leftover
            if x == y_:
                return 1.0
            return float(Cf[x, y_] / p)
        if len(leftover) == 1:
            return 0.0  # odd
        return None  # 4 distinct or 3? 3 leftover is odd -> 0? 3 distinct odd total
        # 3 leftover means 3 odd multiplicities, total degree 4 so one of them
        # was counted wrong. degree 4, leftover odd count even. 0,2,4 leftovers.

    for i in range(n):
        for j in range(n):
            S = {int(i), int(Pi[i]), int(j), int(Pi[j])}
            c = len(S)
            val = float(np.mean(u[:, i] * u[:, j]))
            card_mass[c] += val
            card_count[c] += 1
            w = two_point_moment(i, int(Pi[i]), j, int(Pi[j]))
            if w is None:
                card_wick[c] += float("nan")  # mark
            else:
                card_wick[c] += D[i] * D[j] * w

    print(f"  {'card':>4} {'#(i,j)':>8} {'∑ E[u_i u_j]':>16} {'2-pt pred':>16}")
    tot = 0.0
    tot_w = 0.0
    unk = 0.0
    for c in sorted(card_mass):
        pred = card_wick[c]
        print(
            f"  {c:4d} {card_count[c]:8d} {card_mass[c]:16.6f} "
            f"{pred if pred==pred else float('nan'):16.6f}"
        )
        tot += card_mass[c]
        if pred == pred:
            tot_w += pred
        else:
            unk += card_mass[c]
    print(f"  sum empirical={tot:.6f}  2-pt-known types={tot_w:.6f}  4-distinct mass={unk:.6f}")
    print(f"  4-distinct needed for floor remainder: {tot - tot_w:.6f}")
    return {
        "p": p,
        "Es2": Es2,
        "target": 4 * (n - 1),
        "Gamma": Es2 - 2 * n,
        "card_mass": {str(k): v for k, v in card_mass.items()},
        "two_point_known": tot_w,
        "four_distinct": tot - tot_w,
        "pointwise_constant": len(set(svec.astype(int))) == 1,
    }


def kloosterman_vs_gamma(p: int) -> dict:
    q = p * p
    n = q + 1
    fmul, fadd, fneg, one = field_ops(p)
    C = paley_conference(p)
    finv = inv_table(q, fmul, one)
    gen = primitive_root(q, fmul, one)
    dlog = dlog_table(gen, q, fmul, one)
    sq = {fmul(x, x) for x in range(1, q)}
    chi = np.zeros(q, dtype=np.int8)
    for e in range(1, q):
        chi[e] = 1 if e in sq else -1

    def tr(e):
        a, b = divmod(int(e), p)
        return (2 * a) % p

    w = np.exp(2j * np.pi * np.arange(p) / p)
    Y = MuLab(p, with_deg6=False).Yp.astype(np.int32)
    print(f"\n=== Kloosterman / Bessel vs dilation Γ  p={p} ===", flush=True)

    rows = []
    for t in range(1, q):
        if chi[t] != 1:
            continue
        pi = dilation_pi(t, q, n, fmul)
        d, s = signed_lift(pi, C)
        if s != 1:
            continue
        sval = (Y * (Y[:, pi] * d[None, :])).sum(axis=1).astype(np.float64)
        gam = float((sval * sval).mean() - 2 * n)
        # Kl(1,t)=∑_{x≠0} e(Tr(x + t x^{-1}))
        kl = 0j
        bes = 0j  # ∑ χ(x) e(Tr(t x + x^{-1}))
        bes2 = 0j  # ∑ χ(x(x-1)) e(Tr(t x))
        for x in range(1, q):
            ix = int(finv[x])
            tx = fmul(t, ix)
            kl += w[tr(fadd(x, tx))]
            bes += chi[x] * w[tr(fadd(fmul(t, x), ix))]
            xm1 = fadd(x, fneg(1 if False else one))
            if xm1 != 0:
                bes2 += chi[int(fmul(x, xm1))] * w[tr(fmul(t, x))]
        ord_t = (q - 1) // math.gcd(int(dlog[t]), q - 1)
        rows.append(
            {
                "t": t,
                "dlog": int(dlog[t]),
                "order": ord_t,
                "Gamma": gam,
                "Kl_re": float(kl.real),
                "Kl_im": float(kl.imag),
                "Kl_abs": float(abs(kl)),
                "Bes_re": float(bes.real),
                "Bes_abs": float(abs(bes)),
                "Bes2_re": float(bes2.real),
                "Bes2_abs": float(abs(bes2)),
            }
        )

    # unique Γ vs unique |Kl|
    print(
        f"  {'t':>4} {'ord':>4} {'Γ':>14} {'Re Kl':>10} {'|Kl|':>8} "
        f"{'Re Bes':>10} {'|Bes|':>8} {'|Bes2|':>8}"
    )
    shown = set()
    for r in sorted(rows, key=lambda x: x["dlog"]):
        key = min(r["dlog"], (q - 1 - r["dlog"]) % (q - 1))
        if key in shown:
            continue
        shown.add(key)
        print(
            f"  {r['t']:4d} {r['order']:4d} {r['Gamma']:14.6f} {r['Kl_re']:10.3f} "
            f"{r['Kl_abs']:8.3f} {r['Bes_re']:10.3f} {r['Bes_abs']:8.3f} {r['Bes2_abs']:8.3f}"
        )

    # Does |Kl| (or Re Kl, |Bes|) determine Γ?
    def split_count(field):
        b = defaultdict(set)
        for r in rows:
            b[round(r["Gamma"], 8)].add(round(r[field], 4))
        n_g = len(b)
        n_split = sum(1 for vs in b.values() if len(vs) > 1)
        # reverse: field -> Γ
        br = defaultdict(set)
        for r in rows:
            br[round(r[field], 4)].add(round(r["Gamma"], 8))
        n_fsplit = sum(1 for vs in br.values() if len(vs) > 1)
        return n_g, n_split, len(br), n_fsplit

    print("  invariant → Γ  (nΓ, Γ-split-on-field, nField, field-split-on-Γ)")
    for field in ("Kl_abs", "Kl_re", "Bes_abs", "Bes_re", "Bes2_abs", "Bes2_re"):
        print(f"    {field:10s} {split_count(field)}")

    # linear: Γ = a + b |Kl|²/q + c ReKl /√q
    G = np.array([r["Gamma"] for r in rows])
    qf = float(q)
    feats = {
        "|Kl|^2/q": np.array([r["Kl_abs"] ** 2 / qf for r in rows]),
        "ReKl/p": np.array([r["Kl_re"] / p for r in rows]),
        "|Bes|^2/q": np.array([r["Bes_abs"] ** 2 / qf for r in rows]),
        "ReBes/p": np.array([r["Bes_re"] / p for r in rows]),
        "|Bes2|^2/q": np.array([r["Bes2_abs"] ** 2 / qf for r in rows]),
    }
    print("  lin Γ ~ 1 + feature:")
    for name, col in feats.items():
        X = np.column_stack([np.ones(len(G)), col])
        coef, *_ = np.linalg.lstsq(X, G, rcond=None)
        err = float(np.max(np.abs(X @ coef - G)))
        print(f"    {name:16s} maxerr={err:.4f} coef={np.round(coef, 4)}")
    # 1 + two features
    X = np.column_stack(
        [np.ones(len(G)), feats["|Kl|^2/q"], feats["ReKl/p"], feats["|Bes|^2/q"]]
    )
    coef, *_ = np.linalg.lstsq(X, G, rcond=None)
    err = float(np.max(np.abs(X @ coef - G)))
    print(f"    1+|Kl|²/q+ReKl/p+|Bes|²/q  maxerr={err:.4f} coef={np.round(coef, 4)}")
    return {"p": p, "n_rows": len(rows), "rows": rows}


def main():
    out = {}
    for p in (5, 7):
        out[f"inv_{p}"] = involution_split(p)
        out[f"kl_{p}"] = kloosterman_vs_gamma(p)
    path = ROOT / "evidence" / "involution_kloosterman.json"
    # drop bulky rows from kl for json? keep Gamma and Kl summaries
    slim = {}
    for k, v in out.items():
        if k.startswith("kl_"):
            slim[k] = {
                "p": v["p"],
                "n_rows": v["n_rows"],
                "by_t": [
                    {
                        "t": r["t"],
                        "order": r["order"],
                        "Gamma": r["Gamma"],
                        "Kl_abs": r["Kl_abs"],
                        "Bes_abs": r["Bes_abs"],
                    }
                    for r in v["rows"]
                ],
            }
        else:
            slim[k] = v
    path.write_text(json.dumps(slim, indent=2, default=float))
    print(f"\nwrote {path}", flush=True)


if __name__ == "__main__":
    main()
