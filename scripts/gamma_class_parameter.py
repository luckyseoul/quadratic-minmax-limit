#!/usr/bin/env python3
"""Γ(K) vs PGL class parameter (eigenvalue ratio / τ = tr²/det).

p=5 full Aut(C), then p=7 out-of-sample.  Surviving at both primes:
elliptic Γ=0, involution Γ=2(n-2), unipotent Γ determined by λ_exc,
split Γ a function of τ.  The p=5 O(1)-in-p split formula FAILS at p=7.
See evidence/PLAN_2026-08-22_class_function_route.md Step 3.

Backend: p=5 CPU; p≥7 Γ on the V100 via CuPy (torch in this env is CPU-only).
BFS of Aut is sequential and short.
"""
from __future__ import annotations

import sys
from collections import Counter, defaultdict, deque
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15590 import (  # noqa: E402
    MuLab,
    field_ops,
    paley_conference,
    signed_generators,
)


def frob_fn(p, fmul):
    def frob(e):
        x = e
        for _ in range(p - 1):
            x = fmul(x, e)
        return x

    return frob


def field_inv(p, q, fmul, one):
    finv = [0] * q
    for e in range(1, q):
        finv[e] = next(x for x in range(1, q) if fmul(e, x) == one)
    return finv


def enum_aut(p, C):
    n = C.shape[0]
    gens = signed_generators(p, C)
    auto = [(pi, d) for pi, d, s in gens if s == 1]
    idp = (tuple(range(n)), tuple([1] * n))
    seen = {idp}
    dq = deque([idp])
    elems = [idp]
    while dq:
        pi0, d0 = dq.popleft()
        pa = np.array(pi0)
        da = np.array(d0)
        for pi, d in auto:
            e = (tuple(pi[pa]), tuple(da * d[pa]))
            if e not in seen:
                seen.add(e)
                dq.append(e)
                elems.append(e)
    return elems


def gamma_cpu(Y, elems, n):
    gam = np.empty(len(elems), dtype=np.float64)
    Y32 = Y.astype(np.int32, copy=False)
    for i, (pi, d) in enumerate(elems):
        pa = np.asarray(pi, dtype=np.int64)
        da = np.asarray(d, dtype=np.int32)
        s = (Y32 * (Y32[:, pa] * da[None, :])).sum(axis=1).astype(np.float64)
        gam[i] = (s * s).mean() - 2 * n
    return gam


def gamma_gpu(Y, elems, n, batch=512):
    """Batched Γ on the V100 via CuPy (this env's torch is CPU-only)."""
    import cupy as cp

    Yt = cp.asarray(Y, dtype=cp.float32)
    gam = np.empty(len(elems), dtype=np.float64)
    for lo in range(0, len(elems), batch):
        chunk = elems[lo : lo + batch]
        B = len(chunk)
        pi = cp.asarray([e[0] for e in chunk], dtype=cp.int32)
        d = cp.asarray([e[1] for e in chunk], dtype=cp.float32)
        Yperm = Yt[:, pi] * d[None, :, :]  # (M, B, n)
        s = (Yt[:, None, :] * Yperm).sum(axis=-1)  # (M, B)
        gam[lo : lo + B] = cp.asnumpy((s * s).mean(axis=0)) - 2 * n
        del Yperm, s, pi, d
        cp.get_default_memory_pool().free_all_blocks()
    return gam


def chi_W_of(pi, d, C, p):
    """tr(U|_{V_+}) = (tr U + tr(UC)/p)/2, U_{j, π(j)} = d_j."""
    pa = np.asarray(pi, dtype=np.int64)
    da = np.asarray(d, dtype=np.float64)
    n = len(pa)
    fixed = pa == np.arange(n)
    trU = float((da * fixed).sum())
    trUC = float((da * C[pa, np.arange(n)]).sum())
    return (trU + trUC / p) / 2.0


def classify_perm(pi, p, q, n, fmul, fadd, fneg, one, finv, frob, squares):
    """Recover (possibly Frob-precomposed) Möbius from π; return class data."""
    pa = np.asarray(pi, dtype=np.int64)

    def h_of(idx):
        if idx == 0:
            return (one, 0)  # [1:0] = ∞
        return (idx - 1, one)  # [z:1]

    def apply_A(A, z_idx, use_frob):
        a, b, c, d = A
        if z_idx == 0:
            X, Y = one, 0
        else:
            z = z_idx - 1
            if use_frob:
                z = frob(z)
            X, Y = z, one
        NX = fadd(fmul(a, X), fmul(b, Y))
        NY = fadd(fmul(c, X), fmul(d, Y))
        if NY == 0:
            return 0
        return 1 + fmul(NX, finv[NY])

    def mobius_from_images(img_inf, img_0, img_1):
        u = h_of(img_inf)
        v = h_of(img_0)
        w = h_of(img_1)
        det = fadd(fmul(u[0], v[1]), fneg(fmul(u[1], v[0])))
        if det == 0:
            return None
        # solve α u + β v = w  (2×2)
        # α = det(w,v)/det(u,v), β = det(u,w)/det(u,v)
        inv = finv[det]
        dwv = fadd(fmul(w[0], v[1]), fneg(fmul(w[1], v[0])))
        duw = fadd(fmul(u[0], w[1]), fneg(fmul(u[1], w[0])))
        al = fmul(dwv, inv)
        be = fmul(duw, inv)
        a = fmul(al, u[0])
        c = fmul(al, u[1])
        b = fmul(be, v[0])
        d = fmul(be, v[1])
        return (a, b, c, d)

    A = mobius_from_images(int(pa[0]), int(pa[1]), int(pa[1 + one]))
    if A is None:
        return dict(ok=False, reason="singular3")

    def agrees(use_frob):
        for i in range(n):
            if apply_A(A, i, use_frob) != int(pa[i]):
                return False
        return True

    use_frob = False
    if not agrees(False):
        if agrees(True):
            use_frob = True
        else:
            return dict(ok=False, reason="not_pgl_or_frob", A=A)

    a, b, c, d = A
    tr = fadd(a, d)
    det = fadd(fmul(a, d), fneg(fmul(b, c)))
    if det == 0:
        return dict(ok=False, reason="det0", use_frob=use_frob)
    tau = fmul(fmul(tr, tr), finv[det])  # tr²/det ∈ F_q
    four = 0
    for _ in range(4):
        four = fadd(four, one)
    disc = fadd(fmul(tr, tr), fneg(fmul(four, det)))

    if a == d and b == 0 and c == 0:
        family = "identity"
    elif disc == 0:
        family = "parabolic"
    elif disc in squares:
        family = "split"
    else:
        family = "elliptic"

    # P¹-fixed points of the permutation (not of the Möbius after peeling Frob)
    nfix = int((pa == np.arange(n)).sum())
    return dict(
        ok=True,
        use_frob=use_frob,
        family=family,
        tau=int(tau),
        disc=int(disc),
        tr=int(tr),
        det=int(det),
        nfix=nfix,
        A=(int(a), int(b), int(c), int(d)),
    )


def r4(x):
    return round(float(x), 6)


def predict_psl(p, n, family, tau):
    """p=5-fitted PSL formula, gated at p=7.

    Elliptic 0, identity n(n-2), involution 2(n-2) HOLD at p=7.
    The two split-τ recipes FAIL at p=7 (killed; do not reopen).
    Parabolic is determined by λ_exc, not by this predictor.
    """
    if family == "identity":
        return float(n * (n - 2))
    if family == "elliptic":
        return 0.0
    if family == "parabolic":
        return None
    if family == "split":
        b = tau % p  # encoding e = p*a+b
        if tau == 0:
            return float(2 * (n - 2))
        if b == 0:
            return float(-4 * (n - 2) / n)
        return float(-4 * (n - 2) * p / n)
    return None


def main():
    p = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    q = p * p
    n = q + 1
    print(f"p={p} q={q} n={n}", flush=True)
    fmul, fadd, fneg, one = field_ops(p)
    finv = field_inv(p, q, fmul, one)
    frob = frob_fn(p, fmul)
    squares = {fmul(x, x) for x in range(1, q)}
    C = paley_conference(p)
    lab = MuLab(p, with_deg6=False)
    Y = lab.Yp.astype(np.int8)
    print(f"|Max+|={len(Y)}", flush=True)
    elems = enum_aut(p, C)
    print(f"|G+|={len(elems)}", flush=True)

    if p >= 7:
        print("Γ on GPU (batched gather)", flush=True)
        gam = gamma_gpu(Y, elems, n)
    else:
        gam = gamma_cpu(Y, elems, n)
    print(f"distinct Γ: {len(set(map(r4, gam)))}  Γ(e)={gam[0]:.6f}  n(n-2)={n*(n-2)}", flush=True)

    lbar = 8 * (n - 2) / (n - 6)
    l_exc_census = {5: 176 / 13, 7: 4320 / 409}.get(p)

    rows = []
    nfail = 0
    for i, (pi, d) in enumerate(elems):
        cl = classify_perm(pi, p, q, n, fmul, fadd, fneg, one, finv, frob, squares)
        if not cl.get("ok"):
            nfail += 1
            if nfail <= 5:
                print("classify fail", cl, flush=True)
            continue
        chi = chi_W_of(pi, d, C, p)
        rows.append((cl, float(gam[i]), chi, d))
    print(f"classified {len(rows)}/{len(elems)}  fails={nfail}", flush=True)

    by_key = defaultdict(list)
    for cl, g, chi, d in rows:
        key = (cl["use_frob"], cl["family"], cl["tau"], cl["nfix"])
        by_key[key].append((g, chi))

    if p <= 5:
        print("\n=== Γ by (frob, family, τ, nfix) ===", flush=True)
        print(
            f"{'frob':>4} {'family':>10} {'tau':>6} {'nfix':>4} {'#g':>7} {'nΓ':>3} "
            f"{'Γ uniq':>28} {'χ_W uniq':>24}",
            flush=True,
        )
        for key in sorted(by_key, key=lambda k: (k[0], k[1], k[2], k[3])):
            vs = by_key[key]
            gs = [v[0] for v in vs]
            cs = [v[1] for v in vs]
            gset = sorted(set(map(r4, gs)))
            cset = sorted(set(map(r4, cs)))
            frob, fam, tau, nfix = key
            print(
                f"{int(frob):4d} {fam:>10} {tau:6d} {nfix:4d} {len(vs):7d} {len(gset):3d} "
                f"{str(gset)[:28]:>28} {str(cset)[:24]:>24}",
                flush=True,
            )

    print("\n=== elliptic no-Frob: per-element Γ vs χ_W ===", flush=True)
    ell = [(g, chi, cl) for cl, g, chi, d in rows if (not cl["use_frob"]) and cl["family"] == "elliptic"]
    print(f"  #elliptic no-Frob: {len(ell)}")
    if ell:
        gset = sorted(set(map(r4, [e[0] for e in ell])))
        cset = sorted(set(map(r4, [e[1] for e in ell])))
        print(f"  unique Γ: {gset}")
        print(f"  unique χ_W: {cset}")
        nz = [(g, chi) for g, chi, cl in ell if abs(chi) > 1e-8]
        print(f"  with χ_W≠0: {len(nz)}")
        if nz:
            rats = [g / chi for g, chi in nz]
            print(f"  Γ/χ_W unique rounded: {sorted(set(round(r, 4) for r in rats))}")
            if l_exc_census is not None:
                ok = all(abs(r - l_exc_census) < 1e-3 for r in rats)
                print(f"  census λ_exc={l_exc_census:.6f}  all match: {ok}")
        pairs = Counter((r4(g), r4(chi)) for g, chi, cl in ell)
        print(f"  (Γ, χ_W) pairs: {pairs.most_common(8)}")

    print("\n=== no-Frob split: Γ as function of τ, vs χ_W ===", flush=True)
    spl = defaultdict(list)
    for cl, g, chi, d in rows:
        if (not cl["use_frob"]) and cl["family"] == "split":
            spl[cl["tau"]].append((g, chi))
    nshow = 0
    for tau in sorted(spl):
        vs = spl[tau]
        gset = sorted(set(map(r4, [v[0] for v in vs])))
        cset = sorted(set(map(r4, [v[1] for v in vs])))
        pred = predict_psl(p, n, "split", tau)
        ok = pred is not None and all(abs(x - pred) < 1e-3 for x in gset)
        if p == 5 or nshow < 20 or not ok:
            print(
                f"  τ={tau:6d} inF_p={tau % p == 0} n={len(vs):5d} "
                f"Γ={gset}  pred={pred}  match={ok}  χ_W={cset}"
            )
            nshow += 1
    n_split_tau = len(spl)
    n_split_ok = sum(
        1
        for tau, vs in spl.items()
        if all(abs(v[0] - predict_psl(p, n, "split", tau)) < 1e-3 for v in vs)
    )
    print(f"  split-τ buckets matching prediction: {n_split_ok}/{n_split_tau}")

    print("\n=== PSL (no-Frob) closed-form gate ===", flush=True)
    bad = Counter()
    nchk = 0
    for cl, g, chi, d in rows:
        if cl["use_frob"]:
            continue
        pred = predict_psl(p, n, cl["family"], cl["tau"])
        if pred is None:
            continue
        nchk += 1
        if abs(g - pred) > 1e-3:
            bad[cl["family"]] += 1
    n_psl = sum(1 for cl, g, chi, d in rows if not cl["use_frob"])
    n_parab = sum(1 for cl, g, chi, d in rows if (not cl["use_frob"]) and cl["family"] == "parabolic")
    print(f"  no-Frob elements: {n_psl}  (signed PSL, expect 2|PSL|={q*(q*q-1)})")
    print(f"  checked (id+elliptic+split): {nchk}  mismatches: {dict(bad)}")
    print(f"  GATE {'PASS' if not bad else 'FAIL'}")
    parabs = [g for cl, g, chi, d in rows if (not cl["use_frob"]) and cl["family"] == "parabolic"]
    if parabs:
        print(f"  parabolic n={n_parab} unique Γ: {sorted(set(map(r4, parabs)))}")

    print("\n=== no-Frob parabolic / identity ===", flush=True)
    for fam in ("identity", "parabolic"):
        vs = [(g, chi, cl) for cl, g, chi, d in rows if (not cl["use_frob"]) and cl["family"] == fam]
        gset = sorted(set(map(r4, [v[0] for v in vs])))
        cset = sorted(set(map(r4, [v[1] for v in vs])))
        print(f"  {fam}: n={len(vs)} Γ={gset}  χ_W={cset}")
        if fam == "identity":
            for cl, g, chi, dsgn in rows:
                if (not cl["use_frob"]) and cl["family"] == "identity":
                    pa = None
                    print(f"    Γ={g:.4f} χ_W={chi:.4f} nfix={cl['nfix']} #neg={(np.asarray(dsgn)==-1).sum()}")

    # Does τ determine Γ inside a family (no Frob)?
    print("\n=== does τ determine Γ inside (frob=0, family)? ===", flush=True)
    by_ft = defaultdict(lambda: defaultdict(set))
    for cl, g, chi, d in rows:
        by_ft[(cl["use_frob"], cl["family"])][cl["tau"]].add(r4(g))
    for ft in sorted(by_ft, key=lambda k: (k[0], k[1])):
        sizes = [len(s) for s in by_ft[ft].values()]
        n_tau = len(by_ft[ft])
        print(
            f"  frob={int(ft[0])} {ft[1]:>10}: {n_tau} values of τ, "
            f"Γs-per-τ min={min(sizes)} max={max(sizes)}"
        )

    if p <= 5:
        print("\n=== χ_W ≈ 0 : remaining Γ (principal-series content) ===", flush=True)
        zw = [(cl, g, chi) for cl, g, chi, d in rows if abs(chi) < 1e-6]
        byz = defaultdict(list)
        for cl, g, chi in zw:
            byz[(cl["use_frob"], cl["family"], cl["tau"])].append(g)
        print(f"  #{{g : χ_W=0}} = {len(zw)}")
        for k in sorted(byz)[:30]:
            gs = byz[k]
            print(
                f"  frob={int(k[0])} {k[1]:>10} τ={k[2]:6d}  n={len(gs):6d}  "
                f"Γ={sorted(set(map(r4, gs)))[:4]}"
            )

    print(f"\nlbar={lbar:.6f}", flush=True)


if __name__ == "__main__":
    main()
