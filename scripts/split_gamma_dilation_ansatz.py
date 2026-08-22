#!/usr/bin/env python3
"""Split Γ from square dilations; q-dependent ansatz; p=11 λ_k gate.

Phase 0–3 of the split-Γ route.  No flag flip.  Not an identity file.

For each t∈F_q^× build π: z↦t z on P¹ (∞↦∞, 0↦0), signed lift of Aut(C)
if it exists (s=1 ⇔ t a square).  Γ(t)=E[(y·πy)²]−2n.  Mellin against
ρ_k(t)=2 cos(2π k dlog t /(q−1)) recovers λ_k.  Fit a p-independent
term-count formula with q-dependent coefficients at p=5,7, then predict
the HIP p=11 even-character spectrum.
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


def frac(x, cap=10_000_000):
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


def frob_pow(e, p, fmul, one):
    if e == 0:
        return 0
    r, b, ee = one, e, p
    while ee:
        if ee & 1:
            r = fmul(r, b)
        b = fmul(b, b)
        ee >>= 1
    return r


def signed_lift_pi(pi, C):
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


def dilation_perm(t, q, n, fmul):
    pi = np.zeros(n, dtype=np.int64)
    pi[0] = 0  # ∞
    for e in range(q):
        pi[1 + e] = 1 + fmul(e, t)
    return pi


def analyze_p(p: int) -> dict:
    q = p * p
    n = q + 1
    fmul, fadd, fneg, one = field_ops(p)
    C = paley_conference(p)
    lab = MuLab(p, with_deg6=False)
    Y = lab.Yp.astype(np.int32)
    M = len(Y)
    gen = primitive_root(q, fmul, one)
    dlog = dlog_table(gen, q, fmul, one)
    finv = inv_table(q, fmul, one)
    neg1 = fneg(one)
    sq = {fmul(x, x) for x in range(1, q)}
    chi = np.zeros(q, dtype=np.int8)
    for e in range(1, q):
        chi[e] = 1 if e in sq else -1

    print(f"\n=== p={p} q={q} n={n} |Max+|={M} gen={gen} ===", flush=True)

    # identity / involution / all t
    rows = []
    Yf = Y.astype(np.float64)
    for t in range(1, q):
        pi = dilation_perm(t, q, n, fmul)
        d, s = signed_lift_pi(pi, C)
        if d is None:
            continue
        svec = (Y * (Y[:, pi] * d[None, :])).sum(axis=1).astype(np.float64)
        gam = float((svec * svec).mean() - 2 * n)
        # unsigned multiplicative pairing on F_q^×
        # M = ∑_{e≠0} z_e z_{t e}
        z = Y[:, 1:]
        te = np.array([fmul(e, t) for e in range(q)], dtype=np.int64)
        Mmul = (z * z[:, te]).sum(axis=1).astype(np.float64)  # includes e=0
        # e=0 term is z_0^2=1; want ∑_{e≠0}
        Mstar = Mmul - 1.0
        EM2 = float((Mstar * Mstar).mean())
        # signed overlap s vs Mstar + y_∞ d_∞ y_∞ + z_0 d_0 z_0 = Mstar + d[0] + d[1]
        ord_t = (q - 1) // math.gcd(int(dlog[t]), q - 1)
        tr = fadd(t, frob_pow(t, p, fmul, one))
        nm = fmul(t, frob_pow(t, p, fmul, one))
        ta, tb = divmod(int(tr), p)
        na, nb = divmod(int(nm), p)
        in_fp = (t % p == 0)  # 15590: e=p*a+b, F_p is b=0, e=p*a; one=p so t=one has b=0
        # actually F_p: b=0, a arbitrary. t=p*a, t%p==0. one=p, p%p=0. Good.
        tm1 = fadd(t, fneg(one))
        tp1 = fadd(t, one)
        rows.append(
            {
                "t": int(t),
                "dlog": int(dlog[t]),
                "s_lift": int(s),
                "square": bool(chi[t] == 1),
                "order": int(ord_t),
                "in_Fp": bool(in_fp),
                "chi": int(chi[t]),
                "chi_tm1": 0 if tm1 == 0 else int(chi[tm1]),
                "chi_tp1": 0 if tp1 == 0 else int(chi[tp1]),
                "norm_b": int(nb),
                "tr_b": int(tb),
                "Gamma": gam,
                "Gamma_frac": str(frac(gam, M * n * n)),
                "EM2": EM2,
                "mean_s": float(svec.mean()),
                "mean_Mstar": float(Mstar.mean()),
                "d0": int(d[0]),
                "d_field0": int(d[1]),
                "ndeg": int((np.asarray(d) == -1).sum()),
            }
        )

    auto = [r for r in rows if r["s_lift"] == 1]
    print(f"  lifts: {len(rows)}  Aut s=1: {len(auto)} (expect q-1 squares? {(q-1)//2} plus 1)")
    # identity
    idn = next(r for r in rows if r["t"] == int(one))
    inv = next(r for r in rows if r["t"] == int(neg1))
    print(
        f"  identity t=1  Γ={idn['Gamma_frac']} expect n(n-2)={n*(n-2)}  s={idn['s_lift']}"
    )
    print(
        f"  involution t=-1 Γ={inv['Gamma_frac']} expect 2(n-2)={2*(n-2)}  "
        f"square={inv['square']} s={inv['s_lift']}"
    )

    print(
        f"  {'t':>4} {'dlog':>5} {'ord':>4} {'sq':>2} {'Fp':>2} {'χ':>3} "
        f"{'χ-1':>4} {'χ+1':>4} {'s':>2} {'Γ':>18} {'E M_*^2':>12}"
    )
    # unique Aut Γ by dlog class t ~ t^{-1}
    seen_pair = set()
    split_vals = []
    for r in sorted(auto, key=lambda x: x["dlog"]):
        a, b = r["dlog"], (q - 1 - r["dlog"]) % (q - 1)
        key = (min(a, b),)
        mark = ""
        if key in seen_pair:
            mark = " (inv)"
        else:
            seen_pair.add(key)
            if r["t"] not in (int(one), int(neg1)):
                split_vals.append(r)
        if r["s_lift"] == 1:
            print(
                f"  {r['t']:4d} {r['dlog']:5d} {r['order']:4d} {int(r['square']):2d} "
                f"{int(r['in_Fp']):2d} {r['chi']:3d} {r['chi_tm1']:4d} {r['chi_tp1']:4d} "
                f"{r['s_lift']:2d} {r['Gamma_frac']:>18} {r['EM2']:12.4f}{mark}"
            )

    # Mellin of Aut Γ → λ_k.  ρ_k(t)=2 cos(2π k dlog t /(q-1))
    # Average over Aut dilations (the split torus + {±1}).
    # Proper PSL inner product is class-weighted; for a first check use
    # uniform on {t square} which is the split torus plus identity+involution.
    half = (q - 1) // 2
    ks = list(range(2, half, 2))  # even field k in (0, half); QVAR is (q-1)/4
    # Use square t only (s=1). Class {t,t^{-1}} once.
    reps = []
    used = set()
    for r in auto:
        a = r["dlog"]
        b = (q - 1 - a) % (q - 1)
        key = min(a, b)
        if key in used:
            continue
        used.add(key)
        reps.append(r)
    print(f"  distinct {{t, t^{-1}}} Aut classes: {len(reps)}")

    def mellin(k, use="Gamma"):
        acc = 0.0
        wsum = 0.0
        for r in reps:
            ang = 2 * math.pi * k * r["dlog"] / (q - 1)
            rho = 2.0 * math.cos(ang)
            val = r["Gamma"] if use == "Gamma" else r["EM2"]
            # class size: t=t^{-1} (t^2=1) have smaller classes
            w = 1.0 if r["dlog"] % (q - 1) in (0, half) else 2.0
            # actually we already took one of {t,t^{-1}}, so weight 1 for all
            acc += val * rho
            wsum += 1.0
        return acc / len(reps)

    known = {
        5: {6: 176 / 13, 2: 80 / 13, 4: 144 / 13, 8: 144 / 13, 10: 80 / 13},
        7: {
            12: 4320 / 409,
            2: 3360 / 409,
            4: 4032 / 409,
            6: 3648 / 409,
            8: 3072 / 409,
            10: 3360 / 409,
            14: 3360 / 409,
            16: 3072 / 409,
            18: 3648 / 409,
            20: 4032 / 409,
            22: 3360 / 409,
        },
    }.get(p, {})
    print(f"  Mellin of dilation-Γ vs known even-k λ (NOT the PSL L2 inner product):")
    print(f"  {'k':>4} {'4|k':>4} {'qvar':>4} {'mellin Γ':>12} {'known λ':>12} {'mellin EM2':>12}")
    for k in ks:
        mg = mellin(k, "Gamma")
        mm = mellin(k, "EM2")
        kn = known.get(k)
        print(
            f"  {k:4d} {int(k % 4 == 0):4d} {int(k == (q-1)//4):4d} "
            f"{mg:12.6f} {'' if kn is None else f'{kn:12.6f}'} {mm:12.6f}"
        )

    unique_g = sorted({round(r["Gamma"], 8) for r in auto if r["t"] not in (int(one), int(neg1))})
    print(f"  unique split Aut Γ values: {len(unique_g)}  {unique_g}")

    return {
        "p": p,
        "M": M,
        "identity": idn,
        "involution": inv,
        "involution_ok": abs(inv["Gamma"] - 2 * (n - 2)) < 1e-6,
        "identity_ok": abs(idn["Gamma"] - n * (n - 2)) < 1e-6,
        "n_auto": len(auto),
        "n_split_Gamma": len(unique_g),
        "rows": rows,
        "reps": reps,
    }


def fit_and_gate(recs):
    """Bounded-term q-dependent models of λ(k)-8; gate p=11."""
    npz = np.load(ROOT / "evidence" / "even_char_hip_minmax.npz")
    recs[11] = {
        "p": 11,
        "ks": [int(x) for x in npz["ks"]],
        "lams": [float(x) for x in npz["lams"]],
    }
    # also p=5,7 λ from known / dilation
    known_lams = {
        5: [(2, 80 / 13), (4, 144 / 13), (6, 176 / 13), (8, 144 / 13), (10, 80 / 13)],
        7: [
            (2, 3360 / 409),
            (4, 4032 / 409),
            (6, 3648 / 409),
            (8, 3072 / 409),
            (10, 3360 / 409),
            (12, 4320 / 409),
            (14, 3360 / 409),
            (16, 3072 / 409),
            (18, 3648 / 409),
            (20, 4032 / 409),
            (22, 3360 / 409),
        ],
    }
    for p in (5, 7):
        recs[p]["ks"] = [k for k, _ in known_lams[p]]
        recs[p]["lams"] = [lam for _, lam in known_lams[p]]

    print("\n=== ansatz: λ-8 = a0(q) + a1(q) cos(4πk/(q-1)) + a2(q) cos(8πk/(q-1)) ===")
    print("    a_i linear in (1, q, 1/q, n=q+1) fitted at p=5,7; gate p=11")

    def feats_q(p):
        q = p * p
        n = q + 1
        return np.array([1.0, q, 1.0 / q, n, 1.0 / n, p, 1.0 / p])

    def design_k(p, ks, n_harm):
        q = p * p
        theta = 2 * np.pi * np.asarray(ks, dtype=np.float64) / (q - 1)
        cols = [np.ones(len(ks))]
        for m in range(1, n_harm + 1):
            cols.append(np.cos(2 * m * theta))  # cos(4π m k/(q-1)) for m=1 is 2θ
        return np.column_stack(cols)

    # Model: (λ-8) = sum_{j=0}^{n_harm}  (c_j · feats_q(p))  cos(2π j * 2k/(q-1))
    # Too many params if full feats.  Restrict: each a_j is A + B/q  (2 params).
    for n_harm in (1, 2, 3):
        # params: (n_harm+1) * 2  (const and 1/q for each harmonic including 0)
        Xs, ys = [], []
        for p in (5, 7):
            q = p * p
            ks = recs[p]["ks"]
            y = np.asarray(recs[p]["lams"]) - 8
            Fq = np.array([1.0, 1.0 / q])
            th = 2 * np.pi * np.asarray(ks, dtype=np.float64) / (q - 1)
            blocks = []
            for m in range(n_harm + 1):
                cm = np.ones(len(ks)) if m == 0 else np.cos(2 * m * th)
                blocks.append(cm[:, None] * Fq[None, :])
            X = np.hstack(blocks)
            Xs.append(X)
            ys.append(y)
        X = np.vstack(Xs)
        y = np.concatenate(ys)
        coef, *_ = np.linalg.lstsq(X, y, rcond=None)
        print(f"\n  n_harm={n_harm}  n_param={len(coef)}  coef={np.round(coef, 5)}")
        for p in (5, 7, 11):
            q = p * p
            ks = recs[p]["ks"]
            ytrue = np.asarray(recs[p]["lams"]) - 8
            Fq = np.array([1.0, 1.0 / q])
            th = 2 * np.pi * np.asarray(ks, dtype=np.float64) / (q - 1)
            blocks = []
            for m in range(n_harm + 1):
                cm = np.ones(len(ks)) if m == 0 else np.cos(2 * m * th)
                blocks.append(cm[:, None] * Fq[None, :])
            pred = np.hstack(blocks) @ coef
            err = float(np.max(np.abs(pred - ytrue)))
            print(
                f"    p={p:2d} max|λ-pred|={err:.4f}  "
                f"min_true={ytrue.min()+8:.4f} min_pred={pred.min()+8:.4f}"
            )

    # Jacobi-style: λ-8 = a(q) + b(q) Re[ i^{k} ] wait even.
    # χ-twist: for k,  α_k(-1)=1 (even). Try λ-8 = a + b cos(π k / 2)  useless.
    # Gauss: |G(α_k)| is p for nontrivial. Constant, not useful.
    # J(χ, α_k) for α_k^2 ≠ χ.
    print("\n=== Jacobi-length ansatz: λ-8 = a(q)+b(q) Re J_model(k) ===")
    # J_model(k)= cos(π k / ((q-1)/4)) = cos(4π k /(q-1))  already in cosine.
    # Try a+b χ_4(k) where χ_4(k)= 0 if 4|k else ... skip.

    # Paley-type of t formula for Γ, Mellin to λ — tested via unique Γ vs χ(t),χ(t-1)
    print("\n=== is Aut Γ a function of (χ(t),χ(t-1),χ(t+1),in_Fp,order)? ===")
    for p in (5, 7):
        buckets = defaultdict(set)
        for r in recs[p]["rows"]:
            if r["s_lift"] != 1:
                continue
            if r["t"] in (p,):  # one = p in 15590; skip identity separately
                pass
            key = (
                r["chi"],
                r["chi_tm1"],
                r["chi_tp1"],
                int(r["in_Fp"]),
                r["order"],
            )
            buckets[key].add(round(r["Gamma"], 8))
        nsplit = sum(1 for vs in buckets.values() if len(vs) > 1)
        print(f"  p={p} keys={len(buckets)} split={nsplit}")
        if nsplit:
            shown = 0
            for k, vs in buckets.items():
                if len(vs) > 1:
                    print(f"    split {k} -> {sorted(vs)[:6]}")
                    shown += 1
                    if shown >= 5:
                        break


def sample_p11_dilation_gamma(B=150_000, seed=0):
    """Count unique split dilation Γ at p=11 on a Max+ sample (minmax labels)."""
    p, q, n = 11, 121, 122

    def is_irr(a, b):
        return all((x * x - a * x - b) % p != 0 for x in range(p))

    ia = ib = None
    for a in range(p):
        for b in range(p):
            if is_irr(a, b):
                ia, ib = a, b
                break
        if ia is not None:
            break

    def fmul(u, v):
        c0, c1 = u % p, u // p
        d0, d1 = v % p, v // p
        return ((c0 * d0 + c1 * d1 * ib) % p) + (
            (c0 * d1 + c1 * d0 + c1 * d1 * ia) % p
        ) * p

    def fadd(u, v):
        return ((u % p + v % p) % p) + ((u // p + v // p) % p) * p

    def fneg(u):
        return ((-u) % p) + ((-(u // p)) % p) * p

    one = 1
    C = paley_conference_minmax(p, fmul, fadd, fneg, one)
    path = Path("/mnt/storage/e1work/maxplus_p11/maxplus_p11_eps1.npy")
    Ymm = np.load(path, mmap_mode="r")
    rng = np.random.default_rng(seed)
    idx = np.sort(rng.choice(Ymm.shape[0], size=B, replace=False))
    Y = np.asarray(Ymm[idx], dtype=np.int32)
    print(f"\n=== p=11 sample B={B} irr=({ia},{ib}) ===", flush=True)
    sq = {fmul(x, x) for x in range(1, q)}
    gen = primitive_root(q, fmul, one)
    dlog = dlog_table(gen, q, fmul, one)
    neg1 = fneg(one)
    gammas = []
    for t in range(1, q):
        if t not in sq:
            continue
        pi = dilation_perm(t, q, n, fmul)
        d, s = signed_lift_pi(pi, C)
        if s != 1:
            continue
        svec = (Y * (Y[:, pi] * d[None, :])).sum(axis=1).astype(np.float64)
        gam = float((svec * svec).mean() - 2 * n)
        ord_t = (q - 1) // math.gcd(int(dlog[t]), q - 1)
        tm1 = fadd(t, fneg(one))
        tp1 = fadd(t, one)
        chi = lambda e: 0 if e == 0 else (1 if e in sq else -1)
        gammas.append(
            {
                "t": t,
                "dlog": int(dlog[t]),
                "order": ord_t,
                "chi_tm1": chi(tm1),
                "chi_tp1": chi(tp1),
                "Gamma": gam,
                "is_id": t == one,
                "is_inv": t == neg1,
            }
        )
    split = [g for g in gammas if not g["is_id"] and not g["is_inv"]]
    # cluster unique Γ
    vals = np.array([g["Gamma"] for g in split])
    # round to 0.05 to beat sampling noise
    keys = np.round(vals, 1)
    uniq = sorted(set(keys.tolist()))
    print(f"  Aut square dilations: {len(gammas)}  split: {len(split)}")
    print(f"  unique split Γ (round 0.1): {len(uniq)}  expect (p-1)/2={ (p-1)//2 }")
    print(f"  involution Γ={next(g['Gamma'] for g in gammas if g['is_inv']):.4f} expect {2*(n-2)}")
    print(f"  identity Γ={next(g['Gamma'] for g in gammas if g['is_id']):.4f} expect {n*(n-2)}")
    by = defaultdict(list)
    for g, k in zip(split, keys):
        by[k].append((g["order"], g["chi_tm1"], g["chi_tp1"]))
    for k in uniq:
        ords = sorted({o for o, _, _ in by[k]})
        pal = sorted({(a, b) for _, a, b in by[k]})
        print(f"    Γ≈{k:8.1f}  n={len(by[k]):3d}  orders={ords}  Paley(t±1)={pal}")
    return {"n_unique_split": len(uniq), "involution": next(g["Gamma"] for g in gammas if g["is_inv"])}


def paley_conference_minmax(p, fmul, fadd, fneg, one):
    """C matching paley_conference_prime_power / eps1.npy (e=a+bp)."""
    q = p * p
    n = q + 1
    sq = {fmul(x, x) for x in range(1, q)}
    C = np.zeros((n, n), dtype=np.int64)
    C[0, 1:] = 1
    C[1:, 0] = 1
    for e1 in range(q):
        for e2 in range(q):
            if e1 == e2:
                continue
            d = fadd(e1, fneg(e2))
            C[1 + e1, 1 + e2] = 1 if d in sq else -1
    return C


def main():
    recs = {}
    for p in (5, 7):
        recs[p] = analyze_p(p)
    fit_and_gate(recs)
    p11 = sample_p11_dilation_gamma()
    out = {}
    for p in (5, 7):
        rec = recs[p]
        out[str(p)] = {
            "identity_ok": rec["identity_ok"],
            "involution_ok": rec["involution_ok"],
            "identity": rec["identity"]["Gamma_frac"],
            "involution": rec["involution"]["Gamma_frac"],
            "n_auto": rec["n_auto"],
            "n_split_Gamma": rec["n_split_Gamma"],
            "reps": [
                {
                    "t": r["t"],
                    "dlog": r["dlog"],
                    "order": r["order"],
                    "in_Fp": r["in_Fp"],
                    "chi": r["chi"],
                    "chi_tm1": r["chi_tm1"],
                    "chi_tp1": r["chi_tp1"],
                    "Gamma": r["Gamma_frac"],
                    "EM2": r["EM2"],
                }
                for r in rec["reps"]
            ],
        }
    out["11_sample"] = p11
    path = ROOT / "evidence" / "split_gamma_dilation.json"
    path.write_text(json.dumps(out, indent=2))
    print(f"\nwrote {path}", flush=True)


if __name__ == "__main__":
    main()
