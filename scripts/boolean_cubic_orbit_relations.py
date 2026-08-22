#!/usr/bin/env python3
"""Boolean cubic on Ω as a constraint on Aut-orbit values of Q.

15.279 T, pointwise on every Max+:
    2p ẑ(ξ) + ∑_{t∈R} ẑ(tξ) ẑ((1-t)ξ) = 0    (ξ∈Ω)
R={t≠0,1: χ(t)=χ(1-t)=1}, |R|=(q-5)/4.

Squared and averaged this is a 4-linear of ẑ.  Q(r)=E[u(ξ)u(rξ)]
with u=|ẑ|² is a 4-linear of a special shape.  This script asks:

  (i)  does the (t,s)-Gram Γ_{t,s} of the bilinear live in the
       Aut-orbit span of Q?  (linear relations on leftover Q_O)
  (ii) does Wick of the squared cubic *name* Q(r)?  (closed form)
  (iii) rank of cubic-derived linear forms on the leftover orbits
       after Q(±1)=8q² and the row-sum.

p=5,7 MuLab, vectorized numpy.  No flag flip.  Not an identity file.
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from e1_gmin_m4_prop15590 import MuLab, field_ops  # noqa: E402


def add_mul_neg(q, fadd, fmul, fneg):
    Add = np.empty((q, q), dtype=np.int32)
    Mul = np.empty((q, q), dtype=np.int32)
    Neg = np.empty(q, dtype=np.int32)
    for i in range(q):
        Neg[i] = fneg(i)
        for j in range(q):
            Add[i, j] = fadd(i, j)
            Mul[i, j] = fmul(i, j)
    return Add, Mul, Neg


def tr_table(p, q):
    """15590: e=p*a+b ↔ a+bt, Tr=2a."""
    tr = np.empty(q, dtype=np.int32)
    for e in range(q):
        a, b = divmod(e, p)
        tr[e] = (2 * a) % p
    return tr


def frob_pow(e, p, fmul, one):
    if e == 0:
        return 0
    r, base, ee = one, e, p
    while ee:
        if ee & 1:
            r = fmul(r, base)
        base = fmul(base, base)
        ee >>= 1
    return r


def inv_of(a, q, fmul, one):
    r, base, ee = one, a, q - 2
    while ee:
        if ee & 1:
            r = fmul(r, base)
        base = fmul(base, base)
        ee >>= 1
    return r


def t_rep(r, Mul, neg1):
    return min(int(r), int(Mul[r, neg1]))


def aut_orbits_T(p, squares, Mul, fmul, one, neg1):
    unused = set(t_rep(r, Mul, neg1) for r in squares)
    orbits = []
    while unused:
        start = min(unused)
        stack = [start]
        seen = set()
        while stack:
            x = stack.pop()
            if x in seen:
                continue
            seen.add(x)
            stack.append(t_rep(frob_pow(x, p, fmul, one), Mul, neg1))
            stack.append(t_rep(inv_of(x, p * p, fmul, one), Mul, neg1))
        unused -= seen
        orbits.append(frozenset(seen))
    return orbits


def analyze(p: int) -> dict:
    q = p * p
    fmul, fadd, fneg, one = field_ops(p)
    Add, Mul, Neg = add_mul_neg(q, fadd, fmul, fneg)
    neg1 = int(Neg[one])
    tr = tr_table(p, q)
    chi = np.zeros(q, dtype=np.int8)
    sqset = {fmul(t, t) for t in range(1, q)}
    squares = sorted(sqset)
    for e in range(1, q):
        chi[e] = 1 if e in sqset else -1

    # Ω from χ̂
    Trxi = np.empty((q, q), dtype=np.int32)
    for xi in range(q):
        for x in range(q):
            Trxi[xi, x] = tr[Mul[xi, x]]
    w = np.exp(2j * np.pi * np.arange(p) / p)
    chihat = np.zeros(q, dtype=np.complex128)
    for xi in range(1, q):
        chihat[xi] = (chi[1:].astype(np.float64) * w[Trxi[xi, 1:]]).sum()
    Omega = np.where(np.abs(chihat.real - p) < 0.5)[0]
    Omega = Omega[Omega > 0]
    assert len(Omega) == (q - 1) // 2, (len(Omega), (q - 1) // 2)

    # R
    R = []
    for t in range(1, q):
        if t == one:
            continue
        tm = int(Add[one, Neg[t]])
        if tm != 0 and chi[t] == 1 and chi[tm] == 1:
            R.append(t)
    assert len(R) == (q - 5) // 4, (len(R), (q - 5) // 4)
    R = np.array(R, dtype=np.int32)
    one_minus = np.array([int(Add[one, Neg[t]]) for t in R], dtype=np.int32)

    orbits = aut_orbits_T(p, squares, Mul, fmul, one, neg1)
    oid = {}
    for i, orb in enumerate(orbits):
        for t in orb:
            oid[int(t)] = i
            oid[int(Mul[t, neg1])] = i
    n_orb = len(orbits)
    # identity orbit: t_rep(one)
    id_orb = oid[int(one)]

    lab = MuLab(p, with_deg6=False)
    Z = lab.Yp.astype(np.int8)[:, 1:]
    M = len(Z)
    print(f"\n=== p={p} |Max+|={M} |Ω|={len(Omega)} |R|={len(R)} n_orb={n_orb} ===", flush=True)

    Eom = w[Trxi[Omega]]  # (|Ω|, q)
    zhat = Z.astype(np.complex128) @ Eom.T  # (M, |Ω|)
    u = np.abs(zhat) ** 2
    om_index = {int(xi): j for j, xi in enumerate(Omega)}
    # map field el -> zhat column; 0 if off Ω (should not happen for t*Ω)
    col = np.full(q, -1, dtype=np.int32)
    for j, xi in enumerate(Omega):
        col[int(xi)] = j

    # cubic residual
    max_res = 0.0
    # vectorized over ξ-index: for each t, zhat[:, col[t*Ω]] * zhat[:, col[(1-t)*Ω]]
    Bsum = np.zeros_like(zhat)
    for t, tm in zip(R, one_minus):
        tOm = Mul[t, Omega]
        mOm = Mul[tm, Omega]
        ct = col[tOm]
        cm = col[mOm]
        assert np.all(ct >= 0) and np.all(cm >= 0)
        Bsum += zhat[np.arange(M)[:, None], ct] * zhat[np.arange(M)[:, None], cm]
    z0 = Z.astype(np.complex128).sum(axis=1)  # ẑ(0)=∑z = p y_∞
    acc_2p = 2 * p * zhat + Bsum
    acc_z0 = 2 * z0[:, None] * zhat + Bsum
    yinf = lab.Yp.astype(np.int8)[:, 0]
    print(
        f"  cubic max |2pẑ+∑B|={np.max(np.abs(acc_2p)):.3e}  "
        f"|2ẑ(0)ẑ+∑B|={np.max(np.abs(acc_z0)):.3e}  "
        f"|ẑ(0)| unique {sorted(set(np.round(np.abs(z0), 6)))}",
        flush=True,
    )
    for sign, mask in (("+1", yinf == 1), ("-1", yinf == -1)):
        if not np.any(mask):
            continue
        print(
            f"    y_∞={sign}  |2pẑ+∑B|={np.max(np.abs(acc_2p[mask])):.3e}  "
            f"|2ẑ(0)ẑ+∑B|={np.max(np.abs(acc_z0[mask])):.3e}  "
            f"||∑B|-2p|ẑ||={np.max(np.abs(np.abs(Bsum[mask]) - 2*p*np.abs(zhat[mask]))):.3e}",
            flush=True,
        )
    max_res = float(np.max(np.abs(acc_z0)))
    print(f"  E|ẑ|² on Ω = {u.mean():.4f} expect 2q={2*q}", flush=True)

    # Q(r) = mean_{y,ξ} u(ξ) u(rξ)  for r square, ξ∈Ω
    Q = {}
    for r in squares:
        rOm = Mul[r, Omega]
        cr = col[rOm]
        # if rΩ=Ω, cr all ≥0
        if np.any(cr < 0):
            Q[int(r)] = float("nan")
            continue
        ur = u[np.arange(M)[:, None], cr]
        Q[int(r)] = float(np.mean(u * ur))

    q2 = q * q
    print(f"  Q(1)={Q[int(one)]:.4f} expect 8q²={8*q2}", flush=True)
    print(f"  Q(-1)={Q[int(neg1)]:.4f}", flush=True)
    # orbit means
    print("  Aut-orbit Q:")
    Qorb = []
    for i, orb in enumerate(sorted(orbits, key=lambda o: (len(o), min(o)))):
        reps = sorted(orb)
        r0 = reps[0]
        vals = [Q[t_rep(r, Mul, neg1)] if t_rep(r, Mul, neg1) in Q else Q.get(r, float("nan")) for r in reps]
        # Q is keyed by actual squares; each t-rep is a square
        qs = [Q[r] for r in reps if r in Q]
        # also the -r copies already in squares
        mu = float(np.mean([Q[r] for r in squares if oid.get(r) == oid.get(r0)]))
        # cleaner:
        members = [r for r in squares if oid.get(int(r)) == oid.get(int(r0))]
        mu = float(np.mean([Q[int(r)] for r in members]))
        spread = float(np.ptp([Q[int(r)] for r in members]))
        Qorb.append((i, int(r0), len(orb), mu, spread, r0 == int(one) or r0 == int(neg1)))
        print(
            f"    orb {i} rep={r0:3d} |T|={len(orb)} Q={mu:.4f} "
            f"Wick={8*q2 if (int(r0)==int(one) or int(r0)==int(neg1)) else 4*q2:.0f} "
            f"spread={spread:.3e} leftover={abs(mu-4*q2):.4f}"
        )

    # Γ_{t,s} = mean_{y,ξ} ẑ(tξ)ẑ((1-t)ξ) conj(ẑ(sξ)ẑ((1-s)ξ))
    nR = len(R)
    # For each t, bilinear B_t[y,ξ] = ẑ(tξ)ẑ((1-t)ξ)
    B = np.empty((nR, M, len(Omega)), dtype=np.complex128)
    for k, (t, tm) in enumerate(zip(R, one_minus)):
        ct = col[Mul[t, Omega]]
        cm = col[Mul[tm, Omega]]
        B[k] = zhat[np.arange(M)[:, None], ct] * zhat[np.arange(M)[:, None], cm]
    # Γ[k,l] = mean B[k] conj B[l]
    Gamma = np.empty((nR, nR), dtype=np.complex128)
    for k in range(nR):
        for l in range(nR):
            Gamma[k, l] = np.mean(B[k] * np.conjugate(B[l]))
    # check 4p² E u = ∑_{k,l} Γ  (trace of all-ones sandwich)
    rhs = float(np.real(Gamma.sum()))
    lhs = 4 * p * p * float(u.mean())
    print(f"  squared cubic: 4p² E u={lhs:.4f}  ∑Γ={rhs:.4f}  diff={lhs-rhs:.3e}", flush=True)

    # Wick prediction of Γ: ẑ complex-Gaussian, only pairings
    # E ẑ(a)ẑ(b)conj ẑ(c) conj ẑ(d) with a=tξ, b=(1-t)ξ, c=sξ, d=(1-s)ξ
    # Independent data: E|ẑ|²=2q, ẑ(-η)=conj ẑ(η).
    # Two-pair Wick: pair (t with s) and (1-t with 1-s) if tξ=sξ i.e. t=s,
    # plus pair (t with 1-s) and (1-t with s) if t+(1-s)=0 i.e. s=t? wait
    # a=c ⇒ t=s; a=d ⇒ t=1-s; b=c ⇒ 1-t=s; b=d ⇒ 1-t=1-s ⇒ t=s.
    # Also conjugate pairings through -1: a=-c etc. if -1∈Ω ratios.
    Wick = np.zeros((nR, nR), dtype=np.complex128)
    twoq = 2 * q
    for k, t in enumerate(R):
        tm = one_minus[k]
        for l, s in enumerate(R):
            sm = one_minus[l]
            val = 0.0
            # |ẑ|² pairings (same frequency)
            if t == s:
                val += twoq * twoq  # E|ẑ(tξ)|² E|ẑ((1-t)ξ)|²
            if t == sm and tm == s:
                val += twoq * twoq
            if t == s and tm == sm:
                pass  # already counted if t=s ⇒ 1-t=1-s
            Wick[k, l] = val
    # This naive Wick ignores ẑ(-η)=conj and E ẑ(ξ)².
    print(
        f"  naive |ẑ|²-Wick ∑={Wick.real.sum():.4f} vs ∑Γ={rhs:.4f}",
        flush=True,
    )

    # Ratio of bilinear legs: ρ(t)=(1-t)*inv(t)  (square)
    rho = []
    for t, tm in zip(R, one_minus):
        it = inv_of(int(t), q, fmul, one)
        rho.append(int(Mul[tm, it]))
    rho = np.array(rho, dtype=np.int32)
    # For diagonal Γ_{t,t} = E |ẑ(tξ)ẑ((1-t)ξ)|² = E u(tξ) u((1-t)ξ) = Q(ρ(t))
    print("  diagonal Γ_{t,t} vs Q(ρ(t)):", flush=True)
    diag_ok = True
    rho_orb_count = defaultdict(int)
    for k, (t, r) in enumerate(zip(R, rho)):
        g = float(np.real(Gamma[k, k]))
        qr = Q[int(r)]
        rho_orb_count[oid.get(int(r), -1)] += 1
        if k < 8 or abs(g - qr) > 1.0:
            print(f"    t={int(t):3d} ρ={int(r):3d} orb={oid.get(int(r),-1)} Γtt={g:.4f} Q={qr:.4f} d={g-qr:.3e}")
        if abs(g - qr) > 1e-4 * max(1.0, abs(qr)):
            diag_ok = False
    print(f"  diag Γ=Q(ρ) all t? {diag_ok}  ρ-orb counts {dict(rho_orb_count)}", flush=True)

    # Off-diagonal: is Γ_{t,s} a function of (orb(ρ(t)), orb(ρ(s)), orb(s/t), ...)?
    # Collect unique (orb_ρt, orb_ρs, orb_{s inv t}) -> Γ values
    buckets = defaultdict(list)
    for k, t in enumerate(R):
        for l, s in enumerate(R):
            it = inv_of(int(t), q, fmul, one)
            st = int(Mul[int(s), it])
            key = (
                oid.get(int(rho[k]), -1),
                oid.get(int(rho[l]), -1),
                oid.get(st, -1),
                int(t == s),
                int(t == int(one_minus[l])),
            )
            buckets[key].append(float(np.real(Gamma[k, l])))
    n_split = sum(1 for vs in buckets.values() if (max(vs) - min(vs)) > 1.0)
    print(
        f"  Γ buckets by (orb ρt, ρs, s/t, t=s, t=1-s): {len(buckets)} "
        f"of which split>1.0: {n_split}",
        flush=True,
    )
    # show a few split buckets
    shown = 0
    for key, vs in sorted(buckets.items(), key=lambda kv: -(max(kv[1]) - min(kv[1]))):
        spread = max(vs) - min(vs)
        if spread <= 1.0:
            continue
        print(f"    split key={key} n={len(vs)} mean={np.mean(vs):.2f} spread={spread:.2f}")
        shown += 1
        if shown >= 6:
            break

    # Linear algebra on leftover Q_O.
    # Known: Q on identity orbit = 8q². Row-sum: ∑_{r□} Q(r) = E u(ξ) ∑_{r□} u(rξ).
    # ∑_{r□, rξ∈Ω} u(rξ) = ∑_{η∈Ω} u(η) because squares act on Ω, each η hit
    # |(q-1)/2 / something|... squares have index 2, Ω is a coset of squares
    # so squares act simply transitively? |squares|=(q-1)/2=|Ω|, squares≅Ω
    # as multiplicative sets.  ∑_{r□} u(rξ) = ∑_{η∈Ω} u(η) = q(q-1) constantly.
    # Hence ∑_{r□} Q(r) = E[ u(ξ) * q(q-1) ] = 2q * q(q-1) = 2 q² (q-1).
    members_by_orb = defaultdict(list)
    for r in squares:
        members_by_orb[oid[int(r)]].append(int(r))
    row_sum = sum(Q[r] for r in squares)
    row_pred = 2 * q2 * (q - 1)
    print(f"  row-sum ∑_{{□}} Q = {row_sum:.4f}  2q²(q-1)={row_pred:.4f}", flush=True)

    # Design matrix: leftover orbits (not identity).  Constraints from:
    #   row-sum (1 eq), and mean_ξ 4p² u(ξ) = ∑_{t,s} Γ  (scalar, already used
    #   as E u = 2q).
    # New: for each leftover orbit O, the cubic-derived
    #   E[ u(ξ) * (2p ẑ(ξ) + ∑ ẑẑ)  conj ẑ(rξ) ] = 0  is 0=0.
    # Stronger: E[ |∑_t B_t(ξ)|²  u(rξ) ] = 4p² Q(r)   (cubic tautology).
    # So cubic * u(r) is tautological for every r.  It cannot cut Q-dofs
    # unless Γ_{t,s} is replaced by a Q-expression (Wick).
    #
    # Wick-close test: replace off-diagonal Γ by a linear form of Q(ρ)
    # and see predicted leftover Q vs actual.
    # Model: Γ_{t,s} ≈ A 1_{t=s} Q(ρ(t)) + B 1_{t=1-s} Q(something) + C Q(1)
    # Fit on (t,s) and predict.

    # Features for each (t,s):
    #  Qρt, Qρs, Q_{s/t}, Q_1, 1_{t=s}, 1_{t=1-s}
    Q1 = Q[int(one)]
    feats = []
    ys = []
    for k, t in enumerate(R):
        for l, s in enumerate(R):
            it = inv_of(int(t), q, fmul, one)
            st = int(Mul[int(s), it])
            feats.append(
                [
                    Q[int(rho[k])],
                    Q[int(rho[l])],
                    Q.get(int(st), Q.get(t_rep(st, Mul, neg1), 0.0)),
                    Q1,
                    float(t == s),
                    float(int(t) == int(one_minus[l])),
                ]
            )
            ys.append(float(np.real(Gamma[k, l])))
    Fmat = np.asarray(feats)
    yv = np.asarray(ys)
    coef, resid, rank, sv = np.linalg.lstsq(Fmat, yv, rcond=None)
    pred = Fmat @ coef
    err = float(np.max(np.abs(pred - yv)))
    print(f"  Γ ~ lin(Qρt,Qρs,Q_{{s/t}},Q1,1_{{t=s}},1_{{t=1-s}}): rank={rank} maxerr={err:.4f} coef={np.round(coef, 4)}", flush=True)

    # If maxerr is small, Γ is linear in Q and the squared cubic is
    # a linear relation on Q (likely the already-known E u / row-sum).
    # If maxerr is large, connected 4-point, cubic does not cut Q-dofs.

    # Tautology check: 4p² Q(r) vs mean |∑ B_t|² u(r·)
    taut_err = 0.0
    S = B.sum(axis=0)  # (M, |Ω|)  ∑_t bilinear
    absS2 = np.abs(S) ** 2
    for r in squares[:: max(1, len(squares) // 8)]:
        cr = col[Mul[r, Omega]]
        ur = u[np.arange(M)[:, None], cr]
        lhs_r = 4 * p * p * Q[int(r)]
        rhs_r = float(np.mean(absS2 * ur))
        taut_err = max(taut_err, abs(lhs_r - rhs_r))
    print(f"  tautology 4p² Q(r)=E[|∑B|² u(r·)] maxerr={taut_err:.3e}", flush=True)

    leftover = []
    for r0, orb in ((min(o), o) for o in orbits):
        if oid[int(r0)] == id_orb:
            continue
        members = [r for r in squares if oid[int(r)] == oid[int(r0)]]
        leftover.append(
            {
                "rep": int(r0),
                "size_T": len(orb),
                "Q": float(np.mean([Q[int(r)] for r in members])),
                "Wick": 4 * q2,
                "Q_minus_Wick": float(np.mean([Q[int(r)] for r in members]) - 4 * q2),
            }
        )

    return {
        "p": p,
        "M": M,
        "n_R": int(nR),
        "n_orb": n_orb,
        "cubic_maxres": max_res,
        "Q1": Q[int(one)],
        "Q1_expect": 8 * q2,
        "row_sum": row_sum,
        "row_sum_expect": row_pred,
        "Gamma_sum": rhs,
        "lhs_4p2Eu": lhs,
        "diag_Gamma_is_Qrho": diag_ok,
        "n_Gamma_buckets": len(buckets),
        "n_split_buckets": n_split,
        "Gamma_lin_Q_maxerr": err,
        "Gamma_lin_Q_rank": int(rank),
        "tautology_maxerr": taut_err,
        "leftover": leftover,
        "naive_Wick_sum": float(Wick.real.sum()),
    }


def main():
    out = {}
    for p in (5, 7):
        out[str(p)] = analyze(p)
    path = ROOT / "evidence" / "boolean_cubic_orbit_relations.json"
    # json-friendly
    path.write_text(json.dumps(out, indent=2, default=float))
    print(f"\nwrote {path}", flush=True)
    print("\n=== verdict ===", flush=True)
    for p, rec in out.items():
        print(
            f"  p={p} cubic_res={rec['cubic_maxres']:.1e} "
            f"Γ-lin-Q maxerr={rec['Gamma_lin_Q_maxerr']:.3f} "
            f"split_buckets={rec['n_split_buckets']} "
            f"diag_ok={rec['diag_Gamma_is_Qrho']}",
            flush=True,
        )


if __name__ == "__main__":
    main()
