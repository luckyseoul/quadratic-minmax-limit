#!/usr/bin/env python3
"""
Prop 15.590 — μ/δ eigen-contraction closure; equivariant determination at p=5;
level-4 moment-SDP counter-mechanism for leftover 3.

Does **not** flip type_I / phi_F_ge_6 / residual_ii / e1 / L. Soft-close forbidden.

SETUP (15.189 π=±C/p, 15.191, 15.254, 15.268)
  Paley conference C on P¹(F_q), q=p², n=q+1.  Max±={y∈{±1}ⁿ: Cy=±py}.
  For a four-set S:  m₄^±(S)=E_±[y_S],  μ=½(m₄⁺+m₄⁻),  δ=½(m₄⁺−m₄⁻).
  Six-sets analogously: μ₆, δ₆.  Signed automorphisms: UᵀCU=C (auto) or
  −C (anti, 15.254 nonsquare-scaling); antis swap Max± and twist δ by −1.

PROVED Max-free (only Cy=±py, y_i²=1, E_±[y_iy_j]=±C_ij/p from 15.189)
  A. Inside contraction closes on μ: for all distinct i,j,k
        Σ_{l∉{i,j,k}} C_{kl} μ({i,j,k,l}) = C_ij                     (★μ)
     (the l∈{i,j} boundary terms cancel between the two ensembles), and
        Σ_{l∉{i,j,k}} C_{kl} δ({i,j,k,l}) = −(2/p) C_ik C_jk .       (★δ)
  B. Outside contraction couples μ↔δ: for l∉{i,j,a}
        Σ_{b∉{i,j,a},b≠l} C_lb μ({i,j,a,b}) = p δ({i,j,a,l})
        Σ_{b∉{i,j,a},b≠l} C_lb δ({i,j,a,b}) = p μ({i,j,a,l})
                     − (1/p)(C_li C_ja + C_lj C_ia + C_la C_ij).
  C. Degree-6 versions: for 5 distinct S₅∋k and l∉S₅
        Σ_{l'∉S₅} C_kl' μ₆(S₅∪l') + Σ_{x∈S₅∖k} C_kx μ(S₅∖x) = p δ(S₅∖k)
        Σ_{b∉S₅∪l} C_lb μ₆(S₅∪b) + Σ_{x∈S₅} C_lx μ(S₅∖x) = p δ₆(S₅∪l)
     and the δ-versions with μ↔δ.  (Odd moments vanish by y→−y.)
  D. s-identity: (Cy)_∞=Σ_fin y ⇒ 1ᵀy=(1±p)y_∞ pointwise on Max±, so
     E_±[e₄] is exact: Nh((s²)²−(6n−8)s²+3n²−6n)/24, s²=(p±1)².
  All verified EXACTLY against enumerated Max± at p=5 and p=7 (zero
  violations, exhaustive for (★) at p=5).

THEOREM (p=5, exact elimination; Max-free given the identities)
  E. The complete equivariant system {A,B,C} on signed-orbit variables
     (4 μ-orbits of four-sets, 2 live δ-orbits, 13 live μ₆-orbits,
     19 δ₆-orbits) has RANK = #unknowns = 38: μ, δ, μ₆, δ₆ at p=5 are
     the UNIQUE equivariant solution.  In particular μ on |κ|=1 is
     determined without any census: symmetry + eigen-identities suffice.
  E'. δ is twisted-dead on every |κ|=1 orbit at p=5,7 ⇒ m₄⁺=m₄⁻ on
     |κ|=1 by pure sign-consistency (mechanism independent of 15.268).

CERTIFIED (exact elimination on complete representative-built systems)
  F. Degree-4-only equivariant kernel: dim 1 (p=5), dim 2 (p=7); the
     kernel touches the |κ|=1 μ-coordinates at both primes.
  G. Joint degree-4+6 kernel at p=7: dim 4 (NOT determined; the four
     free directions touch every 4-point coordinate).  Determination
     degree grows with p (conjecturally level p+1; open).

COUNTER-MECHANISM (kills the level-4 moment/SoS route for leftover 3)
  H. The degree-4 relaxation {complete linear theory of A,B + signed
     equivariance + M±=E±[vec(yyᵀ)vec(yyᵀ)ᵀ] ⪰ 0} admits feasible
     points with |μ| > 2/n AND |μ| > L=(p−2)/2p² at p=5 and p=7
     (p=7: feasible |μ_int| 2466 vs L·N=1169 vs true 872).  Any proof
     of |μ|≤2/n or |μ|≤L must use MORE than degree-4 moment closure:
     do not reopen Wick/hull/Ext-style level-4 majorants for leftover 3.

OPEN
  I. Kernel-dimension formula for the complete degree-≤D equivariant
     system (15.589-style character algebra should give it); smallest
     D(p) with kernel 0, if any.  D(5)=6.  If D(p) exists with a
     closed-form solution, leftover 3 reduces to bounding that solution.

Writes evidence/e1_gmin_m4_prop15590.json
"""
from __future__ import annotations

import itertools
import json
import os
import random
import sys
import time
from collections import deque
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]

# ---------------------------------------------------------------- field / C
def _nonresidue(p: int) -> int:
    for r in range(2, p):
        if pow(r, (p - 1) // 2, p) == p - 1:
            return r
    raise ValueError


def field_ops(p: int):
    """F_{p^2}=F_p[t]/(t^2-r); element e=p*a+b <-> a+bt."""
    r = _nonresidue(p)

    def fmul(e1, e2):
        a1, b1 = divmod(e1, p)
        a2, b2 = divmod(e2, p)
        return p * ((a1 * a2 + r * b1 * b2) % p) + ((a1 * b2 + a2 * b1) % p)

    def fadd(e1, e2):
        a1, b1 = divmod(e1, p)
        a2, b2 = divmod(e2, p)
        return p * ((a1 + a2) % p) + ((b1 + b2) % p)

    def fneg(e):
        a, b = divmod(e, p)
        return p * ((-a) % p) + ((-b) % p)

    one = p  # a=1,b=0
    return fmul, fadd, fneg, one


def paley_conference(p: int) -> np.ndarray:
    q = p * p
    n = q + 1
    fmul, fadd, fneg, one = field_ops(p)
    squares = {fmul(x, x) for x in range(1, q)}
    C = np.zeros((n, n), dtype=np.int64)
    C[0, 1:] = 1
    C[1:, 0] = 1
    for e1 in range(q):
        for e2 in range(q):
            if e1 != e2:
                C[1 + e1, 1 + e2] = 1 if fadd(e1, fneg(e2)) in squares else -1
    assert (C == C.T).all() and (C @ C == q * np.eye(n, dtype=np.int64)).all()
    return C


# ---------------------------------------------------- Max± via V±-completion
def _pivoted_basis(Q1: np.ndarray) -> list[int]:
    """Greedy row-pivot selection: rows of Q1 (n x m) spanning R^m."""
    n, m = Q1.shape
    R = Q1.copy()
    picked: list[int] = []
    for _ in range(m):
        norms = (R * R).sum(axis=1)
        norms[picked] = -1.0
        i = int(np.argmax(norms))
        picked.append(i)
        v = R[i] / np.linalg.norm(R[i])
        R = R - np.outer(R @ v, v)
    return picked


def enum_max(p: int, sign: int) -> np.ndarray:
    """All y in {±1}^n with Cy = sign*p*y, by completing over a coordinate
    basis of the eigenspace (dim n/2)."""
    C = paley_conference(p)
    n = C.shape[0]
    A = (p * np.eye(n) + sign * C) / (2 * p)  # projector onto the eigenspace
    w, V = np.linalg.eigh(A)
    Q1 = V[:, w > 0.5]
    m = Q1.shape[1]
    assert m == n // 2
    B = _pivoted_basis(Q1)
    Bc = [i for i in range(n) if i not in B]
    EXT = Q1 @ np.linalg.inv(Q1[B, :])
    EXTc = EXT[Bc, :].astype(np.float64)
    out = []
    CH = 1 << 19
    bits = np.arange(m)
    for lo in range(0, 1 << m, CH):
        idx = np.arange(lo, min(lo + CH, 1 << m), dtype=np.int64)
        YB = (((idx[:, None] >> bits[None, :]) & 1) * 2 - 1).astype(np.float64)
        YR = YB @ EXTc.T
        ok = (np.abs(np.abs(YR) - 1.0) < 1e-3).all(axis=1)
        for jj in np.where(ok)[0]:
            y = np.zeros(n, dtype=np.int64)
            y[B] = YB[jj].astype(np.int64)
            y[Bc] = np.sign(YR[jj]).astype(np.int64)
            if (C @ y == sign * p * y).all():
                out.append(y.copy())
    return np.array(out, dtype=np.int8)


# ------------------------------------------------------------- moments (int)
def four_set_sums(Y: np.ndarray, n: int):
    """Pair-Gram method: exact in float32 because |sums| <= #vectors << 2^24."""
    S4 = np.array(list(itertools.combinations(range(n), 4)))
    pairs = list(itertools.combinations(range(n), 2))
    pidx = {pr: i for i, pr in enumerate(pairs)}
    P = np.array(pairs)
    X2 = (Y[:, P[:, 0]].astype(np.float32) * Y[:, P[:, 1]].astype(np.float32))
    G = X2.T @ X2
    pij = np.array([pidx[(a, b)] for a, b in zip(S4[:, 0], S4[:, 1])])
    pkl = np.array([pidx[(a, b)] for a, b in zip(S4[:, 2], S4[:, 3])])
    return S4, G[pij, pkl].astype(np.int64)


def six_set_sums(Y: np.ndarray, n: int):
    triples = list(itertools.combinations(range(n), 3))
    T3 = np.array(triples)
    X3 = (Y[:, T3[:, 0]].astype(np.int64) * Y[:, T3[:, 1]] * Y[:, T3[:, 2]]).astype(np.float32)
    G3 = X3.T @ X3
    tkeys = (T3[:, 0].astype(np.int64) * n + T3[:, 1]) * n + T3[:, 2]
    order_t = np.argsort(tkeys)
    tk_sorted = tkeys[order_t]
    S6 = np.array(list(itertools.combinations(range(n), 6)), dtype=np.int16)
    a_keys = (S6[:, 0].astype(np.int64) * n + S6[:, 1]) * n + S6[:, 2]
    b_keys = (S6[:, 3].astype(np.int64) * n + S6[:, 4]) * n + S6[:, 5]
    a_idx = order_t[np.searchsorted(tk_sorted, a_keys)]
    b_idx = order_t[np.searchsorted(tk_sorted, b_keys)]
    return S6, G3[a_idx, b_idx].astype(np.int64)


# ------------------------------------------------------- signed automorphisms
def signed_generators(p: int, C: np.ndarray):
    q = p * p
    n = q + 1
    fmul, fadd, fneg, one = field_ops(p)
    finv = [0] * q
    for e in range(1, q):
        finv[e] = next(x for x in range(1, q) if fmul(e, x) == one)

    def order_of(e):
        x, o = e, 1
        while x != one:
            x = fmul(x, e)
            o += 1
        return o

    gen = next(e for e in range(2, q) if order_of(e) == q - 1)
    g2 = fmul(gen, gen)

    def mkperm(fn, inf_to=0):
        pi = np.zeros(n, dtype=np.int64)
        pi[0] = inf_to
        for e in range(q):
            pi[1 + e] = fn(e)
        return pi

    def frob(e):
        x = e
        for _ in range(p - 1):
            x = fmul(x, e)
        return x

    perms = [
        mkperm(lambda e: 1 + fadd(e, one)),        # z -> z+1
        mkperm(lambda e: 1 + fadd(e, 1)),          # z -> z+t
        mkperm(lambda e: 1 + fmul(e, g2)),         # z -> g^2 z  (auto)
        mkperm(lambda e: 1 + fmul(e, gen)),        # z -> g z    (anti)
        mkperm(lambda e: (1 + finv[e]) if e else 0, inf_to=1),  # z -> 1/z
        mkperm(lambda e: 1 + frob(e)),             # Frobenius
    ]

    def signed_lift(pi):
        for s in (1, -1):
            d = np.zeros(n, dtype=np.int64)
            d[0] = 1
            d[1:] = s * C[pi[0], pi[1:]] * C[0, 1:]
            ok = s * d[:, None] * d[None, :] * C == C[np.ix_(pi, pi)]
            np.fill_diagonal(ok, True)
            if ok.all():
                return d, s
        return None, 0

    gens = []
    for pi in perms:
        d, s = signed_lift(pi)
        assert s != 0, "signed lift must exist"
        gens.append((pi, d, s))
        ipi = np.argsort(pi)
        gens.append((ipi, d[ipi], s))
    return gens


def signed_orbits(SETS: np.ndarray, gens, n: int, twist: bool):
    """Label-propagation signed orbits.  Returns (label, sign, dead_labels)."""
    k = SETS.shape[1]
    mult = np.array([n ** e for e in range(k - 1, -1, -1)], dtype=np.int64)
    enc = (SETS.astype(np.int64) * mult[None, :]).sum(axis=1)
    tgts, epss, sss = [], [], []
    for pi, d, s in gens:
        img = np.sort(pi[SETS.astype(np.int64)], axis=1)
        ikey = (img * mult[None, :]).sum(axis=1)
        tgts.append(np.searchsorted(enc, ikey).astype(np.int64))
        epss.append(d[SETS.astype(np.int64)].prod(axis=1).astype(np.int64))
        sss.append(s)
    lab = np.arange(len(SETS), dtype=np.int64)
    sg = np.ones(len(SETS), dtype=np.int64)
    dead = np.zeros(len(SETS), dtype=bool)
    for _ in range(300):
        changed = False
        for tgt, eps, s in zip(tgts, epss, sss):
            e = eps * (s if twist else 1)
            lt = lab[tgt]
            st = sg[tgt] * e
            m = lt < lab
            if m.any():
                lab = np.where(m, lt, lab)
                sg = np.where(m, st, sg)
                changed = True
            m2 = (lt == lab) & (~m) & (st != sg)
            if m2.any():
                dead |= m2
        l2 = lab[lab]
        s2 = sg[lab] * sg
        if (l2 != lab).any():
            lab, sg = l2, s2
            changed = True
        if not changed:
            break
    dead_labels = set(int(x) for x in np.unique(lab[dead]))
    return lab, sg, dead_labels


# ------------------------------------------------------------- the lab per p
class MuLab:
    def __init__(self, p: int, with_deg6: bool = True):
        self.p = p
        self.q = p * p
        self.n = self.q + 1
        self.C = paley_conference(p)
        self.Yp = enum_max(p, +1)
        self.Ym = enum_max(p, -1)
        assert len(self.Yp) == len(self.Ym)
        self.N = 2 * len(self.Yp)
        n = self.n
        self.S4, Sp = four_set_sums(self.Yp, n)
        _, Sm = four_set_sums(self.Ym, n)
        self.mu4 = Sp + Sm
        self.d4 = Sp - Sm
        i_, j_, k_, l_ = (self.S4[:, 0], self.S4[:, 1], self.S4[:, 2], self.S4[:, 3])
        self.kap = (self.C[i_, j_] * self.C[k_, l_] + self.C[i_, k_] * self.C[j_, l_]
                    + self.C[i_, l_] * self.C[j_, k_])
        key4 = (self.S4[:, 0] * n ** 3 + self.S4[:, 1] * n ** 2 + self.S4[:, 2] * n + self.S4[:, 3])
        self.key4idx = {int(k): i for i, k in enumerate(key4)}
        self.gens = signed_generators(p, self.C)
        self.labM4, self.sgM4, self.deadM4 = signed_orbits(self.S4, self.gens, n, twist=False)
        self.labD4, self.sgD4, self.deadD4 = signed_orbits(self.S4, self.gens, n, twist=True)
        if with_deg6:
            self.S6, Sp6 = six_set_sums(self.Yp, n)
            _, Sm6 = six_set_sums(self.Ym, n)
            self.mu6 = Sp6 + Sm6
            self.d6 = Sp6 - Sm6
            mult6 = np.array([n ** e for e in range(5, -1, -1)], dtype=np.int64)
            self.enc6 = (self.S6.astype(np.int64) * mult6[None, :]).sum(axis=1)
            self.mult6 = mult6
            self.labM6, self.sgM6, self.deadM6 = signed_orbits(self.S6, self.gens, n, twist=False)
            self.labD6, self.sgD6, self.deadD6 = signed_orbits(self.S6, self.gens, n, twist=True)

    def i4(self, *args) -> int:
        x = sorted(args)
        n = self.n
        return self.key4idx[x[0] * n ** 3 + x[1] * n ** 2 + x[2] * n + x[3]]

    def i6(self, *args) -> int:
        x = sorted(args)
        return int(np.searchsorted(self.enc6, sum(int(v) * int(m) for v, m in zip(x, self.mult6))))

    # ---- identity checks (exact integers; return violation counts)
    def check_star_mu(self, exhaustive: bool = False, samples: int = 400, seed: int = 0):
        n, C, N, p = self.n, self.C, self.N, self.p
        rng = random.Random(seed)
        bad = tested = 0
        if exhaustive:
            configs = ((k, i, j) for k in range(n)
                       for i, j in itertools.combinations([x for x in range(n) if x != k], 2))
        else:
            def gen():
                for _ in range(samples):
                    k, i, j = rng.sample(range(n), 3)
                    yield k, i, j
            configs = gen()
        for k, i, j in configs:
            lhs = sum(int(C[k, l]) * int(self.mu4[self.i4(i, j, k, l)])
                      for l in range(n) if l not in (i, j, k) and C[k, l])
            tested += 1
            if lhs != N * int(C[i, j]):
                bad += 1
        return tested, bad

    def check_out_mu(self, samples: int = 400, seed: int = 1):
        n, C, p = self.n, self.C, self.p
        rng = random.Random(seed)
        bad = 0
        for _ in range(samples):
            i, j, a, l = rng.sample(range(n), 4)
            lhs = sum(int(C[l, b]) * int(self.mu4[self.i4(i, j, a, b)])
                      for b in range(n) if b not in (i, j, a, l) and C[l, b])
            if lhs != p * int(self.d4[self.i4(i, j, a, l)]):
                bad += 1
        return samples, bad

    def check_deg6(self, samples: int = 200, seed: int = 2):
        n, C, p = self.n, self.C, self.p
        rng = random.Random(seed)
        bad = [0, 0, 0, 0]
        for _ in range(samples):
            S5 = rng.sample(range(n), 5)
            k = S5[0]
            S4r = S5[1:]
            lhs = sum(int(C[k, l]) * int(self.mu6[self.i6(*S5, l)])
                      for l in range(n) if l not in S5 and C[k, l])
            lhs += sum(int(C[k, x]) * int(self.mu4[self.i4(*[y for y in S5 if y != x])])
                       for x in S4r if C[k, x])
            if lhs != p * int(self.d4[self.i4(*S4r)]):
                bad[0] += 1
            lhs = sum(int(C[k, l]) * int(self.d6[self.i6(*S5, l)])
                      for l in range(n) if l not in S5 and C[k, l])
            lhs += sum(int(C[k, x]) * int(self.d4[self.i4(*[y for y in S5 if y != x])])
                       for x in S4r if C[k, x])
            if lhs != p * int(self.mu4[self.i4(*S4r)]):
                bad[1] += 1
            l = rng.choice([x for x in range(n) if x not in S5])
            lhs = sum(int(C[l, b]) * int(self.mu6[self.i6(*S5, b)])
                      for b in range(n) if b not in S5 and b != l and C[l, b])
            lhs += sum(int(C[l, x]) * int(self.mu4[self.i4(*[y for y in S5 if y != x])])
                       for x in S5 if C[l, x])
            if lhs != p * int(self.d6[self.i6(*S5, l)]):
                bad[2] += 1
            lhs = sum(int(C[l, b]) * int(self.d6[self.i6(*S5, b)])
                      for b in range(n) if b not in S5 and b != l and C[l, b])
            lhs += sum(int(C[l, x]) * int(self.d4[self.i4(*[y for y in S5 if y != x])])
                       for x in S5 if C[l, x])
            if lhs != p * int(self.mu6[self.i6(*S5, l)]):
                bad[3] += 1
        return samples, bad

    def s_identity(self):
        sP = set(self.Yp.astype(np.int64).sum(axis=1).tolist())
        sM = set(self.Ym.astype(np.int64).sum(axis=1).tolist())
        return sP == {1 + self.p, -(1 + self.p)} and sM == {1 - self.p, self.p - 1}

    # ---- unknown indexing over live orbits
    def _index_maps(self, deg6: bool):
        ulM4 = sorted(set(np.unique(self.labM4).tolist()) - self.deadM4)
        ulD4 = sorted(set(np.unique(self.labD4).tolist()) - self.deadD4)
        maps = {"M4": {l: i for i, l in enumerate(ulM4)}}
        off = len(ulM4)
        maps["D4"] = {l: off + i for i, l in enumerate(ulD4)}
        off += len(ulD4)
        if deg6:
            ulM6 = sorted(set(np.unique(self.labM6).tolist()) - self.deadM6)
            ulD6 = sorted(set(np.unique(self.labD6).tolist()) - self.deadD6)
            maps["M6"] = {l: off + i for i, l in enumerate(ulM6)}
            off += len(ulM6)
            maps["D6"] = {l: off + i for i, l in enumerate(ulD6)}
            off += len(ulD6)
        return maps, off

    def _c(self, cvec, kind, maps, si, v):
        lab, sg, dead = {
            "M4": (self.labM4, self.sgM4, self.deadM4),
            "D4": (self.labD4, self.sgD4, self.deadD4),
            "M6": (self.labM6, self.sgM6, self.deadM6) if hasattr(self, "labM6") else (None,) * 3,
            "D6": (self.labD6, self.sgD6, self.deadD6) if hasattr(self, "labD6") else (None,) * 3,
        }[kind]
        l = int(lab[si])
        if l not in dead:
            cvec[maps[kind][l]] += v * int(sg[si])

    def equivariant_system(self, deg6: bool):
        """Complete system via orbit representatives.  Returns dict with rank,
        kernel vectors, and the unique solution when kernel is 0."""
        n, C, p, N = self.n, self.C, self.p, self.N
        maps, NV = self._index_maps(deg6)
        rows = set()
        # (★μ),(★δ): representative-complete = all (k,{i,j}) is cheap only for
        # p=5; use marked-triple representatives via 3-set orbits for larger p.
        S3 = np.array(list(itertools.combinations(range(n), 3)))
        lab3, _, _ = signed_orbits(S3, self.gens, n, twist=False)
        reps3 = [S3[int(np.where(lab3 == l)[0][0])] for l in np.unique(lab3)]
        for T in reps3:
            for k in map(int, T):
                i, j = [int(x) for x in T if x != k]
                cM = [0] * NV
                cD = [0] * NV
                for l in range(n):
                    if l in (i, j, k) or C[k, l] == 0:
                        continue
                    si = self.i4(i, j, k, l)
                    self._c(cM, "M4", maps, si, int(C[k, l]))
                    self._c(cD, "D4", maps, si, int(C[k, l]))
                rows.add((tuple(cM), N * int(C[i, j])))
                rows.add((tuple(cD), -(2 * N // p) * int(C[i, k]) * int(C[j, k])))
        # (out): marked-4-set representatives
        for l4 in np.unique(self.labM4):
            T = [int(x) for x in self.S4[int(np.where(self.labM4 == l4)[0][0])]]
            for l_ in T:
                tri = [x for x in T if x != l_]
                st = self.i4(*T)
                cM = [0] * NV
                for x in range(n):
                    if x in T or C[l_, x] == 0:
                        continue
                    self._c(cM, "M4", maps, self.i4(*tri, x), int(C[l_, x]))
                self._c(cM, "D4", maps, st, -p)
                rows.add((tuple(cM), 0))
                cD = [0] * NV
                for x in range(n):
                    if x in T or C[l_, x] == 0:
                        continue
                    self._c(cD, "D4", maps, self.i4(*tri, x), int(C[l_, x]))
                self._c(cD, "M4", maps, st, -p)
                a, b, c_ = tri
                corr = (int(C[l_, a]) * int(C[b, c_]) + int(C[l_, b]) * int(C[a, c_])
                        + int(C[l_, c_]) * int(C[a, b]))
                rows.add((tuple(cD), -(N // p) * corr))
        if deg6:
            S5 = np.array(list(itertools.combinations(range(n), 5)), dtype=np.int16)
            lab5, _, _ = signed_orbits(S5, self.gens, n, twist=False)
            reps5 = [S5[int(np.where(lab5 == l)[0][0])] for l in np.unique(lab5)]
            for R5 in reps5:
                S5l = [int(x) for x in R5]
                for k in S5l:
                    S4r = [y for y in S5l if y != k]
                    for KIND6, KIND4a, KIND4b in (("M6", "M4", "D4"), ("D6", "D4", "M4")):
                        cv = [0] * NV
                        for l in range(n):
                            if l in S5l or C[k, l] == 0:
                                continue
                            self._c(cv, KIND6, maps, self.i6(*S5l, l), int(C[k, l]))
                        for x in S4r:
                            if C[k, x]:
                                self._c(cv, KIND4a, maps,
                                        self.i4(*[y for y in S5l if y != x]), int(C[k, x]))
                        self._c(cv, KIND4b, maps, self.i4(*S4r), -p)
                        rows.add((tuple(cv), 0))
            for l6 in np.unique(self.labM6):
                T = [int(x) for x in self.S6[int(np.where(self.labM6 == l6)[0][0])]]
                for l_ in T:
                    S5l = [y for y in T if y != l_]
                    st6 = self.i6(*T)
                    for KIND6, KIND4 in (("M6", "M4"), ("D6", "D4")):
                        cv = [0] * NV
                        for b in range(n):
                            if b in S5l or b == l_ or C[l_, b] == 0:
                                continue
                            self._c(cv, KIND6, maps, self.i6(*S5l, b), int(C[l_, b]))
                        for x in S5l:
                            if C[l_, x]:
                                self._c(cv, KIND4, maps,
                                        self.i4(*[y for y in S5l if y != x]), int(C[l_, x]))
                        self._c(cv, "D6" if KIND6 == "M6" else "M6", maps, st6, -p)
                        rows.add((tuple(cv), 0))
        # exact elimination
        M = [[Fraction(x) for x in r] + [Fraction(rhs)] for r, rhs in rows]
        r = 0
        piv_cols = []
        for c in range(NV):
            piv = next((i for i in range(r, len(M)) if M[i][c] != 0), None)
            if piv is None:
                continue
            M[r], M[piv] = M[piv], M[r]
            M[r] = [x / M[r][c] for x in M[r]]
            for i in range(len(M)):
                if i != r and M[i][c] != 0:
                    M[i] = [u - M[i][c] * v for u, v in zip(M[i], M[r])]
            piv_cols.append(c)
            r += 1
        incons = any(all(x == 0 for x in row[:-1]) and row[-1] != 0 for row in M)
        free_cols = [c for c in range(NV) if c not in piv_cols]
        kvs = []
        for fc in free_cols:
            kv = [Fraction(0)] * NV
            kv[fc] = Fraction(1)
            for ri, c in enumerate(piv_cols):
                kv[c] = -M[ri][fc]
            kvs.append(kv)
        sol = None
        particular = None
        if not incons:
            particular = [Fraction(0)] * NV
            for ri, c in enumerate(piv_cols):
                particular[c] = M[ri][NV]
            if not free_cols:
                sol = particular
        return dict(maps=maps, NV=NV, rank=r, kernel_dim=NV - r,
                    inconsistent=incons, kernel=kvs, solution=sol,
                    particular=particular, free_cols=free_cols, nrows=len(rows))

    def solution_matches_data(self, sysres) -> bool:
        """Reconstruct per-set μ from the unique solution and compare EXACTLY."""
        if sysres["solution"] is None:
            return False
        maps = sysres["maps"]
        sol = sysres["solution"]
        for si in range(len(self.S4)):
            l = int(self.labM4[si])
            v = Fraction(0) if l in self.deadM4 else sol[maps["M4"][l]] * int(self.sgM4[si])
            if v != Fraction(int(self.mu4[si])):
                return False
        return True

    def kappa1_delta_dead(self) -> bool:
        """δ twisted-dead on every |κ|=1 orbit, and d4 data vanishes there."""
        k1 = np.abs(self.kap) == 1
        labs = set(int(x) for x in np.unique(self.labD4[k1]))
        return labs <= self.deadD4 and (self.d4[k1] == 0).all()


# ------------------------------------------------- level-4 PSD counter-mech
def psd_kill(lab: MuLab, deg4sys) -> dict:
    """Exhibit a PSD-feasible point of the degree-4 relaxation whose |μ|
    exceeds both 2/n and L thresholds.  Returns the certificate data."""
    p, n, C, N = lab.p, lab.n, lab.C.astype(float), lab.N
    Nh = N // 2
    maps = deg4sys["maps"]
    kvs = deg4sys["kernel"]
    if not kvs:
        return {"killed": False, "reason": "no kernel"}
    # kernel set-functions
    def setfuncs(kv):
        fM = np.zeros(len(lab.S4))
        fD = np.zeros(len(lab.S4))
        for si in range(len(lab.S4)):
            l = int(lab.labM4[si])
            if l not in lab.deadM4:
                fM[si] = float(kv[maps["M4"][l]]) * lab.sgM4[si]
            l = int(lab.labD4[si])
            if l not in lab.deadD4:
                fD[si] = float(kv[maps["D4"][l]]) * lab.sgD4[si]
        return fM, fD
    KF = [setfuncs(kv) for kv in kvs]

    def T4_build(fvals, pair):
        T = np.zeros((n, n, n, n))
        idx = lab.S4
        for perm in itertools.permutations(range(4)):
            T[idx[:, perm[0]], idx[:, perm[1]], idx[:, perm[2]], idx[:, perm[3]]] = fvals
        for x in range(n):  # one equal pair (least degenerate repeats first)
            T[x, x, :, :] = pair
            T[x, :, x, :] = pair
            T[x, :, :, x] = pair
            T[:, x, x, :] = pair
            T[:, x, :, x] = pair
            T[:, :, x, x] = pair
        for x in range(n):  # two pairs
            T[x, x, :, :][np.arange(n), np.arange(n)] = 1.0
            T[x, :, x, :][np.arange(n), np.arange(n)] = 1.0
            T[x, :, :, x][np.arange(n), np.arange(n)] = 1.0
        for x in range(n):  # triple + singleton
            T[:, x, x, x] = pair[:, x]
            T[x, :, x, x] = pair[x, :]
            T[x, x, :, x] = pair[x, :]
            T[x, x, x, :] = pair[x, :]
        for x in range(n):
            T[x, x, x, x] = 1.0
        return T

    def compress(T, Q):
        m = Q.shape[1]
        X = np.einsum('ijkl,ia->ajkl', T, Q)
        X = np.einsum('ajkl,jb->abkl', X, Q)
        X = np.einsum('abkl,kc->abcl', X, Q)
        X = np.einsum('abcl,ld->abcd', X, Q)
        pr = [(a, b) for a in range(m) for b in range(a, m)]
        R = np.zeros((len(pr), len(pr)))
        for u, (a, b) in enumerate(pr):
            fu = np.sqrt(2) if a != b else 1.0
            for v, (c, d) in enumerate(pr):
                fv = np.sqrt(2) if c != d else 1.0
                R[u, v] = X[a, b, c, d] * fu * fv
        return R

    A = (p * np.eye(n) + C) / (2 * p)
    w, V = np.linalg.eigh(A)
    Qp = V[:, w > 0.5]
    w2, V2 = np.linalg.eigh((p * np.eye(n) - C) / (2 * p))
    Qm = V2[:, w2 > 0.5]
    R0p = compress(T4_build((lab.mu4 + lab.d4) / N, C / p), Qp)
    R0m = compress(T4_build((lab.mu4 - lab.d4) / N, -C / p), Qm)
    RKp = [compress(T4_build((fM + fD) / N, np.zeros((n, n))), Qp) for fM, fD in KF]
    RKm = [compress(T4_build((fM - fD) / N, np.zeros((n, n))), Qm) for fM, fD in KF]
    # subtract the spurious repeated-pattern fills for kernel tensors
    # (T4_build with pair=0 still sets (2,2)->1 and (4)->1)
    corr = compress(T4_build(np.zeros(len(lab.S4)), np.zeros((n, n))), Qp)
    corrm = compress(T4_build(np.zeros(len(lab.S4)), np.zeros((n, n))), Qm)
    RKp = [R - corr for R in RKp]
    RKm = [R - corrm for R in RKm]

    def lmin_at(tvec):
        Mp = R0p + sum(t * R for t, R in zip(tvec, RKp))
        Mm = R0m + sum(t * R for t, R in zip(tvec, RKm))
        return min(np.linalg.eigvalsh(Mp)[0], np.linalg.eigvalsh(Mm)[0])

    thr_2n = N * 2.0 / n
    thr_L = N * (p - 2) / (2.0 * p * p)
    k1 = np.abs(lab.kap) == 1
    k1labs = sorted(set(int(x) for x in np.unique(lab.labM4[k1])) - lab.deadM4)
    dim = len(kvs)
    best = {"val": 0.0, "t": None, "orb": None}
    rng = np.linspace(0, 2 * np.pi, 181)[:-1]
    dirs = ([(np.cos(a), np.sin(a)) for a in rng] if dim == 2
            else [(1.0,), (-1.0,)] if dim == 1 else None)
    assert dirs is not None, f"unhandled kernel dim {dim}"
    for d in dirs:
        lo, hi = 0.0, 1.0
        while lmin_at([hi * x for x in d]) > -1e-9 and hi < 1e6:
            hi *= 2
        for _ in range(60):
            mid = (lo + hi) / 2
            if lmin_at([mid * x for x in d]) > -1e-9:
                lo = mid
            else:
                hi = mid
        t_edge = [0.98 * lo * x for x in d]
        for l4 in k1labs:
            si = int(np.where(lab.labM4 == l4)[0][0])
            pert = sum(t * float(kv[maps["M4"][l4]]) for t, kv in zip(t_edge, kvs))
            val = abs(float(lab.mu4[si]) + int(lab.sgM4[si]) * pert)
            if val > best["val"]:
                best = {"val": val, "t": t_edge, "orb": int(l4)}
    feas = lmin_at(best["t"]) > -1e-8
    return {
        "killed": bool(feas and best["val"] > thr_2n and best["val"] > thr_L),
        "feasible_mu_int": best["val"],
        "threshold_2n_int": thr_2n,
        "threshold_L_int": thr_L,
        "true_max_mu_int": float(np.abs(lab.mu4[k1]).max()),
        "lambda_min_at_point": float(lmin_at(best["t"])),
    }


# ------------------------------------------------------------------- main
def main():
    t0 = time.time()
    out = {"prop": "15.590",
           "title": "mu/delta contraction closure; determination at p=5; level-4 SDP kill"}
    full = os.environ.get("PROP15590_FULL", "") == "1"

    lab5 = MuLab(5, with_deg6=True)
    tested, bad = lab5.check_star_mu(exhaustive=True)
    out["star_mu_p5"] = {"tested": tested, "violations": bad}
    _, bado = lab5.check_out_mu(samples=800)
    _, bad6 = lab5.check_deg6(samples=400)
    out["out_mu_p5_violations"] = bado
    out["deg6_p5_violations"] = bad6
    out["s_identity_p5"] = lab5.s_identity()
    out["kappa1_delta_dead_p5"] = lab5.kappa1_delta_dead()
    d4sys5 = lab5.equivariant_system(deg6=False)
    out["deg4_kernel_p5"] = d4sys5["kernel_dim"]
    j5 = lab5.equivariant_system(deg6=True)
    out["joint_kernel_p5"] = j5["kernel_dim"]
    out["determined_p5"] = lab5.solution_matches_data(j5)
    out["psd_kill_p5"] = psd_kill(lab5, d4sys5)

    lab7 = MuLab(7, with_deg6=full)
    t7, b7 = lab7.check_star_mu(samples=300)
    out["star_mu_p7"] = {"tested": t7, "violations": b7}
    out["s_identity_p7"] = lab7.s_identity()
    out["kappa1_delta_dead_p7"] = lab7.kappa1_delta_dead()
    d4sys7 = lab7.equivariant_system(deg6=False)
    out["deg4_kernel_p7"] = d4sys7["kernel_dim"]
    out["psd_kill_p7"] = psd_kill(lab7, d4sys7)
    if full:
        _, bad67 = lab7.check_deg6(samples=200)
        out["deg6_p7_violations"] = bad67
        j7 = lab7.equivariant_system(deg6=True)
        out["joint_kernel_p7"] = j7["kernel_dim"]

    out["flags_not_flipped"] = ["type_I", "phi_F_ge_6", "residual_ii", "e1", "L"]
    out["L_status"] = "OPEN"
    out["seconds"] = round(time.time() - t0, 1)
    path = ROOT / "evidence" / "e1_gmin_m4_prop15590.json"
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print("Prop 15.590  contraction closure; p=5 determination; level-4 SDP kill")
    for k in ("star_mu_p5", "out_mu_p5_violations", "deg6_p5_violations",
              "deg4_kernel_p5", "joint_kernel_p5", "determined_p5",
              "deg4_kernel_p7", "joint_kernel_p7" if full else "star_mu_p7"):
        if k in out:
            print(f"  {k} = {out[k]}")
    print(f"  psd_kill_p5: killed={out['psd_kill_p5']['killed']}")
    print(f"  psd_kill_p7: killed={out['psd_kill_p7']['killed']}")
    print(f"  ({out['seconds']}s)")
    return out


if __name__ == "__main__":
    main()
