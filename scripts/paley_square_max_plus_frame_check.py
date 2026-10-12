#!/usr/bin/env python3
"""Checks for the Max_+ frame identity at Paley orders n = p^2 + 1 (NOTE_2026-10-11_MAX_PLUS_FRAME_IDENTITY_ALL_P).

For q = p^2 the script builds GF(q) as F_p[w]/(w^2 - a w - b), the symmetric Paley conference matrix C of order q+1
(indices GF(q) and infinity), and the explicit Boolean eigenvector y = (h(Tr(lambda x))_{x in GF(q)}, sum_t h(t)),
with h : F_p -> {+1,-1} taking the value +1 exactly (p+1)/2 times (or (p-1)/2 times) and chi(lambda) = sign of the
quadratic Gauss sum.  It then checks, in exact integer arithmetic where possible:
  (1) C is symmetric, zero diagonal, C^2 = q I;
  (2) C y = p y with y in {+1,-1}^n;
  (3) the average of y y^T over the orbit of y under H = {x -> s x + b : s a nonzero square, b in GF(q)} equals I + C/p;
  (4) for random A in S_n, the orbit average of Q_A(y) equals <A,C>/p (sum over i<j), hence Phi(A) >= |<A,C>|/p;
  (5) optional (--full): all Boolean vectors in the +p eigenspace, enumerated through a pivot basis (p <= 7),
      have average yy^T equal to I + C/p.
Usage: python3 paley_square_max_plus_frame_check.py [--full] p1 p2 ...   (default 3 5 7 11 13)
Requires numpy only."""
import sys, itertools
import numpy as np

def build_field(p):
    for a in range(p):
        for b in range(1, p):
            if all((x * x - a * x - b) % p for x in range(p)):
                break
        else:
            continue
        break
    q = p * p
    def mul(u, v):
        c0, c1, d0, d1 = u % p, u // p, v % p, v // p
        return (c0 * d0 + c1 * d1 * b) % p + ((c0 * d1 + c1 * d0 + c1 * d1 * a) % p) * p
    def add(u, v):
        return (u % p + v % p) % p + ((u // p + v // p) % p) * p
    MUL = np.array([[mul(u, v) for v in range(q)] for u in range(q)])
    ADD = np.array([[add(u, v) for v in range(q)] for u in range(q)])
    NEG = np.array([next(w for w in range(q) if ADD[u][w] == 0) for u in range(q)])
    for g in range(2, q):
        x, seen = 1, set()
        for _ in range(q - 1):
            seen.add(x); x = MUL[x][g]
        if len(seen) == q - 1:
            break
    log = np.full(q, -1); x = 1
    for k in range(q - 1):
        log[x] = k; x = MUL[x][g]
    chi = np.array([0] + [1 if log[u] % 2 == 0 else -1 for u in range(1, q)])
    return q, MUL, ADD, NEG, chi

def trace(u, p, MUL, ADD):
    up = 1
    for _ in range(p):
        up = MUL[up][u]
    s = ADD[u][up] if u else 0
    assert s // p == 0
    return int(s % p)

def check(p, full=False):
    q, MUL, ADD, NEG, chi = build_field(p)
    n = q + 1
    C = np.zeros((n, n), dtype=np.int64)
    for x in range(q):
        for y in range(q):
            if x != y:
                C[x, y] = chi[ADD[x][NEG[y]]]
    C[q, :q] = 1; C[:q, q] = 1
    assert (C == C.T).all() and (np.diag(C) == 0).all() and (C @ C == q * np.eye(n, dtype=np.int64)).all()
    # explicit Boolean eigenvector (search over lambda and the two admissible h)
    found = None
    for lam in range(1, q):
        for plus in ((p + 1) // 2, (p - 1) // 2):
            h = [1 if t < plus else -1 for t in range(p)]
            f = np.array([h[trace(int(MUL[lam][x]), p, MUL, ADD)] for x in range(q)])
            y = np.append(f, sum(h))
            if np.array_equal(C @ y, p * y):
                found = (lam, y); break
        if found: break
    assert found, "no explicit Boolean eigenvector found"
    lam, y = found
    assert set(np.abs(y)) == {1}
    squares = [u for u in range(1, q) if chi[u] == 1]
    M = np.zeros((n, n)); orbit = set()
    for s in squares:
        for b in range(q):
            perm = np.array([ADD[MUL[s][x]][b] for x in range(q)] + [q])
            z = np.empty(n, dtype=np.int64); z[perm] = y
            assert np.array_equal(C @ z, p * z)
            orbit.add(tuple(z.tolist()))
    Z = np.array(sorted(orbit))
    M = Z.T @ Z / len(Z)
    err = np.abs(M - (np.eye(n) + C / p)).max()
    assert err < 1e-12, err
    # (4) averaged quadratic form
    rng = np.random.default_rng(p)
    for _ in range(3):
        A = np.triu(rng.choice([-1, 1], size=(n, n)), 1); A = A + A.T
        avgQ = np.mean([(z @ A @ z) // 2 for z in Z])
        assert abs(avgQ - np.triu(A * C, 1).sum() / p) < 1e-9
    msg = f"p={p} n={n} orbit size={len(Z)} (p(p^2-1)/4={p*(p*p-1)//4}) max|M-(I+C/p)|={err:.1e} y_inf={int(y[-1])} chi(lambda)={int(chi[lam])}"
    if full and p <= 7:
        P = (np.eye(n) + C / p) / 2
        w, V = np.linalg.eigh(P); B = V[:, w > 0.5]; d = B.shape[1]
        # greedy pivot rows
        rows = []; R = B.copy()
        Rm = B.copy()
        for _ in range(d):
            i = int(np.argmax(np.abs(Rm).sum(axis=1) * np.array([0 if r in rows else 1 for r in range(n)])))
            rows.append(i)
            # Gram-Schmidt on columns to keep pivoting stable
            col = int(np.argmax(np.abs(Rm[i])))
            Rm = Rm - np.outer(Rm[:, col] / Rm[i, col], Rm[i]); Rm[:, col] = 0
        Binv = np.linalg.inv(B[rows])
        cnt = 0; S = np.zeros((n, n))
        chunk = 1 << 16
        for start in range(0, 1 << d, chunk):
            idx = np.arange(start, min(start + chunk, 1 << d))
            sg = 1 - 2 * ((idx[:, None] >> np.arange(d)) & 1)
            Y = (sg @ Binv.T) @ B.T
            ok = np.all(np.abs(np.abs(Y) - 1) < 1e-6, axis=1)
            Yok = np.rint(Y[ok]); cnt += len(Yok); S += Yok.T @ Yok
        S /= cnt
        msg += f" | all Boolean vectors in V_+: {cnt} (incl. sign pairs), max|avg yy^T-(I+C/p)|={np.abs(S-(np.eye(n)+C/p)).max():.1e}"
    return msg

if __name__ == "__main__":
    args = sys.argv[1:]
    full = "--full" in args
    ps = [int(a) for a in args if a != "--full"] or [3, 5, 7, 11, 13]
    for p in ps:
        print(check(p, full), flush=True)
