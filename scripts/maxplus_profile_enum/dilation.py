"""Dilation/Frobenius gauge: the 120-element group (square multiplications x
Frobenius) acting on solutions, directions, and outer coefficient tuples.

Element rep: (perm_fin, dirperm, alphas)
  perm_fin: q-array, x -> gamma(x) on the finite plane (no signs)
  dirperm:  m-array, direction j -> gamma(j)
  alphas:   m-array in F_p^*, with t_{dirperm[j]}(gamma x) = alphas[j]*t_j(x).
Action on profiles: rho'_{dirperm[j]}(t) = rho_j(alphas[j]^{-1} t), so a
level-d coefficient c_{j,d} maps to position dirperm[j] with value
c_{j,d} * alphas[j]^{-d}.
"""
import numpy as np
import sys
sys.path.insert(0,'/tmp/e1work')
from kgen import field_ctx


def build_group(p):
    q, mul, chi, tr = field_ctx(p)
    m = (p + 1) // 2

    def powm(u, e):
        r, base = 1, u
        while e:
            if e & 1:
                r = mul(r, base)
            base = mul(base, base)
            e >>= 1
        return r

    # square directions in the same order as kgen.square_coords
    dirs = []
    coords = []
    seen = set()
    for g in range(1, q):
        if g in seen:
            continue
        line = [mul(t, g) for t in range(1, p)]
        seen.update(line)
        if chi(g) == 1:
            cj = next(c for c in range(1, q) if tr(mul(c, g)) == 0)
            dirs.append(g)
            coords.append(np.array([tr(mul(cj, x)) for x in range(q)],
                                   dtype=np.int64))
    assert len(dirs) == m
    sq_elems = [c for c in range(1, q) if chi(c) == 1]
    group = []
    for c in sq_elems:
        for frob in (False, True):
            if frob:
                perm = np.array([powm(mul(c, x), p) for x in range(q)],
                                dtype=np.int64)
            else:
                perm = np.array([mul(c, x) for x in range(q)], dtype=np.int64)
            # direction permutation + alphas
            dirperm = np.zeros(m, dtype=np.int64)
            alphas = np.zeros(m, dtype=np.int64)
            ok = True
            for j in range(m):
                tj = coords[j]
                tg = tj[np.argsort(perm)]  # t_j(gamma^{-1} x) -- not needed
                # t_{j'}(gamma x) as a function of x:
                comp = None
                for j2 in range(m):
                    v = coords[j2][perm]     # t_{j2}(gamma x)
                    # is v = alpha * t_j(x)? test on two points where t_j != 0
                    nz = np.where(tj != 0)[0]
                    a = (v[nz[0]] * pow(int(tj[nz[0]]), p - 2, p)) % p
                    if a != 0 and ((v - a * tj) % p == 0).all():
                        comp = (j2, int(a))
                        break
                if comp is None:
                    ok = False
                    break
                dirperm[j], alphas[j] = comp
            assert ok
            group.append((perm, dirperm, alphas))
    assert len(group) == (q - 1)  # (q-1)/2 dilations x 2
    return group, coords


def outer_key(subset, coeffs, p):
    """Canonical hashable key of an outer: subset (sorted tuple of global dir
    indices) + per-level coefficient tuples aligned to the sorted subset."""
    sub = tuple(sorted(subset))
    parts = [sub]
    for d in sorted(coeffs):
        vec = coeffs[d]
        parts.append((d, tuple(int(x) % p for x in vec)))
    return tuple(parts)


def act_outer(gamma, subset, coeffs, p):
    """Apply gamma to an outer state. subset: sorted tuple of global dir
    indices; coeffs: {level d: array aligned to subset}. Returns new
    (subset', coeffs') with subset' sorted and coeffs aligned."""
    perm, dirperm, alphas = gamma
    sub2 = [int(dirperm[j]) for j in subset]
    order = np.argsort(sub2)
    new_sub = tuple(int(sub2[i]) for i in order)
    out = {}
    for d, vec in coeffs.items():
        newv = np.zeros(len(vec), dtype=np.int64)
        for pos, j in enumerate(subset):
            a_inv = pow(int(alphas[j]), p - 2, p)
            val = (int(vec[pos]) * pow(a_inv, d, p)) % p
            newv[list(order).index(pos)] = val
        out[d] = newv
    return new_sub, out


def orbits(states, group, p):
    """states: list of (subset, coeffs). Returns list of
    (rep_state, transversal) where transversal is a list of group elements
    gamma with gamma . rep covering the whole orbit exactly once."""
    keyof = {}
    for i, (sub, cf) in enumerate(states):
        keyof[outer_key(sub, cf, p)] = i
    unvisited = set(keyof)
    out = []
    for k0 in list(keyof):
        if k0 not in unvisited:
            continue
        sub0, cf0 = states[keyof[k0]]
        orb = {}
        for gamma in group:
            s2, c2 = act_outer(gamma, sub0, cf0, p)
            k2 = outer_key(s2, c2, p)
            if k2 not in orb:
                orb[k2] = gamma
        for k2 in orb:
            if k2 in unvisited:
                unvisited.remove(k2)
            else:
                if k2 != k0 and k2 not in keyof:
                    raise RuntimeError("orbit left the state space")
        out.append(((sub0, cf0), list(orb.values())))
    return out
