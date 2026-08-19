"""Numba flip solver: per-candidate propagation + iterative DFS."""
import numpy as np
from numba import njit


@njit(cache=True)
def flip_batch(p, k, q, SIG, FV, Tm, LIDX, thi, tlo, ybuf, ycount, node_cap):
    """SIG: (B,k,p) std sigma tables (+p at class p-1 slots).
    FV: (B,k) required flip counts. LIDX: (k,p,p) point ids per line.
    Appends solutions (int8 q-vectors) into ybuf, count in ycount[0]."""
    B = SIG.shape[0]
    two_p = 2 * p
    maxbits = k * p
    bj = np.zeros(maxbits, np.int64)
    bs = np.zeros(maxbits, np.int64)
    state = np.zeros(maxbits, np.int8)      # -1 free, 0/1 fixed
    A = np.zeros(q, np.int64)
    g = np.zeros(q, np.int64)
    fix1 = np.zeros(q, np.int64)            # sum of fixed-1 bits covering x
    free = np.zeros(q, np.int64)            # free bits covering x
    ones_j = np.zeros(k, np.int64)
    free_j = np.zeros(k, np.int64)
    stack_bit = np.zeros(maxbits + 1, np.int64)
    stack_val = np.zeros(maxbits + 1, np.int64)
    for r in range(B):
        # build A, g
        for x in range(q):
            s0 = 0
            for j in range(k):
                s0 += SIG[r, j, Tm[j, x]]
            A[x] = s0
        okc = True
        for x in range(q):
            d = A[x] - thi
            if d % two_p != 0:
                okc = False
                break
            g[x] = d // two_p
        if not okc:
            continue
        # bits
        nb = 0
        for j in range(k):
            for s1 in range(p):
                if SIG[r, j, s1] == p:
                    bj[nb] = j
                    bs[nb] = s1
                    nb += 1
        # coverage / init
        for x in range(q):
            fix1[x] = 0
            free[x] = 0
        for b in range(nb):
            j = bj[b]
            s1 = bs[b]
            for t in range(p):
                free[LIDX[j, s1, t]] += 1
        feas = True
        for x in range(q):
            if g[x] > free[x] or g[x] + 1 < 0:
                feas = False
                break
        if not feas:
            continue
        for b in range(nb):
            state[b] = -1
        for j in range(k):
            ones_j[j] = 0
            free_j[j] = 0
        for b in range(nb):
            free_j[bj[b]] += 1
        wok = True
        for j in range(k):
            if FV[r, j] > free_j[j]:
                wok = False
                break
        if not wok:
            continue
        # propagation
        changed = True
        dead = False
        while changed and not dead:
            changed = False
            for b in range(nb):
                if state[b] != -1:
                    continue
                j = bj[b]
                s1 = bs[b]
                can0 = True
                can1 = True
                # try val: check the p points of this line
                for t in range(p):
                    x = LIDX[j, s1, t]
                    # val 0: fix1 same, free-1
                    if fix1[x] > g[x] + 1 or fix1[x] + free[x] - 1 < g[x]:
                        can0 = False
                    if fix1[x] + 1 > g[x] + 1 or fix1[x] + 1 + free[x] - 1 < g[x]:
                        can1 = False
                    if not can0 and not can1:
                        break
                # weight feasibility
                if can1 and ones_j[j] + 1 > FV[r, j]:
                    can1 = False
                if can0 and ones_j[j] + (free_j[j] - 1) < FV[r, j]:
                    can0 = False
                if not can0 and not can1:
                    dead = True
                    break
                if can0 != can1:
                    val = 1 if can1 else 0
                    state[b] = val
                    ones_j[j] += val
                    free_j[j] -= 1
                    for t in range(p):
                        x = LIDX[j, s1, t]
                        fix1[x] += val
                        free[x] -= 1
                    changed = True
        if dead:
            continue
        # DFS over remaining free bits
        order = np.empty(nb, np.int64)
        no = 0
        for b in range(nb):
            if state[b] == -1:
                order[no] = b
                no += 1
        nodes = 0
        sp = 0
        stack_bit[0] = 0
        stack_val[0] = 0        # next value to try at this depth: 0 then 1
        # iterative DFS: depth = sp; at each depth try val 0 then 1
        depth = 0
        tryval = 0
        while True:
            nodes += 1
            if nodes > node_cap:
                raise RuntimeError("flip node cap")
            if depth == no:
                # leaf: verify and emit
                good = True
                for x in range(q):
                    xi = fix1[x]
                    if xi != g[x] and xi != g[x] + 1:
                        good = False
                        break
                if good:
                    m = ycount[0]
                    if m < ybuf.shape[0]:
                        for x in range(q):
                            ybuf[m, x] = 1 if fix1[x] == g[x] else -1
                    ycount[0] = m + 1
                # backtrack
                depth -= 1
                if depth < 0:
                    break
                b = order[depth]
                j = bj[b]
                s1 = bs[b]
                val = stack_val[depth]
                ones_j[j] -= val
                free_j[j] += 1
                for t in range(p):
                    x = LIDX[j, s1, t]
                    fix1[x] -= val
                    free[x] += 1
                state[b] = -1
                tryval = val + 1
                continue
            b = order[depth]
            j = bj[b]
            s1 = bs[b]
            advanced = False
            while tryval <= 1:
                val = tryval
                okv = True
                if val == 1 and ones_j[j] + 1 > FV[r, j]:
                    okv = False
                if val == 0 and ones_j[j] + (free_j[j] - 1) < FV[r, j]:
                    okv = False
                if okv:
                    for t in range(p):
                        x = LIDX[j, s1, t]
                        nf1 = fix1[x] + val
                        nfr = free[x] - 1
                        if nf1 > g[x] + 1 or nf1 + nfr < g[x]:
                            okv = False
                            break
                if okv:
                    # apply
                    state[b] = val
                    ones_j[j] += val
                    free_j[j] -= 1
                    for t in range(p):
                        x = LIDX[j, s1, t]
                        fix1[x] += val
                        free[x] -= 1
                    stack_val[depth] = val
                    depth += 1
                    tryval = 0
                    advanced = True
                    break
                tryval += 1
            if advanced:
                continue
            # exhausted values at this depth: backtrack
            depth -= 1
            if depth < 0:
                break
            b = order[depth]
            j = bj[b]
            s1 = bs[b]
            val = stack_val[depth]
            ones_j[j] -= val
            free_j[j] += 1
            for t in range(p):
                x = LIDX[j, s1, t]
                fix1[x] -= val
                free[x] += 1
            state[b] = -1
            tryval = val + 1
