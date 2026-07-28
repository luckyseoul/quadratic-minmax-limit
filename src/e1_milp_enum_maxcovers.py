#!/usr/bin/env python3
"""
Enumerate Max-cover perfect matchings at p=5 via MILP + no-good cuts.
Verify each with maxR, MITM Phi, clique-flip, op-norm.
Writes evidence/e1_milp_enum_maxcovers.json
"""
from __future__ import annotations

import json
import os
import sys
import time
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_k] = "1"

from minmax_quadratic import (  # noqa: E402
    maxR_matching_level,
    paley_conference_prime_power,
    phi_mitm,
)


def has_clique_flip(C, M, Maxf, p=5) -> bool:
    pairs = list(M)
    for y in Maxf:
        chis = []
        s = 0.0
        for i, j in pairs:
            ch = C[i, j] * y[i] * y[j]
            chis.append((ch, i, j))
            s += ch
        if abs(s - 1) > 1e-9:
            continue
        pos = [t for t, ch in enumerate(chis) if ch[0] > 0]
        neg = [t for t, ch in enumerate(chis) if ch[0] < 0]
        if len(pos) < 4 or len(neg) < 1:
            continue
        for Psel in combinations(pos, 4):
            for n1 in neg:
                edges = list(Psel) + [n1]
                for bits in range(32):
                    F = []
                    for t, ei in enumerate(edges):
                        _ch, i, j = chis[ei]
                        F.append(i if ((bits >> t) & 1) == 0 else j)
                    ok = True
                    for a in range(5):
                        for b in range(a + 1, 5):
                            if C[F[a], F[b]] * y[F[a]] * y[F[b]] < 0.5:
                                ok = False
                                break
                        if not ok:
                            break
                    if ok:
                        return True
    return False


def main():
    p = 5
    C = paley_conference_prime_power(p).astype(np.float64)
    n = C.shape[0]
    m = n // 2
    Phi = 0.5 * n * p
    Max = np.load("/tmp/maxplus_p5.npy").astype(np.float64)
    Maxf = np.vstack([Max, -Max])
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    E = len(edges)

    Chi = np.zeros((len(Maxf), E))
    for yi, y in enumerate(Maxf):
        for k, (i, j) in enumerate(edges):
            Chi[yi, k] = C[i, j] * y[i] * y[j]

    A_eq = np.zeros((n, E))
    for k, (i, j) in enumerate(edges):
        A_eq[i, k] = 1.0
        A_eq[j, k] = 1.0

    base_cons = [
        LinearConstraint(A_eq, np.ones(n), np.ones(n)),
        LinearConstraint(Chi, np.ones(len(Maxf)), np.full(len(Maxf), np.inf)),
    ]
    bounds = Bounds(0, 1)
    integrality = np.ones(E)

    found = {}
    rng = np.random.default_rng(0)
    t_start = time.time()
    max_trials = 100
    time_limit_feas = 90.0
    time_limit_each = 60.0

    # First feasibility
    res0 = milp(
        np.zeros(E),
        integrality=integrality,
        bounds=bounds,
        constraints=base_cons,
        options={"time_limit": time_limit_feas, "disp": False},
    )
    print(f"feas: success={res0.success} msg={res0.message}", flush=True)
    if not res0.success:
        out = {"error": "no Max-cover found", "message": str(res0.message)}
        (ROOT / "evidence" / "e1_milp_enum_maxcovers.json").write_text(
            json.dumps(out, indent=2)
        )
        return

    def add_solution(x, trial_tag):
        x = np.rint(x).astype(int)
        M = [edges[k] for k in range(E) if x[k] == 1]
        assert len(M) == m
        key = tuple(sorted(tuple(sorted(e)) for e in M))
        if key in found:
            print(f"{trial_tag}: rediscover", flush=True)
            return False
        # quick minS check
        s = np.zeros(len(Maxf))
        for i, j in M:
            s += C[i, j] * Maxf[:, i] * Maxf[:, j]
        minS = float(s.min())
        assert minS >= 1 - 1e-9, minS
        print(f"{trial_tag}: NEW #{len(found)+1} minS={minS} verifying...", flush=True)
        mr = int(maxR_matching_level(C, M, p))
        A = C.copy()
        for i, j in M:
            A[i, j] *= -1
            A[j, i] *= -1
        ph = float(phi_mitm(A))
        op = float(np.linalg.norm(A, ord=2))
        flip = has_clique_flip(C, M, Maxf, p)
        rec = {
            "matching": [list(e) for e in M],
            "minS": minS,
            "maxR": mr,
            "phi": ph,
            "op": op,
            "op_sqrt41": bool(abs(op - np.sqrt(41)) < 1e-8),
            "clique_flip": bool(flip),
            "undercut": bool(ph < Phi - 0.5),
            "criterion": bool(mr >= 60),
        }
        found[key] = rec
        print(
            f"  maxR={mr} phi={ph} op={op:.4f} flip={flip} undercut={rec['undercut']}",
            flush=True,
        )
        return True

    add_solution(res0.x, "seed0")

    for trial in range(1, max_trials + 1):
        if time.time() - t_start > 3600:
            print("time budget exceeded", flush=True)
            break
        cons = list(base_cons)
        for key in found:
            row = np.zeros(E)
            eset = set(key)
            for k, e in enumerate(edges):
                if e in eset:
                    row[k] = 1.0
            cons.append(LinearConstraint(row, -np.inf, m - 1))

        c = rng.normal(size=E)
        res = milp(
            c,
            integrality=integrality,
            bounds=bounds,
            constraints=cons,
            options={"time_limit": time_limit_each, "disp": False},
        )
        if res.success:
            add_solution(res.x, f"trial{trial}")
            continue

        # Check pure feasibility with no-goods
        res2 = milp(
            np.zeros(E),
            integrality=integrality,
            bounds=bounds,
            constraints=cons,
            options={"time_limit": time_limit_feas, "disp": False},
        )
        print(
            f"trial{trial}: opt fail; pure feas={res2.success} ({res2.message})",
            flush=True,
        )
        if res2.success:
            add_solution(res2.x, f"trial{trial}-feas")
        else:
            print(
                f"ENUMERATION COMPLETE at {len(found)} unique Max-covers "
                f"(no further feasible cover under no-goods)",
                flush=True,
            )
            break

    covers = list(found.values())
    out = {
        "n_unique": len(covers),
        "all_no_undercut": all(not r["undercut"] for r in covers) if covers else True,
        "all_criterion": all(r["criterion"] for r in covers) if covers else True,
        "all_clique_flip": all(r["clique_flip"] for r in covers) if covers else True,
        "all_op_sqrt41": all(r["op_sqrt41"] for r in covers) if covers else True,
        "all_phi_eq_Phi": all(abs(r["phi"] - Phi) < 1e-9 for r in covers)
        if covers
        else True,
        "enumeration_complete": False,  # set True only if no-goods made infeasible
        "covers": covers,
        "status": (
            "OPEN — MILP no-good enumeration of Max-covers; "
            "not a mathematical forall proof. Existence of lim alpha_n remains OPEN."
        ),
        "elapsed_sec": time.time() - t_start,
    }
    # Detect completion: last message pattern — if we broke on infeas
    # Heuristic: if trials stopped early due to infeas, mark complete
    path = ROOT / "evidence" / "e1_milp_enum_maxcovers.json"
    path.write_text(json.dumps(out, indent=2))
    print("wrote", path, "n=", len(covers), flush=True)
    print(
        json.dumps({k: out[k] for k in out if k != "covers"}, indent=2),
        flush=True,
    )


if __name__ == "__main__":
    main()
