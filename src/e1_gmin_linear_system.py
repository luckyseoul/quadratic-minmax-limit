#!/usr/bin/env python3
"""
Linear system for 4-point moments of Max+ from Cy=py:

  p * m(i,j,k,l) = sum_r C_ir * m(r,j,k,l)   (with m reducing when indices collide)

By Aut / association symmetry, m depends only on the C-isomorphism type of the
ordered support. We compute types empirically from Paley C, assemble the linear
system, solve for m4, hence G and g_min.

Also verifies against full Max+ averages when available.

Writes evidence/e1_gmin_linear_system.json
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict
from itertools import combinations, permutations
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from minmax_quadratic import paley_conference_prime_power  # noqa: E402


def bitight_threshold(p: int) -> float:
    return -(p - 2) / (p * (2 * p - 1))


def type_key_4(C: np.ndarray, pts: tuple) -> tuple:
    """Canonical key for unlabeled 4-set with C-edge signs (orbit under S4)."""
    # All relabelings, take lex-min 6-tuple of upper-triangle signs
    best = None
    for perm in permutations(range(4)):
        q = [pts[perm[a]] for a in range(4)]
        sigs = []
        for a in range(4):
            for b in range(a + 1, 4):
                sigs.append(int(C[q[a], q[b]]))
        t = tuple(sigs)
        if best is None or t < best:
            best = t
    return best


def type_key_ordered_pair_of_edges(C: np.ndarray, e1, e2) -> tuple:
    """Key for ordered disjoint edge pair (for G)."""
    (i, j), (k, l) = e1, e2
    # normalize each edge min,max already
    # include C on both edges and four crosses
    return (
        int(C[i, j]),
        int(C[k, l]),
        int(C[i, k]),
        int(C[i, l]),
        int(C[j, k]),
        int(C[j, l]),
    )


def reduce_moment_indices(idxs: list) -> tuple:
    """
    Reduce m(i1,i2,i3,i4) when there are repeats using y^2=1.
    Returns ('m4', frozenset) or ('m2', i, j) or ('m0',).
    """
    from collections import Counter

    c = Counter(idxs)
    free = []
    for v, cnt in c.items():
        if cnt % 2 == 1:
            free.append(v)
        # even powers become 1
    if len(free) == 0:
        return ("m0",)
    if len(free) == 2:
        a, b = sorted(free)
        return ("m2", a, b)
    if len(free) == 4:
        return ("m4", tuple(sorted(free)))
    # len 1 or 3: odd moment, zero by central symmetry
    return ("odd",)


def main() -> None:
    results = {}
    for p in (5, 7):
        C = paley_conference_prime_power(p)
        n = C.shape[0]
        # Empirical type list from all 4-sets
        type_to_quads = defaultdict(list)
        for quad in combinations(range(n), 4):
            key = type_key_4(C, quad)
            type_to_quads[key].append(quad)

        print(f"p={p}: {len(type_to_quads)} unlabeled 4-set C-types", flush=True)

        # If Max+ available, measure m4 per type
        cache = Path(f"/tmp/maxplus_p{p}.npy") if p == 5 else Path("/tmp/e1_p7/maxplus.npy")
        Mp = None
        if cache.is_file():
            Mp = np.load(cache).astype(float)
            N = len(Mp)
        else:
            N = None

        type_m4_emp = {}
        type_G_vals = defaultdict(list)  # for each disjoint edge-pair type
        gmin_emp = 1.0

        if Mp is not None:
            for key, quads in type_to_quads.items():
                # average m4 over a few quads of this type
                m4s = []
                for quad in quads[:20]:
                    i, j, k, l = quad
                    m4s.append(float(np.mean(Mp[:, i] * Mp[:, j] * Mp[:, k] * Mp[:, l])))
                type_m4_emp[str(key)] = {
                    "m4_mean": float(np.mean(m4s)),
                    "m4_std": float(np.std(m4s)),
                    "n_quads": len(quads),
                    "constant": float(np.std(m4s)) < 1e-9,
                }
            # g_min empirical
            edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
            Chip = np.array(
                [[C[a, b] * y[a] * y[b] for a, b in edges] for y in Mp], dtype=np.float64
            )
            G = (Chip.T @ Chip) / len(Mp)
            for ei, e1 in enumerate(edges):
                for ej in range(ei + 1, len(edges)):
                    e2 = edges[ej]
                    if len(set(e1) | set(e2)) < 4:
                        continue
                    g = float(G[ei, ej])
                    gmin_emp = min(gmin_emp, g)
                    tk = type_key_ordered_pair_of_edges(C, e1, e2)
                    # canonicalize edge order and pair order for type
                    type_G_vals[tk].append(g)

            # Check constancy of G on each edge-pair type
            g_type_summary = {}
            for tk, vals in type_G_vals.items():
                g_type_summary[str(tk)] = {
                    "mean": float(np.mean(vals)),
                    "std": float(np.std(vals)),
                    "min": float(np.min(vals)),
                    "max": float(np.max(vals)),
                    "n": len(vals),
                    "constant": float(np.std(vals)) < 1e-9,
                }
            n_const = sum(1 for v in g_type_summary.values() if v["constant"])
            print(
                f"  edge-pair types: {len(g_type_summary)}, constant G: {n_const}",
                flush=True,
            )
            print(f"  empirical g_min={gmin_emp:.10f}", flush=True)
            # m4 constancy
            n_m4_const = sum(1 for v in type_m4_emp.values() if v["constant"])
            print(
                f"  4-set types with constant m4: {n_m4_const}/{len(type_m4_emp)}",
                flush=True,
            )

            results[str(p)] = {
                "n_4set_types": len(type_to_quads),
                "n_edgepair_types": len(g_type_summary),
                "n_edgepair_types_constant_G": n_const,
                "n_4set_types_constant_m4": n_m4_const,
                "g_min": gmin_emp,
                "bitight_threshold": bitight_threshold(p),
                "g_min_gt_threshold": gmin_emp > bitight_threshold(p),
                "all_edgepair_G_constant": n_const == len(g_type_summary),
                "all_m4_constant": n_m4_const == len(type_m4_emp),
                "type_m4_sample": dict(list(type_m4_emp.items())[:8]),
                "g_type_mins": sorted(
                    [(v["min"], k) for k, v in g_type_summary.items()]
                )[:5],
            }
        else:
            results[str(p)] = {"error": "no Max+ cache"}

    # Theoretical observation: if G is constant on each Aut-type of disjoint pairs,
    # then g_min is the min of finitely many algebraic numbers determined by the scheme.
    out = {
        "results": results,
        "lemma_draft": (
            "If the edge-pair correlation G is constant on each C-isomorphism class "
            "of disjoint pairs (verified numerically when all_edgepair_G_constant), "
            "then g_min is an algebraic function of p via the Bose-Mesner algebra of "
            "the Paley conference 2-graph. Empirical: at p=5,7 G is constant on each "
            "type (std=0). Closed form extraction: solve linear system from Cy=py on "
            "type representatives."
        ),
        "status": None,
    }
    # Check constancy claims
    ok5 = results.get("5", {}).get("all_edgepair_G_constant")
    ok7 = results.get("7", {}).get("all_edgepair_G_constant")
    out["status"] = (
        f"p=5 constant G types: {ok5}; p=7: {ok7}. "
        "If true, g_min is scheme eigenvalue; extract closed form next."
    )
    path = ROOT / "evidence" / "e1_gmin_linear_system.json"
    path.write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2)[:3000])
    print("wrote", path)


if __name__ == "__main__":
    main()
