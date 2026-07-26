#!/usr/bin/env python3
"""
n=10 structural follow-up: perfect-matching flips of Paley C_10.

Campaign result: min Hamming distance from Paley to an m_10-optimum is 5,
and the first witness is a perfect matching of K_10. This script certifies:

  1. Every (or which) perfect matching flip yields Φ=13.
  2. Degree multiset of all 5-edge undercutters (not just matchings).
  3. Spectral fingerprints of matching-optima.
"""
from __future__ import annotations

import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from minmax_quadratic import paley_conference_prime_power, phi, is_conference_matrix
from n10_structure import (
    _flip_edge,
    _edges,
    spectral_fingerprint,
    all_switchings_of,
    min_switch_hamming,
)

for _k in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(_k, "1")


def perfect_matchings(n: int = 10) -> List[List[Tuple[int, int]]]:
    """All perfect matchings of K_n (n even). Count = n!/(2^{n/2}(n/2)!)."""
    assert n % 2 == 0
    verts = list(range(n))
    out: List[List[Tuple[int, int]]] = []

    def rec(remaining: List[int], partial: List[Tuple[int, int]]) -> None:
        if not remaining:
            out.append(list(partial))
            return
        a = remaining[0]
        for i in range(1, len(remaining)):
            b = remaining[i]
            edge = (a, b) if a < b else (b, a)
            nxt = remaining[1:i] + remaining[i + 1 :]
            partial.append(edge)
            rec(nxt, partial)
            partial.pop()

    rec(verts, [])
    return out


def _eval_matching(args):
    for _k in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[_k] = "1"
    C_list, matchings = args
    C = np.array(C_list, dtype=np.float64)
    rows = []
    for m in matchings:
        A = C.copy()
        for i, j in m:
            _flip_edge(A, i, j)
        ph = phi(A)
        rows.append({"matching": m, "phi": float(ph)})
    return rows


def _five_flip_shard(args):
    for _k in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[_k] = "1"
    C_list, edges, combos = args
    C = np.array(C_list, dtype=np.float64)
    hits = []
    for combo in combos:
        A = C.copy()
        verts = []
        for t in combo:
            i, j = edges[t]
            _flip_edge(A, i, j)
            verts += [i, j]
        ph = phi(A)
        if ph <= 13.0 + 1e-9:
            from collections import Counter

            deg = sorted(Counter(verts).values(), reverse=True)
            is_matching = deg == [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
            hits.append(
                {
                    "edges": [edges[t] for t in combo],
                    "degree_multiset": deg,
                    "is_perfect_matching": is_matching,
                    "phi": float(ph),
                }
            )
    return hits


def main():
    workers = int(os.environ.get("WORKERS", max(1, (os.cpu_count() or 4) - 2)))
    out_dir = Path(os.environ.get("OUT", str(ROOT / "evidence")))
    out_dir.mkdir(parents=True, exist_ok=True)

    C = paley_conference_prime_power(3)
    assert is_conference_matrix(C)
    report: Dict = {"n": 10, "workers": workers, "t0": time.time()}

    # --- perfect matchings ---
    print("[1] enumerate perfect matchings of K_10 ...", flush=True)
    pms = perfect_matchings(10)
    report["n_perfect_matchings"] = len(pms)
    print(f"  count={len(pms)} (expect 945)", flush=True)

    print("[2] exact Φ after each perfect-matching flip of Paley ...", flush=True)
    shard = max(1, (len(pms) + workers - 1) // workers)
    jobs = []
    C_list = C.tolist()
    for s in range(0, len(pms), shard):
        jobs.append((C_list, pms[s : s + shard]))
    rows = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for part in ex.map(_eval_matching, jobs):
            rows.extend(part)
    phis = [r["phi"] for r in rows]
    from collections import Counter

    hist = Counter(phis)
    n13 = sum(1 for p in phis if abs(p - 13) < 1e-9)
    n15 = sum(1 for p in phis if abs(p - 15) < 1e-9)
    n_other = len(phis) - n13 - n15
    report["matching_flip"] = {
        "phi_histogram": {str(k): int(v) for k, v in sorted(hist.items())},
        "n_phi13": n13,
        "n_phi15": n15,
        "n_other": n_other,
        "all_matchings_give_m_or_paley": n_other == 0 and n13 + n15 == len(pms),
        "fraction_optimal": n13 / len(pms),
    }
    print(
        f"  Φ hist={dict(hist)}  n13={n13}/{len(pms)} "
        f"frac={n13/len(pms):.4f}",
        flush=True,
    )

    # fingerprint a few φ=13 matchings
    sc = all_switchings_of(C)
    fps = []
    for r in rows:
        if abs(r["phi"] - 13) < 1e-9 and len(fps) < 5:
            A = C.copy()
            for i, j in r["matching"]:
                _flip_edge(A, i, j)
            fp = spectral_fingerprint(A)
            fp["matching"] = r["matching"]
            fp["hamming_to_switch_class"] = min_switch_hamming(A, sc)
            fps.append(fp)
    report["matching_optima_fingerprints"] = fps

    # --- all 5-edge undercutters: are they all perfect matchings? ---
    print("[3] classify all 5-edge undercutters of Paley ...", flush=True)
    edges = _edges(10)
    combos = list(combinations(range(len(edges)), 5))
    shard2 = max(1, (len(combos) + workers * 4 - 1) // (workers * 4))
    jobs2 = []
    for s in range(0, len(combos), shard2):
        jobs2.append((C_list, edges, combos[s : s + shard2]))
    hits = []
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(_five_flip_shard, j) for j in jobs2]
        for fut in as_completed(futs):
            hits.extend(fut.result())
    n_match_hits = sum(1 for h in hits if h["is_perfect_matching"])
    deg_hist = Counter(tuple(h["degree_multiset"]) for h in hits)
    report["five_edge_undercutters"] = {
        "n_combos": len(combos),
        "n_phi13": len(hits),
        "n_perfect_matching": n_match_hits,
        "n_non_matching": len(hits) - n_match_hits,
        "all_are_perfect_matchings": n_match_hits == len(hits) and len(hits) > 0,
        "degree_multiset_histogram": {
            str(k): int(v) for k, v in sorted(deg_hist.items(), key=lambda x: -x[1])
        },
        "time_s": time.time() - t0,
    }
    print(
        f"  5-edge Φ=13 hits={len(hits)}  matchings={n_match_hits}  "
        f"non-matching={len(hits)-n_match_hits}  t={time.time()-t0:.1f}s",
        flush=True,
    )
    print(f"  degree hist: {report['five_edge_undercutters']['degree_multiset_histogram']}", flush=True)

    # Cross-check: n_match_hits should equal n13 from matching enumeration
    report["consistency"] = {
        "matching_enum_n13": n13,
        "five_flip_matching_hits": n_match_hits,
        "agree": n13 == n_match_hits,
    }

    report["total_time_s"] = time.time() - report.pop("t0")
    report["verdict"] = {
        "min_hamming_paley_to_m10": 5,
        "first_undercutters_are_perfect_matchings": report["five_edge_undercutters"][
            "all_are_perfect_matchings"
        ],
        "n_optimal_matching_flips": n13,
        "n_perfect_matchings_total": len(pms),
        "relative_gap": (15 - 13) / (10 ** 1.5),
        "absolute_gap": 2,
        "r_at_optima": fps[0]["r"] if fps else None,
        "note": (
            "Exact optima at Hamming-5 from Paley C_10 are precisely the "
            "perfect-matching flips with Φ=13 (subset of all 945 perfect matchings). "
            "SA finds more distant optima (Hamming 11–16 to switch class) with the "
            "same r=13/15 and often δ=320."
        ),
    }

    out = out_dir / "n10_matching_optima.json"
    out.write_text(json.dumps(report, indent=2, default=str))
    print(f"wrote {out}", flush=True)
    print(json.dumps(report["verdict"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
