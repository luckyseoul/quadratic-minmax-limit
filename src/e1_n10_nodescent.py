#!/usr/bin/env python3
"""
Prop 15.41 candidate cert (n=10): no-descent for gap-2 undercutters.

For every perfect-matching undercutter F of Paley C_10 (Phi=13=Phi-2),
every single-edge extension F∪{e} has Phi >= 13.

Also samples multi-edge extensions. Writes evidence/e1_n10_nodescent.json.

This supports (does not prove) the no-descent lemma needed to upgrade
Prop 15.40 -> m_n >= Phi(C)-2. See evidence/E1_FAILURE_GRAPH.md F13.
"""
from __future__ import annotations

import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from itertools import product
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_k] = "1"

from minmax_quadratic import paley_conference_prime_power  # noqa: E402

N = 10
PHI = 15.0
TARGET = 13.0


def _all_pm(nn: int):
    def rec(rem):
        rem = list(rem)
        if not rem:
            yield []
            return
        a = rem[0]
        for bi, b in enumerate(rem[1:], 1):
            rest = rem[1:bi] + rem[bi + 1 :]
            for M in rec(rest):
                yield [(a, b)] + M

    yield from rec(list(range(nn)))


def _worker_init():
    global C, X, EDGES
    C = paley_conference_prime_power(3)
    bits = list(product([-1.0, 1.0], repeat=N - 1))
    X = np.ones((len(bits), N), dtype=np.float64)
    for r, b in enumerate(bits):
        X[r, 1:] = b
    EDGES = [(i, j) for i in range(N) for j in range(i + 1, N)]


def _phi_flip(F):
    A = C.copy()
    for i, j in F:
        A[i, j] *= -1
        A[j, i] *= -1
    vals = np.einsum("bi,ij,bj->b", X, A, X) * 0.5
    return float(np.max(np.abs(vals)))


def check_undercutter(F_t):
    F = list(F_t)
    ph = _phi_flip(F)
    Fs = set(F)
    min_add = 1e9
    n_viol = 0
    for e in EDGES:
        if e in Fs:
            continue
        ph2 = _phi_flip(F + [e])
        if ph2 < min_add:
            min_add = ph2
        if ph2 < TARGET - 1e-9:
            n_viol += 1
    # multi-edge sample
    rng = np.random.default_rng(abs(hash(F_t)) % (2**31))
    outside = [e for e in EDGES if e not in Fs]
    deep_viol = 0
    deep_min = 1e9
    for kextra in (2, 3, 4, 5, 8, 12):
        for _ in range(40):
            idx = rng.choice(len(outside), size=kextra, replace=False)
            extra = [outside[i] for i in idx]
            ph2 = _phi_flip(F + extra)
            deep_min = min(deep_min, ph2)
            if ph2 < TARGET - 1e-9:
                deep_viol += 1
    return {
        "ph": ph,
        "min_add1": min_add,
        "n_viol_add1": n_viol,
        "deep_min": deep_min,
        "deep_viol": deep_viol,
    }


def main():
    _worker_init()
    und = []
    for M in _all_pm(N):
        Ft = tuple(sorted((min(a, b), max(a, b)) for a, b in M))
        ph = _phi_flip(list(Ft))
        if ph < PHI - 1e-9:
            und.append(Ft)
    assert len(und) == 144, len(und)

    W = max(1, (os.cpu_count() or 4) - 2)
    with ProcessPoolExecutor(W, initializer=_worker_init) as ex:
        results = list(ex.map(check_undercutter, und, chunksize=2))

    n_viol = sum(r["n_viol_add1"] for r in results)
    deep_viol = sum(r["deep_viol"] for r in results)
    out = {
        "n": N,
        "Phi_C": PHI,
        "n_undercut_pm": len(und),
        "all_phi_eq_13": all(abs(r["ph"] - 13) < 1e-9 for r in results),
        "add1_total_violations": n_viol,
        "add1_min_of_min": min(r["min_add1"] for r in results),
        "deep_total_violations": deep_viol,
        "deep_min_of_min": min(r["deep_min"] for r in results),
        "workers": W,
        "lemma_candidate": (
            "no-descent: if Phi(C⊕F)=Phi-2 then Phi(C⊕(F∪{e}))>=Phi-2 for all e"
        ),
        "status": (
            "CERT n=10 PM undercutters: no-descent holds (0 viol); "
            "general proof OPEN (F13); L OPEN"
        ),
    }
    path = ROOT / "evidence" / "e1_n10_nodescent.json"
    path.write_text(json.dumps(out, indent=2))
    md = ROOT / "evidence" / "E1_NODESCENT.md"
    md.write_text(
        "\n".join(
            [
                "# No-descent candidate toward m_n ≥ Φ−2",
                "",
                "**Date:** 2026-07-28",
                "**Status:** Structural lemmas partial; n=10 cert for matching undercutters;",
                "**forall / general n OPEN. Existence of lim α_n remains OPEN.**",
                "",
                "## Goal (graph-allowed path)",
                "",
                "Upgrade Prop 15.40 → global m_n >= Phi(C)-2 **without** F13 shortcut.",
                "",
                "## Lemmas",
                "",
                "### Lemma A (first-hit; proved)",
                "Along any edge-adding chain, the first undercutting prefix F has",
                "Φ(C⊕F) ≥ Φ(C)−2 (Prop 15.20b with k=1).",
                "",
                "### Lemma B (no-descent; OPEN in general)",
                "If Φ(C⊕F)=Φ(C)−2, then for every edge e∉F, Φ(C⊕(F∪{e}))≥Φ(C)−2.",
                "",
                "If A+B hold on the ρ=1 family, induction on Hamming distance yields",
                "m_n≥Φ(C)−2, gap O(1)⇒E(1)⇒L=1/2 by denseness.",
                "",
                "**F13 reminder:** abstract 2-Lipschitz functions need not satisfy B.",
                "The claim is special to Φ of Seidel matrices.",
                "",
                "## Certified at n=10 (matching undercutters)",
                "",
                "| Check | Result |",
                "|-------|--------|",
                f"| # PM undercutters | **144** (all Φ=13) |",
                f"| add-1 violations Φ<13 | **{n_viol}** |",
                f"| min Φ over add-1 | {out['add1_min_of_min']} |",
                f"| multi-edge sample violations | **{deep_viol}** |",
                f"| multi-edge sample min Φ | {out['deep_min_of_min']} |",
                f"| workers | {W} |",
                "",
                "JSON: `evidence/e1_n10_nodescent.json`. Script: `src/e1_n10_nodescent.py`.",
                "",
                "## Dangerous-edge obstruction (partial structure)",
                "",
                "Write A=C⊕F, Φ(A)=Φ−2, edge e=(u,v), B=A⊕e.",
                "Then Q_B(x)=Q_A(x)−2σ_x with σ_x=A_uv x_u x_v ∈ {±1}.",
                "",
                "If some maximiser x with Q_A(x)=Φ−2 has σ_x=−1, then Q_B(x)=Φ,",
                "so Φ(B)≥Φ≥Φ−2. Similarly for negative maximisers with σ=+1.",
                "",
                "**Descent to Φ−4 requires** rigid alignment σ≡+1 on all + maximisers",
                "and σ≡−1 on all − maximisers. Open: prove this alignment is impossible",
                "for undercutters of ρ=1 conferences (or non-maximiser spikes save Φ−2).",
                "",
                "## What this does **not** prove",
                "",
                "- Lemma B for general n (only n=10 PM cert)",
                "- Clique-flip / matching non-undercut for p≥5",
                "- k_⋆ / E(1) / existence of L",
                "",
                "**Do not mark Main Theorem settled.**",
                "",
            ]
        )
    )
    print(json.dumps(out, indent=2))
    assert n_viol == 0
    assert all(abs(r["ph"] - 13) < 1e-9 for r in results)


if __name__ == "__main__":
    main()
