#!/usr/bin/env python3
"""
Prop 15.525 — Type-index even-character Gram is named and does
not have live δ as an eigenvector; D=N(a+bi) has no catalog
p-law for (a,b).  Does not force ⟨δ,ψ⟩≤2.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.524 flags.

============================================================================
Setup.  T=□ ⊂ F_q^×, q=p².  Types τ are 15.290 Paley×norm classes
on T\\{±1}, merged as in 15.298 B (++, N*, --N1).  Even k∉{0,χ},
J_τ(ψ_k)=∑_{r∈τ} ψ_k(r).  Type-index Gram
    G_τσ = ∑_k J_τ(ψ_k) conj(J_σ(ψ_k)).
δ_τ=4−Q_τ/q².  D=|H_+|/(2p) (13, 409).  13=2²+3², 409=3²+20².

============================================================================
Theorem A — PROVED (character orthogonality; certified p=5..19).
  ∑_{k even} ψ_k(r) conj(ψ_k(s)) = ((q−1)/2) 1_{r=±s}.
  k=0 and k=χ each contribute 1 on T, so
      G_τσ = ((q−1)/2)(|τ∩σ|+|τ∩(−σ)|) − 2|τ||σ|.
  Fail: drop the minus intersection (p=5 named/live max-err > 1).  ∎

Theorem B — PROVED (A + live type values).
  Merged G at p=5 is [[72,−48],[−48,64]] on (++ , N*).
  Live δ ∝ (4,20) has Gδ=(−672,1088), wedge 17792≠0.
  At p=7, (92,260,140) on (++, --N1, N*) is not an eigenvector
  (second singular value of [δ | Gδ] > 0).  Fail: claim the
  p=5 wedge vanishes.  ∎

Theorem C — PROVED as a negative (named p-law catalog).
  No pair (f,g) in the catalog satisfies f(p)²+g(p)²=D(p) at
  both p=5 and p=7.  (p−1)/2 is the unique catalog law through
  (2,3).  No catalog law through (3,20).  Fail: claim
  D=((p−1)/2)²+((p+1)/2)² (that is dim V_+=(p²+1)/2; 25≠409).  ∎

Theorem D — PROVED (Gauss/Jacobi magnitudes).
  |G(ψ)|²=q and |J(χ,ψ)|²=q on nontrivial even ψ of F_q.
  Neither equals D (13≠25, 409≠49).  Fail: D=q at p=7.  ∎

Theorem E — OPEN.  Q_τ unnamed.  Pairing ≤2 not forced.
  phi_F not imported.

============================================================================
Backend: field character sums (ProcessPool) + Fraction catalog.
Writes evidence/e1_gmin_m4_prop15525.json
"""
from __future__ import annotations

import json
import os
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import _kernel_field, _psi_vec
from e1_gmin_m4_prop15298 import coarse_Q_class
from e1_gmin_m4_prop15301 import jacobi_chi_psi
from e1_gmin_m4_prop15522 import D_named

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "evidence" / "e1_gmin_m4_prop15525.json"

CERT = (5, 7, 11, 13, 17, 19)
W = 86

LIVE_DELTA_NUM = {
    5: {"++": 4, "N*": 20},
    7: {"++": 92, "--N1": 260, "N*": 140},
}


def even_ks_ex(q: int) -> list[int]:
    chi_k = (q - 1) // 2
    return [k for k in range(2, q - 1, 2) if k != chi_k]


def merged_buckets(F) -> dict[str, list[int]]:
    buckets: dict[str, list[int]] = defaultdict(list)
    for r in F["squares"]:
        r = int(r)
        if r in (1, int(F["neg1"])):
            continue
        c = coarse_Q_class(F, r)
        if c.startswith("mixed"):
            c = "++"
        buckets[c].append(r)
    return dict(buckets)


def named_G(F, buckets: dict[str, list[int]], drop_minus: bool = False):
    q = F["q"]
    keys = sorted(buckets)
    sets = {k: set(buckets[k]) for k in keys}
    n = {k: len(buckets[k]) for k in keys}
    G = np.zeros((len(keys), len(keys)), dtype=np.int64)
    for i, a in enumerate(keys):
        na = {int(F["neg"][x]) for x in sets[a]}
        for j, b in enumerate(keys):
            inter = len(sets[a] & sets[b])
            inter_m = 0 if drop_minus else len(na & sets[b])
            G[i, j] = (q - 1) // 2 * (inter + inter_m) - 2 * n[a] * n[b]
    return keys, G, n


def live_G(F, buckets: dict[str, list[int]]):
    q = F["q"]
    ks = even_ks_ex(q)
    keys = sorted(buckets)
    J = np.zeros((len(keys), len(ks)), dtype=np.complex128)
    for j, k in enumerate(ks):
        psi = _psi_vec(F, k)
        for i, key in enumerate(keys):
            J[i, j] = sum(psi[r] for r in buckets[key])
    return J @ J.conj().T


def _gram_job(p: int) -> dict:
    F = _kernel_field(p)
    buckets = merged_buckets(F)
    keys, Gn, n = named_G(F, buckets)
    Gl = live_G(F, buckets)
    err = float(np.max(np.abs(Gl - Gn)))
    _, Gdrop, _ = named_G(F, buckets, drop_minus=True)
    drop_err = float(np.max(np.abs(Gl - Gdrop)))
    return {
        "p": p,
        "keys": keys,
        "n": n,
        "G": Gn.tolist(),
        "named_live_err": err,
        "drop_minus_err": drop_err,
    }


def dim_Vplus(p: int) -> int:
    return (p * p + 1) // 2


def D_as_half_sumsq(p: int) -> int:
    """Fail-when-wrong: ((p-1)/2)²+((p+1)/2)² = dim V_+."""
    a = (p - 1) // 2
    b = (p + 1) // 2
    return a * a + b * b


def plaws(p: int) -> dict[str, int]:
    q = p * p
    return {
        "0": 0,
        "1": 1,
        "2": 2,
        "3": 3,
        "4": 4,
        "p": p,
        "p-1": p - 1,
        "p+1": p + 1,
        "p-2": p - 2,
        "p+2": p + 2,
        "p-3": p - 3,
        "p+3": p + 3,
        "2p": 2 * p,
        "4p": 4 * p,
        "(p-1)/2": (p - 1) // 2,
        "(p+1)/2": (p + 1) // 2,
        "(p-3)/2": (p - 3) // 2,
        "(p+3)/2": (p + 3) // 2,
        "2(p-2)": 2 * (p - 2),
        "4(p-2)": 4 * (p - 2),
        "(p2+1)/2": (q + 1) // 2,
        "(p2-1)/2": (q - 1) // 2,
        "(p2-1)/4": (q - 1) // 4,
        "(p2-1)/8": (q - 1) // 8,
        "(p2-9)/2": (q - 9) // 2,
        "(p2-5)/4": (q - 5) // 4,
        "n_nstar": (p - 1) * (p - 3) // 2,
        "k": p * (p - 1) // 2,
        "m": (p + 1) // 2,
        "p2-p": q - p,
        "3p-1": 3 * p - 1,
        "3p+1": 3 * p + 1,
    }


def catalog_hits_D() -> list[tuple[str, str]]:
    L5, L7 = plaws(5), plaws(7)
    names = sorted(set(L5) & set(L7))
    hits = []
    for fa in names:
        for gb in names:
            n5 = L5[fa] ** 2 + L5[gb] ** 2
            n7 = L7[fa] ** 2 + L7[gb] ** 2
            if n5 == D_named(5) and n7 == D_named(7):
                hits.append((fa, gb))
    return hits


def laws_through(vals5: int, vals7: int) -> list[str]:
    L5, L7 = plaws(5), plaws(7)
    return [
        n
        for n in set(L5) & set(L7)
        if abs(L5[n]) == vals5 and abs(L7[n]) == vals7
    ]


def prove_A() -> dict:
    with ProcessPoolExecutor(max_workers=W) as ex:
        rows = list(ex.map(_gram_job, CERT))
    ok = all(r["named_live_err"] < 1e-6 for r in rows)
    p5 = next(r for r in rows if r["p"] == 5)
    if p5["drop_minus_err"] < 1.0:
        ok = False
    return {
        "proved": bool(ok),
        "rows": {str(r["p"]): r for r in rows},
        "theorem": (
            "G_τσ=((q−1)/2)(|τ∩σ|+|τ∩(−σ)|)−2|τ||σ|. "
            "Fail drop minus (p=5 err>1)."
        ),
    }


def prove_B() -> dict:
    F5 = _kernel_field(5)
    k5, G5, n5 = named_G(F5, merged_buckets(F5))
    d5 = np.array([LIVE_DELTA_NUM[5][k] for k in k5], dtype=np.int64)
    gd5 = G5 @ d5
    wedge5 = int(d5[0] * gd5[1] - d5[1] * gd5[0])
    F7 = _kernel_field(7)
    k7, G7, n7 = named_G(F7, merged_buckets(F7))
    d7 = np.array([LIVE_DELTA_NUM[7][k] for k in k7], dtype=np.int64)
    gd7 = G7 @ d7
    M = np.vstack([d7.astype(float), gd7.astype(float)])
    svals = np.linalg.svd(M, compute_uv=False)
    ok = (
        k5 == ["++", "N*"]
        and G5.tolist() == [[72, -48], [-48, 64]]
        and wedge5 != 0
        and wedge5 == 17792
        and float(svals[1]) > 1.0
    )
    return {
        "proved": bool(ok),
        "p5": {
            "keys": k5,
            "n": n5,
            "G": G5.tolist(),
            "delta": d5.tolist(),
            "Gdelta": gd5.tolist(),
            "wedge": wedge5,
        },
        "p7": {
            "keys": k7,
            "n": n7,
            "G": G7.tolist(),
            "delta": d7.tolist(),
            "Gdelta": gd7.tolist(),
            "svals": [float(x) for x in svals],
        },
        "theorem": (
            "Live δ is not an eigenvector of the type-index Gram. "
            "Fail p=5 wedge=0."
        ),
    }


def prove_C() -> dict:
    hits = catalog_hits_D()
    a_laws = laws_through(2, 3)
    b_laws = laws_through(3, 20)
    wrong7 = D_as_half_sumsq(7)
    ok = (
        hits == []
        and a_laws == ["(p-1)/2"]
        and b_laws == []
        and wrong7 == 25
        and wrong7 != D_named(7)
        and D_as_half_sumsq(5) == D_named(5)
    )
    return {
        "proved": bool(ok),
        "catalog_hits": hits,
        "laws_2_3": a_laws,
        "laws_3_20": b_laws,
        "dimVplus_7": wrong7,
        "D7": D_named(7),
        "theorem": (
            "No catalog (f,g) names D=N(f+gi) at both p=5,7. "
            "Fail D=dim V_+ at p=7."
        ),
    }


def prove_D() -> dict:
    F5 = _kernel_field(5)
    F7 = _kernel_field(7)
    j5 = jacobi_chi_psi(F5, 2)
    j7 = jacobi_chi_psi(F7, 2)
    n5 = int(round(abs(j5) ** 2))
    n7 = int(round(abs(j7) ** 2))
    ok = (
        n5 == 25
        and n7 == 49
        and D_named(5) != 25
        and D_named(7) != 49
        and D_named(7) != 7
    )
    return {
        "proved": bool(ok),
        "Jchi_k2_N": {"5": n5, "7": n7},
        "q": {"5": 25, "7": 49},
        "theorem": "|J(χ,ψ)|²=q ≠ D. Fail D=q at p=7.",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "Q_tau_named_in_p": False,
        "D_named_in_p": False,
        "pairing_forced": False,
        "phi_F_imported": False,
        "note": (
            "Type-index Gram does not pin δ. Catalog p-laws miss "
            "(3,20). Q_τ unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("15.525 type-Gram is named and misses δ; D has no catalog p-law", flush=True)
    A, B, C, D, E = prove_A(), prove_B(), prove_C(), prove_D(), prove_open()
    out = {
        "prop": "15.525",
        "title": "Type-index Gram named and misses live δ; no catalog D=N(a+bi)",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "C": C,
        "D": D,
        "E": E,
        "e1_closed_general": bool(e1_closed_general()),
        "gsum_disj_lb": bool(gsum_disj_lb_proved_general()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "flags_not_flipped": [
            "phi_F_ge_6",
            "e1",
            "L",
            "Aut-Schur",
            "Gsum",
            "pairing",
        ],
        "L_status": "OPEN",
    }
    EV.parent.mkdir(parents=True, exist_ok=True)
    EV.write_text(json.dumps(out, indent=2) + "\n")
    print(
        f"A={A['proved']} B={B['proved']} C={C['proved']} D={D['proved']} "
        f"phi_F={out['phi_F_ge_6']}",
        flush=True,
    )
    return out


if __name__ == "__main__":
    main()
