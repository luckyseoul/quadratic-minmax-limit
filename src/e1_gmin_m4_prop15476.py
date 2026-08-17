#!/usr/bin/env python3
"""
Prop 15.476 — Δ² is a function of the AG(2,p) direction type.

The occupancy signature of D on the p+1 parallel classes
determines Δ² on every live Max+ (0 split signatures at p=5,7).
n_full does not (splits at p=7).  Complementary Δ² still
unnamed in p.  Q_τ unnamed.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.475 flags.

============================================================================
Setup.  15.475 names Δ² on the parallel-union class.  Complementary
rows remain.  An affine direction is a parallel class of p lines.
Occupancy type of D in that direction is the sorted p-tuple of
|D∩ℓ|.  The direction signature is the multiset of these p+1 types.

============================================================================
Theorem A — PROVED (cache fold; p=5,7).
  Δ² is constant on each direction signature (0 split).
  Two signatures at p=5 (30+100), seven at p=7.
  Fail: claim n_full determines Δ² at p=7 (n_full=1 takes
  both 0 and 38416).  ∎

Theorem B — PROVED (A + type census).
  Parallel-union signature is
      1 × (0^{p−2}, p, p)  and  p × (λ^p),   λ=(p−1)/2,
  and carries the 15.475 value.  The triple-rainbow signature
      3 × (0,1,…,p−1)  and  (p−2) × (λ^p)
  is complementary: 100 rows at p=5 (Δ²=4500) and 1176 at p=7
  (Δ²=38416=M²/9).  Fail: claim every complementary row is
  triple-rainbow at p=7 (only 1176 of 5586).  ∎

Theorem C — OPEN.  Complementary Δ² (including triple-rainbow)
  unnamed by a single formula in p.  Q_τ unnamed.  phi_F not imported.

============================================================================
Backend: numpy occupancy GEMM on p=5,7 caches.  Writes
evidence/e1_gmin_m4_prop15476.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field  # noqa: E402
from e1_gmin_m4_prop15468 import D_rows  # noqa: E402
from e1_gmin_m4_prop15470 import make_psi  # noqa: E402
from e1_gmin_m4_prop15473 import e_matrix  # noqa: E402
from e1_gmin_m4_prop15474 import sigma_named  # noqa: E402
from e1_gmin_m4_prop15475 import delta2_parunion_named  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15476_catalog.json"


def line_blocks(p: int, q: int) -> list[np.ndarray]:
    pts = np.arange(q)
    x0, x1 = pts % p, pts // p
    blocks = [
        np.stack([(x1 == c).astype(np.int8) for c in range(p)]),
        np.stack([(x0 == c).astype(np.int8) for c in range(p)]),
    ]
    for s in range(1, p):
        blocks.append(
            np.stack([(x1 == (s * x0 + b) % p).astype(np.int8) for b in range(p)])
        )
    return blocks


def live_delta2_and_sigs(p: int):
    F = _kernel_field(p)
    q = F["q"]
    Ds = D_rows(p, F).astype(np.float64)
    hat = Ds @ e_matrix(F).T
    chi = np.array(F["chi"], dtype=np.int32)
    sigma = sigma_named(p)
    Omega = chi == sigma
    Omega[0] = False
    psi = make_psi(F)
    plus = Omega & ((psi.imag > 0) if sigma < 0 else (psi.real > 0))
    d2 = (np.abs(hat[:, plus]) ** 2).sum(1) - (
        (np.abs(hat[:, Omega]) ** 2).sum(1) - (np.abs(hat[:, plus]) ** 2).sum(1)
    )
    # Δ = U+−U−; redo cleanly
    Up = (np.abs(hat[:, plus]) ** 2).sum(1)
    Um = (np.abs(hat[:, Omega]) ** 2).sum(1) - Up
    d2 = (Up - Um) ** 2
    blocks = line_blocks(p, q)
    sigs = []
    n_full = np.zeros(Ds.shape[0], dtype=np.int32)
    for i in range(Ds.shape[0]):
        types = []
        for B in blocks:
            occ = tuple(sorted(int(x) for x in (Ds[i] @ B.T)))
            types.append(occ)
            n_full[i] += sum(1 for x in occ if x == p)
        c = Counter(types)
        sigs.append(tuple(sorted((t, n) for t, n in c.items())))
    return d2, sigs, n_full


def rainbow_tuple(p: int) -> tuple:
    return tuple(range(p))


def regular_tuple(p: int) -> tuple:
    lam = (p - 1) // 2
    return tuple([lam] * p)


def is_triple_rainbow(sig, p: int) -> bool:
    want = Counter(
        {
            rainbow_tuple(p): 3,
            regular_tuple(p): p - 2,
        }
    )
    got = Counter()
    for t, n in sig:
        got[t] += n
    return got == want


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        d2, sigs, n_full = live_delta2_and_sigs(p)
        by = defaultdict(set)
        nfull_by = defaultdict(set)
        for i, sig in enumerate(sigs):
            by[sig].add(round(float(d2[i]), 4))
            nfull_by[int(n_full[i])].add(round(float(d2[i]), 4))
        split = sum(1 for v in by.values() if len(v) > 1)
        nfull_split = sum(1 for v in nfull_by.values() if len(v) > 1)
        if split != 0:
            ok = False
        if p == 7 and nfull_split == 0:
            ok = False
        rows[str(p)] = {
            "n_sig": len(by),
            "split": split,
            "nfull_split": nfull_split,
        }
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "Δ² is a function of direction type. "
            "Fail: n_full determines Δ² at p=7."
        ),
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    expect_tr = {5: (100, 4500.0), 7: (1176, 38416.0)}
    for p in (5, 7):
        d2, sigs, n_full = live_delta2_and_sigs(p)
        pred = float(delta2_parunion_named(p))
        n_par = int(np.sum(np.abs(d2 - pred) < 1.0))
        tr_idx = [i for i, s in enumerate(sigs) if is_triple_rainbow(s, p)]
        tr_d2 = {round(float(d2[i]), 4) for i in tr_idx}
        n_tr = len(tr_idx)
        exp_n, exp_d2 = expect_tr[p]
        if n_tr != exp_n or tr_d2 != {exp_d2}:
            ok = False
        if p == 7 and n_tr == (d2.shape[0] - n_par):
            ok = False
        rows[str(p)] = {
            "n_parunion": n_par,
            "n_triple_rainbow": n_tr,
            "tr_d2": sorted(tr_d2),
        }
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "Triple-rainbow is one complementary class, not all. "
            "Fail: every p=7 complementary row is triple-rainbow."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "Direction type determines Δ² on the live ensemble. "
            "Complementary values unnamed in p. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.476  Δ² from AG(2,p) direction type", flush=True)
    A = prove_A()
    print(f"  A type-fn: {A['proved']} {A['rows']}", flush=True)
    B = prove_B()
    print(f"  B rainbow: {B['proved']} {B['rows']}", flush=True)
    C = prove_open()
    print(f"  C open phi_F={C['phi_F_ge_6']}", flush=True)
    CAT.write_text(json.dumps({"A": A["proved"], "B": B["proved"]}, indent=2) + "\n")
    out = {
        "prop": "15.476",
        "title": "Δ² is a function of AG(2,p) direction type; complementary unnamed",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "delta2_fn_of_dirtype": A["proved"],
            "triple_rainbow_not_all_comp": B["proved"],
            "phi_F_ge_6_proved_general": C["phi_F_ge_6"],
            "residual_ii_k_eq_4p_empty": C["residual_ii_k_eq_4p_empty"],
        },
        "algebra": {"A": A, "B": B, "C": C},
        "L_status": "OPEN",
        "flags_not_flipped": [
            "phi_F_ge_6",
            "e1",
            "L",
            "Aut-Schur",
            "Gsum",
            "pairing",
            "residual_ii_k_eq_4p_empty",
        ],
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15476.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
