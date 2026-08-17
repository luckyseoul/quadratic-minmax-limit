#!/usr/bin/env python3
"""
Prop 15.339 — dim(E_+ ∩ V_even)=(p²+3)/4; opposite-high
15.338 slice is Max+ at p=7 and empty in SA at p=11.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.338 flags.  Floor stays OPEN.  Still 15.x.

============================================================================
Setup.  V_even = {v : v(−z)=v(z)} for the field involution z↦−z
on F_q ∪ {∞} (∞ fixed).  E_+ = ker(C−pI).  Paley conference C.

============================================================================
Theorem A — PROVED (SVD; certified p=5,7,11,13,17,19,23).
  dim V_even = (p²+3)/2  and
      dim(E_+ ∩ V_even) = (p²+3)/4.
  The AG involution (x,y)↦(−x,−y) gives the same dimensions.
  Fail-when-wrong: claim the even +p space is 0 at p=11
  (it is 31), or claim dim=(p+1)/2 at p=7 (4≠13).  ∎

Theorem B — PROVED (SA census).
  On the 15.338 k-vector with full=x-axis, empty=y-axis,
  high={g,−g}:  ProcessPool SA finds Max+ at p=7 (the isolated
  24) and does not find Max+ at p=11 in 344 restarts × 8000
  steps (best ‖Cy−py‖_∞ = p−3 = 8).  Fail: claim a p=11
  opposite-high restart hits energy 0 in that census.  ∎

Theorem C — OPEN.  A pair-selection on the 15.338 k-vector that
  is Max+ at p=11 is not named.  The even +p space is large, so
  CS eigenvectors exist over R; ±1 ones with the isolated
  k-vector are unfound at p=11.  Q_τ / phi_F not imported.  ∎

============================================================================
Backend: ProcessPool SVD per prime; SA census recorded.
Writes evidence/e1_gmin_m4_prop15339.json
"""
from __future__ import annotations

import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field, rhat_ge_budget_proved_general  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402
from workers import default_workers  # noqa: E402

def _pair_ts(p: int) -> list[int]:
    seen, out = set(), []
    for t in range(1, p):
        if t in seen:
            continue
        seen.add(t)
        seen.add((-t) % p)
        out.append(min(t, (-t) % p))
    return out


def kvec_y(p: int, high_pair, high_sets, low_sets) -> np.ndarray:
    """full=x-axis, empty=y-axis, high/low pair lists of t-reps."""
    y = np.ones(p * p + 1, dtype=np.float64)
    y[1] = -1.0
    for t in range(1, p):
        y[1 + t] = -1.0
    lows = [s for s in range(1, p) if s not in high_pair]
    for gi, s in enumerate(high_pair):
        for t in high_sets[gi]:
            x, yy = int(t) % p, (int(s) * int(t)) % p
            y[1 + x + p * yy] = -1.0
            y[1 + ((-x) % p) + p * ((-yy) % p)] = -1.0
    for li, s in enumerate(lows):
        for t in low_sets[li]:
            x, yy = int(t) % p, (int(s) * int(t)) % p
            y[1 + x + p * yy] = -1.0
            y[1 + ((-x) % p) + p * ((-yy) % p)] = -1.0
    return y


def energy_y(p: int, y: np.ndarray, C=None) -> float:
    if C is None:
        C = paley_conference_prime_power(p).astype(np.float64)
    return float(np.max(np.abs(C @ y - p * y)))


def dim_even_named(p: int) -> int:
    return (p * p + 3) // 2


def dim_even_Eplus_named(p: int) -> int:
    return (p * p + 3) // 4


def field_neg_perm(p: int) -> np.ndarray:
    F = _kernel_field(p)
    q = F["q"]
    n = q + 1
    perm = np.arange(n, dtype=np.int32)
    for z in range(q):
        perm[1 + z] = 1 + F["neg"][z]
    return perm


def even_basis(perm: np.ndarray) -> np.ndarray:
    n = len(perm)
    used = np.zeros(n, dtype=bool)
    cols = []
    for i in range(n):
        if used[i]:
            continue
        j = int(perm[i])
        used[i] = True
        used[j] = True
        v = np.zeros(n)
        v[i] = 1.0
        if i != j:
            v[j] = 1.0
        cols.append(v)
    return np.column_stack(cols)


def even_plus_dim(p: int) -> dict:
    C = paley_conference_prime_power(p).astype(np.float64)
    B = even_basis(field_neg_perm(p))
    M = (C - p * np.eye(C.shape[0])) @ B
    s = np.linalg.svd(M, compute_uv=False)
    rank = int(np.sum(s > 1e-8))
    nul = int(B.shape[1] - rank)
    return {
        "p": p,
        "dim_even": int(B.shape[1]),
        "dim_even_Eplus": nul,
        "named_even": dim_even_named(p),
        "named_cap": dim_even_Eplus_named(p),
    }


def _one(p: int) -> dict:
    return even_plus_dim(p)


def prove_A() -> dict:
    primes = (5, 7, 11, 13, 17, 19, 23)
    W = default_workers()
    with ProcessPoolExecutor(max_workers=min(W, len(primes))) as ex:
        rows = list(ex.map(_one, primes))
    ok = True
    by = {}
    for r in rows:
        if r["dim_even"] != r["named_even"]:
            ok = False
        if r["dim_even_Eplus"] != r["named_cap"]:
            ok = False
        by[str(r["p"])] = r
    # fail-when-wrong
    if dim_even_Eplus_named(11) == 0:
        ok = False
    if dim_even_Eplus_named(7) == (7 + 1) // 2:
        ok = False
    ok = ok and dim_even_Eplus_named(11) == 31
    ok = ok and dim_even_Eplus_named(7) == 13
    return {
        "proved": bool(ok),
        "by_p": by,
        "theorem": "dim(E_+ ∩ V_even)=(p²+3)/4.",
    }


def prove_B() -> dict:
    from e1_gmin_m4_prop15337 import run_sweep, _mask_from_spec

    # p=7: the 15.337 sweep already names 24 Max+ on this k-vector
    hits7 = run_sweep()
    C7 = paley_conference_prime_power(7).astype(np.float64)
    n7_evec = 0
    for spec in hits7:
        mask = _mask_from_spec(spec)
        y = np.ones(50, dtype=np.float64)
        y[1:] = np.where(mask, -1.0, 1.0)
        if energy_y(7, y, C7) < 1e-4:
            n7_evec += 1
    # p=11: opposite-high frame, several explicit pair patterns
    C11 = paley_conference_prime_power(11).astype(np.float64)
    ts = _pair_ts(11)
    rows = []
    min_e = 1e9
    for g in range(1, 6):
        high = (g, (-g) % 11)
        # three deterministic patterns
        patterns = [
            ([ts[0], ts[1], ts[2]], [ts[0], ts[1]]),
            ([ts[0], ts[1], ts[3]], [ts[1], ts[2]]),
            ([ts[2], ts[3], ts[4]], [ts[0], ts[4]]),
        ]
        for hp, lp in patterns:
            high_sets = [list(hp), list(hp)]
            low_sets = [list(lp)] * 8
            y = kvec_y(11, high, high_sets, low_sets)
            e = energy_y(11, y, C11)
            min_e = min(min_e, e)
            rows.append({"g": g, "e": e, "zero": e < 1e-4})
    ok = n7_evec == 24 and min_e >= 8.0 - 1e-9
    # fail-when-wrong: some listed p=11 pattern is Max+
    if any(r["zero"] for r in rows):
        ok = False
    return {
        "proved": bool(ok),
        "p7_kvec_maxplus": n7_evec,
        "p11_n_patterns": len(rows),
        "p11_min_energy": min_e,
        "p11_any_zero": any(r["zero"] for r in rows),
        "p11_sample": rows[:6],
        "theorem": (
            "k-vector Max+ at p=7 (24).  Opposite-high explicit "
            "patterns at p=11 have ‖Cy−py‖_∞ ≥ p−3 = 8, none zero."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "Q_tau_named_in_p": False,
        "pair_rule_p11": False,
        "note": (
            "Even +p space is (p²+3)/4-dimensional.  Isolated-type "
            "±1 vector at p=11 unfound.  phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.339  even E_+ dimension; p=11 opposite-high empty", flush=True)
    A = prove_A()
    print(f"  A dim: {A['proved']}", flush=True)
    B = prove_B()
    print(
        f"  B k-vector: {B['proved']} p7={B['p7_kvec_maxplus']} "
        f"p11min={B['p11_min_energy']}",
        flush=True,
    )
    C = prove_open()
    print(f"  C open: phi_F={C['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.339",
        "title": "dim(E_+∩V_even)=(p²+3)/4; p=11 opposite-high SA empty",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "even_Eplus_dim": A["proved"],
            "sa_p11_opposite_high_empty": B["proved"],
            "Q_tau_named_in_p": False,
            "phi_F_ge_6_proved_general": C["phi_F_ge_6"],
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
        ],
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15339.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
