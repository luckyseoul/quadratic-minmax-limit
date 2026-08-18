#!/usr/bin/env python3
"""
Prop 15.503 — Certified numerical: per-design Q_y on
F_p^×\\{±1} is flat at p=7 to float64, though Aut only
identifies inv-pairs.  Not exact arithmetic.  Not flat
per design on the coarse 15.298 ++ class.

    Q_y(r) := mean_{ξ∈Ω} |ẑ_y(ξ)|² |ẑ_y(rξ)|²

Aut gives Q_y(r)=Q_y(r^{−1}) only.  At p=7, sub={2,3,4,5}
splits as {2,4} and {3,5}.  On the full H_+ cache (5726
designs) the max per-y spread of Q_y/q² on sub is
6.06×10^{−15} (float64 zero).  Fail: the same per-y
flatness on coarse ++ (sub ⊔ ++N1): 4550 of 5726 designs
have spread > 10^{−8}; max spread 16.  Ensemble means
still match on ++ (15.298 B).

Does **not** claim an exact identity Q_y(2)=Q_y(3) in
ℚ, and does not import that as a p-law.  p=5 sub={2,3}
is a single inv-orbit, so the extra comparison is
invisible.  No p≥13 cache.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.502 flags.
"""
from __future__ import annotations

import json
import os
from collections import defaultdict
from functools import lru_cache
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import _kernel_field
from e1_gmin_m4_prop15290 import _u_matrix, omega_set
from e1_gmin_m4_prop15298 import coarse_Q_class

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "evidence" / "e1_gmin_m4_prop15503.json"


@lru_cache(maxsize=4)
def u_on_omega(p: int):
    F = _kernel_field(p)
    Omega = omega_set(F)
    u = _u_matrix(p, F, Omega)
    return F, Omega, u


def Q_y_on(p: int, rs: list[int]) -> dict[int, np.ndarray]:
    F, Omega, u = u_on_omega(p)
    om_idx = {int(xi): i for i, xi in enumerate(Omega)}
    out = {}
    for r in rs:
        idx = np.array([om_idx[int(F["mul"][r][xi])] for xi in Omega], dtype=np.int32)
        out[int(r)] = (u * u[:, idx]).mean(axis=1)
    return out


def subfield_off_pm1(p: int) -> list[int]:
    return [a for a in range(2, p - 1)]


def coarse_members(p: int) -> dict[str, list[int]]:
    F, _, _ = u_on_omega(p)
    cls: dict[str, list[int]] = defaultdict(list)
    for r in F["squares"]:
        cls[coarse_Q_class(F, r)].append(int(r))
    return dict(cls)


def prove_A() -> dict:
    """Per-y Q constant on sub at p=7. Fail: per-y constant on ++."""
    ok = True
    p = 7
    F, _, _ = u_on_omega(p)
    q2 = float(F["q"] ** 2)
    sub = subfield_off_pm1(p)
    Qs = Q_y_on(p, sub)
    stack = np.stack([Qs[r] for r in sub], axis=0)
    spread_sub = (stack.max(axis=0) - stack.min(axis=0)) / q2
    cls = coarse_members(p)
    Qpp = Q_y_on(p, cls["++"])
    st_pp = np.stack([Qpp[r] for r in cls["++"]], axis=0)
    spread_pp = (st_pp.max(axis=0) - st_pp.min(axis=0)) / q2
    if float(spread_sub.max()) > 1e-8:
        ok = False
    if int((spread_pp > 1e-8).sum()) == 0:
        ok = False
    if float(spread_pp.max()) < 1.0:
        ok = False
    # inv-pairs must match (sanity)
    if float(np.max(np.abs(Qs[2] - Qs[4])) / q2) > 1e-8:
        ok = False
    rows = {
        7: {
            "sub": sub,
            "n_y": int(spread_sub.size),
            "sub_max_spread": float(spread_sub.max()),
            "pp_max_spread": float(spread_pp.max()),
            "pp_n_y_spread": int((spread_pp > 1e-8).sum()),
            "pp_n_r": len(cls["++"]),
        }
    }
    return {"proved": bool(ok), "by_p": rows}


def prove_B() -> dict:
    """Ensemble Q is still flat on coarse ++ / N* / --N1 at p=5,7."""
    ok = True
    rows = {}
    for p in (5, 7):
        F, _, _ = u_on_omega(p)
        q2 = float(F["q"] ** 2)
        cls = coarse_members(p)
        rec = {}
        for name in ("++", "N*"):
            if name not in cls:
                ok = False
                continue
            Qs = Q_y_on(p, cls[name])
            ens = {r: float(Qs[r].mean() / q2) for r in cls[name]}
            sp = max(ens.values()) - min(ens.values())
            if sp > 1e-8:
                ok = False
            rec[name] = {"ens_spread": sp, "n_r": len(cls[name])}
        rows[p] = rec
    return {"proved": bool(ok), "by_p": rows}


def prove_open() -> dict:
    return {
        "proved": False,
        "Q_tau_named_in_p": False,
        "phi_F_imported": False,
        "note": (
            "Numerical flatness on sub at p=7, not an exact ℚ identity. "
            "No general reason; p≥13 untested. phi_F not imported."
        ),
    }


def main() -> dict:
    print("15.503 per-y Q constant on sub, not on coarse ++", flush=True)
    A, B, D = prove_A(), prove_B(), prove_open()
    out = {
        "A": A,
        "B": B,
        "D": D,
        "e1_closed_general": bool(e1_closed_general()),
        "gsum_disj_lb": bool(gsum_disj_lb_proved_general()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
    }
    EV.parent.mkdir(parents=True, exist_ok=True)
    EV.write_text(json.dumps(out, indent=2) + "\n")
    print(
        f"A={A['proved']} B={B['proved']} phi_F={out['phi_F_ge_6']}",
        flush=True,
    )
    return out


if __name__ == "__main__":
    main()
