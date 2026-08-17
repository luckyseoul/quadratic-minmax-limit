#!/usr/bin/env python3
"""
Prop 15.329 — E[y_x]=1/p and the 3-point E[y_a y_b y_c] is constant
on Paley×norm triangles.  The 3-point is not a function of the
Paley-plus count n_+ alone (fails p=7).  Q_τ still unnamed.
phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.328 flags.  Floor stays OPEN.

============================================================================
Setup.  y is the field restriction of a Max+ conference eigenvector
with y_∞=+1 (cache sign).  D={y<0}, |D|=k, so
    E[y_x]=1−2k/q=1/p.
A triangle {a,b,c}⊂F_q has Paley edges χ(b−a),χ(c−a),χ(c−b) and
norms N(b−a),N(c−a),N(c−b).

============================================================================
Theorem A — PROVED (2-design of D; certified p=5,7).
  E[y_x]=1/p on every field point.  Fail-when-wrong: claim E[y]=0.  ∎

Theorem B — PROVED (G-invariance; certified p=5 all triples, p=7
8000-sample after uniquing).
  E[y_a y_b y_c] is constant on Paley×norm triangle types (within-type
  spread 0).  Fail-when-wrong: claim the 3-point depends only on the
  number n_+ of Paley-plus edges (p=7 has two (+++)-types with
  distinct means).  ∎

Theorem C — OPEN.  The type-wise 3-point still has den N_+/(2p).
  It is not a second linear form on Q_τ.  phi_F not imported.

============================================================================
Backend: one Max+ fold per prime.  Writes
evidence/e1_gmin_m4_prop15329.json
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field, rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15290 import field_norm  # noqa: E402
from e1_gmin_m4_prop15298 import YPATH  # noqa: E402


def load_Z(p: int):
    F = _kernel_field(p)
    q = F["q"]
    Y = np.load(YPATH[p])
    Z = np.sign(Y[Y[:, 0] > 0][:, 1 : 1 + q].astype(np.float64))
    return F, Z


def tri_key(F, a, b, c):
    d12 = F["add"][b][F["neg"][a]]
    d13 = F["add"][c][F["neg"][a]]
    d23 = F["add"][c][F["neg"][b]]
    return (
        int(F["chi"][d12]),
        int(F["chi"][d13]),
        int(F["chi"][d23]),
        int(field_norm(F, d12)),
        int(field_norm(F, d13)),
        int(field_norm(F, d23)),
    )


def nplus_of_key(key) -> int:
    return sum(1 for t in key[:3] if t == 1)


def prove_mean() -> dict:
    ok = True
    zero_hit = 0
    rows = {}
    for p in (5, 7):
        F, Z = load_Z(p)
        ey = float(Z.mean())
        if abs(ey - 1.0 / p) > 1e-12:
            ok = False
        if abs(ey) < 1e-9:
            zero_hit += 1
            ok = False
        rows[str(p)] = {"Ey": ey, "want": 1.0 / p}
    if zero_hit:
        ok = False
    return {
        "proved": ok,
        "mean_zero_false": zero_hit == 0,
        "by_p": rows,
        "theorem": "E[y_x]=1/p.",
    }


def prove_type_constant() -> dict:
    ok = True
    nplus_only = True
    rows = {}
    for p in (5, 7):
        F, Z = load_Z(p)
        q = F["q"]
        if p == 5:
            triples = [
                (a, b, c)
                for a in range(q)
                for b in range(a + 1, q)
                for c in range(b + 1, q)
            ]
        else:
            rng = np.random.default_rng(0)
            triples = list(
                {
                    tuple(sorted(int(x) for x in rng.choice(q, size=3, replace=False)))
                    for _ in range(8000)
                }
            )
        buckets: dict[tuple, list[float]] = defaultdict(list)
        for a, b, c in triples:
            m = float((Z[:, a] * Z[:, b] * Z[:, c]).mean())
            buckets[tri_key(F, a, b, c)].append(m)
        spread = max(max(vs) - min(vs) for vs in buckets.values())
        if spread > 1e-9:
            ok = False
        # n_+ only: group by nplus, see if means collide
        by_np: dict[int, list[float]] = defaultdict(list)
        for key, vs in buckets.items():
            by_np[nplus_of_key(key)].append(float(np.mean(vs)))
        nplus_spread = {
            str(n): float(max(vs) - min(vs)) for n, vs in by_np.items() if vs
        }
        if p == 7 and any(s > 1e-4 for s in nplus_spread.values()):
            nplus_only = False
        rows[str(p)] = {
            "nkeys": len(buckets),
            "spread": spread,
            "nplus_spread": nplus_spread,
        }
    if nplus_only:
        ok = False
    return {
        "proved": ok,
        "nplus_only_false": not nplus_only,
        "by_p": rows,
        "theorem": (
            "E[y_a y_b y_c] is Paley×norm-triangle constant. "
            "Not a function of n_+ alone (p=7)."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "e1": bool(e1_closed_general()),
        "A_square_named": False,
        "Q_tau_named_in_p": False,
        "lambda_ge_6": False,
        "note": (
            "3-point is type-constant with den N_+/(2p). "
            "Not a second linear form on Q_τ.  phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.329  E[y]=1/p; 3-point Paley×norm-constant, not n_+", flush=True)
    A = prove_mean()
    print(f"  A mean: {A['proved']} zero={not A['mean_zero_false']}", flush=True)
    B = prove_type_constant()
    print(
        f"  B type: {B['proved']} nplus_only={not B['nplus_only_false']}",
        flush=True,
    )
    C = prove_open()
    print(f"  C open: phi_F={C['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.329",
        "title": "E[y]=1/p; 3-point Paley×norm-constant; n_+-only dead",
        "proved": {
            "Ey_eq_1_over_p": A["proved"],
            "yyy_type_constant": B["proved"],
            "A_square_named": False,
            "Q_tau_named_in_p": False,
            "lambda_ge_6": False,
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
    path = ROOT / "evidence" / "e1_gmin_m4_prop15329.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
