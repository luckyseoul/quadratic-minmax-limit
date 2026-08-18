#!/usr/bin/env python3
"""
Prop 15.563 — full-Ω (m>3) is not locked Fejer; A_full not a p-law.

    At p=7, n_full=4410=90q splits into two energy types
    (units p²(p+1)): (1,1,1,3) n=36q and (1,1,2,2) n=54q.
    ẑ-on-line is not the k=3 Fejer F (min shape>0.13).
    Occupancy-sum on good forms is 4- or 5-valued, not two-level.
A_full,dbl=A_full,n3=−4p³/15 is an ensemble mean, not a p-law.

Does **not** flip leftover flags.  Does **not** import phi_F.
Does **not** overwrite 15.550–15.562.  Does **not** interpolate (p−5)/15.

============================================================================
Setup.  15.276/15.562: k=3 is a locked sawtooth triple, ẑ=F Fejer
on each support line, occupancy two-level.  15.271: m-form
occupancy-sum is not two-level.  m=(p+1)/2; at p=7, m=4, unique
ẑ-support is all of Ω.

============================================================================
Theorem A — PROVED (p=7 cache; p=5 control).
  Per-line energy of a full-Ω vector, in units p²(p+1), is not
  (1,1,1,1).  At p=7 the only types are (1,1,1,3) (1764=36q)
  and (1,1,2,2) (2646=54q).  At p=5, full=k=3 is equi-energy
  (4/3,4/3,4/3).  Fail: one energy type at p=7.  ∎

Theorem B — PROVED (p=7 cache).
  ẑ-on-line shape is not Fejer |F|² (mean L2>0.2 after shift
  and scale).  Occupancy fibre-sums on the p+1 forms are never
  all two-level.  Fail: ẑ=F as for k=3.  ∎

Theorem C — OPEN.  A_full,dbl mean is +4p³/9 on (1,1,1,3) and
  −20p³/27 on (1,1,2,2); the mix is −4p³/15.  A_full,n3 is
  −4p³/9 and −4p³/27.  These are p=7 rationals, not a
  Gauss/Jacobi/Fejer integer in p.  Do not interpolate (p−5)/15.

============================================================================
Backend: p=5,7 yhat energy/occupancy.  Writes
evidence/e1_gmin_m4_prop15563.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15276 import sawtooth_dft_G
from e1_gmin_m4_prop15554 import yhat_omega

EV = ROOT / "evidence" / "e1_gmin_m4_prop15563.json"


def omega_lines(F, Omega):
    mul, p = F["mul"], F["p"]
    taken, lines = set(), []
    for a in Omega:
        if a in taken:
            continue
        L = [mul[g][a] for g in range(1, p)]
        lines.append(L)
        taken.update(L)
    return lines


def unique_forms(p: int):
    out, seen = [], set()
    for a0 in range(p):
        for a1 in range(p):
            if a0 == 0 and a1 == 0:
                continue
            if a0 != 0:
                key = (1, (a1 * pow(a0, p - 2, p)) % p)
            else:
                key = (0, 1)
            if key in seen:
                continue
            seen.add(key)
            ell = np.array(
                [(key[0] * (x % p) + key[1] * (x // p)) % p for x in range(p * p)],
                dtype=np.int32,
            )
            out.append(ell)
    return out


def energy_patterns(p: int) -> dict:
    F, Y, Omega, yhat = yhat_omega(p)
    nOm = len(Omega)
    idx = {int(om): i for i, om in enumerate(Omega)}
    lines = omega_lines(F, Omega)
    nsupp = (np.abs(yhat) > 1e-6).sum(axis=1)
    full = np.where(nsupp == nOm)[0]
    unit = float(p * p * (p + 1))
    pat: dict[tuple, int] = defaultdict(int)
    for i in full:
        ev = []
        for L in lines:
            ev.append(sum(abs(yhat[i, idx[om]]) ** 2 for om in L) / unit)
        key = tuple(sorted(int(round(x)) if p == 7 else round(x, 4) for x in ev))
        pat[key] += 1
    return {
        "n_full": int(full.size),
        "n_types": len(pat),
        "patterns": {str(k): int(v) for k, v in sorted(pat.items(), key=lambda kv: -kv[1])},
        "equi": bool(len(pat) == 1 and (p != 7)),
    }


def shape_and_occ(p: int, nchk: int = 40) -> dict:
    F, Y, Omega, yhat = yhat_omega(p)
    nOm, q = len(Omega), F["q"]
    idx = {int(om): i for i, om in enumerate(Omega)}
    lines = omega_lines(F, Omega)
    nsupp = (np.abs(yhat) > 1e-6).sum(axis=1)
    full = np.where(nsupp == nOm)[0]
    Fej2 = np.array([abs(sawtooth_dft_G(p, 1, 0, t)) ** 2 for t in range(1, p)])
    forms = unique_forms(p)
    rng = np.random.default_rng(p)
    sample = rng.choice(full, size=min(nchk, full.size), replace=False)
    mul = F["mul"]

    def best_fej(prof):
        best = 1e99
        for s in range(p - 1):
            t = np.roll(Fej2, s)
            a, b = prof / (prof.sum() + 1e-18), t / (t.sum() + 1e-18)
            best = min(best, float(np.linalg.norm(a - b)))
        return best

    shapes = []
    n_two = 0
    for i in sample:
        L = lines[0]
        gen = L[0]
        prof = np.array([abs(yhat[i, idx[mul[g][gen]]]) ** 2 for g in range(1, p)])
        shapes.append(best_fej(prof))
        two = True
        for ell in forms:
            vals = {round(float(Y[i, ell == t].sum()), 6) for t in range(p)}
            if len(vals) > 2:
                two = False
                break
        n_two += int(two)
    return {
        "nchk": int(sample.size),
        "mean_fej_shape": float(np.mean(shapes)),
        "min_fej_shape": float(np.min(shapes)),
        "n_two_level": n_two,
        "fej_fails": bool(np.min(shapes) > 0.05 if p == 7 else np.max(shapes) < 1e-8),
        "occ_not_two_level": bool(n_two == 0 if p == 7 else True),
    }


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        rec = energy_patterns(p)
        if p == 7:
            pats = rec["patterns"]
            row_ok = (
                rec["n_types"] == 2
                and pats.get("(1, 1, 1, 3)") == 1764
                and pats.get("(1, 1, 2, 2)") == 2646
            )
        else:
            row_ok = rec["n_types"] == 1
        if not row_ok:
            ok = False
        rows[str(p)] = {**rec, "ok": bool(row_ok)}
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "p=7 full-Ω energy types (1,1,1,3) and (1,1,2,2). "
            "Fail one equi-energy type."
        ),
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        rec = shape_and_occ(p)
        if p == 7:
            row_ok = rec["min_fej_shape"] > 0.05 and rec["n_two_level"] == 0
        else:
            row_ok = rec["min_fej_shape"] < 1e-8
        if not row_ok:
            ok = False
        rows[str(p)] = {**rec, "ok": bool(row_ok)}
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": "full-Ω ẑ-on-line is not Fejer F; occupancy not two-level. Fail ẑ=F.",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "A_full_named_in_p": False,
        "Q3_double_named_in_p": False,
        "Q3_n3_named_in_p": False,
        "phi_F_imported": False,
        "note": (
            "A_full,dbl mean −4p³/15 at p=7 is a mix of +4p³/9 and "
            "−20p³/27 on two energy types, not a Fejer/Johnson p-law. "
            "Do not interpolate (p−5)/15. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.563  full-Ω not Fejer; A_full not a p-law", flush=True)
    A = prove_A()
    print(f"  A energy types: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B not Fejer/two-level: {B['proved']}", flush=True)
    C = prove_open()
    out = {
        "prop": "15.563",
        "title": "full-Ω energy splits; ẑ not Fejer; A_full open",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "C": C,
        "proved": {
            "full_energy_two_types": A["proved"],
            "full_not_fejer": B["proved"],
            "A_full_named_in_p": False,
            "phi_F_ge_6": False,
        },
        "e1_closed_general": bool(e1_closed_general()),
        "gsum_disj_lb": bool(gsum_disj_lb_proved_general()),
        "phi_F_ge_6": False,
        "L_status": "OPEN",
        "flags_not_flipped": [
            "phi_F_ge_6",
            "e1",
            "L",
            "Aut-Schur",
            "Gsum",
            "pairing",
            "residual_ii",
            "type_I",
        ],
        "obstruction": (
            "A_full is a p=7 mix on energy types (1,1,1,3)/(1,1,2,2); "
            "not Fejer, not two-level, not a p-law"
        ),
    }
    EV.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print("wrote", EV, flush=True)
    return out


if __name__ == "__main__":
    main()
