#!/usr/bin/env python3
"""
Prop 15.347 — Field norms N_{q/p} of the six (resp. three) differences
are a p=7-stable type law: n_4(2+2) and n_3(1+2) are constant on
(K_{2,2}, sorted norms, λ-orbit) at p=5 and p=7.  χ_p alone still
splits at p=7 (15.346 C).  E[|ĥ(α)|⁴]=p⁴/2 for α≠0 on a square
class.  E[n_i² n_j²] still unnamed in p.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.346 flags.  Floor stays OPEN.

============================================================================
Setup.  Two parallel square F_p-lines in F_q.  Differences live in
F_q; N(x)=x^{p+1}∈F_p.  Paley χ_q=χ_p∘N, so N strictly refines χ.
ĥ(α)=∑_t n_t ω^{αt} is the cyclic DFT of the intersection vector.

============================================================================
Theorem A — PROVED (certified p=5,7).
  n_4 of a 2+2 set and n_3 of a 1+2 set are constant on the type
      (K_{2,2} matching data, sorted N of sides and transversals,
       {λ,−λ} of the F_p-parameter ratio).
  p=5: 9+6 types, 0 splits.  p=7: 36+12 types, 0 splits.
  Fail-when-wrong: drop N (use only χ_p / 15.346 C) and claim
  n_4 is constant at p=7.  ∎

Theorem B — CERTIFIED p=5,7 (cache DFT).  Not claimed for general p.
  For α∈F_p^×, E[|ĥ(α)|⁴] = p⁴/2 = q²/2 at p=5,7
  (625/2 and 2401/2).  Fail: claim that value equals
  (E[|ĥ|²])²=p⁴/4 at p=5 (that would be constant |ĥ|²).  ∎

Theorem C — OPEN.  The type law names the *support* of the 4-point,
  not the values.  n_4/|H_+| still has den |H_+|/(2p) at p=7
  (e.g. 154/5726=11/409).  Deg≤3/deg≤2 rationals and
  {m2,enn,m1}-polynomials miss 735/26 and 101585/818.
  phi_F not imported.

============================================================================
Backend: field_norm + cache n_4; numpy FFT of n-vectors.
Writes evidence/e1_gmin_m4_prop15347.json
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict
from fractions import Fraction
from itertools import combinations
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15290 import field_norm  # noqa: E402
from e1_gmin_m4_prop15307 import load_D, square_classes  # noqa: E402
from e1_gmin_m4_prop15321 import Qpp_floor_ub  # noqa: E402
from e1_gmin_m4_prop15326 import Q_lo_named  # noqa: E402
from e1_gmin_m4_prop15330 import LIVE_QPP  # noqa: E402
from e1_gmin_m4_prop15344 import HPLUS  # noqa: E402
from e1_gmin_m4_prop15346 import (  # noqa: E402
    live_cross,
    two_two_types,
)


def n_on(D, idxs):
    return int(D[:, list(idxs)].min(axis=1).sum())


def diff(F, x, y):
    return F["add"][int(y)][F["neg"][int(x)]]


def chi_q(F, x, y):
    d = diff(F, x, y)
    return 0 if d == 0 else int(F["chi"][d])


def Nq(F, x, y):
    d = diff(F, x, y)
    return 0 if d == 0 else int(field_norm(F, d))


def k22(F, a, b, c, d):
    e = [chi_q(F, a, c), chi_q(F, a, d), chi_q(F, b, c), chi_q(F, b, d)]
    nm = sum(1 for v in e if v == -1)
    if nm != 2:
        return nm, ""
    if (e[0] == -1 and e[3] == -1) or (e[1] == -1 and e[2] == -1):
        return 2, "match"
    return 2, "path"


def norm_22_types(p: int) -> dict:
    D, F = load_D(p)
    lines = square_classes(F)[0]
    L0 = [int(x) for x in lines[0]]
    L1 = [int(x) for x in lines[1]]
    by = defaultdict(list)
    for i, j in combinations(range(p), 2):
        for k, m in combinations(range(p), 2):
            a, b, c, d = L0[i], L0[j], L1[k], L1[m]
            nm, kind = k22(F, a, b, c, d)
            sides = tuple(sorted((Nq(F, a, b), Nq(F, c, d))))
            trans = tuple(
                sorted((Nq(F, a, c), Nq(F, a, d), Nq(F, b, c), Nq(F, b, d)))
            )
            den = (m - k) % p
            lam = ((j - i) * pow(den, p - 2, p)) % p if den else 0
            lam_key = tuple(sorted({lam, (-lam) % p}))
            by[(nm, kind, sides, trans, lam_key)].append(n_on(D, (a, b, c, d)))
    n_const = 0
    types = {}
    for key, vs in by.items():
        u = sorted(set(vs))
        const = len(u) == 1
        n_const += int(const)
        types[str(key)] = {"unique": u, "n": len(vs), "const": const}
    return {
        "n_types": len(types),
        "n_const": n_const,
        "n_split": len(types) - n_const,
        "types": types,
    }


def norm_12_types(p: int) -> dict:
    D, F = load_D(p)
    lines = square_classes(F)[0]
    L0 = [int(x) for x in lines[0]]
    L1 = [int(x) for x in lines[1]]
    by = defaultdict(list)
    for i in range(p):
        for k, m in combinations(range(p), 2):
            a, c, d = L0[i], L1[k], L1[m]
            chs = tuple(sorted((chi_q(F, a, c), chi_q(F, a, d), chi_q(F, c, d))))
            ns = tuple(sorted((Nq(F, a, c), Nq(F, a, d), Nq(F, c, d))))
            den = (m - k) % p
            sig = ((i - k) * pow(den, p - 2, p)) % p if den else 0
            sig_key = tuple(sorted({sig, (1 - sig) % p}))
            by[(chs, ns, sig_key)].append(n_on(D, (a, c, d)))
    n_const = 0
    types = {}
    for key, vs in by.items():
        u = sorted(set(vs))
        const = len(u) == 1
        n_const += int(const)
        types[str(key)] = {"unique": u, "n": len(vs), "const": const}
    return {
        "n_types": len(types),
        "n_const": n_const,
        "n_split": len(types) - n_const,
        "types": types,
    }


def hat_nz4(p: int) -> float:
    D, F = load_D(p)
    acc = n = 0
    for lines in square_classes(F):
        ns = np.stack([D[:, ix].sum(axis=1) for ix in lines], axis=1).astype(
            np.float64
        )
        hats = np.fft.fft(ns, axis=1)
        acc += float((np.abs(hats[:, 1:]) ** 4).sum())
        n += ns.shape[0] * (p - 1)
    return acc / n


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        t22 = norm_22_types(p)
        t12 = norm_12_types(p)
        if t22["n_split"] != 0 or t12["n_split"] != 0:
            ok = False
        rows[str(p)] = {
            "22": {k: t22[k] for k in ("n_types", "n_const", "n_split")},
            "12": {k: t12[k] for k in ("n_types", "n_const", "n_split")},
        }
    chi_only = two_two_types(7)
    if chi_only["n_split"] == 0:
        ok = False
    rows["p7_chi_only_splits"] = chi_only["n_split"]
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "N_{q/p} types make 2+2 and 1+2 n_t constant at p=5,7. "
            "χ_p alone splits at p=7."
        ),
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        live = hat_nz4(p)
        named = (p ** 4) / 2.0
        if abs(live - named) > 1e-6:
            ok = False
        if abs(live - (p ** 4) / 4.0) < 1e-4:
            ok = False
        rows[str(p)] = {"live": live, "named_p4_over_2": named}
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": "E[|ĥ(α)|⁴]=p⁴/2 at p=5,7 (not p⁴/4). Not general-p.",
    }


def prove_open() -> dict:
    c5 = live_cross(5)
    c7 = live_cross(7)
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "Q_tau_named_in_p": False,
        "interval_closed": False,
        "cross": {"5": str(c5["E_n2n2"]), "7": str(c7["E_n2n2"])},
        "Qlo5": str(Q_lo_named(5)),
        "Qub5": str(Qpp_floor_ub(5)),
        "live_Qpp": {str(p): str(LIVE_QPP[p]) for p in LIVE_QPP},
        "note": (
            "Norm types classify the 4-point support.  Values of "
            "n_4/|H_+| and E[n_i² n_j²] still carry den |H_+|/p. "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.347  N_{q/p} type law; hat 4th moment; cross unnamed", flush=True)
    A = prove_A()
    print(
        f"  A norms: {A['proved']} p7_22={A['by_p']['7']['22']} "
        f"chi_splits={A['by_p']['p7_chi_only_splits']}",
        flush=True,
    )
    B = prove_B()
    print(f"  B hat4: {B['proved']} p5={B['by_p']['5']['live']}", flush=True)
    C = prove_open()
    print(f"  C open: phi_F={C['phi_F_ge_6']} cross={C['cross']}", flush=True)
    out = {
        "prop": "15.347",
        "title": "N_{q/p} kills p=7 2+2/1+2 splits; E[|ĥ|⁴]=p⁴/2; Q_τ open",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "norm_types_constant": A["proved"],
            "hat_nz4_is_p4_over_2": B["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15347.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
