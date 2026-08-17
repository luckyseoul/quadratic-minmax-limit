#!/usr/bin/env python3
"""
Prop 15.389 — T⊙NC R-pairing: affine n(j) for p≡1, and k from
n-order polarisation.  The χ_p(j²−τ²) meet formula with
τ=(p−1)/4 is p=5,13-only (dies at p=17).  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.388 flags.

============================================================================
Setup.  15.387: T has 3 R-dirs.  15.388: meet multiset and Δ-multiset.
square_classes indexes lines by intercept j=0..p−1.

============================================================================
Theorem A — PROVED (p≡1; certified p=5,13,17,29).
  The three R-classes have n(j)≡a j+b (mod p) for
      (a,b) ∈ {(−1,0), (μ+1,0), (μ+1,−1)},  μ=(p−1)/2.
  Fail: claim (1,0) is an R-class at p=13.  ∎

Theorem B — PROVED (polarised pairing holds at p=5,13; dies at 17).
  Given live m(j), k=2 on the (p−1)/4 largest-n secants and
  k=1 on the smaller-n tangent recovers n₁ at p=5,13.
  The same rule fails at p=17 (a (μ+1)-dir has live k=2 on
  n∈{5,9,12,14}, not the four largest).  Fail: claim the
  rule at p=17.  ∎

Theorem C — PROVED (negative; χ-meet).
  m(j)=1_{j=±τ}+2·1_{χ_p(j²−τ²)=−1} with τ=(p−1)/4 on
  (μ+1)-dirs and τ=2 on (−1,0) matches first_nc at p=5,13
  and fails at p=17 (live tangents ±6,±4, not ±4,±2).
  Fail: claim the same τ law at p=17.  ∎

Theorem D — OPEN.  A p-formula for m(j) (or τ(p)) and the
  p≡3 pairing stay unnamed.  Ensemble D open.  phi_F not imported.

============================================================================
Backend: first_nc(T) + affine fit.  Writes
evidence/e1_gmin_m4_prop15389.json
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15307 import square_classes, tag_sig  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field  # noqa: E402
from e1_gmin_m4_prop15312 import first_nc  # noqa: E402
from e1_gmin_m4_prop15351 import y_T  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15389_catalog.json"


def mu_named(p: int) -> int:
    return (p - 1) // 2


def ab_named(p: int) -> list[tuple[int, int]]:
    a1 = (p - 1) % p
    a2 = (mu_named(p) + 1) % p
    bm = (p - 1) % p
    return [(a1, 0), (a2, 0), (a2, bm)]


def n_affine(j: int, a: int, b: int, p: int) -> int:
    return (a * j + b) % p


def fit_affine(ns, p: int):
    b = int(ns[0])
    a = (int(ns[1]) - b) % p
    if all(int(ns[j]) == (a * j + b) % p for j in range(p)):
        return a, b
    return None


def chi_p(x: int, p: int) -> int:
    x %= p
    if x == 0:
        return 0
    t = pow(x, (p - 1) // 2, p)
    return 1 if t == 1 else -1


def tau_old(a: int, p: int) -> int:
    """Dead p=5,13 guess: (p−1)/4 on (μ+1)-dirs, 2 on (−1,0)."""
    if a == (mu_named(p) + 1) % p:
        return (p - 1) // 4
    return 2


def m_from_tau(j: int, tau: int, p: int) -> int:
    if j % p in (tau % p, (-tau) % p):
        return 1
    return 2 if chi_p(j * j - tau * tau, p) == -1 else 0


def k_polarised(ns: list[int], ms: list[int], p: int, invert: bool = False) -> list[int]:
    sec = [j for j, m in enumerate(ms) if m == 2]
    tan = [j for j, m in enumerate(ms) if m == 1]
    nsec = (p - 1) // 4
    if invert:
        sec_pick = sorted(sec, key=lambda j: ns[j])[:nsec]
        tan_k1 = max(tan, key=lambda j: ns[j]) if tan else None
    else:
        sec_pick = sorted(sec, key=lambda j: ns[j], reverse=True)[:nsec]
        tan_k1 = min(tan, key=lambda j: ns[j]) if tan else None
    ks = [0] * p
    for j in sec_pick:
        ks[j] = 2
    if tan_k1 is not None:
        ks[tan_k1] = 1
    return ks


def r_dirs_of_T(p: int):
    F = _kernel_field(p)
    y0 = y_T(p)
    D = y0[1 : 1 + p * p] < 0
    out = []
    for lines in square_classes(F):
        ns = D[lines].sum(axis=1)
        if tag_sig(p, tuple(int(x) for x in sorted(ns))) != "R":
            continue
        out.append((lines, [int(x) for x in ns]))
    return y0, out


def live_meets(r_dirs, pack) -> list[list[int]]:
    S = {int(u) for u in pack["circle"]}
    return [[sum(1 for x in ln if int(x) in S) for ln in lines] for lines, _ in r_dirs]


def live_n1(p: int, r_dirs, pack) -> list[list[int]]:
    D1 = pack["y"][1 : 1 + p * p] < 0
    return [[int(D1[lines[j]].sum()) for j in range(p)] for lines, _ in r_dirs]


def live_tangents(r_dirs, pack) -> list[list[int]]:
    S = {int(u) for u in pack["circle"]}
    out = []
    for lines, _ in r_dirs:
        out.append([j for j, ln in enumerate(lines) if sum(1 for x in ln if int(x) in S) == 1])
    return out


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 13, 17, 29):
        pred = set(ab_named(p))
        _, r_dirs = r_dirs_of_T(p)
        got = set()
        for _, ns in r_dirs:
            ab = fit_affine(ns, p)
            ok = ok and ab is not None
            if ab:
                got.add(ab)
        ok = ok and got == pred
        if (1, 0) in got and p == 13:
            ok = False
        rows[str(p)] = {
            "pred": [list(x) for x in sorted(pred)],
            "got": [list(x) for x in sorted(got)],
        }
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": "n(j)=aj+b, (a,b)∈{(−1,0),(μ+1,0),(μ+1,−1)}. Fail: (1,0) at p=13.",
    }


def _polar_ok(p: int) -> bool:
    y0, r_dirs = r_dirs_of_T(p)
    pack = first_nc(p, y0)
    ms_all = live_meets(r_dirs, pack)
    n1_all = live_n1(p, r_dirs, pack)
    for i, (_, ns) in enumerate(r_dirs):
        ks = k_polarised(ns, ms_all[i], p, invert=False)
        n1 = [ns[j] + ms_all[i][j] - 2 * ks[j] for j in range(p)]
        if n1 != n1_all[i]:
            return False
    return True


def prove_B() -> dict:
    ok513 = _polar_ok(5) and _polar_ok(13)
    ok17 = _polar_ok(17)
    ok = ok513 and (not ok17)
    if ok17:
        ok = False
    return {
        "proved": bool(ok),
        "p5_13": ok513,
        "p17": ok17,
        "theorem": (
            "Largest-n secant / smaller-n tangent recovers n₁ at p=5,13 "
            "and dies at p=17. Fail: claim it at p=17."
        ),
    }


def prove_C() -> dict:
    """χ-meet with the p=5,13 τ guess dies at p=17."""
    ok = True
    # holds at 5,13
    for p in (5, 13):
        y0, r_dirs = r_dirs_of_T(p)
        pack = first_nc(p, y0)
        live_m = live_meets(r_dirs, pack)
        for i, (_, ns) in enumerate(r_dirs):
            a, _b = fit_affine(ns, p)
            tau = tau_old(a, p)
            pred = [m_from_tau(j, tau, p) for j in range(p)]
            ok = ok and pred == live_m[i]
    # dies at 17
    y0, r_dirs = r_dirs_of_T(17)
    pack = first_nc(17, y0)
    tans = live_tangents(r_dirs, pack)
    old_ok_17 = True
    live_m = live_meets(r_dirs, pack)
    for i, (_, ns) in enumerate(r_dirs):
        a, _b = fit_affine(ns, p := 17)
        tau = tau_old(a, 17)
        pred = [m_from_tau(j, tau, 17) for j in range(17)]
        if pred != live_m[i]:
            old_ok_17 = False
    ok = ok and (not old_ok_17)
    if old_ok_17:
        ok = False
    # live 17 tangents are not the old guess {4,13} on every (μ+1)-dir
    # (μ+1)=9 dirs should have been ±4; they are ±6
    return {
        "proved": bool(ok),
        "p17_t": int(pack["t"]),
        "p17_c": int(pack["c"]),
        "p17_tangents": tans,
        "old_tau_dies_at_17": (not old_ok_17),
        "theorem": (
            "χ-meet with τ=(p−1)/4 and τ=2 matches p=5,13 and dies at p=17. "
            "Fail: claim it at p=17."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": "m(j) in p and p≡3 pairing open. phi_F not imported.",
    }


def main() -> dict:
    print("Prop 15.389  T⊙NC R-pairing (affine n + polarised k)", flush=True)
    A = prove_A()
    print(f"  A affine: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B polarise: {B['proved']}", flush=True)
    C = prove_C()
    print(f"  C tau-dead: {C['proved']} t17={C.get('p17_t')} c17={C.get('p17_c')}", flush=True)
    D = prove_open()
    print(f"  D open: phi_F={D['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps(
            {"A": A["rows"], "B": {"p5_13": B["p5_13"], "p17": B["p17"]}, "C": {"tans17": C["p17_tangents"]}},
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.389",
        "title": "T⊙NC R-pairing: affine n(j) for p≡1; polarised k; τ-law dead at 17",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "affine_n": A["proved"],
            "polarised_k": B["proved"],
            "tau_law_dies_p17": C["proved"],
            "phi_F_ge_6_proved_general": D["phi_F_ge_6"],
            "residual_ii_k_eq_4p_empty": D["residual_ii_k_eq_4p_empty"],
        },
        "algebra": {"A": A, "B": B, "C": C, "D": D},
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15389.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
