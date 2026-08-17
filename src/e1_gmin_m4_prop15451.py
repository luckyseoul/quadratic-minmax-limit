#!/usr/bin/env python3
"""
Prop 15.451 — Two 19-free constraints on c, neither a pin
for p≥7.  (1) Q++≥Q_N* on S0 ⇔ c≤c_eq(p).  (2) At p=5
every μ-half-net of the nonsquare Desarguesian r-net is
Max+ (count 130=|H+|).  Type− is not a half-net.
c still unnamed.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.450 flags.

============================================================================
Setup.  15.450: p≡1 S0 even-S floor ⇔ c∈c_window (Q_lo,Q^{ub}).
15.398: Q++≥Q_N* ⇔ Q++≥Q_eq on S0.  15.305 C: every Max+
meets each nonsquare-direction line in μ=(p−1)/2.  r=(p+1)/2.
15.331: Type− ker is a nonsquare line.

============================================================================
Theorem A — PROVED (15.398 + 15.409 lattice; p≡1).
  Q++≥Q_eq ⇔ 8A/D≥Q_eq ⇔ c≤c_eq(p), the Hoffman cut
  of u_window built with Q_eq in place of Q_lo.  Never
  divides by 19.  Combined with 15.450,
      c ∈ [c_min(p), c_eq(p)].
  At p=5 this is still {1}.  At p=13 it is 67 points
  [551,617], not a singleton.  Fail: use Q_lo as the
  Q++≥Q_N* cut (then c_max=642≠617 at p=13).  ∎

Theorem B — PROVED (full enum; p=5).
  The number of k-subsets of F_25 meeting every
  nonsquare-direction line in μ=2 points is 130=|H+|.
  Combined with 15.305 C (Max+⊆half-nets), every such
  half-net is Max+.  Fail: claim the count is n_1d=30.  ∎

Theorem C — PROVED (15.331 ker; certified p=5,7).
  A Type− lift is monochromatic on its nonsquare ker line,
  so n_ℓ∈{0,p}≠μ.  Hence Type− is not a half-net and
  not Max+ (0 of n_1d Type− lifts lie in the p=5,7 caches).
  Fail: claim Type− ⊂ H+.  ∎

Theorem D — OPEN.  c_eq does not pin c for p≥13.
  Half-net count unnamed in p for p≥7.  Q_τ unnamed.
  phi_F not imported.

============================================================================
Backend: Fraction window + p=5 combination enum.  Writes
evidence/e1_gmin_m4_prop15451.json
"""
from __future__ import annotations

import json
import os
import sys
from itertools import combinations
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import (  # noqa: E402
    _kernel_field,
    _linear_forms,
    is_type_plus_form,
    lift_sigma_vector,
)
from e1_gmin_m4_prop15292 import n_1d  # noqa: E402
from e1_gmin_m4_prop15298 import YPATH  # noqa: E402
from e1_gmin_m4_prop15321 import Qpp_floor_ub  # noqa: E402
from e1_gmin_m4_prop15326 import Q_lo_named  # noqa: E402
from e1_gmin_m4_prop15330 import HPLUS  # noqa: E402
from e1_gmin_m4_prop15367 import A_num  # noqa: E402
from e1_gmin_m4_prop15398 import Q_eq_named  # noqa: E402
from e1_gmin_m4_prop15409 import D_1d, _ceil_frac, _floor_frac  # noqa: E402
from e1_gmin_m4_prop15450 import c_window_ends, in_c_window, u_from_c  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15451_catalog.json"


def u_window_from_Q(p: int, qlo, qhi) -> tuple[int, int]:
    A = A_num(p)
    d1 = D_1d(p)
    dmin = (8 * A) / qhi
    dmax = (8 * A) / qlo
    return _ceil_frac((dmin - d1) / p), _floor_frac((dmax - d1) / p)


def c_cut_from_u(p: int, ulo: int, uhi: int) -> tuple[int, int, int]:
    step = (p - 1) // 2
    u0 = ((ulo + step - 1) // step) * step
    if u0 > uhi:
        return (0, -1, 0)
    u1 = (uhi // step) * step
    return (2 * u0 // (p - 1), 2 * u1 // (p - 1), (u1 - u0) // step + 1)


def c_eq_ends(p: int) -> tuple[int, int, int]:
    """Hoffman ∩ [Q_eq, Q^{ub}] window."""
    ulo, uhi = u_window_from_Q(p, Q_eq_named(p), Qpp_floor_ub(p))
    return c_cut_from_u(p, ulo, uhi)


def ns_parallel_classes(F, p: int) -> list[list[set[int]]]:
    used: set[int] = set()
    gens: list[list[int]] = []
    for v in range(1, F["q"]):
        if v in used or int(F["chi"][v]) != -1:
            continue
        line = [0]
        x = 0
        for _ in range(1, p):
            x = F["add"][x][v]
            line.append(x)
            used.add(x)
        gens.append(line)
    classes: list[list[set[int]]] = []
    for gen in gens:
        covered: set[int] = set()
        cl: list[set[int]] = []
        for x in range(F["q"]):
            if x in covered:
                continue
            ell = {F["add"][x][s] for s in gen}
            cl.append(ell)
            covered |= ell
        classes.append(cl)
    return classes


def count_halfnets_p5() -> int:
    p = 5
    F = _kernel_field(p)
    classes = ns_parallel_classes(F, p)
    mu = (p - 1) // 2
    k = p * (p - 1) // 2
    line_sets = [[set(ell) for ell in cl] for cl in classes]
    n_ok = 0
    for comb in combinations(range(F["q"]), k):
        D = set(comb)
        good = True
        for cl in line_sets:
            for ell in cl:
                if len(D & ell) != mu:
                    good = False
                    break
            if not good:
                break
        if good:
            n_ok += 1
    return n_ok


def type_minus_halfnet_count(p: int) -> int:
    F = _kernel_field(p)
    classes = ns_parallel_classes(F, p)
    mu = (p - 1) // 2
    m = (p + 1) // 2
    forms = [ab for ab in _linear_forms(p) if not is_type_plus_form(p, *ab)]
    plus = list(combinations(range(p), m))
    n_half = 0
    for a, b in forms:
        for A in plus:
            y = lift_sigma_vector(p, a, b, set(A))
            D = {x for x in range(F["q"]) if y[1 + x] < 0}
            if all(len(D & ell) == mu for cl in classes for ell in cl):
                n_half += 1
    return n_half


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 13, 17, 29):
        cmin, cmax, nlo = c_window_ends(p)
        emin, emax, neq = c_eq_ends(p)
        rows[str(p)] = {
            "c_min": cmin,
            "c_max_lo": cmax,
            "n_lo": nlo,
            "c_eq": emax,
            "n_eq": neq,
        }
        if emin != cmin:
            ok = False
        if emax > cmax:
            ok = False
        if p == 5 and (neq, emin, emax) != (1, 1, 1):
            ok = False
        if p == 13:
            if neq != 67 or emax != 617 or cmax != 642:
                ok = False
            if emax == cmax:
                ok = False
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "Q++≥Q_eq ⇔ c≤c_eq. Fail: Q_lo cut at p=13 (642≠617)."
        ),
    }


def prove_B() -> dict:
    n = count_halfnets_p5()
    ok = n == HPLUS[5] == 130
    if n == n_1d(5):
        ok = False
    return {
        "proved": bool(ok),
        "n_halfnet": n,
        "Hplus": HPLUS[5],
        "n_1d": n_1d(5),
        "theorem": (
            "p=5 half-net count=130=|H+|. Fail: count=n_1d=30."
        ),
    }


def prove_C() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        nh = type_minus_halfnet_count(p)
        rows[str(p)] = {"n_type_minus": n_1d(p), "n_halfnet": nh}
        if nh != 0:
            ok = False
        if p == 5:
            Y = np.load(YPATH[p])
            nplus = int((Y[:, 0] > 0).sum())
            if nplus != HPLUS[p]:
                ok = False
    if type_minus_halfnet_count(5) == n_1d(5):
        ok = False
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": "Type− is not a half-net. Fail: Type− ⊂ H+.",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "n_eq_13": c_eq_ends(13)[2],
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "c_eq leaves 67 lattice points at p=13. Half-net count "
            "unnamed for p≥7. Q_τ unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.451  Q_eq cut of c-window; p=5 half-net IFF", flush=True)
    A = prove_A()
    print(f"  A c_eq: {A['proved']} p13={A['by_p']['13']}", flush=True)
    B = prove_B()
    print(f"  B p5 enum: {B['proved']} n={B['n_halfnet']}", flush=True)
    C = prove_C()
    print(f"  C Type−: {C['proved']} {C['by_p']}", flush=True)
    D = prove_open()
    print(f"  D open n13={D['n_eq_13']} phi_F={D['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps(
            {
                "A": A["proved"],
                "B": B["proved"],
                "C": C["proved"],
                "p13": A["by_p"]["13"],
                "n_halfnet5": B["n_halfnet"],
            },
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.451",
        "title": (
            "Q++≥Q_eq cuts c_window at c_eq (67 points at p=13); "
            "p=5 half-net count equals |H+|; Type− is not a half-net"
        ),
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "c_eq_cut": A["proved"],
            "p5_halfnet_iff": B["proved"],
            "type_minus_not_halfnet": C["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15451.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
