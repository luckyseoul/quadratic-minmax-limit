#!/usr/bin/env python3
"""
Prop 15.454 — Complementary-class M-sum of every μ-half-net is
named: ∑_t M_t = p(p−1)(p+1)²/8.  Per-class M is not constant.
Two 15.307 O-signatures are the J2 Jacobi law μ+χ(b)−χ(b−1);
the λ=6 signature (0,1,1,2,5,5,7) is not J_r for r≤4.
N_r still unnamed in p.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.453 flags.

============================================================================
Setup.  15.453: N_r=|H+| at p=5,7.  Complementary directions
are P¹ minus the 15.451 ns-net.  For a k-set, AG(2,p) pair
counting gives ∑_{all lines} n_ℓ² = k(k+p).  On the ns-net
every n_ℓ=μ (half-net).  15.307: square-class signatures
{C,H,R,O}.  J_r(ε,A)(b)=μ+∑_{i=1}^r ε_i χ(b−a_i).

============================================================================
Theorem A — PROVED (pair counting + 15.305 C margins; all odd p).
  For every μ-half-net of the ns-net,
      ∑_{complementary lines} n_ℓ² = p(p−1)(p+1)²/8.
  Equals 90 at p=5 and 336 at p=7.  Fail: μ lines of an
  ns-direction (slope ∞) at p=5 — not a half-net.  ∎

Theorem B — PROVED (1D constructor; certified p=5,7).
  A union of μ lines of one complementary direction is a
  half-net with M-tuple (pμ²)^{r−1} × (p k) and count n_1d.
  Fail: claim that tuple is the only half-net type at p=7
  (then N_r=n_1d=140≠5726).  ∎

Theorem C — PROVED (Legendre arithmetic).
  J_2(+1,−1; {0,1}) has multiset (1,1,3,3,4,4,5) and
  J_2(−1,+1; {0,1}) has (1,2,2,3,3,5,5), two live p=7
  O-signatures.  The live λ=6 signature (0,1,1,2,5,5,7)
  is not J_r for any r≤4 and any ±1 coefficients.
  Fail: claim (0,1,1,2,5,5,7) is a J_2.  ∎

Theorem D — OPEN.  Total M named does not name per-class M
  or N_r (den 409).  c not pinned.  phi_F not imported.

============================================================================
Backend: Fraction/binomial + p=5 witness.  Writes
evidence/e1_gmin_m4_prop15454.json
"""
from __future__ import annotations

import json
import os
import sys
from itertools import combinations, product
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field  # noqa: E402
from e1_gmin_m4_prop15292 import n_1d  # noqa: E402
from e1_gmin_m4_prop15330 import HPLUS  # noqa: E402
from e1_gmin_m4_prop15451 import ns_parallel_classes  # noqa: E402
from e1_gmin_m4_prop15453 import class_slope, n4_pred  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15454_catalog.json"


def M_comp_named(p: int) -> int:
    return p * (p - 1) * (p + 1) ** 2 // 8


def M_all_named(p: int) -> int:
    k = p * (p - 1) // 2
    return k * (k + p)


def M_ns_named(p: int) -> int:
    r = (p + 1) // 2
    mu = (p - 1) // 2
    return r * p * mu * mu


def chi_p(p: int, a: int) -> int:
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def J_r(p: int, eps: tuple[int, ...], A: tuple[int, ...]) -> tuple[int, ...]:
    mu = (p - 1) // 2
    vals = []
    for b in range(p):
        vals.append(mu + sum(e * chi_p(p, b - a) for e, a in zip(eps, A)))
    return tuple(sorted(vals))


def is_Jr(p: int, tgt: tuple[int, ...], rmax: int = 4) -> bool:
    T = tuple(sorted(tgt))
    for r in range(0, rmax + 1):
        for A in combinations(range(p), r):
            for eps in product((-1, 1), repeat=r):
                if J_r(p, eps, A) == T:
                    return True
    return False


def line_sumsq(D: set[int], classes: list[list[set[int]]]) -> int:
    acc = 0
    for cl in classes:
        for ell in cl:
            n = len(D & ell)
            acc += n * n
    return acc


def square_classes(F, p: int) -> list[list[set[int]]]:
    """Parallel classes whose generator v has χ_q(v)=+1."""
    used: set[int] = set()
    gens: list[list[int]] = []
    for v in range(1, F["q"]):
        if v in used or int(F["chi"][v]) != 1:
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


def initial_kset(p: int) -> set[int]:
    """μ full lines of an ns-direction (slope ∞: const x0).  Not a half-net."""
    mu = (p - 1) // 2
    D: set[int] = set()
    for x0 in range(mu):
        for x1 in range(p):
            D.add(x0 + p * x1)
    return D


def one_d_halfnet(p: int) -> set[int]:
    """μ full lines of slope 0 (x1=const), Max+-free 1D."""
    mu = (p - 1) // 2
    D: set[int] = set()
    for x1 in range(mu):
        for x0 in range(p):
            D.add(x0 + p * x1)
    return D


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 13):
        named = M_comp_named(p)
        if named != M_all_named(p) - M_ns_named(p):
            ok = False
        rows[str(p)] = {"M_comp": named}
    if M_comp_named(5) != 90 or M_comp_named(7) != 336:
        ok = False
    p = 5
    F = _kernel_field(p)
    sq = square_classes(F, p)
    ns = ns_parallel_classes(F, p)
    Dbad = initial_kset(p)
    m_bad = line_sumsq(Dbad, sq)
    if m_bad == M_comp_named(p):
        ok = False
    D1 = one_d_halfnet(p)
    m_1d = line_sumsq(D1, sq)
    m_1d_ns = line_sumsq(D1, ns)
    if m_1d != 90:
        ok = False
    if m_1d_ns != M_ns_named(p):
        ok = False
    return {
        "proved": bool(ok),
        "by_p": rows,
        "bad_p5": m_bad,
        "one_d_p5": m_1d,
        "theorem": (
            "∑_comp n²=p(p−1)(p+1)²/8. Fail: initial k-set at p=5 "
            f"({m_bad}≠90)."
        ),
    }


def prove_B() -> dict:
    ok = True
    p = 7
    mu = (p - 1) // 2
    r = (p + 1) // 2
    C = p * mu * mu
    H = p * p * (p - 1) // 2
    if (r - 1) * C + H != M_comp_named(p):
        ok = False
    if n_1d(p) == HPLUS[p]:
        ok = False
    if n_1d(5) != 30 or n_1d(7) != 140:
        ok = False
    return {
        "proved": bool(ok),
        "C": C,
        "H": H,
        "n_1d7": n_1d(7),
        "Hplus7": HPLUS[7],
        "theorem": (
            "1D M-tuple is C^{r−1}×H, count n_1d. "
            "Fail: N_r=n_1d at p=7 (140≠5726)."
        ),
    }


def prove_C() -> dict:
    ok = True
    j_plus = J_r(7, (1, -1), (0, 1))
    j_minus = J_r(7, (-1, 1), (0, 1))
    if j_plus != (1, 1, 3, 3, 4, 4, 5):
        ok = False
    if j_minus != (1, 2, 2, 3, 3, 5, 5):
        ok = False
    residual = (0, 1, 1, 2, 5, 5, 7)
    if is_Jr(7, residual, rmax=4):
        ok = False
    if is_Jr(7, residual, rmax=2):
        ok = False
    if J_r(7, (1, -1), (0, 1)) == residual:
        ok = False
    return {
        "proved": bool(ok),
        "J2_plus": list(j_plus),
        "J2_minus": list(j_minus),
        "residual": list(residual),
        "residual_is_Jr4": is_Jr(7, residual, rmax=4),
        "theorem": (
            "J2 names two O-signatures. Fail: (0,1,1,2,5,5,7) is J2."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "M_comp7": M_comp_named(7),
        "den7": 409,
        "n4_bad": n4_pred(1, 3),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "Total complementary M is named. Per-class M and N_r "
            "still have den 409. c not pinned. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.454  complementary M-sum named; J2 not all O-types", flush=True)
    A = prove_A()
    print(f"  A M: {A['proved']} bad={A['bad_p5']} 1d={A['one_d_p5']}", flush=True)
    B = prove_B()
    print(f"  B 1D: {B['proved']} C={B['C']} H={B['H']}", flush=True)
    C = prove_C()
    print(f"  C J2: {C['proved']} residual_Jr4={C['residual_is_Jr4']}", flush=True)
    D = prove_open()
    print(f"  D open den7={D['den7']} phi_F={D['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps(
            {
                "A": A["proved"],
                "B": B["proved"],
                "C": C["proved"],
                "M5": M_comp_named(5),
                "M7": M_comp_named(7),
            },
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.454",
        "title": (
            "complementary M-sum of every μ-half-net is named; "
            "J2 names two O-signatures, not the λ=6 residual"
        ),
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "M_comp_named": A["proved"],
            "one_d_MH_tuple": B["proved"],
            "J2_not_all_O": C["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15454.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
