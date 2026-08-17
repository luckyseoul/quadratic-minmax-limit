#!/usr/bin/env python3
"""
Prop 15.453 — The 15.451 ns-net μ-half-net count equals |H+|
at p=7, not only p=5.  Every 3-net of order 5 has N_3=130.
A 4-net of order 7 has N_4∈{4844,5726} by the PGL type of
its direction 4-set; the ns-net is the 5726 type.
N_r still unnamed in p.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.452 flags.

============================================================================
Setup.  15.305 C / 15.451: Max+ ⊆ μ-half-nets of the
nonsquare-direction r-net in AG(2,p), r=(p+1)/2, μ=(p−1)/2.
15.451 B: at p=5 the converse is an enum (N_r=130=|H+|).
PGL(2,p) is 3-transitive on directions, not 4-transitive.
N_k = # of (0,1) p×p matrices with all line sums μ on a
fixed k-net (Max+-free).  N_2(5)=2040, N_2(7)=68938800.
N_4({0,∞,s,t}) at p=7 is the certified 15-pair table of
cpu_N3_profile (rows+cols+two extra slopes).

============================================================================
Theorem A — PROVED (full recurse; p=5).
  N_3(s)=130 for every extra slope s∈F_5^×, and N_2=2040.
  PGL 3-transitive ⇒ every 3-net has 130 μ-half-nets.
  The ns-net is a 3-net, so N_r(5)=130=|H+|.
  Fail: claim N_3=n_1d=30.  ∎

Theorem B — PROVED (15-pair table).
  At p=7, N_4({0,∞,s,t})=5726 if χ(t/s)=+1 or t/s=−1,
  and 4844 otherwise.  Fail: claim every 4-net has N_4=5726
  (six pairs have 4844).  Difference 882=2k², k=21.  ∎

Theorem C — PROVED (15.451 ns-classes + Möbius + B).
  The p=7 ns-net has slopes {1,2,5,6}.  For every ordered
  pair (a,b) from this set, the Möbius x↦(x−a)/(x−b) sends
  the other two slopes to a pair (s,t) of type 5726.
  Hence N_r(7)=5726=|H+|.  Fail: claim the ns-net is
  PGL-equivalent to {0,∞,1,3} (N_4=4844).  ∎

Theorem D — OPEN.  N_r=|H+| at p=5,7 is a census IFF, not
  a formula in p (5726=2·7·409).  c not pinned.  Q_τ
  unnamed.  phi_F not imported.

============================================================================
Backend: p=5 3-net recurse + p=7 Möbius on a certified N_4
table.  Writes evidence/e1_gmin_m4_prop15453.json
"""
from __future__ import annotations

import json
import os
import sys
from functools import lru_cache
from itertools import combinations, permutations
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

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15453_catalog.json"

# Certified N_4({0,∞,s,t}) at p=7 (cpu_N3_profile; Max+-free).
N4_P7: dict[tuple[int, int], int] = {
    (1, 2): 5726,
    (1, 3): 4844,
    (1, 4): 5726,
    (1, 5): 4844,
    (1, 6): 5726,
    (2, 3): 4844,
    (2, 4): 5726,
    (2, 5): 5726,
    (2, 6): 4844,
    (3, 4): 5726,
    (3, 5): 5726,
    (3, 6): 5726,
    (4, 5): 4844,
    (4, 6): 4844,
    (5, 6): 5726,
}

NS_NET7 = frozenset({1, 2, 5, 6})
BAD_NET7 = frozenset({0})  # {0,∞,1,3} in the axis frame


def chi_p(p: int, a: int) -> int:
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def n4_pred(s: int, t: int) -> int:
    """Named type of N_4({0,∞,s,t}) at p=7."""
    r = (t * pow(s, 5, 7)) % 7
    if chi_p(7, r) == 1 or r == 6:
        return 5726
    return 4844


def mobius_pair(a: int, b: int, x: int, p: int) -> int:
    """Send a→0, b→∞."""
    return ((x - a) * pow((x - b) % p, p - 2, p)) % p


def point_xy(x: int, p: int) -> tuple[int, int]:
    return x % p, x // p


def class_slope(cl: list[set[int]], p: int) -> str:
    ell = next(iter(cl))
    pts = list(ell)
    x0, y0 = point_xy(pts[0], p)
    x1, y1 = point_xy(pts[1], p)
    dx = (x1 - x0) % p
    dy = (y1 - y0) % p
    if dx == 0:
        return "inf"
    s = (dy * pow(dx, p - 2, p)) % p
    return str(s)


def ns_net_slopes(p: int) -> list[str]:
    F = _kernel_field(p)
    classes = ns_parallel_classes(F, p)
    return [class_slope(cl, p) for cl in classes]


@lru_cache(maxsize=1)
def n3_by_slope_p5() -> dict[int, int]:
    p, mu = 5, 2
    rows = list(combinations(range(p), mu))
    nch = len(rows)
    n3 = {s: 0 for s in range(1, p)}
    n2 = 0
    col = [0] * p
    stack = [0] * p

    def rec(d: int) -> None:
        nonlocal n2
        if d == p:
            if all(c == mu for c in col):
                n2 += 1
                mat = [[0] * p for _ in range(p)]
                for x in range(p):
                    for y in rows[stack[x]]:
                        mat[x][y] = 1
                for s in range(1, p):
                    ok = True
                    for b in range(p):
                        sm = 0
                        for x in range(p):
                            sm += mat[x][(s * x + b) % p]
                        if sm != mu:
                            ok = False
                            break
                    if ok:
                        n3[s] += 1
            return
        for i in range(nch):
            ch = rows[i]
            if any(col[y] + 1 > mu for y in ch):
                continue
            stack[d] = i
            for y in ch:
                col[y] += 1
            rec(d + 1)
            for y in ch:
                col[y] -= 1

    rec(0)
    n3[0] = n2  # stash N_2
    return n3


def prove_A() -> dict:
    n3 = dict(n3_by_slope_p5())
    n2 = n3.pop(0)  # copy: cached dict is immutable to callers
    ok = n2 == 2040
    for s, v in n3.items():
        if v != HPLUS[5]:
            ok = False
    if any(v == n_1d(5) for v in n3.values()):
        ok = False
    sl = ns_net_slopes(5)
    if len(sl) != 3:
        ok = False
    return {
        "proved": bool(ok),
        "n2": n2,
        "n3": {str(s): v for s, v in n3.items()},
        "ns_slopes": sl,
        "n_1d": n_1d(5),
        "theorem": "N_3(5)=130=|H+| every slope. Fail: N_3=n_1d=30.",
    }


def prove_B() -> dict:
    ok = True
    n5726 = 0
    n4844 = 0
    rows = {}
    for (s, t), v in N4_P7.items():
        pred = n4_pred(s, t)
        rows[f"{s},{t}"] = {"N4": v, "pred": pred}
        if pred != v:
            ok = False
        if v == 5726:
            n5726 += 1
        elif v == 4844:
            n4844 += 1
        else:
            ok = False
    if n4844 != 6 or n5726 != 9:
        ok = False
    if 5726 - 4844 != 2 * 21 * 21:
        ok = False
    if n4_pred(1, 3) == 5726:
        ok = False
    return {
        "proved": bool(ok),
        "n5726": n5726,
        "n4844": n4844,
        "diff": 5726 - 4844,
        "bad_pair_13": n4_pred(1, 3),
        "theorem": (
            "N_4=5726 iff χ(t/s)=+1 or t/s=−1. Fail: all 4-nets 5726."
        ),
    }


def prove_C() -> dict:
    sl = ns_net_slopes(7)
    finite = {int(s) for s in sl if s != "inf"}
    ok = finite == set(NS_NET7) and "inf" not in sl and 0 not in finite
    images = []
    for a, b in permutations(NS_NET7, 2):
        rest = [mobius_pair(a, b, x, 7) for x in NS_NET7 if x not in (a, b)]
        s, t = sorted(rest)
        v = N4_P7[(s, t)]
        images.append({"ab": [a, b], "st": [s, t], "N4": v})
        if v != 5726:
            ok = False
    if n4_pred(1, 3) != 4844:
        ok = False
    if len(images) != 12:
        ok = False
    if HPLUS[7] != 5726:
        ok = False
    return {
        "proved": bool(ok),
        "ns_slopes": sl,
        "n_images": len(images),
        "all_5726": all(r["N4"] == 5726 for r in images),
        "bad_type": n4_pred(1, 3),
        "Hplus": HPLUS[7],
        "theorem": (
            "p=7 ns-net is 5726-type; N_r=5726=|H+|. "
            "Fail: equivalent to {0,∞,1,3} (4844)."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "Nr5": 130,
        "Nr7": 5726,
        "den7": 409,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "N_r=|H+| at p=5,7 is a type-IFF, not a formula in p. "
            "c not pinned. Q_τ unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.453  ns-net N_r = |H+| at p=7 (5726-type 4-net)", flush=True)
    A = prove_A()
    print(f"  A N3: {A['proved']} n2={A['n2']} n3={A['n3']}", flush=True)
    B = prove_B()
    print(f"  B N4 type: {B['proved']} 5726={B['n5726']} 4844={B['n4844']}", flush=True)
    C = prove_C()
    print(f"  C ns-net: {C['proved']} slopes={C['ns_slopes']}", flush=True)
    D = prove_open()
    print(f"  D open den7={D['den7']} phi_F={D['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps(
            {
                "A": A["proved"],
                "B": B["proved"],
                "C": C["proved"],
                "ns5": A["ns_slopes"],
                "ns7": C["ns_slopes"],
            },
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.453",
        "title": (
            "ns-net μ-half-net count equals |H+| at p=7 "
            "(5726-type 4-net), not only p=5"
        ),
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "n3_p5_is_Hplus": A["proved"],
            "n4_p7_type_split": B["proved"],
            "nr7_is_Hplus": C["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15453.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
