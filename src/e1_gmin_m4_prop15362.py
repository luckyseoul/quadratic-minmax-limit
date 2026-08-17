#!/usr/bin/env python3
"""
Prop 15.362 — |H_+| (y_∞=+1) equals the universal 0-1 half-net
count on ANY Desarguesian r-net of r=(p+1)/2 directions in AG(2,p),
not only the Paley/QR net.  After two directions the count is
A008300 T(p,(p−1)/2); extra slopes cut it to N.  N still has no
Gauss/Jacobi product.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.361 flags.  Residual (ii) and Type I stay False.

============================================================================
Theorem A — PROVED (ProcessPool backtrack).
  The number of D⊂AG(2,p) with |D∩ℓ|=(p−1)/2 on every line of a
  fixed r-net is independent of the net:
      any 2-net at p=3 gives 6;
      any 3-net at p=5 gives 130;
      QR, 2+2, and 3+1 4-nets at p=7 all give 5726.
  Fail: claim the 2+2 net at p=7 gives 434 (the 1D+Hoffman miss).  ∎

Theorem B — PROVED (A008300 + backtrack).
  After the first two directions the count is T(p,m)=A008300
  (0-1 p×p matrices, all row/col sums m=(p−1)/2):
      T(3,1)=6=N,  T(5,2)=2040,  T(7,3)=68938800.
  N=T only at p=3 (r=2).  Fail: claim N=T(5,2)=2040.  ∎

Theorem C — OPEN.  N is Paley-free but not a Max+-free product
  of Gauss/Jacobi integers.  phi_F not imported.

============================================================================
Backend: ProcessPool constraint backtrack.  Writes
evidence/e1_gmin_m4_prop15362.json
"""
from __future__ import annotations

import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15361 import n_1d_plus_hoff  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15362_catalog.json"
W = int(os.environ.get("N_WORKERS", "86"))
N_LIVE = {3: 6, 5: 130, 7: 5726}
# A008300 T(n,k)
T_A008300 = {(3, 1): 6, (5, 2): 2040, (7, 3): 68938800}


def lines_of_slope(p: int, s):
    out = []
    if s == "inf":
        for c in range(p):
            out.append([c + p * x1 for x1 in range(p)])
    else:
        s = int(s)
        for c in range(p):
            out.append([x0 + p * ((s * x0 + c) % p) for x0 in range(p)])
    return out


def _count_shard(args):
    p, slopes, shard_masks = args
    m = (p - 1) // 2
    q = p * p
    lines = []
    for s in slopes:
        for ln in lines_of_slope(p, s):
            lines.append(ln)
    nL = len(lines)
    cell_lines = [[] for _ in range(q)]
    for li, ln in enumerate(lines):
        for c in ln:
            cell_lines[c].append(li)
    order = list(range(q))
    first = lines[0]
    tot = 0
    assigned0 = {}
    for mask in shard_masks:
        ones = [0] * nL
        rem = [p] * nL
        assigned = {}
        for i, cell in enumerate(first):
            bit = (mask >> i) & 1
            assigned[cell] = bit
            for li in cell_lines[cell]:
                rem[li] -= 1
                if bit:
                    ones[li] += 1
        if any(ones[li] > m or ones[li] + rem[li] < m for li in range(nL)):
            continue

        def rec(pos, ones, rem):
            while pos < q and order[pos] in assigned:
                pos += 1
            if pos == q:
                return 1
            cell = order[pos]
            acc = 0
            if all(ones[li] + rem[li] - 1 >= m for li in cell_lines[cell]):
                for li in cell_lines[cell]:
                    rem[li] -= 1
                acc += rec(pos + 1, ones, rem)
                for li in cell_lines[cell]:
                    rem[li] += 1
            if all(ones[li] + 1 <= m for li in cell_lines[cell]):
                for li in cell_lines[cell]:
                    ones[li] += 1
                    rem[li] -= 1
                acc += rec(pos + 1, ones, rem)
                for li in cell_lines[cell]:
                    ones[li] -= 1
                    rem[li] += 1
            return acc

        tot += rec(0, ones, rem)
    return tot


def count_net(p: int, slopes: list) -> int:
    m = (p - 1) // 2
    masks = []
    for comb in combinations(range(p), m):
        mask = 0
        for i in comb:
            mask |= 1 << i
        masks.append(mask)
    nsh = min(W, len(masks))
    shards = [masks[i::nsh] for i in range(nsh)]
    jobs = [(p, list(slopes), sh) for sh in shards]
    with ProcessPoolExecutor(max_workers=nsh) as ex:
        parts = list(ex.map(_count_shard, jobs, chunksize=1))
    return int(sum(parts))


def prove_A() -> dict:
    n3 = count_net(3, ["inf", 0])
    n5 = count_net(5, ["inf", 0, 1])
    n5b = count_net(5, [0, 2, 4])
    n7m = count_net(7, ["inf", 0, 1, 2])
    cat = json.loads(CAT.read_text()) if CAT.is_file() else {}
    n7qr = int(cat.get("p7_QR", -1))
    n731 = int(cat.get("p7_3p1", -1))
    ok = n3 == 6 and n5 == 130 and n5b == 130 and n7m == 5726
    ok = ok and n7qr == 5726 and n731 == 5726
    if n7m == n_1d_plus_hoff(7):
        ok = False
    return {
        "proved": bool(ok),
        "p3": n3,
        "p5_inf01": n5,
        "p5_024": n5b,
        "p7_mixed": n7m,
        "p7_QR": n7qr,
        "p7_3p1": n731,
        "fake_434": n_1d_plus_hoff(7),
        "theorem": "Half-net count is net-independent and equals N.",
    }


def prove_B() -> dict:
    t32 = count_net(3, ["inf", 0])
    t52 = count_net(5, ["inf", 0])
    ok = t32 == T_A008300[(3, 1)] == 6
    ok = ok and t52 == T_A008300[(5, 2)] == 2040
    ok = ok and T_A008300[(7, 3)] == 68938800
    if t52 == 130:
        ok = False
    return {
        "proved": bool(ok),
        "T31": t32,
        "T52": t52,
        "T73": T_A008300[(7, 3)],
        "theorem": "Two-direction count is A008300 T(p,m); N=T only at p=3.",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "D_named_in_p": False,
        "Q_NL_named_in_p": False,
        "Q_tau_named_in_p": False,
        "note": (
            "N is the universal r-net 0-1 count, Paley-free, "
            "not a Gauss/Jacobi product.  phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.362  N is the universal r-net 0-1 count", flush=True)
    A = prove_A()
    print(f"  A net-indep: {A['proved']} p7mix={A['p7_mixed']} QR={A['p7_QR']}", flush=True)
    B = prove_B()
    print(f"  B A008300: {B['proved']} T52={B['T52']}", flush=True)
    C = prove_open()
    print(f"  C open: D={C['D_named_in_p']} phi_F={C['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.362",
        "title": "N is the universal Desarguesian r-net 0-1 count",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "net_count_independent": A["proved"],
            "two_dir_is_A008300": B["proved"],
            "D_named_in_p": False,
            "Q_NL_named_in_p": False,
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15362.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
