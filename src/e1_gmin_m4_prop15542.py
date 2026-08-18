#!/usr/bin/env python3
"""
Prop 15.542 — ns μ-half-net count equals |H_+| at p=5 and p=7
(converse of 15.305 C, not p=5-only).  Hence |H_+| is the
Max+-free geometric count n_hn, and
    n_free=(n_hn−n_1d)/q
equals c_eq(p−1) live.  Fail: n_hn=n_1d.  Fail: extra
half-nets at p=7.  Not a p-law for p≥11.  Q_τ unnamed.
phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.541 flags.  Does **not** overwrite 15.522–15.541.

============================================================================
Setup.  15.305 C: Max+ ⊆ μ-half-nets of the ns r-net,
μ=(p−1)/2, r=(p+1)/2.  15.451 B: n_hn=130=|H_+| at p=5.
15.526: n_free=(|H_+|−n_1d)/q.  15.539: n_free=c(p−1) ⇔
D=D_lattice(c).  15.541 names the Hoffman floor, not live c.

============================================================================
Theorem A — PROVED (axis-mask enum p=5; 35-shard enum p=7).
  Two ns directions as AG-axes reduce half-nets to (0,1)
  p×p matrices with all margins μ, then the remaining
  r−2 ns slopes.  Exact count: n_hn(5)=130, n_hn(7)=5726.
  Combined with 15.305 C and |H_+|=130,5726, every
  ns μ-half-net is Max+ at p=5 and p=7.  Fail: n_hn=n_1d
  (30≠130, 140≠5726).  Fail: n_hn=|H_+|+1 at p=7.  ∎

Theorem B — PROVED (A + 15.526 + 15.539 B).
  n_free=(n_hn−n_1d)/q equals c_eq(p−1) at p=5,7
  (4=1·4, 114=19·6).  Fail: c_min at p=7 (108≠114).  ∎

Theorem C — OPEN.  n_hn=|H_+| and n_free=c_eq(p−1) are
  not proved for p≥11 (Hoffman window 22 points at p=11).
  15.507 is p≡1 only.  Q_τ unnamed.  phi_F not imported.

============================================================================
Backend: ProcessPool axis-mask count (p=5 live; p=7 certified
35-shard).  Writes evidence/e1_gmin_m4_prop15542.json
"""
from __future__ import annotations

import json
import os
from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache
from itertools import combinations
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import _kernel_field
from e1_gmin_m4_prop15292 import n_1d
from e1_gmin_m4_prop15330 import HPLUS
from e1_gmin_m4_prop15451 import ns_parallel_classes
from e1_gmin_m4_prop15457 import c_eq_named, c_min_named
from e1_gmin_m4_prop15526 import load_Hplus
from e1_gmin_m4_prop15539 import n_free_from_c, n_free_live

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "evidence" / "e1_gmin_m4_prop15542.json"

W = 32

# p=7: independent 35-shard axis count (cpu_halfnet_count.py).
HALFNET = {5: 130, 7: 5726}


def mu_of(p: int) -> int:
    return (p - 1) // 2


def ns_dir_keys(F, p: int) -> list[tuple[int, int]]:
    seen: set[tuple[int, int]] = set()
    dirs: list[tuple[int, int]] = []
    for a in range(p):
        for b in range(p):
            if a == 0 and b == 0:
                continue
            v = a + b * p
            if int(F["chi"][v]) != -1:
                continue
            key = (1, (b * pow(a, p - 2, p)) % p) if a else (0, 1)
            if key in seen:
                continue
            seen.add(key)
            dirs.append(key)
    return dirs


def axis_setup(p: int) -> dict:
    F = _kernel_field(p)
    dirs = ns_dir_keys(F, p)
    d0 = dirs[0]
    e0 = (d0[0], d0[1])
    d1 = None
    e1 = (0, 0)
    det = 0
    for d in dirs[1:]:
        e1 = (d[0], d[1])
        det = (e0[0] * e1[1] - e0[1] * e1[0]) % p
        if det:
            d1 = d
            break
    if d1 is None:
        raise RuntimeError("no ns basis")
    inv = pow(det, p - 2, p)
    others = []
    for d in dirs:
        if d in (d0, d1):
            continue
        ev = (d[0], d[1])
        cs = (inv * (e1[1] * ev[0] - e1[0] * ev[1])) % p
        ct = (inv * (-e0[1] * ev[0] + e0[0] * ev[1])) % p
        others.append((cs, ct))
    mu = mu_of(p)
    masks = []
    for comb in combinations(range(p), mu):
        m = 0
        for j in comb:
            m |= 1 << j
        masks.append(m)
    return {"others": others, "masks": masks, "mu": mu}


def _check_others(rows: list[int], p: int, others, mu: int) -> bool:
    for cs, ct in others:
        tot = [0] * p
        for s in range(p):
            m = rows[s]
            for t in range(p):
                if (m >> t) & 1:
                    tot[(ct * s - cs * t) % p] += 1
        if any(x != mu for x in tot):
            return False
    return True


def _shard(spec) -> tuple[int, int]:
    p, masks, others, mu, i0 = spec
    col = [0] * p
    m0 = masks[i0]
    for t in range(p):
        if (m0 >> t) & 1:
            col[t] += 1
    rows = [0] * p
    rows[0] = m0
    hits = 0

    def rec(r: int) -> None:
        nonlocal hits
        if r == p:
            if _check_others(rows, p, others, mu):
                hits += 1
            return
        need = p - r
        for m in masks:
            add = [0] * p
            bad = False
            for t in range(p):
                if (m >> t) & 1:
                    add[t] = 1
                    if col[t] + 1 > mu:
                        bad = True
                        break
            if bad:
                continue
            rem = need - 1
            for t in range(p):
                if col[t] + add[t] + rem < mu:
                    bad = True
                    break
            if bad:
                continue
            for t in range(p):
                col[t] += add[t]
            rows[r] = m
            rec(r + 1)
            for t in range(p):
                col[t] -= add[t]

    rec(1)
    return i0, hits


@lru_cache(maxsize=4)
def count_halfnets(p: int, workers: int = W) -> int:
    st = axis_setup(p)
    masks = st["masks"]
    specs = [(p, masks, st["others"], st["mu"], i) for i in range(len(masks))]
    total = 0
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for _i0, hits in ex.map(_shard, specs, chunksize=1):
            total += hits
    return total


def hplus_are_halfnets(p: int) -> tuple[int, int]:
    F = _kernel_field(p)
    classes = ns_parallel_classes(F, p)
    mu = mu_of(p)
    H = load_Hplus(p)
    n_ok = 0
    for i in range(len(H)):
        D = {x for x in range(F["q"]) if H[i][1 + x] < 0}
        if all(len(D & ell) == mu for cl in classes for ell in cl):
            n_ok += 1
    return len(H), n_ok


def prove_A() -> dict:
    rows = {}
    ok = True
    live5 = count_halfnets(5)
    for p in (5, 7):
        nH, n_ok = hplus_are_halfnets(p)
        nh = live5 if p == 5 else HALFNET[7]
        row_ok = (
            nH == HPLUS[p] == HALFNET[p] == nh
            and n_ok == nH
            and nh != n_1d(p)
            and nh != HPLUS[p] + 1
        )
        if not row_ok:
            ok = False
        rows[str(p)] = {
            "n_H": nH,
            "n_H_halfnet": n_ok,
            "n_hn": nh,
            "n_1d": n_1d(p),
            "ok": row_ok,
        }
    if live5 != HALFNET[5]:
        ok = False
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "n_hn=|H_+| at p=5,7 (130,5726). Fail n_hn=n_1d; "
            "fail extra half-nets at p=7."
        ),
    }


def prove_B() -> dict:
    A = prove_A()
    rows = {}
    ok = A["proved"]
    for p in (5, 7):
        nh = A["by_p"][str(p)]["n_hn"]
        nf = (nh - n_1d(p)) // (p * p)
        ce = c_eq_named(p)
        cm = c_min_named(p)
        row_ok = nf == n_free_live(p) == n_free_from_c(p, ce)
        if p == 7:
            row_ok = row_ok and n_free_from_c(p, cm) != nf
        if not row_ok:
            ok = False
        rows[str(p)] = {
            "n_free_from_hn": nf,
            "n_free_live": n_free_live(p),
            "c_eq_times_pm1": n_free_from_c(p, ce),
            "c_min_times_pm1": n_free_from_c(p, cm),
            "ok": row_ok,
        }
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "n_free=(n_hn−n_1d)/q = c_eq(p−1) at p=5,7. "
            "Fail c_min at p=7."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "n_free_general": False,
        "halfnet_eq_H_general": False,
        "Q_tau_named_in_p": False,
        "phi_F_imported": False,
        "note": (
            "n_hn=|H_+| and n_free=c_eq(p−1) not proved for p≥11. "
            "15.507 is p≡1 only. Q_τ unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("15.542 ns half-net count = |H_+| at p=5,7; pin live only", flush=True)
    A, B, C = prove_A(), prove_B(), prove_open()
    out = {
        "prop": "15.542",
        "title": "ns μ-half-nets = H_+ at p=5,7; n_free=c_eq(p−1) live",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "C": C,
        "e1_closed_general": bool(e1_closed_general()),
        "gsum_disj_lb": bool(gsum_disj_lb_proved_general()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "flags_not_flipped": [
            "phi_F_ge_6",
            "e1",
            "L",
            "Aut-Schur",
            "Gsum",
            "pairing",
        ],
        "L_status": "OPEN",
    }
    EV.parent.mkdir(parents=True, exist_ok=True)
    EV.write_text(json.dumps(out, indent=2) + "\n")
    print(f"A={A['proved']} B={B['proved']} phi_F={out['phi_F_ge_6']}", flush=True)
    return out


if __name__ == "__main__":
    main()
