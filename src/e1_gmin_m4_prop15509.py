#!/usr/bin/env python3
"""
Prop 15.509 — Burnside involution class u↦−u+b in Aut_∞.

  |class|=q.  Cycle type (1,1,2^{(q−1)/2}).  On the p=3,5,7 H_+
  caches every such g has
      Fix(g) = 2 · 3^{C((p−1)/2, 2)}
  (2, 6, 54).  Fail: |class|=2q (the 15.508 ω-size);
  Fix = n_1d/p (hits 2,6 and dies at p=7: 20≠54).
  Does not name |H_+| or Q_τ.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.508 flags.

============================================================================
Setup.  G={u↦a u^{p^k}+b : a∈□, k∈{0,1}} as in 15.308, |G|=q(q−1).
χ_q(−1)=(−1)^{(q−1)/2}.  For odd p, (q−1)/2=(p²−1)/2 is even, so
−1 is a square in F_q and every
    g_b: u ↦ −u + b
lies in G (k=0).  There are q such maps.

Fixed points of g_b on PG(1): ∞ and the unique u with 2u=b
(char ≠ 2).  The remaining q−1 points form (q−1)/2 transpositions.

============================================================================
Theorem A — PROVED (field; odd p).
  |{g_b : b∈F_q}|=q and χ_q(−1)=1.  Fail: |class|=2q.  ∎

Theorem B — PROVED (field; odd p).
  Every g_b has cycle type (1,1,2^{(q−1)/2}).
  Fail: the 15.508 3-cycle type.  ∎

Theorem C — PROVED (cache; certified p=3,5,7).
  Every g_b has
      Fix_{H_+}(g_b) = 2 · 3^{C((p−1)/2, 2)}
  constantly (2, 6, 54).  Fail: n_1d/p (2, 6, 20).  ∎

Theorem D — OPEN.  The Fix formula is a 3-point cache identity
  (Claude direction: C((p−1)/2,2) counts pairs of F_p-squares).
  Factorization Fix ≅ {±}×{0,+,−}^{pairs} is not exhibited.
  Predicted p=11 value 118098 is untested (no H_+ cache).
  n_orbits and |H_+| / D / Q_τ stay unnamed.  phi_F not imported.  ∎

============================================================================
Backend: apply_map + cache Fix (ProcessPool over the q maps).
Writes evidence/e1_gmin_m4_prop15509.json
"""
from __future__ import annotations

import json
import os
from concurrent.futures import ProcessPoolExecutor
from math import comb
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import _kernel_field
from e1_gmin_m4_prop15292 import n_1d
from e1_gmin_m4_prop15307 import load_D
from e1_gmin_m4_prop15308 import pack_y
from e1_gmin_m4_prop15508 import apply_map

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "evidence" / "e1_gmin_m4_prop15509.json"

CERT = (3, 5, 7)
W = 86


def class_size_named(p: int) -> int:
    return p * p


def class_size_wrong(p: int) -> int:
    """Fail-when-wrong: 15.508 ω-class size 2q."""
    return 2 * p * p


def named_cyc(p: int) -> tuple[int, ...]:
    q = p * p
    return tuple([1, 1] + [2] * ((q - 1) // 2))


def fix_named(p: int) -> int:
    return 2 * (3 ** comb((p - 1) // 2, 2))


def fix_wrong(p: int) -> int:
    """Fail-when-wrong: n_1d/p (hits p=3,5; dies at p=7)."""
    return n_1d(p) // p


def chi_q_neg1(p: int) -> int:
    q = p * p
    return 1 if ((q - 1) // 2) % 2 == 0 else -1


def cycle_type(perm: np.ndarray) -> tuple[int, ...]:
    seen = np.zeros(perm.shape[0], dtype=np.uint8)
    cyc = []
    for i in range(perm.shape[0]):
        if seen[i]:
            continue
        L = 0
        j = i
        while not seen[j]:
            seen[j] = 1
            j = int(perm[j])
            L += 1
        cyc.append(L)
    return tuple(sorted(cyc))


def _fix_job(spec):
    p, bs = spec
    D, F = load_D(p)
    Y = np.stack([pack_y(D[i]) for i in range(D.shape[0])])
    neg = int(F["neg1"])
    want = named_cyc(p)
    out = []
    for b in bs:
        g = apply_map(F, neg, int(b))
        Yg = Y[:, g]
        fx = int(np.all(Yg == Y, axis=1).sum())
        out.append((fx, cycle_type(g) == want))
    return out


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (3, 5, 7, 11, 13, 17, 19):
        F = _kernel_field(p)
        chi = chi_q_neg1(p)
        n = F["q"]
        rows[str(p)] = {
            "n": n,
            "named": class_size_named(p),
            "chi_q_neg1": chi,
            "neg1_in_squares": int(F["neg1"]) in set(int(s) for s in F["squares"]),
        }
        if n != class_size_named(p) or chi != 1:
            ok = False
        if n == class_size_wrong(p):
            ok = False
        if int(F["neg1"]) not in set(int(s) for s in F["squares"]):
            ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": "|−id-class|=q and χ_q(−1)=1. Fail |class|=2q.",
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in (3, 5, 7, 11, 13):
        F = _kernel_field(p)
        neg = int(F["neg1"])
        want = named_cyc(p)
        n_ok = 0
        for b in range(F["q"]):
            g = apply_map(F, neg, b)
            if cycle_type(g) == want:
                n_ok += 1
        rows[str(p)] = {"n_ok": n_ok, "q": F["q"], "want_n2": (F["q"] - 1) // 2}
        if n_ok != F["q"]:
            ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": "cycle type (1,1,2^{(q−1)/2}). Fail 3-cycles.",
    }


def prove_C() -> dict:
    ok = True
    rows = {}
    for p in CERT:
        q = p * p
        sw = max(1, (q + W - 1) // W)
        jobs = [(p, list(range(lo, min(q, lo + sw)))) for lo in range(0, q, sw)]
        fixes = []
        cyc_ok = True
        with ProcessPoolExecutor(max_workers=min(W, len(jobs))) as ex:
            for part in ex.map(_fix_job, jobs):
                for fx, cok in part:
                    fixes.append(fx)
                    cyc_ok = cyc_ok and cok
        named = fix_named(p)
        wrong = fix_wrong(p)
        if not fixes or min(fixes) != named or max(fixes) != named:
            ok = False
        # n_1d/p interpolates p=3,5 and dies at p=7 (20≠54).
        if p == 7 and named == wrong:
            ok = False
        if p == 7 and any(f == wrong for f in fixes):
            ok = False
        if not cyc_ok:
            ok = False
        rows[str(p)] = {
            "n": len(fixes),
            "fix": named,
            "wrong": wrong,
            "live_min": min(fixes) if fixes else None,
            "live_max": max(fixes) if fixes else None,
            "cyc_ok": cyc_ok,
        }
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "Fix(g_b)=2·3^{C((p−1)/2,2)} constantly on H_+. "
            "Fail n_1d/p."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "Q_tau_named_in_p": False,
        "D_named_in_p": False,
        "phi_F_imported": False,
        "p11_predicted": fix_named(11),
        "p13_predicted": fix_named(13),
        "note": (
            "Involution class size and cycle type are field-named. "
            "Fix formula is a p=3,5,7 cache identity. Factorization "
            "and p≥11 untested. |H+| / D / Q_τ unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("15.509 Burnside −id class |class|=q, Fix=2·3^{C((p-1)/2,2)}", flush=True)
    A, B, C, D = prove_A(), prove_B(), prove_C(), prove_open()
    out = {
        "prop": "15.509",
        "title": "Burnside −id class size q, Fix=2·3^{C((p−1)/2,2)}; D unnamed",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "C": C,
        "D": D,
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
    print(
        f"A={A['proved']} B={B['proved']} C={C['proved']} phi_F={out['phi_F_ge_6']}",
        flush=True,
    )
    return out


if __name__ == "__main__":
    main()
