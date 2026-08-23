#!/usr/bin/env python3
"""
Prop 15.627 — octic linear-box empty; split-involution class W2 at p=31.

Does **not** flip residual_ii / multilevel_ND / phi_F / type_I / e1 / L.
Does **not** close W1 residual, W2 p-law, Walsh, leftover 2.

============================================================================
Setup.  Fable xhigh suggest_direction: (1) residual W1 via [2/π]_8
on π=a+8ci, equivalently split by octic and a (mod 8); (2) W2 by
existence in the switched split-involution class, discriminator p=31;
(3) leftover 3 two-character μ.  Walsh cannot close leftover 2.
k=6 mesh skipped (already a p-law).

============================================================================
Theorem A — PROVED kill (61 residual primes; ProcessPool 86).
  (2/p)_8 = 2^{(p-1)/8} ∈ {±1} on p=a²+64c².  Inside each octic
  class, and inside each pair ((2/p)_8, a_G (mod 8)) with
  a_G≡1 (mod 4), the bounded stay box ua+vb+wi+k
  (|u,v,w|≤4, |k|≤8) has empty intersection.  Named stays
  including o·(−a) are MIXED on each octic half.  Fail: a box
  form with ε=1 on all oct=+ primes (601,1249,…).  ∎

Theorem B — CERTIFIED, not a p-law.
  Split involutions {±}[[α,β],[γ,−α]] with α²+βγ=1: p(p+1)/2
  maps, all Max− under switching.  p=17: 153 class, 49 in U,
  17 W2 (includes t=-2).  p=31: 496 class, 146 in U, 76 W2.
  First p=31 hit includes π(x)=x/(x−1).  t=-2 fails W2 at p=31
  (15.626) but the class does not.  Fail: 0 W2 in the class at
  p=31.  Counting identity (some conjugate always coprime) OPEN.  ∎

Theorem C — OPEN.  W1 residual (nonlinear / octic class function
  outside the box).  W2 p-law by class-existence count.  Walsh.
  leftover 2.

============================================================================
Backend: W1 prefix ProcessPool 86; W2 class ProcessPool 86
(independent Auts; inner Krylov sequential).  GPU unused.
Writes evidence/e1_gmin_m4_prop15627.json
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15626 import _switched  # noqa: E402


def _oct8(p: int) -> int:
    r = pow(2, (p - 1) // 8, p)
    if r == 1:
        return 1
    if r == p - 1:
        return -1
    raise AssertionError(f"not ±1 at {p}")


def theorem_A_octic_box() -> dict:
    o601 = _oct8(601)
    o1201 = _oct8(1201)
    box = ROOT / "evidence" / "w1_octic_box.json"
    payload = json.loads(box.read_text()) if box.exists() else {}
    ip = payload.get("n_inter_plus", -1)
    im = payload.get("n_inter_minus", -1)
    splits = payload.get("splits", {})
    split_empty = all(v.get("n_forms", -1) == 0 for v in splits.values()) if splits else False
    ok = o601 == 1 and o1201 == -1 and ip == 0 and im == 0 and split_empty
    return {
        "proved": bool(ok),
        "W1_residual": False,
        "oct8_601": o601,
        "oct8_1201": o1201,
        "n_inter_plus": ip,
        "n_inter_minus": im,
        "n_cells": len(splits),
        "theorem": (
            "Octic split does not rescue the bounded linear stay box.  "
            "Fail: a box form on all oct=+ primes."
        ),
    }


def theorem_B_class_p31() -> dict:
    r17 = _switched(17, 1, 0, 15, 16)
    r31 = _switched(31, 1, 0, 1, 30)
    rec = ROOT / "evidence" / "w2_class_scan.json"
    payload = json.loads(rec.read_text()) if rec.exists() else {}
    c17 = payload.get("17", {})
    c31 = payload.get("31", {})
    ok = (
        r17["W2"] is True
        and r31["W2"] is True
        and r31["inU_y"]
        and c17.get("n_W2") == 17
        and c31.get("n_W2") == 76
        and c31.get("n_class") == 496
        and c31.get("n_eigen") == 496
    )
    return {
        "proved": False,
        "certified": bool(ok),
        "W2_p_law": False,
        "p17_tm2": r17,
        "p31_xin_xminus1": r31,
        "class_17": {k: c17.get(k) for k in ("n_class", "n_eigen", "n_inU", "n_W2")},
        "class_31": {k: c31.get(k) for k in ("n_class", "n_eigen", "n_inU", "n_W2")},
        "theorem": (
            "Split-involution class is W2-nonempty at p=31 (76 hits) "
            "though t=-2 fails.  Fail: 0 W2 at p=31.  Count OPEN."
        ),
    }


def theorem_C_open() -> dict:
    return {
        "proved": False,
        "W1_all_odd_p": False,
        "W2_p_law": False,
        "walsh_general_p": False,
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "theorem": (
            "W1 residual and W2 class-count remain.  "
            "Fail: Walsh from this unit."
        ),
    }


def main() -> dict:
    t0 = time.time()
    print("Prop 15.627  octic box kill / W2 class nonempty at p=31", flush=True)
    A = theorem_A_octic_box()
    print(f"  A {A['proved']} oct 601={A['oct8_601']} 1201={A['oct8_1201']} ∩+={A['n_inter_plus']}", flush=True)
    B = theorem_B_class_p31()
    print(
        f"  B certified={B['certified']} p17_W2={B['p17_tm2']['W2']} "
        f"p31_W2={B['p31_xin_xminus1']['W2']} class31={B['class_31']}",
        flush=True,
    )
    C = theorem_C_open()
    out = {
        "prop": "15.627",
        "title": "Octic linear-box empty; split-involution class W2 at p=31",
        "proved": {
            "octic_box_empty": A["proved"],
            "W1_residual": False,
            "W2_class_nonempty_p31": bool(B["certified"]),
            "W2_p_law": False,
            "walsh_general_p": False,
        },
        "A": A,
        "B": B,
        "C": C,
        "flags_not_flipped": [
            "residual_ii_k_eq_4p_empty",
            "phi_F_ge_6_proved_general",
            "e1",
            "L",
        ],
        "L_status": "OPEN",
        "walsh_15_406_E": "OPEN",
        "backend": "W1 prefix ProcessPool 86; W2 class ProcessPool 86; GPU unused",
        "fable": (
            "suggest_direction rank2: class existence, discriminator p=31.  "
            "Rank1 octic linear class function is empty in the tested box."
        ),
        "openai_referee": "math_review PASS; do_not_branch",
        "seconds": round(time.time() - t0, 3),
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15627.json"
    dest.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"  wrote {dest}", flush=True)
    return out


if __name__ == "__main__":
    main()
