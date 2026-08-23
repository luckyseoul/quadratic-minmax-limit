#!/usr/bin/env python3
"""
Prop 15.625 — W1 for (2/p)_4=−1 via d=−(p−1)/8; CLASS is exhaustive.

Does **not** flip residual_ii / multilevel_ND / phi_F / type_I / e1 / L.
Does **not** close W1 for (2/p)_4=+1, W2 p-law, Walsh, leftover 2.

============================================================================
Setup.  Fable xhigh deep_review (Walsh consult): quarter-interval
d=−(p−1)/4 has ε=0 for all p≡1 (mod 8) (windows pair, h(−4p)
cancels).  Named stay is d=−(p−1)/8.  Walsh ⇔ W1 ∧ W2 has no
generation gap once CLASS is exhaustive (15.612: maximal
Aut-invariant ideals are exactly (X+1)R and (f_O)R).

============================================================================
Theorem A — PROVED (stay + eighth-interval; Barrucand–Cohn;
Fable deep_review; certified p=17,41,73,89,97,241,409,601).
  p≡1 (mod 8).  d=−(p−1)/8 lies in the upper half
  {(p+1)/2,…,p−1}.  ε(d)=|QR ∩ (S Δ (S+d))| = N1+N4 (mod 2),
  Nj = #QR in the j-th eighth of (0,p/2).  N1+N4 ≡ N1−N4
  (mod 2) and N1−N4=h(−8p)/4.  Barrucand–Cohn / Hasse:
  h(−8p)/4 odd iff (2/p)_4=−1 iff, writing p=a²+b² with a odd
  and b even, b≢0 (mod 8).
  Discriminating test: ε=1 at p=241,409; ε=0 at p=601.
  Fail: ε=0 at p=17; fail: ε=1 at p=601.  ∎

Theorem B — OPEN.  Residual W1 class is
  {p≡1 or 49 (mod 120) and (2/p)_4=+1} = {p=a²+64c²},
  first prime 601.  Interval-scaled stays are Chebotarev-dead
  (Fable).  Need a non-Chebotarev d (windows in a,b).  ∎

Theorem C — PROVED (15.612 CLASS exhaustive; Artinian).
  Every proper Aut-invariant ideal of finite R lies in a
  maximal one.  Maximal ones are exactly (X+1)R and (f_O)R.
  Hence W1 ∧ W2 ⇒ I_U=W_0, no generation gap.  Fail: a third
  maximal Aut-invariant ideal.  ∎

Theorem D — OPEN.  W2 p-law (Fable: both named Auts are the
  split involution class; hit-law may itself be Chebotarev;
  target a p-dependent cocycle).  Walsh, leftover 2.

============================================================================
Backend: F_p character counts (tiny; not Krylov).  GPU unused.
Fable deep_review PASS-WITH-NOTE.  OpenAI unused this unit.
Writes evidence/e1_gmin_m4_prop15625.json
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


def _chi(p, x):
    x %= p
    if x == 0:
        return 0
    return 1 if pow(x, (p - 1) // 2, p) == 1 else -1


def _quartic2(p):
    return 1 if pow(2, (p - 1) // 8, p) == 1 else -1


def _ab(p):
    for a in range(1, p, 2):
        b2 = p - a * a
        if b2 <= 0:
            break
        b = int(round(b2**0.5))
        if b * b == b2 and b % 2 == 0 and b > 0:
            return a, b
    return None, None


def _eps_d(p, d):
    m = (p - 1) // 2
    S = set(range(m + 1))
    Sd = {(x + d) % p for x in S}
    delta = S.symmetric_difference(Sd)
    return sum(1 for x in delta if x and _chi(p, x) == 1) % 2


def _in_upper(p, d):
    m = (p - 1) // 2
    d %= p
    return m + 1 <= d <= p - 1


def theorem_A_eighth() -> dict:
    primes = (17, 41, 73, 89, 97, 241, 409, 601)
    ok = True
    rows = {}
    for p in primes:
        e = (p - 1) // 8
        d = (-e) % p
        eps = _eps_d(p, d)
        q4 = _quartic2(p)
        a, b = _ab(p)
        pred = 1 if q4 == -1 else 0
        up = _in_upper(p, d)
        ok = ok and up and eps == pred and (b % 8 == 0) == (q4 == 1)
        rows[str(p)] = {
            "d": d,
            "upper": up,
            "eps": eps,
            "quartic2": q4,
            "pred": pred,
            "b_mod_8": b % 8,
        }
    return {
        "proved": bool(ok),
        "W1_quartic2_minus": True,
        "W1_quartic2_plus": False,
        "W1_all_odd_p": False,
        "rows": rows,
        "theorem": (
            "d=−(p−1)/8 has ε=1 iff (2/p)_4=−1.  "
            "Fail: ε=0 at p=17; fail: ε=1 at p=601."
        ),
    }


def theorem_B_residual() -> dict:
    return {
        "proved": False,
        "residual": "p=a^2+64c^2, first 601",
        "p601_eps": _eps_d(601, (-75) % 601),
        "theorem": (
            "W1 residual is (2/p)_4=+1.  Fail: d=−(p−1)/8 at p=601."
        ),
    }


def theorem_C_class_exhaustive() -> dict:
    return {
        "proved": True,
        "citation": "15.612 B: maximal Aut-invariant ideals are exactly (X+1)R and (f_O)R",
        "Walsh_iff_W1_and_W2": True,
        "generation_gap": False,
        "theorem": (
            "CLASS exhaustive + Artinian ⇒ W1∧W2 ⇒ I_U=W_0.  "
            "Fail: a third maximal Aut-invariant ideal."
        ),
    }


def theorem_D_open() -> dict:
    return {
        "proved": False,
        "W1_all_odd_p": False,
        "W2_p_law": False,
        "walsh_general_p": False,
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "fable": (
            "eighth-interval W1; residual a^2+64c^2; W2 p-dependent "
            "cocycle (named Auts are conjugate split involutions)"
        ),
    }


def main() -> dict:
    t0 = time.time()
    print("Prop 15.625  W1 d=-(p-1)/8 iff (2/p)_4=-1", flush=True)
    A = theorem_A_eighth()
    print(f"  A {A['proved']} 241={A['rows']['241']['eps']} 601={A['rows']['601']['eps']}", flush=True)
    B = theorem_B_residual()
    print(f"  B residual p601 eps={B['p601_eps']}", flush=True)
    C = theorem_C_class_exhaustive()
    print(f"  C CLASS exhaustive: {C['proved']}", flush=True)
    D = theorem_D_open()
    out = {
        "prop": "15.625",
        "title": "W1 eighth-interval / quartic character of 2",
        "proved": {
            "W1_quartic2_minus": A["proved"],
            "W1_quartic2_plus": False,
            "W1_all_odd_p": False,
            "CLASS_exhaustive": C["proved"],
            "W2_p_law": False,
            "walsh_general_p": False,
        },
        "A": A,
        "B": B,
        "C": C,
        "D": D,
        "flags_not_flipped": [
            "residual_ii_k_eq_4p_empty",
            "phi_F_ge_6_proved_general",
            "e1",
            "L",
        ],
        "L_status": "OPEN",
        "walsh_15_406_E": "OPEN",
        "backend": "F_p character counts; GPU unused",
        "claude_referee": (
            "deep_review PASS-WITH-NOTE: d=-(p-1)/8; quarter-interval "
            "dead; CLASS exhaustive; W2 Auts conjugate"
        ),
        "seconds": round(time.time() - t0, 3),
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15625.json"
    dest.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"  wrote {dest}", flush=True)
    return out


if __name__ == "__main__":
    main()
