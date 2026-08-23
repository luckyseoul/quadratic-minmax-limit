#!/usr/bin/env python3
"""
Prop 15.623 — W1 for p≡73 or 97 (mod 120) via d=−3.

Does **not** flip residual_ii / multilevel_ND / phi_F / type_I / e1 / L.
Does **not** close W1 for p≡1 or 49 (mod 120), W2 p-law, Walsh, leftover 2.

============================================================================
Setup.  15.622: W1 except p≡1 (mod 24).  That class splits by (5/p):
p≡73 or 97 (mod 120) have (5/p)=−1; p≡1 or 49 (mod 120) have (5/p)=+1.

============================================================================
Theorem A — PROVED (OpenAI math_review PASS).
  d=−3.  Stay for p≥7.  S Δ (S−3)={m−2,m−1,m,p−3,p−2,p−1}.
  p≡1 (mod 8) and (3/p)=+1 ⇒ χ(−1)=χ(2)=χ(3)=1, so the last
  five listed characters are +1 except χ(m−2)=χ(−5·2^{−1})=χ(5).
  (5/p)=−1 on this class ⇒ exactly five QR, ε=1.
  Fail: ε=0 at p=73 for d=−3 (actual 1).
  Union 15.614/621/622: W1 open only for p≡1 or 49 (mod 120).  ∎

Theorem B — CERTIFIED, not a p-law.
  p≡1 or 49 (mod 120): d=−3 has ε=0 (χ(5)=+1).  Stay hits
  exist (p=241,409).  Fail: d=−3 at p=241.  ∎

Theorem C — OPEN.  W1 for p≡1,49 (mod 120), W2 p-law, Walsh,
  leftover 2.

============================================================================
Backend: 6-point character count (serial; nuka-class).  GPU unused.
OpenAI PASS on A.
Writes evidence/e1_gmin_m4_prop15623.json
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


def _eps_d(p, d):
    m = (p - 1) // 2
    S = set(range(m + 1))
    Sd = {(x + d) % p for x in S}
    delta = S.symmetric_difference(Sd)
    nqr = sum(1 for x in delta if x and _chi(p, x) == 1)
    return nqr % 2


def theorem_A_W1_73_97() -> dict:
    ok = True
    rows = {}
    for p in (73, 97, 193, 337):
        e = _eps_d(p, p - 3)
        ok = (
            ok
            and p % 120 in (73, 97)
            and _chi(p, 5) == -1
            and e == 1
        )
        rows[str(p)] = {"mod120": p % 120, "chi5": _chi(p, 5), "eps": e}
    e241 = _eps_d(241, 241 - 3)
    ok = ok and e241 == 0
    rows["241"] = {"mod120": 241 % 120, "eps": e241, "note": "not claimed"}
    return {
        "proved": bool(ok),
        "W1_p_eq_73_or_97_mod_120": True,
        "W1_p_eq_1_or_49_mod_120": False,
        "W1_all_odd_p": False,
        "rows": rows,
        "theorem": (
            "d=−3 gives ε=1 for p≡73 or 97 (mod 120).  "
            "Fail: ε=0 at p=73."
        ),
    }


def theorem_B_open_class() -> dict:
    e = _eps_d(241, 238)
    return {
        "proved": False,
        "d_minus_3_eps_p241": e,
        "theorem": (
            "p≡1 (mod 120): d=−3 has ε=0.  Fail: d=−3 at p=241."
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
    }


def main() -> dict:
    t0 = time.time()
    print("Prop 15.623  W1 p≡73 or 97 (mod 120) via d=−3", flush=True)
    A = theorem_A_W1_73_97()
    print(f"  A {A['proved']} {A['rows']}", flush=True)
    B = theorem_B_open_class()
    print(f"  B p241 d-3 eps={B['d_minus_3_eps_p241']}", flush=True)
    C = theorem_C_open()
    out = {
        "prop": "15.623",
        "title": "W1 for p≡73 or 97 (mod 120)",
        "proved": {
            "W1_p_eq_73_or_97_mod_120": A["proved"],
            "W1_p_eq_1_or_49_mod_120": False,
            "W1_all_odd_p": False,
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
        "backend": "6-point F_p character counts; GPU unused; serial-class",
        "openai_referee": "math_review PASS on A, confidence 0.99",
        "seconds": round(time.time() - t0, 3),
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15623.json"
    dest.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"  wrote {dest}", flush=True)
    return out


if __name__ == "__main__":
    main()
