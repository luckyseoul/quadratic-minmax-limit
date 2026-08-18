#!/usr/bin/env python3
"""
Prop 15.532 — T_v acts freely on H_+ iff χ_q(v)=−1.
Fix(T_v)=C(p,m) on square lines (Type+ 1D with that ker)
and 0 on nonsquare lines (Type− ⊄ H_+, 15.451 C).
Hence |H_+|=p·n_orb(T_ns).  D=n_orb/2 still unnamed.
Q_τ unnamed.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.531 flags.  Does **not** overwrite 15.522–15.531.

============================================================================
Setup.  T_v is translation by a generator of the F_p-line ⟨v⟩⊂F_q.
y is T_v-fixed iff y(x+v)=y(x), iff y is a 1D lift with ker=⟨v⟩.
Those lifts are Type+ ⇔ ⟨v⟩ is a square line ⇔ χ_q(v)=+1, and
every Type+ 1D lies in H_+ (15.292 / 15.526).  Type− ⊄ H_+
(15.451 C).  Per square line there are C(p,m) such lifts,
m=(p+1)/2, and n_1d=m C(p,m) with m square lines.

============================================================================
Theorem A — PROVED (15.292+15.451 C+15.526; certified p=5,7).
  Fix_H(T_v)=C(p,m) if χ(v)=1, else 0.  Fail: T_1 is free
  (χ(1)=1 and Fix=10,35).  ∎

Theorem B — PROVED (orbit-stabilizer on free T_ns).
  For χ(v)=−1, T_v acts freely, so |H_+|=p n_orb(T_v).
  Live n_orb=26,818=2D.  Fail: n_orb=D.  ∎

Theorem C — OPEN.  The pairing n_orb=2D is not yet a named
  involution with 0 invariant orbits (Frob/neg both fix some
  T_ns-orbits).  D and Q_τ unnamed.  phi_F not imported.

============================================================================
Backend: H_+ caches for Fix counts.  Writes
evidence/e1_gmin_m4_prop15532.json
"""
from __future__ import annotations

import json
from math import comb
from pathlib import Path

import numpy as np

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import _kernel_field, is_type_plus_form
from e1_gmin_m4_prop15292 import n_1d
from e1_gmin_m4_prop15330 import HPLUS
from e1_gmin_m4_prop15522 import D_named
from e1_gmin_m4_prop15526 import load_Hplus

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "evidence" / "e1_gmin_m4_prop15532.json"

# One square and one nonsquare direction (field element a+bω).
DIRS = {
    5: {"sq": (1, 0), "ns": (1, 2)},
    7: {"sq": (1, 0), "ns": (1, 1)},
}


def m_of(p: int) -> int:
    return (p + 1) // 2


def Fix_named(p: int, chi_v: int) -> int:
    return comb(p, m_of(p)) if chi_v == 1 else 0


def perm_T(p: int, da: int, db: int) -> list[int]:
    q = p * p
    return [((x % p + da) % p) + ((x // p + db) % p) * p for x in range(q)]


def count_fix(H: np.ndarray, perm: list[int]) -> int:
    q = len(perm)
    n = 0
    for i in range(len(H)):
        y = H[i]
        ok = True
        for x in range(q):
            if y[1 + perm[x]] != y[1 + x]:
                ok = False
                break
        if ok:
            n += 1
    return n


def chi_dir(F, a: int, b: int) -> int:
    return int(F["chi"][a + b * F["p"]])


def prove_A() -> dict:
    rows = {}
    ok = True
    for p in (5, 7):
        F = _kernel_field(p)
        H = load_Hplus(p)
        rec = {"n_1d": n_1d(p), "Cpm": comb(p, m_of(p))}
        for name, (a, b) in DIRS[p].items():
            ch = chi_dir(F, a, b)
            live = count_fix(H, perm_T(p, a, b))
            named = Fix_named(p, ch)
            rec[name] = {
                "dir": [a, b],
                "chi": ch,
                "type_plus": bool(is_type_plus_form(p, a, b)),
                "live_fix": live,
                "named_fix": named,
                "ok": live == named,
            }
            if live != named:
                ok = False
        # fail T_1 free
        if rec["sq"]["live_fix"] == 0 or rec["ns"]["live_fix"] != 0:
            ok = False
        if rec["sq"]["chi"] != 1 or rec["ns"]["chi"] != -1:
            ok = False
        if n_1d(p) != m_of(p) * rec["Cpm"]:
            ok = False
        rows[str(p)] = rec
        rows[str(p)]["ok"] = ok
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "Fix(T_v)=C(p,m) on squares, 0 on nonsquares. "
            "Fail T_1 free."
        ),
    }


def prove_B() -> dict:
    A = prove_A()
    rows = {}
    ok = A["proved"]
    for p in (5, 7):
        nH = HPLUS[p]
        n_orb = nH // p
        rec = {
            "nH": nH,
            "n_orb": n_orb,
            "D": D_named(p),
            "two_D": 2 * D_named(p),
            "ok": n_orb == 2 * D_named(p) and n_orb != D_named(p),
        }
        if not rec["ok"]:
            ok = False
        rows[str(p)] = rec
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": "|H_+|=p n_orb(T_ns) and n_orb=2D live. Fail n_orb=D.",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "D_named_as_orbits": False,
        "Q_tau_named_in_p": False,
        "phi_F_imported": False,
        "note": (
            "n_orb=2D is live, not a named involution pairing. "
            "Q_τ unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("15.532 T_v free on H_+ iff chi(v)=-1", flush=True)
    A, B, C = prove_A(), prove_B(), prove_open()
    out = {
        "prop": "15.532",
        "title": "T_v free on H_+ iff χ(v)=−1; |H_+|=p n_orb(T_ns)",
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
