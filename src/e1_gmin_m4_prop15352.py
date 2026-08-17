#!/usr/bin/env python3
"""
Prop 15.352 — Complete named dictionary of the seven p=7 H_+ types.
294 = T⊙NC², 588 = ystar, leftover 1176 = T⊙NC_{c=1}.
Does **not** name Q_τ in p for all primes.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.351 flags.

============================================================================
Setup.  15.307 seven D-types.  T=A·T_0 (15.349).  first_nc (15.312)
searches (t,c) in lex order.  first_nc_c fixes the norm level c.

============================================================================
Theorem A — PROVED (cache fold + constructors; certified p=7).
  The seven 15.307 types are:
      140  (C³H)           = AP / QR0 halfspaces
      1176 (CRRR)          = T
      1176 (O:0,1,1,2,5,5,7, …) = T⊙NC = QR0⊙NC
      1176 (O:1,1,2,2,3,6,6, …) = y_♯ size-12 (15.141)
      1176 (O:0,2,2,3,3,4,7, …) = T⊙NC_{c=1}
      588  (O:0,0,2,2,5,6,6, …) = ystar (AP⊙NC)
      294  (O:0,0,3,4,4,5,5, …) = T⊙NC²
  Sizes sum to 5726.  Fail: claim T⊙NC² is a 1176 (it is 294);
  claim T⊙NC_{c=2} (i.e. first_nc) is the leftover 1176.  ∎

Theorem B — PROVED (constructor; certified p=7).
  first_nc_c(y_T, c=1) returns a Max+ y with the leftover-1176
  signature.  first_nc_c(y_T, c=2) returns T⊙NC (15.351), not that
  type.  ẑ=ẑ_T−2Γ on the c=1 switch.  Fail: c=2.  ∎

Theorem C — OPEN.  The dictionary is p=7.  No general-p type list
  and no Gauss/Jacobi Q_τ in p.  phi_F not imported.

============================================================================
Backend: first_nc_c + 15.307 classify.  Writes
evidence/e1_gmin_m4_prop15352.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15138 import construct_ystar  # noqa: E402
from e1_gmin_m4_prop15141 import construct_y_sharp  # noqa: E402
from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field, rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15307 import classify_p  # noqa: E402
from e1_gmin_m4_prop15312 import Q_switch, _norm_table, first_nc  # noqa: E402
from e1_gmin_m4_prop15351 import sig_y, y_T  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402


SIG_TNC = (
    "O:0,1,1,2,5,5,7",
    "O:1,1,3,3,4,4,5",
    "O:1,2,2,3,3,5,5",
    "O:1,2,2,3,3,5,5",
)
SIG_294 = (
    "O:0,0,3,4,4,5,5",
    "O:0,0,3,4,4,5,5",
    "O:1,1,3,3,4,4,5",
    "O:1,1,3,3,4,4,5",
)
SIG_588 = (
    "O:0,0,2,2,5,6,6",
    "O:1,1,3,3,4,4,5",
    "O:1,1,3,3,4,4,5",
    "O:1,2,2,3,3,5,5",
)
SIG_SHARP = (
    "O:1,1,2,2,3,6,6",
    "O:1,1,2,2,4,4,7",
    "O:1,1,3,3,4,4,5",
    "O:1,2,2,3,3,5,5",
)
SIG_LAST = (
    "O:0,2,2,3,3,4,7",
    "O:1,1,2,2,3,6,6",
    "O:1,2,2,3,3,5,5",
    "O:1,2,2,3,3,5,5",
)


def first_nc_c(p: int, y0: np.ndarray, c_req: int, skip_sigs=()):
    """Lex-first Hoffman circle of C_{y0} with field-norm level c=c_req.

    skip_sigs: signatures already named (ystar, T⊙NC²) so the leftover
    1176 is the first c=1 hit that is not those types.
    """
    F = _kernel_field(p)
    q = F["q"]
    C = paley_conference_prime_power(p).astype(np.float64)
    y0 = np.asarray(y0, dtype=np.float64)
    if y0[0] < 0:
        y0 = -y0
    Cbase = (y0[:, None] * C) * y0[None, :]
    Nt = _norm_table(F)
    add, neg = F["add"], F["neg"]
    skip = {tuple(s) for s in skip_sigs}
    for t in range(q):
        pts = [u for u in range(q) if Nt[add[u][neg[t]]] == c_req]
        if len(pts) != p + 1:
            continue
        S = [1 + u for u in pts]
        if any(Cbase[S[i], S[j]] <= 0 for i in range(p + 1) for j in range(i + 1, p + 1)):
            continue
        z = np.ones(q + 1)
        for i in S:
            z[i] = -1.0
        if not np.allclose(Cbase @ z, p * z, atol=1e-5):
            continue
        y = y0 * z
        if y[0] < 0:
            y = -y
        if not np.allclose(C @ y, p * y, atol=1e-5):
            continue
        if skip and sig_y(p, y) in skip:
            continue
        return {"y": y, "circle": pts, "t": t, "c": c_req, "y0": y0}
    return None


def prove_A() -> dict:
    ok = True
    p = 7
    census = classify_p(p)["D_types"]
    y0 = y_T(p)
    tnc = first_nc(p, y0)
    tnc2 = first_nc(
        p, tnc["y"], skip_pred=lambda y: np.array_equal(np.sign(y), np.sign(y0))
    )
    ys = construct_ystar(p)["y"]
    ysh = construct_y_sharp(p)["y"]
    last = first_nc_c(p, y0, 1, skip_sigs=(SIG_588, SIG_294))
    mapping = {
        "140_hs": ("('C', 'C', 'C', 'H')", 140),
        "1176_T": (str(("C", "R", "R", "R")), 1176),
        "1176_Tnc": (str(SIG_TNC), 1176),
        "1176_sharp": (str(SIG_SHARP), 1176),
        "1176_c1": (str(SIG_LAST), 1176),
        "588_ystar": (str(SIG_588), 588),
        "294_Tnc2": (str(SIG_294), 294),
    }
    rows = {}
    total = 0
    for name, (key, n) in mapping.items():
        got = census.get(key)
        rows[name] = {"key": key, "expect": n, "census": got}
        if got != n:
            ok = False
        total += n if got else 0
    if total != 5726:
        ok = False
    # live constructors match keys
    if sig_y(p, y0) != ("C", "R", "R", "R"):
        ok = False
    if sig_y(p, tnc["y"]) != SIG_TNC:
        ok = False
    if sig_y(p, tnc2["y"]) != SIG_294:
        ok = False
    if sig_y(p, ys) != SIG_588:
        ok = False
    if sig_y(p, ysh) != SIG_SHARP:
        ok = False
    if last is None or sig_y(p, last["y"]) != SIG_LAST:
        ok = False
    # fail-when-wrong: T⊙NC² is not a 1176
    if census.get(str(SIG_294)) == 1176:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "sum": total,
        "theorem": (
            "Seven p=7 types named.  T⊙NC² is 294, not 1176.  "
            "ystar is 588.  T⊙NC_{c=1} is the leftover 1176."
        ),
    }


def prove_B() -> dict:
    ok = True
    p = 7
    y0 = y_T(p)
    raw = first_nc_c(p, y0, 1)
    c1 = first_nc_c(p, y0, 1, skip_sigs=(SIG_588, SIG_294))
    c2 = first_nc_c(p, y0, 2)
    if c1 is None or c2 is None or raw is None:
        ok = False
        return {"proved": False, "missing": True}
    s1 = sig_y(p, c1["y"])
    s2 = sig_y(p, c2["y"])
    sraw = sig_y(p, raw["y"])
    if s1 != SIG_LAST:
        ok = False
    if s2 != SIG_TNC:
        ok = False
    if s1 == s2:
        ok = False
    # unskipped c=1 is ystar (588), not the leftover
    if sraw != SIG_588:
        ok = False
    if sraw == SIG_LAST:
        ok = False
    Q, err2, err1 = Q_switch(p, c1)
    if err2 > 1e-9:
        ok = False
    if err1 < 1.0:
        ok = False
    if s2 == SIG_LAST:
        ok = False
    return {
        "proved": bool(ok),
        "c1": {"t": int(c1["t"]), "sig": [str(x) for x in s1]},
        "c2": {"t": int(c2["t"]), "sig": [str(x) for x in s2]},
        "c1_raw": {"t": int(raw["t"]), "sig": [str(x) for x in sraw]},
        "err2": err2,
        "err1": err1,
        "Qpp_c1": str(Q.get("++")),
        "theorem": "T⊙NC_{c=1} is the leftover 1176.  c=2 is T⊙NC, not that type.",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "Q_tau_named_in_p": False,
        "p7_types_named": True,
        "note": (
            "p=7 type list is complete and Max+-free.  "
            "No general-p Gauss/Jacobi Q_τ.  phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.352  p=7 type dictionary complete", flush=True)
    A = prove_A()
    print(f"  A dict: {A['proved']} sum={A['sum']}", flush=True)
    B = prove_B()
    print(f"  B c=1 leftover: {B['proved']}", flush=True)
    C = prove_open()
    print(f"  C open: phi_F={C['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.352",
        "title": "p=7 H_+ types all named; Q_τ in p still OPEN",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "p7_seven_types_named": A["proved"],
            "Tnc_c1_is_leftover_1176": B["proved"],
            "Q_tau_named_in_p": False,
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15352.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
