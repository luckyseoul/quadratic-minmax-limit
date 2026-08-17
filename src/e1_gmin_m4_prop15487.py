#!/usr/bin/env python3
"""
Prop 15.487 — Paley-type masses on T are CM-named; the mixed
Paley type-indicator Fourier is not a name of Q++ (or of Q_N*).

    n_{+-}=n_{-+}=(p²−2p−3+a_E²)/8,
    n_{--}=((p+1)²−a_E²)/8,
    n_N1=p−1,   n_sub=p−3,   n_N*=(p−1)(p−3)/2.

I_{+-} pairing equals live Q++ at p=5 and Paley-N* at p=7.
phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.486 flags.

============================================================================
Setup.  T=□ ⊂ F_q^×, |T|=(q−1)/2.  15.478 named n_++ of Paley
(1,1).  a_E is the trace of E:y²=x³−x (only a_E² is used).
χ=χ_p∘N so N(r)=1 ⇒ χ(r)=1.  F_p^× ⊂ T.

============================================================================
Theorem A — PROVED (15.478 Jacobi + checksum; certified p=3..31).
  n_{+-}=n_{-+}=(p²−2p−3+a_E²)/8,
  n_{--}=((p+1)²−a_E²)/8,
  and 2+n_+++2 n_{+-}+n_{--}=(q−1)/2.
  Fail: drop a_E² in n_{+-} at p=5 (3/2 ≠ 2).  ∎

Theorem B — PROVED (χ=χ_p∘N + F_p^×⊂T; certified p=3..31).
  n_N1=p−1, n_sub=p−3, n_N*=(p−1)(p−3)/2.
  Fail: n_N1=(p+1)/2 at p=5 (3 ≠ 4);
  n_N*=(p−1)(p−3)/4 at p=7 (6 ≠ 12).  ∎

Theorem C — PROVED (type-indicator pairing; p=5,7 caches).
  I_{+-}(r)=1_{χ(r)=1, χ(r−1)=−1, χ(r+1)=+1}.
  n_{+-}^{−1} E[u(ξ) ∑_r I_{+-}(r) u(rξ)] equals 48/13 at p=5
  and 1496/409 at p=7 (Paley-N*, not Q++=1544/409).
  Fail: claim the pairing equals 48/13 at p=7;
  fail: claim it equals 32/13 at p=5.  ∎

Theorem D — OPEN.  Mixed Paley Fourier is not a name of Q_τ.
  Q++ / Q_N* unnamed as a rational in p.  phi_F not imported.

============================================================================
Backend: field counts + 15.478-style pairing.  Writes
evidence/e1_gmin_m4_prop15487.json
"""
from __future__ import annotations

import json
import os
import sys
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _e_tr, _kernel_field  # noqa: E402
from e1_gmin_m4_prop15290 import paley_pair, norm_type  # noqa: E402
from e1_gmin_m4_prop15298 import load_N  # noqa: E402
from e1_gmin_m4_prop15330 import LIVE_QPP  # noqa: E402
from e1_gmin_m4_prop15468 import D_rows  # noqa: E402
from e1_gmin_m4_prop15473 import e_matrix  # noqa: E402
from e1_gmin_m4_prop15474 import sigma_named  # noqa: E402
from e1_gmin_m4_prop15478 import a_E_sq, n_pp_named  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15487_catalog.json"

LIVE_QNSTAR = {5: Fraction(32, 13), 7: Fraction(1496, 409)}
NPP_CERT = (3, 5, 7, 11, 13, 17, 19, 23, 29, 31)


def n_pm_named(p: int) -> int:
    return (p * p - 2 * p - 3 + a_E_sq(p)) // 8


def n_pm_drop_aE(p: int) -> Fraction:
    return Fraction(p * p - 2 * p - 3, 8)


def n_mm_named(p: int) -> int:
    return ((p + 1) ** 2 - a_E_sq(p)) // 8


def n_N1_named(p: int) -> int:
    return p - 1


def n_sub_named(p: int) -> int:
    return p - 3


def n_Nstar_named(p: int) -> int:
    return (p - 1) * (p - 3) // 2


def n_Nstar_half(p: int) -> int:
    return (p - 1) * (p - 3) // 4


def I_pm(F, r: int) -> float:
    if r == 0 or F["chi"][r] != 1:
        return 0.0
    rm = F["add"][r][F["neg"][1]]
    rp = F["add"][r][1]
    if rm == 0 or rp == 0:
        return 0.0
    return 1.0 if (F["chi"][rm] == -1 and F["chi"][rp] == 1) else 0.0


def count_paley(F) -> dict[str, int]:
    pm = mm = pp = mp = p1 = 0
    for r in range(1, F["q"]):
        if F["chi"][r] != 1:
            continue
        if r in (1, F["neg1"]):
            p1 += 1
            continue
        a, b = paley_pair(F, r)
        if a == 1 and b == 1:
            pp += 1
        elif a == -1 and b == 1:
            pm += 1
        elif a == 1 and b == -1:
            mp += 1
        else:
            mm += 1
    return {"pm1": p1, "++": pp, "+-": pm, "-+": mp, "--": mm}


def count_norm(F) -> dict[str, int]:
    out = {"pm1": 0, "sub": 0, "N1": 0, "N*": 0}
    for r in range(1, F["q"]):
        if F["chi"][r] != 1:
            continue
        out[norm_type(F, r)] += 1
    return out


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in NPP_CERT:
        F = _kernel_field(p)
        c = count_paley(F)
        npm, nmm = n_pm_named(p), n_mm_named(p)
        npp = n_pp_named(p)
        if c["+-"] != npm or c["-+"] != npm or c["--"] != nmm:
            ok = False
        if c["++"] != npp:
            ok = False
        tot = c["pm1"] + c["++"] + c["+-"] + c["-+"] + c["--"]
        if tot != (p * p - 1) // 2:
            ok = False
        rows[str(p)] = {"+-": c["+-"], "named": npm, "--": c["--"]}
    if n_pm_drop_aE(5) == 2:
        ok = False
    if n_pm_drop_aE(5) != Fraction(3, 2):
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "n_{+-}=(p²−2p−3+a_E²)/8, n_{--}=((p+1)²−a_E²)/8. "
            "Fail: drop a_E² at p=5."
        ),
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in NPP_CERT:
        F = _kernel_field(p)
        c = count_norm(F)
        if c["N1"] != n_N1_named(p) or c["sub"] != n_sub_named(p):
            ok = False
        if c["N*"] != n_Nstar_named(p):
            ok = False
        rows[str(p)] = {"N1": c["N1"], "sub": c["sub"], "N*": c["N*"]}
    if n_N1_named(5) == (5 + 1) // 2:
        ok = False
    if n_Nstar_half(7) == n_Nstar_named(7):
        ok = False
    if n_Nstar_half(7) != 6 or n_Nstar_named(7) != 12:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "n_N1=p−1, n_sub=p−3, n_N*=(p−1)(p−3)/2. "
            "Fail: n_N1=(p+1)/2 at p=5; half n_N* at p=7."
        ),
    }


@lru_cache(maxsize=4)
def fold_pm(p: int):
    F = _kernel_field(p)
    q = F["q"]
    I = np.array([I_pm(F, r) for r in range(q)], dtype=np.float64)
    n = int(round(I.sum()))
    Ds = D_rows(p, F).astype(np.float64)
    u = 4.0 * np.abs(Ds @ e_matrix(F).T) ** 2
    sigma = sigma_named(p)
    Omega = np.array(F["chi"]) == sigma
    Omega[0] = False
    xi0 = int(np.where(Omega)[0][0])
    N = load_N(p, F).astype(np.float64)
    e_tab = np.array([_e_tr(F, t) for t in range(q)], dtype=np.complex128)
    mul = F["mul"]
    K = np.array([(I * e_tab[mul[a]]).sum() for a in range(q)], dtype=np.complex128)
    Kd = np.array([K[mul[xi0][d]] for d in range(q)])
    Iu = 4.0 * (N * Kd[None, :]).sum(axis=1)
    pred = float(np.mean((u[:, xi0] * Iu).real / n) / (q * q))
    return {"n": n, "pred": pred, "named_n": n_pm_named(p)}


def prove_C() -> dict:
    ok = True
    rows = {}
    rec5, rec7 = fold_pm(5), fold_pm(7)
    p5 = Fraction(rec5["pred"]).limit_denominator(20000)
    p7 = Fraction(rec7["pred"]).limit_denominator(20000)
    if p5 != LIVE_QPP[5]:
        ok = False
    if p7 != LIVE_QNSTAR[7]:
        ok = False
    if p7 == LIVE_QPP[7]:
        ok = False
    if p5 == LIVE_QNSTAR[5]:
        ok = False
    if rec5["n"] != 2 or rec7["n"] != 4:
        ok = False
    rows["5"] = {"pred": str(p5), "Qpp": str(LIVE_QPP[5]), "Qn": str(LIVE_QNSTAR[5])}
    rows["7"] = {"pred": str(p7), "Qpp": str(LIVE_QPP[7]), "Qn": str(LIVE_QNSTAR[7])}
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "I_{+-} pairing is 48/13 at p=5 and 1496/409 at p=7. "
            "Fail: =48/13 at p=7; =32/13 at p=5."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "Paley masses named. Mixed Paley Fourier is not Q++. "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.487  Paley masses + mixed type-indicator", flush=True)
    A = prove_A()
    print(f"  A masses: {A['proved']} p5+−={A['rows']['5']}", flush=True)
    B = prove_B()
    print(f"  B norm: {B['proved']} p7 N*={B['rows']['7']['N*']}", flush=True)
    C = prove_C()
    print(f"  C pairing: {C['proved']} {C['rows']}", flush=True)
    D = prove_open()
    print(f"  D open phi_F={D['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps({"A": A["proved"], "B": B["proved"], "C": C["proved"]}, indent=2)
        + "\n"
    )
    out = {
        "prop": "15.487",
        "title": (
            "Paley-type masses CM-named; mixed type-indicator "
            "is not a name of Q++"
        ),
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "paley_masses_CM": A["proved"],
            "norm_masses_named": B["proved"],
            "mixed_I_not_Qpp": C["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15487.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
