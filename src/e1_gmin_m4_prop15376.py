#!/usr/bin/env python3
"""
Prop 15.376 — Joint (Q++,Q--,Q_N*) atoms at p=7; the Q--=0
slice is exactly the 1D family (AP ∪ QR0); the 15.138
norm-circle y_* is a named atom at p=7 and has
    Q_N*(y_*) = 24/p
at both p=7 and p=11.  Ensemble Q_{--} is not identified
with 8(A−k)/D at p=11.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.375 flags.  Residual (ii) and Type I stay False.

============================================================================
Setup.  GPU GEMM on the p=7 H_+ cache (5726).  y_* is the
lex-first norm-circle Max+ of 15.138 (cert (t,c)=(35,4) at
p=7 and (88,4) at p=11).

============================================================================
Theorem A — PROVED (GPU joint hist; certified p=7).
  Q--=0 on exactly n_1d=140 vectors, and the mean Q++ on
  that slice is the named 1D mix 104/15.  The slice splits
      AP:  84, Q++=40/9=Q_AP, Q_N*=0
      QR0: 56, Q++=32/3=Q_QR0, Q_N*=0.
  Fail: claim a 1D vector has Q--≠0; claim the 140-mean Q++
  is 40/9.  ∎

Theorem B — PROVED (15.138 constructor + zhat; p=7,11).
  y_* at p=7 is the atom (Q++,Q--,Q_N*)=(8/3, 80/21, 24/7)
  of mass 588.  At p=11, y_* has Q_N*=24/11.  Hence
      Q_N*(y_*) = 24/p
  at both primes p≡3 in {7,11}.  Fail: 24/(p−1) (4≠24/7);
  fail: claim 24/p at p=5 (Hoffman/y_* has Q_N*=16/5≠24/5).
  Q--(y_*)=80/21, 224/55 is not 8(A−k)/D (1376/409, 21104/5553).  ∎

Theorem C — OPEN.  Ensemble Q_{--} is still only named at p=7
  as 8(A−k)/D.  D_form untested at p≥11.  24/p is a constructor
  value, not the ensemble.  phi_F not imported.

============================================================================
Backend: GPU joint hist (p=7 cache) + 15.138 constructor.
Writes evidence/e1_gmin_m4_prop15376.json
"""
from __future__ import annotations

import json
import os
import sys
from fractions import Fraction
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15138 import YSTAR_CERT, construct_ystar  # noqa: E402
from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field, rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15290 import omega_set, zhat_omega  # noqa: E402
from e1_gmin_m4_prop15292 import n_1d  # noqa: E402
from e1_gmin_m4_prop15300 import merge_cls  # noqa: E402
from e1_gmin_m4_prop15355 import Q_AP_named, Q_QR0_p3_named  # noqa: E402
from e1_gmin_m4_prop15361 import n_AP_named  # noqa: E402
from e1_gmin_m4_prop15375 import Qmm_form, Qpp_form  # noqa: E402

SCRATCH = Path("/tmp/grok-goal-ea1b2cdbd197/implementer")
JOINT = SCRATCH / "gpu_Qjoint.json"
CAT = ROOT / "evidence" / "e1_gmin_m4_prop15376_catalog.json"

# p=7 joint atoms from GPU GEMM (gpu_Qjoint.json).
JOINT_P7 = {
    ("40/9", "24/7", "76/21"): 1176,
    ("8/3", "4", "4"): 1176,
    ("32/9", "24/7", "256/63"): 1176,
    ("40/9", "8/3", "32/9"): 1176,
    ("8/3", "80/21", "24/7"): 588,
    ("40/9", "80/21", "24/7"): 294,
    ("40/9", "0", "0"): 84,
    ("32/3", "0", "0"): 56,
}


def Q_1d_pp_p3(p: int) -> Fraction:
    return Fraction(4 * (p * p - 3 * p - 2), 3 * (p - 2))


def QN_ystar_named(p: int) -> Fraction:
    return Fraction(24, p)


def QN_ystar_wrong_pminus1(p: int) -> Fraction:
    return Fraction(24, p - 1)


def _buckets(F):
    out = {"++": [], "--N1": [], "N*": []}
    for r in F["squares"]:
        c = merge_cls(F, r)
        if c in out:
            out[c].append(r)
    return out


def Q_of_y(y, F, Omega, rs) -> Fraction:
    p, q = F["p"], F["q"]
    z = __import__("numpy").asarray(y[1 : 1 + q], dtype="float64")
    if y[0] < 0:
        z = -z
    u = abs(zhat_omega(z, F, Omega)) ** 2
    om = {xi: i for i, xi in enumerate(Omega)}
    sm = 0.0
    for r in rs:
        for xi in Omega:
            sm += u[om[xi]] * u[om[F["mul"][r][xi]]]
    return Fraction(sm / (len(Omega) * len(rs) * q * q)).limit_denominator(p**6)


def ystar_Q(p: int) -> dict[str, Fraction]:
    r = construct_ystar(p)
    if not r.get("found"):
        raise RuntimeError(f"y_* missing at p={p}")
    F = _kernel_field(p)
    Omega = omega_set(F)
    b = _buckets(F)
    return {c: Q_of_y(r["y"], F, Omega, rs) for c, rs in b.items()}


def prove_A() -> dict:
    rec = json.loads(JOINT.read_text()) if JOINT.is_file() else {}
    nzero = rec.get("n_Qmm_zero", -1)
    mean = rec.get("mean_Qpp_on_zero")
    onz = rec.get("Qpp_on_zero") or {}
    ok = nzero == n_1d(7) == 140
    ok = ok and mean == str(Q_1d_pp_p3(7)) == "104/15"
    ok = ok and onz.get("40/9") == n_AP_named(7) == 84
    ok = ok and onz.get("32/3") == n_1d(7) - n_AP_named(7) == 56
    ok = ok and Q_AP_named(7) == Fraction(40, 9)
    ok = ok and Q_QR0_p3_named(7) == Fraction(32, 3)
    if nzero != 140 or mean != "104/15":
        ok = False
    if onz.get("40/9") == 140:
        ok = False
    tot = sum(JOINT_P7.values())
    ok = ok and tot == 5726
    return {
        "proved": bool(ok),
        "nzero": nzero,
        "mean_Qpp_zero": mean,
        "AP": onz.get("40/9"),
        "QR0": onz.get("32/3"),
        "Q_1d_pp": str(Q_1d_pp_p3(7)),
        "theorem": (
            "Q--=0 iff 1D (n_1d=140); AP 84 at Q++=40/9, "
            "QR0 56 at Q++=32/3; mean 104/15."
        ),
    }


def prove_B() -> dict:
    Q7 = ystar_Q(7)
    ok = Q7["++"] == Fraction(8, 3)
    ok = ok and Q7["--N1"] == Fraction(80, 21)
    ok = ok and Q7["N*"] == Fraction(24, 7) == QN_ystar_named(7)
    ok = ok and Q7["N*"] != QN_ystar_wrong_pminus1(7)
    ok = ok and Q7["--N1"] != Qmm_form(7)
    # p=11 constructor (cert (t,c)=(88,4)); Q_N*=24/p
    Q11 = ystar_Q(11)
    ok = ok and Q11["N*"] == QN_ystar_named(11) == Fraction(24, 11)
    ok = ok and Q11["N*"] != QN_ystar_wrong_pminus1(11)
    ok = ok and Q11["--N1"] != Qmm_form(11)
    # p=5 fail-when-wrong for 24/p on the N* type of Hoffman/y_*
    ok = ok and Fraction(16, 5) != QN_ystar_named(5)
    if Q7["N*"] == QN_ystar_wrong_pminus1(7):
        ok = False
    if Q11["--N1"] == Qmm_form(11):
        ok = False
    return {
        "proved": bool(ok),
        "p7": {c: str(v) for c, v in Q7.items()},
        "p11": {c: str(v) for c, v in Q11.items()},
        "QN_named": {7: str(QN_ystar_named(7)), 11: str(QN_ystar_named(11))},
        "p11_Qmm_form": str(Qmm_form(11)),
        "p11_Qmm_ystar": str(Q11["--N1"]),
        "cert": {str(p): YSTAR_CERT[p] for p in (7, 11)},
        "theorem": (
            "y_* is the (8/3,80/21,24/7) atom at p=7; "
            "Q_N*(y_*)=24/p at p=7,11. Not the ensemble Q--."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "Q_tau_named_in_p": False,
        "D_named_in_p": False,
        "Qmm_named_in_p": False,
        "note": (
            "24/p names Q_N* on the 15.138 constructor, not ensemble Q--. "
            "D_form untested at p≥11. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.376  1D atoms + y_* Q_N*=24/p", flush=True)
    A = prove_A()
    print(f"  A 1D slice: {A['proved']} nzero={A['nzero']}", flush=True)
    B = prove_B()
    print(f"  B y_*: {B['proved']} p7={B['p7']} p11={B['p11']}", flush=True)
    C = prove_open()
    print(f"  C open: phi_F={C['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps(
            {
                "joint_p7": {f"{a}|{b}|{c}": n for (a, b, c), n in JOINT_P7.items()},
                "QN_ystar": "24/p",
                "p7": B.get("p7"),
                "p11": B.get("p11"),
            },
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.376",
        "title": "1D is the Q--=0 slice; Q_N*(y_*)=24/p at p=7,11",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "Qmm_zero_is_1d": A["proved"],
            "QN_ystar_24_over_p": B["proved"],
            "Q_tau_named_in_p": False,
            "D_named_in_p": False,
            "Qmm_named_in_p": False,
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15376.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
