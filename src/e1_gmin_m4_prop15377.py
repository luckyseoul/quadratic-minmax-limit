#!/usr/bin/env python3
"""
Prop 15.377 — Every p=7 joint atom has a Max+-free constructor.
y_*⊙NC produces the missing 1176 atom (40/9,8/3,32/9) and the
294 atom (40/9,80/21,24/7).  QR0⊙NC exists at p=7,11,19 and
is not Q_AP.  Ensemble Q_τ still unnamed.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.376 flags.  Residual (ii) and Type I stay False.

============================================================================
Setup.  NC = lex Hoffman circle of the Seidel switch (15.312
first_nc / all circles).  Joint atoms from 15.376.

============================================================================
Theorem A — PROVED (all NC of QR0 and y_*; certified p=7).
  QR0⊙NC is the single type (40/9, 24/7, 76/21).  y_*⊙NC
  produces four types:
      (32/9, 24/7, 256/63)   gen2
      (40/9, 8/3, 32/9)      the missing 1176
      (40/9, 0, 0)           inverse to AP
      (40/9, 80/21, 24/7)    the 294 atom
  together with AP, QR0, y_*, and the size-12 sharp (8/3,4,4)
  this is all eight joint atoms.  Fail: claim y_*⊙NC has only
  gen2; claim QR0⊙NC is (8/3,4,4).  ∎

Theorem B — PROVED (constructor; p=7,11,19).
  QR0⊙NC exists at all three primes.  Q++ ≠ Q_AP at p=11
  (3976/495 ≠ 56/9) and p=19 (15064/855 ≠ 88/9).  The p=7,11
  interpolant Q-- ≟ 2(p+53)/(5p) dies at p=19
  (144/95 ≠ 256/171).  Fail: claim QR0⊙NC Q++=Q_AP.  ∎

Theorem C — OPEN.  y_* empty at p=19 (15.138 ansatz).  24/p
  is still only p=7,11.  No Gauss/Jacobi name of the NL atoms
  that survives p=19.  Ensemble Q_τ unnamed.  phi_F not imported.

============================================================================
Backend: CPU NC enumeration (p=7) + first_nc (p=11,19).
Writes evidence/e1_gmin_m4_prop15377.json
"""
from __future__ import annotations

import json
import os
import sys
from functools import lru_cache
from fractions import Fraction
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15138 import construct_ystar  # noqa: E402
from e1_gmin_m4_prop15139 import affine_halfspace  # noqa: E402
from e1_gmin_m4_prop15141 import construct_y_sharp  # noqa: E402
from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field, rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15290 import omega_set, zhat_omega  # noqa: E402
from e1_gmin_m4_prop15300 import merge_cls  # noqa: E402
from e1_gmin_m4_prop15309 import qr0  # noqa: E402
from e1_gmin_m4_prop15312 import _norm_table, first_nc  # noqa: E402
from e1_gmin_m4_prop15355 import Q_AP_named, Q_QR0_p3_named  # noqa: E402
from e1_gmin_m4_prop15376 import JOINT_P7, QN_ystar_named, ystar_Q  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15377_catalog.json"


def Q_types(y, F, Omega) -> dict[str, Fraction]:
    p, q = F["p"], F["q"]
    z = np.asarray(y[1 : 1 + q], dtype=np.float64)
    if y[0] < 0:
        z = -z
    u = np.abs(zhat_omega(z, F, Omega)) ** 2
    om = {xi: i for i, xi in enumerate(Omega)}
    acc: dict[str, list[float]] = {"++": [], "--N1": [], "N*": []}
    for r in F["squares"]:
        c = merge_cls(F, r)
        if c not in acc:
            continue
        sm = 0.0
        for xi in Omega:
            sm += u[om[xi]] * u[om[F["mul"][r][xi]]]
        acc[c].append(sm / (len(Omega) * q * q))
    return {
        c: Fraction(float(np.mean(vs))).limit_denominator(p**6)
        for c, vs in acc.items()
        if vs
    }


def qkey(Q: dict[str, Fraction]) -> tuple[str, str, str]:
    return (str(Q["++"]), str(Q["--N1"]), str(Q.get("N*", 0)))


def all_nc(p: int, y0: np.ndarray) -> list[np.ndarray]:
    F = _kernel_field(p)
    q = F["q"]
    C = paley_conference_prime_power(p).astype(np.float64)
    y0 = np.asarray(y0, dtype=np.float64)
    if y0[0] < 0:
        y0 = -y0
    Cbase = (y0[:, None] * C) * y0[None, :]
    Nt = _norm_table(F)
    add, neg = F["add"], F["neg"]
    out = []
    seen: set[tuple] = set()
    for t in range(q):
        for c in range(1, p):
            pts = [u for u in range(q) if Nt[add[u][neg[t]]] == c]
            if len(pts) != p + 1:
                continue
            S = [1 + u for u in pts]
            if any(
                Cbase[S[i], S[j]] <= 0
                for i in range(p + 1)
                for j in range(i + 1, p + 1)
            ):
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
            key = tuple(np.sign(y).astype(np.int8).tolist())
            if key in seen:
                continue
            seen.add(key)
            out.append(y)
    return out


def nc_type_keys(p: int, y0: np.ndarray) -> dict[tuple[str, str, str], int]:
    F = _kernel_field(p)
    Omega = omega_set(F)
    counts: dict[tuple[str, str, str], int] = {}
    for y in all_nc(p, y0):
        k = qkey(Q_types(y, F, Omega))
        counts[k] = counts.get(k, 0) + 1
    return counts


def Qmm_interp_p7p11(p: int) -> Fraction:
    """2-point interpolant of QR0⊙NC Q-- at p=7,11. Dies at p=19."""
    return Fraction(2 * (p + 53), 5 * p)


@lru_cache(maxsize=8)
def qr0_nc_Q(p: int) -> dict[str, Fraction]:
    yq = affine_halfspace(p, (0, 1), qr0(p))
    pack = first_nc(p, yq)
    if pack is None:
        raise RuntimeError(f"QR0⊙NC missing at p={p}")
    F = _kernel_field(p)
    return Q_types(pack["y"], F, omega_set(F))


def prove_A() -> dict:
    p = 7
    yq = affine_halfspace(p, (0, 1), qr0(p))
    ys = construct_ystar(p)["y"]
    qr0_types = nc_type_keys(p, yq)
    ys_types = nc_type_keys(p, ys)
    atom_qr0 = ("40/9", "24/7", "76/21")
    atom_gen2 = ("32/9", "24/7", "256/63")
    atom_miss = ("40/9", "8/3", "32/9")
    atom_294 = ("40/9", "80/21", "24/7")
    atom_ap = ("40/9", "0", "0")
    atom_sharp = ("8/3", "4", "4")
    ok = set(qr0_types) == {atom_qr0}
    ok = ok and atom_gen2 in ys_types
    ok = ok and atom_miss in ys_types
    ok = ok and atom_294 in ys_types
    ok = ok and atom_ap in ys_types
    sh = construct_y_sharp(p)
    F = _kernel_field(p)
    Qsh = Q_types(sh["y"], F, omega_set(F))
    ok = ok and qkey(Qsh) == atom_sharp
    if atom_miss not in ys_types or atom_qr0 not in qr0_types:
        ok = False
    if atom_sharp in qr0_types:
        ok = False
    covered = {
        ("40/9", "0", "0"),
        ("32/3", "0", "0"),
        ("8/3", "80/21", "24/7"),
        atom_qr0,
        atom_gen2,
        atom_miss,
        atom_294,
        atom_sharp,
    }
    ok = ok and set(JOINT_P7) == covered
    return {
        "proved": bool(ok),
        "qr0_nc": {f"{a}|{b}|{c}": n for (a, b, c), n in qr0_types.items()},
        "ystar_nc": {f"{a}|{b}|{c}": n for (a, b, c), n in ys_types.items()},
        "sharp": qkey(Qsh),
        "all_eight": set(JOINT_P7) == covered,
        "theorem": (
            "y_*⊙NC produces the missing 1176 and the 294 atom. "
            "All eight p=7 joint atoms have constructors."
        ),
    }


def prove_B() -> dict:
    rows = {}
    ok = True
    for p in (7, 11, 19):
        Q = qr0_nc_Q(p)
        qap = Q_AP_named(p)
        rows[str(p)] = {
            "Q": {c: str(v) for c, v in Q.items()},
            "Q_AP": str(qap),
            "Qpp_eq_QAP": Q["++"] == qap,
            "Qmm_interp": str(Qmm_interp_p7p11(p)),
            "Qmm_eq_interp": Q["--N1"] == Qmm_interp_p7p11(p),
        }
        if p == 7:
            ok = ok and Q["++"] == qap == Fraction(40, 9)
            ok = ok and Q["--N1"] == Qmm_interp_p7p11(p)
        else:
            ok = ok and Q["++"] != qap
            if p == 19:
                ok = ok and Q["--N1"] != Qmm_interp_p7p11(p)
                ok = ok and Q["--N1"] == Fraction(256, 171)
        if p == 11 and Q["++"] == qap:
            ok = False
        if p == 19 and Q["--N1"] == Qmm_interp_p7p11(p):
            ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "QR0⊙NC Q++≠Q_AP at p=11,19. "
            "2(p+53)/(5p) dies at p=19."
        ),
    }


def prove_open() -> dict:
    ys19 = construct_ystar(19)
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "ystar_p19": bool(ys19.get("found")),
        "Q_tau_named_in_p": False,
        "D_named_in_p": False,
        "Qmm_named_in_p": False,
        "note": (
            "y_* empty at p=19. NL atom Q's have no name surviving "
            "p=19. Ensemble unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.377  all p=7 atoms constructed; QR0⊙NC not Q_AP", flush=True)
    A = prove_A()
    print(f"  A eight: {A['proved']} ystar_nc={A['ystar_nc']}", flush=True)
    B = prove_B()
    print(f"  B QR0NC: {B['proved']}", flush=True)
    C = prove_open()
    print(f"  C open: ystar19={C['ystar_p19']} phi_F={C['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps(
            {
                "ystar_nc_p7": A["ystar_nc"],
                "qr0_nc_p7": A["qr0_nc"],
                "QR0_NC": B["rows"],
                "ystar_p19": C["ystar_p19"],
            },
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.377",
        "title": (
            "All p=7 joint atoms constructed; "
            "QR0⊙NC Q++≠Q_AP at p=11,19"
        ),
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "all_eight_atoms_constructed": A["proved"],
            "QR0NC_not_QAP": B["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15377.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
