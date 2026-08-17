#!/usr/bin/env python3
"""
Prop 15.469 — Cov(N(d),N(sd)) on coarse ++ is 10/39 and
2793/1636.  Named L/Johnson/μ atoms and their ratios miss
both primes.  Not a 2-point interpolant.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.468 flags.

============================================================================
Setup.  15.467: coarse L++(+)=985/39, 91581/818.  15.468: the
2-point of 1_D is named and pointwise-constant; L++=μ+²+Cov is
an ensemble 4-point.  Aut_∞ conference+row-sum still has
nullity 1 (p=5) and 5 (p=7) (15.330).

============================================================================
Theorem A — PROVED (15.467 coarse L − μ+²; p=5,7).
  Cov++ := L++(+) − μ+² equals 10/39 at p=5 and 2793/1636 at p=7.
  Fail: claim Cov++=0 (i.e. L++=μ+²).  ∎

Theorem B — PROVED (same Cov vs named catalog).
  Lpar_1d, Lpar(R), Lspl_1d, L1d_pp, μ+², μ+μ−, S_χ, χ+, k,
  and all pairwise ratios miss Cov++ at both primes.
  Fail: Lpar_1d − μ+² equals 10/39 at p=5.  ∎

Theorem C — OPEN.  Cov++ / L++ still unnamed in p.  Aut_∞
  Jacobi pairings that cut the 4-point kernel still carry den
  |H+|/(2p).  phi_F not imported.

============================================================================
Backend: 15.467 coarse_L + named atoms.  Writes
evidence/e1_gmin_m4_prop15469.json
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

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15298 import mu_minus, mu_plus  # noqa: E402
from e1_gmin_m4_prop15414 import Lpar_1d_named  # noqa: E402
from e1_gmin_m4_prop15415 import L1d_pp_named, Lspl_1d_named  # noqa: E402
from e1_gmin_m4_prop15430 import Lpar_R_named  # noqa: E402
from e1_gmin_m4_prop15467 import coarse_L  # noqa: E402
from e1_gmin_m4_prop15468 import S_chi_named, chi_plus_named, k_named  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15469_catalog.json"

LIVE_COV = {5: Fraction(10, 39), 7: Fraction(2793, 1636)}
ATOMS = ("Lpar_1d", "Lpar_R", "Lspl_1d", "L1d_pp", "mup_sq", "mup_mum", "S_chi", "chi+", "k")


def cov_pp(p: int) -> Fraction:
    return coarse_L(p)["++"][1] - mu_plus(p) ** 2


def atom(name: str, p: int) -> Fraction:
    if name == "Lpar_1d":
        return Lpar_1d_named(p)
    if name == "Lpar_R":
        return Lpar_R_named(p)
    if name == "Lspl_1d":
        return Lspl_1d_named(p)
    if name == "L1d_pp":
        return L1d_pp_named(p)
    if name == "mup_sq":
        return mu_plus(p) ** 2
    if name == "mup_mum":
        return mu_plus(p) * mu_minus(p)
    if name == "S_chi":
        return Fraction(S_chi_named(p))
    if name == "chi+":
        return chi_plus_named(p)
    if name == "k":
        return Fraction(k_named(p))
    raise KeyError(name)


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        c = cov_pp(p)
        rows[str(p)] = str(c)
        if c != LIVE_COV[p]:
            ok = False
        if c == 0:
            ok = False
    if cov_pp(5) == 0:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": "Cov++=10/39, 2793/1636. Fail: Cov++=0.",
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    t = {p: cov_pp(p) for p in (5, 7)}
    for p in (5, 7):
        rec = {}
        for name in ATOMS:
            v = atom(name, p)
            rec[name] = str(v)
            if v == t[p]:
                ok = False
        # pairwise ratios
        for i, ni in enumerate(ATOMS):
            ai = atom(ni, p)
            if ai == 0:
                continue
            for nj in ATOMS[i + 1 :]:
                aj = atom(nj, p)
                if aj == 0:
                    continue
                if ai / aj == t[p] or aj / ai == t[p]:
                    ok = False
        rows[str(p)] = rec
    if atom("Lpar_1d", 5) - mu_plus(5) ** 2 == Fraction(10, 39):
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "Named atoms/ratios miss Cov++ at both primes. "
            "Fail: Lpar_1d−μ+²=10/39."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "Cov++ isolated as 10/39, 2793/1636. "
            "Jac/edge6 pairings that cut Aut_∞ ker still have den D. "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.469  Cov++ isolated; named atoms miss", flush=True)
    A = prove_A()
    print(f"  A cov: {A['proved']} {A['rows']}", flush=True)
    B = prove_B()
    print(f"  B atoms: {B['proved']}", flush=True)
    C = prove_open()
    print(f"  C open phi_F={C['phi_F_ge_6']}", flush=True)
    CAT.write_text(json.dumps({"A": A["proved"], "B": B["proved"]}, indent=2) + "\n")
    out = {
        "prop": "15.469",
        "title": "Cov++ is 10/39 and 2793/1636; named atoms miss both",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "cov_values": A["proved"],
            "named_atoms_miss_cov": B["proved"],
            "phi_F_ge_6_proved_general": C["phi_F_ge_6"],
            "residual_ii_k_eq_4p_empty": C["residual_ii_k_eq_4p_empty"],
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
            "residual_ii_k_eq_4p_empty",
        ],
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15469.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
