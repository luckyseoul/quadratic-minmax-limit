#!/usr/bin/env python3
"""
Prop 15.444 — The 15.301 p=5,7 pattern
    J(χ,ψ_2)=-(p-2)±2i√(p-1)
is not a name: it fails at every checked prime p≥11.
The half-Gauss slot C_τ(r)=∑_{s∈τ} S_χ(1+rs) is not
constant on coarse 15.298 types, so Q is not L_τ C_τ.
phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.443 flags.  Residual (ii) and Type I stay False.

============================================================================
Setup.  J(χ,ψ_k)=G(χ)G(ψ_k)/G(χψ_k) (15.301 A).  At p=5,7 one
has Re J(χ,ψ_2)=-(p-2) and |Im|=2√(p-1), and |J|=p.  C_ns is
the named constant −(q−1)/4 (15.429 A).  Coarse types as in
15.298 B: ++ / N* / --N1.

============================================================================
Theorem A — PROVED (Gauss product; certified p=5,7,11,13,17,19,23).
  J(χ,ψ_2)=-(p-2)±2i√(p-1) holds at p=5,7 (|J|=p) and fails
  at p=11 (Re≈−10.708 ≠ −9) and every later checked prime.
  Fail: claim the formula at p=11.  ∎

Theorem B — PROVED (half-Gauss on type buckets; certified 5,7,11).
  C_++(r) takes two values at p=5 already ({9,14}), so it is
  not a function of the coarse type.  Same at p=7,11.
  Fail: claim C_++ is constant at p=5.  ∎

Theorem C — OPEN.  Q_τ still unnamed.  phi_F not imported.

============================================================================
Backend: field Gauss/Jacobi (Max+-free).  Writes
evidence/e1_gmin_m4_prop15444.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _e_tr, _kernel_field  # noqa: E402
from e1_gmin_m4_prop15290 import norm_type, omega_set, paley_pair  # noqa: E402
from e1_gmin_m4_prop15298 import half_gauss  # noqa: E402
from e1_gmin_m4_prop15301 import jacobi_chi_psi  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15444_catalog.json"
CERT_J = (5, 7, 11, 13, 17, 19, 23)
CERT_C = (5, 7, 11)


def formula_re(p: int) -> int:
    return -(p - 2)


def formula_im(p: int) -> float:
    return 2.0 * (p - 1) ** 0.5


def j_matches_formula(p: int, J: complex, tol: float = 1e-6) -> bool:
    re, im = float(J.real), float(J.imag)
    return abs(re - formula_re(p)) < tol and (
        abs(im - formula_im(p)) < tol or abs(im + formula_im(p)) < tol
    )


def coarse_type(F, r: int) -> str:
    if r in (1, F["neg1"]):
        return "pm1"
    a, b = paley_pair(F, r)
    nt = norm_type(F, r)
    if a == 1 and b == 1:
        return "++"
    if nt == "N*":
        return "N*"
    if a == -1 and b == -1 and nt == "N1":
        return "--N1"
    return f"mixed:{a}:{b}:{nt}"


def C_values_on_type(F, tau: str, n_r: int = 4) -> list[str]:
    q = F["q"]
    buckets: dict[str, list[int]] = defaultdict(list)
    for r in range(1, q):
        if F["chi"][r] == 1:
            buckets[coarse_type(F, r)].append(r)
    rs_tau = buckets.get(tau, [])
    if not rs_tau:
        return []
    xi = omega_set(F)[0]
    G = float(
        sum(
            F["chi"][x] * _e_tr(F, F["mul"][xi][x]) for x in range(1, q)
        ).real
    )
    vals = set()
    for r in rs_tau[:n_r]:
        for cd in (1, -1):
            acc = 0.0
            for s in rs_tau:
                beta = F["add"][1][F["mul"][r][s]]
                acc += half_gauss(F, beta, cd, G, q)
            vals.add(str(Fraction(acc).limit_denominator(2000)))
    return sorted(vals)


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in CERT_J:
        F = _kernel_field(p)
        J = jacobi_chi_psi(F, 2)
        hit = j_matches_formula(p, J) and abs(abs(J) - p) < 1e-6
        rows[str(p)] = {
            "re": float(J.real),
            "im": float(J.imag),
            "want_re": formula_re(p),
            "hit": hit,
        }
        if p in (5, 7):
            if not hit:
                ok = False
        else:
            if hit:
                ok = False
    if j_matches_formula(11, jacobi_chi_psi(_kernel_field(11), 2)):
        ok = False
    if abs(rows["11"]["re"] + 9.0) < 1e-3:
        ok = False
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "J(χ,ψ_2)=-(p-2)±2i√(p-1) at p=5,7 and not at p≥11. "
            "Fail: claim the formula at p=11."
        ),
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in CERT_C:
        F = _kernel_field(p)
        vals = C_values_on_type(F, "++")
        const = len(vals) == 1
        rows[str(p)] = {"C_pp": vals, "const": const}
        if const:
            ok = False
    if len(C_values_on_type(_kernel_field(5), "++")) < 2:
        ok = False
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "C_++ takes ≥2 values at p=5,7,11. "
            "Fail: claim C_++ constant at p=5."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "J(χ,ψ_2) has no 0-parameter name from the p=5,7 "
            "pattern. C_τ is not type-constant. Q_τ unnamed. "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.444  J(χ,ψ_2) 2-point formula dies at p≥11", flush=True)
    A = prove_A()
    print(
        f"  A J-formula: {A['proved']} p11 Re={A['by_p']['11']['re']:.4f} "
        f"want={A['by_p']['11']['want_re']}",
        flush=True,
    )
    B = prove_B()
    print(f"  B C_++: {B['proved']} p5={B['by_p']['5']['C_pp']}", flush=True)
    C = prove_open()
    print(f"  C open phi_F={C['phi_F_ge_6']} e1={C['e1']}", flush=True)
    CAT.write_text(
        json.dumps(
            {
                "A": A["proved"],
                "B": B["proved"],
                "p11_re": A["by_p"]["11"]["re"],
                "p5_Cpp": B["by_p"]["5"]["C_pp"],
            },
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.444",
        "title": "J(χ,ψ_2)=-(p-2)±2i√(p-1) dies at p≥11; C_τ not type-constant",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "J_formula_dies": A["proved"],
            "C_tau_not_const": B["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15444.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
