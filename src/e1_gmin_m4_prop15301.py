#!/usr/bin/env python3
"""
Prop 15.301 — S_k formula 16−16(p−Re J(χ,ψ_k))/d (d=dim V₊).

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.300 flags.  Floor stays OPEN unless inversion recovers
live Q at p=5 and p=7 (it does not).

============================================================================
Setup.  ψ_k = F_q^× character of index k (even).  J(χ,ψ)=G(χ)G(ψ)/G(χψ)
for nontrivial ψ, |J|=p.  d=dim V_+=(q+1)/2.  S_k/q²=∑_{r∈T} ψ_k(r) Q(r)/q²
= λ(ψ_k) on F for k≠0.  S_0=2(q−1) named (15.300).

============================================================================
Theorem A — PROVED (Max+-free Gauss sums).
  J(χ,ψ_k) is the Jacobi sum G(χ)G(ψ_k)/G(χψ_k).  |J|=p on the
  15.300 basis for k≠0.  Fail-when-wrong: replace J by G(χ)G(ψ)/G(χψ̄)
  (inverted character); |J−J_inv|>1 at p=5 k=2.  ∎

Theorem B — PROVED (p=5 identity; p=7 counterexample).
  The unique formula linear in Re J that matches both p=5 basis
  characters k=2,6 (S=80/13 and 176/13=16(d−2)/d) is
      S_k/q² = 16 − 16(p − Re J(χ,ψ_k)) / d.
  It holds at p=5 (err < 1e-12) and fails at p=7 (k=12: 16 vs
  live 4320/409; max err > 1).  Fail-when-wrong: claim the
  formula at p=7; then J^{-1} S_formula misses live Q_++.  ∎

Theorem C — OPEN.  S_{k≠0} has no named-den formula that works
  at both p=5 and p=7.  Q=J^{-1}S therefore still has den N/(4p).
  ⟨δ,ψ⟩≤2 / phi_F stay OPEN.  ∎

============================================================================
Backend: CPU Gauss sums (no Max+ for A–B except the live-Q invert).
Writes evidence/e1_gmin_m4_prop15301.json
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import (  # noqa: E402
    _gauss,
    _kernel_field,
    _psi_vec,
    rhat_ge_budget_proved_general,
)
from e1_gmin_m4_prop15298 import live_Q_any  # noqa: E402
from e1_gmin_m4_prop15300 import classes_of, merge_cls  # noqa: E402

BASES = {5: (0, 2, 6), 7: (0, 2, 8, 12)}


def dim_Vplus(p: int) -> int:
    return (p * p + 1) // 2


def jacobi_chi_psi(F, k: int) -> complex:
    """J(χ,ψ_k)=G(χ)G(ψ)/G(χψ).  0 if any factor is trivial."""
    psi = _psi_vec(F, k)
    Gpsi = _gauss(F, psi)
    Gchi = _gauss(F, [complex(F["chi"][x]) for x in range(F["q"])])
    Gchipsi = _gauss(F, [F["chi"][x] * psi[x] for x in range(F["q"])])
    if abs(Gpsi) < 1e-8 or abs(Gchipsi) < 1e-8:
        return 0j
    return Gchi * Gpsi / Gchipsi


def jacobi_inverted(F, k: int) -> complex:
    """Wrong: G(χ)G(ψ)/G(χψ̄)."""
    psi = _psi_vec(F, k)
    psibar = [z.conjugate() for z in psi]
    Gpsi = _gauss(F, psi)
    Gchi = _gauss(F, [complex(F["chi"][x]) for x in range(F["q"])])
    Gchipsi_bar = _gauss(F, [F["chi"][x] * psibar[x] for x in range(F["q"])])
    if abs(Gpsi) < 1e-8 or abs(Gchipsi_bar) < 1e-8:
        return 0j
    return Gchi * Gpsi / Gchipsi_bar


def S_formula(p: int, k: int) -> Fraction:
    """16 − 16(p − Re J)/d.  For k=0 return the named mass 2(q−1)."""
    if k == 0:
        return Fraction(2 * (p * p - 1))
    F = _kernel_field(p)
    J = jacobi_chi_psi(F, k)
    d = dim_Vplus(p)
    return Fraction(16) - Fraction(16 * (p - int(round(J.real))), d)


def live_S_vector(p: int) -> tuple[list[str], np.ndarray, np.ndarray, np.ndarray]:
    F = _kernel_field(p)
    Q = live_Q_any(p)
    q2 = float(F["q"] ** 2)
    C = classes_of(F)
    names = list(C)
    Qv = np.array([float(np.mean([Q[r] / q2 for r in C[c]])) for c in names])
    ks = BASES[p]
    Jmat = []
    S_live = []
    for k in ks:
        psi = _psi_vec(F, k)
        row = []
        s = 0j
        for c in names:
            acc = 0j
            for r in C[c]:
                acc += psi[r]
                s += psi[r] * Q[r] / q2
            row.append(acc)
        Jmat.append(row)
        S_live.append(s.real)
    return names, np.real(np.array(Jmat)), np.array(S_live), Qv


def prove_jacobi_named() -> dict:
    ok = True
    invert_hit = 0
    invert_tot = 0
    rows = {}
    for p in (5, 7):
        F = _kernel_field(p)
        rec = {}
        for k in BASES[p]:
            if k == 0:
                continue
            J = jacobi_chi_psi(F, k)
            Jinv = jacobi_inverted(F, k)
            invert_tot += 1
            if abs(J) < p - 0.5 or abs(abs(J) - p) > 1e-6:
                ok = False
            if abs(J - Jinv) > 1.0:
                invert_hit += 1
            rec[str(k)] = {
                "J_re": float(J.real),
                "J_im": float(J.imag),
                "absJ": float(abs(J)),
                "absJ_eq_p": abs(abs(J) - p) < 1e-6,
                "inv_diff": float(abs(J - Jinv)),
            }
        rows[str(p)] = rec
    if invert_hit == 0:
        ok = False
    return {
        "proved": ok,
        "invert_hit": invert_hit,
        "invert_tot": invert_tot,
        "by_p": rows,
        "theorem": (
            "J(χ,ψ_k)=G(χ)G(ψ)/G(χψ), |J|=p on the 15.300 basis k≠0. "
            "Inverted G(χψ̄) differs."
        ),
    }


def prove_formula_p5_only() -> dict:
    ok = True
    p5_errs = {}
    p7_errs = {}
    for k in BASES[5]:
        if k == 0:
            continue
        names, Jmat, S_live, Qv = live_S_vector(5)
        # compute formula S
        Sf = float(S_formula(5, k))
        # find this k's live S
        idx = list(BASES[5]).index(k)
        err = abs(Sf - S_live[idx])
        p5_errs[str(k)] = err
        if err > 1e-8:
            ok = False
    # p=7 must fail
    p7_max = 0.0
    for k in BASES[7]:
        if k == 0:
            continue
        _, _, S_live, _ = live_S_vector(7)
        idx = list(BASES[7]).index(k)
        err = abs(float(S_formula(7, k)) - S_live[idx])
        p7_errs[str(k)] = err
        p7_max = max(p7_max, err)
    if p7_max < 0.5:
        ok = False
    # inversion with formula S recovers Q at p=5, misses at p=7
    recov = {}
    for p in (5, 7):
        names, Jmat, S_live, Qv = live_S_vector(p)
        Sf = np.array([float(S_formula(p, k)) for k in BASES[p]])
        Qf = np.linalg.solve(Jmat, Sf)
        jpp = names.index("++")
        recov[str(p)] = {
            "Qpp_live": float(Qv[jpp]),
            "Qpp_formula": float(Qf[jpp]),
            "err": float(abs(Qf[jpp] - Qv[jpp])),
        }
    if recov["5"]["err"] > 1e-6:
        ok = False
    if recov["7"]["err"] < 0.05:
        ok = False
    return {
        "proved": ok,
        "p5_S_err": p5_errs,
        "p7_S_err": p7_errs,
        "p7_max_S_err": p7_max,
        "inversion": recov,
        "formula": "16 - 16(p - Re J(χ,ψ_k))/dim V_+",
        "theorem": (
            "S=16−16(p−Re J)/d matches both p=5 basis characters and "
            "recovers Q via J^{-1}S.  At p=7 the same formula misses "
            "live S (max err>1) and live Q_++."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "e1": bool(e1_closed_general()),
        "Q_tau_named_in_p": False,
        "S_k_named_in_p": False,
        "note": (
            "The named-den formula 16−16(p−Re J)/d is p=5-only. "
            "S_{k≠0} still unnamed for general p. Floor OPEN."
        ),
    }


def main() -> dict:
    print("Prop 15.301  S=16−16(p−Re J)/d is p=5-only", flush=True)
    A = prove_jacobi_named()
    print(f"  A Jacobi |J|=p: {A['proved']} invert={A['invert_hit']}/{A['invert_tot']}", flush=True)
    B = prove_formula_p5_only()
    print(
        f"  B p=5-only formula: {B['proved']} p7_max_err={B['p7_max_S_err']:.3f} "
        f"Qpp_err7={B['inversion']['7']['err']:.3f}",
        flush=True,
    )
    C = prove_open()
    print(f"  C floor: {C['proved']} phi_F={C['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.301",
        "title": "S_k=16−16(p−Re J)/d holds at p=5 and fails at p=7",
        "proved": {
            "jacobi_abs_p": A["proved"],
            "formula_p5_only": B["proved"],
            "phi_F_ge_6_proved_general": C["phi_F_ge_6"],
            "Q_tau_named_in_p": False,
            "S_k_named_in_p": False,
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
    path = ROOT / "evidence" / "e1_gmin_m4_prop15301.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
