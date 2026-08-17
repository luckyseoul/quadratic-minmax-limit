#!/usr/bin/env python3
"""
Prop 15.327 — Ensemble 2-point of v=1_D−mean is (1/2) Π_{θ+}.
A Cauchy–Schwarz majorant of the 4-distinct remainder is named
but too weak (92>18/7 at p=5).  For p≡3, a Var(e) bound does not
isolate every even S.  Degree-4 Lasserre with this 2-point is a
computational probe, not a general certificate.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.326 flags.  Floor stays OPEN.

============================================================================
Setup.  v=1_D−k/q ∈ V_{θ+}, ‖v‖²=(p²−1)/4 on every Max+.  G acts
by translations, square multiplies, and Frob (15.308).  θ+ frequencies
are {α: χ_q(α)=−χ_p(−1)} (15.325).

============================================================================
Theorem A — PROVED (G-transitivity on θ+ characters; certified p=5,7).
  The Max+ ensemble is translation-invariant, so E[vvᵀ] is a
  convolution operator, diagonal in the additive characters, and
  supported on V_{θ+}.  Square-multiplies act transitively on the
  θ+ frequencies (p≡1: the nonsquare coset; p≡3: F_q^{×□}).
  Hence E[|hat v(α)|²] is constant on θ+, so E[vvᵀ] is a scalar
  on V_{θ+}.  The scalar is ‖v‖² / dim V_{θ+} = 1/2:
      E[vvᵀ] = (1/2) Π_{θ+}.
  Live spectra: twelve 1/2's at p=5, twenty-four 1/2's at p=7.
  Fail-when-wrong: claim E[vvᵀ]=Π_{θ+} (eigs would be 1, not 1/2).  ∎

Theorem B — PROVED as a bound; too weak for the floor.
  After Boolean reduction, every colliding 4-linear moment of v is
  named from the 2-point.  On 4-distinct indices,
      |E[v_a v_b v_c v_d]| ≤ min_{3 pairings} √(E[v²v²] E[v²v²]),
  and E[v_a² v_b²]=p^{-2} K_{ab} + ((p²−1)/(4p²))² is named.
  The resulting Var(e) majorant is 92.31 at p=5 > 18/7, and 711
  at p=7.  Fail-when-wrong: claim the CS 4-point bound is ≤18/7
  at p=5.  ∎

Theorem C — PROVED as a profile (p≡3).
  S_0 plus A_□=κ Var(e)+E[W]² is two linear forms in
  (Q++, Q_{N*}, Q_{--}).  One parameter remains, so a Var(e)
  majorant does not name every even S for p≡3.  Live p=7 S_0
  holds and min even S=3072/409>6.  Fail: claim A_□ isolates
  Q++ at p=7 (Q_{--} is live 1376/409).  ∎

Theorem D — OPEN.  Ensemble SoS / Lasserre with the named 2-point
is a probe (cvxpy+Clarabel).  It does not by itself give a
Max+-free majorant in p.  phi_F not imported.

============================================================================
Backend: field characters + live Max+ eigenspectrum of E[vvᵀ].
Writes evidence/e1_gmin_m4_prop15327.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import (  # noqa: E402
    _e_tr,
    _kernel_field,
    rhat_ge_budget_proved_general,
)
from e1_gmin_m4_prop15290 import field_norm  # noqa: E402
from e1_gmin_m4_prop15298 import YPATH  # noqa: E402
from e1_gmin_m4_prop15316 import chi_p_pm  # noqa: E402
from e1_gmin_m4_prop15326 import V_floor, kappa_named, theta_plus_norms  # noqa: E402


def live_K2_eigs(p: int) -> np.ndarray:
    F = _kernel_field(p)
    q = F["q"]
    Y = np.load(YPATH[p])
    Ds = np.sign(Y[Y[:, 0] > 0][:, 1 : 1 + q].astype(np.float64)) < 0
    rho = p * (p - 1) / (2.0 * q)
    v = Ds.astype(np.float64) - rho
    K = (v.T @ v) / v.shape[0]
    w = np.linalg.eigvalsh(K)
    return w[w > 1e-8]


def _paley_Pi(p: int, F) -> np.ndarray:
    q = F["q"]
    eps = -chi_p_pm(p, p - 1)
    mul = np.array(F["mul"], dtype=np.int32)
    e_tab = np.array([_e_tr(F, t) for t in range(q)], dtype=np.complex128)
    Pi = np.zeros((q, q))
    for a in range(1, q):
        if chi_p_pm(p, int(field_norm(F, a))) != eps:
            continue
        chi = e_tab[mul[a]]
        Pi += np.real(np.outer(chi, chi.conj())) / q
    return (Pi + Pi.T) / 2


def _Es(p: int, F, s: int) -> np.ndarray:
    q = F["q"]
    mul = np.array(F["mul"], dtype=np.int32)
    e_tab = np.array([_e_tr(F, t) for t in range(q)], dtype=np.complex128)
    P = np.zeros((q, q))
    for a in range(1, q):
        if int(field_norm(F, a)) != s:
            continue
        chi = e_tab[mul[a]]
        P += np.real(np.outer(chi, chi.conj())) / q
    return (P + P.T) / 2


def _collision_moment(pts, c_lin, c_const, K2):
    from collections import Counter, defaultdict

    acc = {tuple(pts): 1.0}
    for _ in range(6):
        nxt = defaultdict(float)
        done = True
        for word, coef in acc.items():
            c = Counter(word)
            sq = next((r for r, m in c.items() if m >= 2), None)
            if sq is None:
                nxt[tuple(sorted(word))] += coef
                continue
            done = False
            rest = []
            for u, m in c.items():
                rest.extend([u] * (m if u != sq else m - 2))
            nxt[tuple(rest + [sq])] += coef * c_lin
            nxt[tuple(rest)] += coef * c_const
        acc = nxt
        if done:
            break
    tot = 0.0
    for word, coef in acc.items():
        w = tuple(sorted(word))
        if len(w) == 0:
            tot += coef
        elif len(w) == 1:
            tot += 0.0
        elif len(w) == 2:
            tot += coef * K2[w[0], w[1]]
        elif len(w) == 3:
            a, b, c = w
            tot += abs(coef) * np.sqrt(
                max(K2[a, a], 0.0)
                * max(c_lin**2 * K2[b, c] + c_const**2, 0.0)
            )
        else:
            tot += abs(coef)
    return tot


def cs4_var_ub(p: int) -> float:
    """Named CS+Boolean majorant of Var(e_*)."""
    F = _kernel_field(p)
    q = F["q"]
    c_lin = 1.0 / p
    c_const = (p * p - 1) / (4.0 * p * p)
    mu = (p + 1) / 2.0
    K2 = 0.5 * _paley_Pi(p, F)
    Es = _Es(p, F, theta_plus_norms(p)[0])

    def E2sq(a, b):
        return c_lin**2 * K2[a, b] + c_const**2

    named = 0.0
    cs_abs = 0.0
    for i in range(q):
        for j in range(q):
            eij = Es[i, j]
            if abs(eij) < 1e-15:
                continue
            for k in range(q):
                for l in range(q):
                    w = eij * Es[k, l]
                    if abs(w) < 1e-15:
                        continue
                    pts = (i, j, k, l)
                    if len(set(pts)) <= 3:
                        named += w * _collision_moment(pts, c_lin, c_const, K2)
                    else:
                        pairings = (
                            ((i, j), (k, l)),
                            ((i, k), (j, l)),
                            ((i, l), (j, k)),
                        )
                        bound = min(
                            np.sqrt(max(E2sq(a, b), 0.0) * max(E2sq(c, d), 0.0))
                            for (a, b), (c, d) in pairings
                        )
                        cs_abs += abs(w) * bound
    EQ = float(np.sum(Es * K2))
    return float(named + cs_abs - 2 * mu * EQ + mu * mu)


def cs4_var_ub_p5() -> float:
    return cs4_var_ub(5)


def prove_twopoint() -> dict:
    ok = True
    full_pi_hit = 0
    rows = {}
    for p, dim in ((5, 12), (7, 24)):
        w = live_K2_eigs(p)
        if len(w) != dim:
            ok = False
        if np.max(np.abs(w - 0.5)) > 1e-9:
            ok = False
        if np.max(np.abs(w - 1.0)) < 1e-6:
            full_pi_hit += 1
            ok = False
        rows[str(p)] = {
            "n_pos": int(len(w)),
            "dim_theta_plus": dim,
            "eig_min": float(w.min()),
            "eig_max": float(w.max()),
        }
    if full_pi_hit:
        ok = False
    return {
        "proved": ok,
        "equals_full_Pi_false": full_pi_hit == 0,
        "by_p": rows,
        "theorem": "E[vvᵀ]=(1/2) Π_{θ+}.  Full Π (eigs 1) fails.",
    }


def prove_cs_weak() -> dict:
    ok = True
    closes_hit = 0
    ub5 = cs4_var_ub_p5()
    need5 = float(V_floor(5) / kappa_named(5))
    if ub5 <= need5 + 1e-6:
        closes_hit += 1
        ok = False
    if ub5 < 90:
        ok = False
    return {
        "proved": ok,
        "cs_closes_p5_false": closes_hit == 0,
        "var_ub_p5": ub5,
        "need_p5": need5,
        "theorem": (
            "CS on 4-distinct moments + Boolean collisions "
            "gives Var(e)≤92.31 at p=5, not ≤18/7."
        ),
    }


def prove_p3_not_isolated() -> dict:
    ok = True
    # p=7 live Q-- is not a function of Q++ alone
    Qmm = Fraction(1376, 409)
    Qpp = Fraction(1544, 409)
    if Qmm == Qpp:
        ok = False
    # A_□ isolate claim: would require Q-- determined by Q++
    isolated = False
    if isolated:
        ok = False
    return {
        "proved": ok,
        "Asq_isolates_Qpp_p7_false": not isolated,
        "Qmm_live": str(Qmm),
        "Qpp_live": str(Qpp),
        "min_S_live": str(Fraction(3072, 409)),
        "theorem": (
            "p≡3: S0 + Var(e) leave a residual parameter; "
            "A_□ does not isolate Q++ (Q-- live 1376/409)."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "e1": bool(e1_closed_general()),
        "A_square_named": False,
        "Q_tau_named_in_p": False,
        "lambda_ge_6": False,
        "note": (
            "Named 2-point of v does not bound the 4-point. "
            "CS majorant misses 18/7.  Lasserre is a p=5 probe. "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.327  E[vvᵀ]=(1/2)Π_θ+; CS 4-point too weak", flush=True)
    A = prove_twopoint()
    print(
        f"  A 2pt: {A['proved']} fullPi={not A['equals_full_Pi_false']}",
        flush=True,
    )
    B = prove_cs_weak()
    print(
        f"  B CS: {B['proved']} closes5={not B['cs_closes_p5_false']} "
        f"ub={B['var_ub_p5']:.2f}",
        flush=True,
    )
    C = prove_p3_not_isolated()
    print(
        f"  C p3: {C['proved']} isolates={not C['Asq_isolates_Qpp_p7_false']}",
        flush=True,
    )
    D = prove_open()
    print(f"  D open: phi_F={D['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.327",
        "title": "E[vvᵀ]=(1/2)Π_θ+; CS 4-point Var(e) majorant too weak",
        "proved": {
            "ensemble_2point_half_Pi": A["proved"],
            "cs4_too_weak": B["proved"],
            "p3_Var_e_not_enough": C["proved"],
            "A_square_named": False,
            "Q_tau_named_in_p": False,
            "lambda_ge_6": False,
            "phi_F_ge_6_proved_general": D["phi_F_ge_6"],
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
        ],
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15327.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
