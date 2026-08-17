#!/usr/bin/env python3
"""
Prop 15.305 — hatN variance; F_p-orbit equicorrelation for Q_++;
nonsquare-direction lines meet every Max+ in (p−1)/2 points.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.304 flags.  Does **not** prove λ_j≥6.  Floor stays OPEN.

============================================================================
Setup.  D⊂F_q is {z=−1} for a Max+ with y_∞=+1, |D|=k=p(p−1)/2.
hatN(ξ)=∑_δ N(δ)e(ξδ)=|1̂_D(ξ)|²=w(ξ)/4 on Ω.  S=∑_{a∈F_p^×} hatN(aξ).
A direction ⟨v⟩ is nonsquare if χ(v)=−1.  Those lines are Paley cocliques
of Hoffman size p.

============================================================================
Theorem A — PROVED (Max+-free; Q(±1)=8q² and E[hatN]=q/2).
  Var(hatN(ξ))=q²/4 on Ω.  Fail-when-wrong: claim Var=q²/2.  ∎

Theorem B — PROVED (equicorrelation on F_p^×·ξ ⊂ Ω).
  For a≠±a' in F_p^×, E[hatN(aξ) hatN(a'ξ)]=Q_{++} q²/16; antipodes
  give Q(−1)/16=q²/2.  Hence
      E[S²]=(p−1)q² + (p−1)(p−3) q² Q_{++}/16.
  Fail-when-wrong: drop antipodes (use (p-1)(p-2) instead of (p-3)).  ∎

Theorem C — PROVED (regular-set identity + F_p quadratic sum).
  For every nonsquare direction and every line ℓ in that parallel class,
      |D ∩ ℓ| = (p−1)/2.
  Proof.  Cy=py  ⇒  ∑_{δ≠0} χ(δ) 1_D(x+δ) = p 1_D(x) − (p−1)/2.
  Sum over x∈ℓ.  LHS=∑_δ χ(δ)|D∩(ℓ−δ)|.  On ⟨v⟩^×, χ=−1 and ℓ−δ=ℓ,
  giving −(p−1)n_ℓ.  On any other parallel line the coset sum
  ∑_{t∈F_p} χ(w+tv) equals −χ(v)=+1 (disc≠0 iff w∉⟨v⟩; standard
  ∑ χ(At²+Bt+C)=−χ(A)).  So LHS=k−p n_ℓ.  RHS=p n_ℓ−k.  Thus
  n_ℓ=k/p=(p−1)/2.  Fail-when-wrong: claim a square-direction line
  is constant (spread p at p=5,7).  ∎

Theorem D — OPEN.  Square directions give the F_p-orbit of Ω and
  carry Q_{++} via (B).  E[(∑ n_ℓ²)²] on a square class still has
  den N/(4p).  The 15.304 comparison E|⟨ΔN,G_U⟩|²≥(3/4)(p+1)q²
  is not implied.  phi_F not imported.  ∎

============================================================================
Backend: CPU field sums for C; N-cache only to check A/B and the
square-direction spread.  Writes evidence/e1_gmin_m4_prop15305.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field, rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15290 import in_subfield, live_Q_on_T, omega_set  # noqa: E402
from e1_gmin_m4_prop15298 import load_N  # noqa: E402
from e1_gmin_m4_prop15300 import merge_cls  # noqa: E402
from e1_gmin_m4_prop15304 import e_vec  # noqa: E402

YPATH = {5: "/tmp/maxplus_p5.npy", 7: "/tmp/maxplus_p7.npy"}


def chi_line_sum(F, w: int, v: int) -> int:
    """∑_{t=0}^{p-1} χ(w+t v), with χ(0)=0."""
    acc = 0
    for t in range(F["p"]):
        x = F["add"][w][F["mul"][t][v]]
        acc += 0 if x == 0 else int(F["chi"][x])
    return acc


def in_span(F, w: int, v: int) -> bool:
    return any(F["mul"][t][v] == w for t in range(F["p"]))


def prove_variance() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        F = _kernel_field(p)
        q = F["q"]
        N = load_N(p, F).astype(np.float64)
        xi = omega_set(F)[0]
        hat = N @ e_vec(F, xi)
        mean = complex(np.mean(hat))
        var = float(np.mean(np.abs(hat) ** 2) - abs(mean) ** 2)
        expect_mean = q / 2
        expect_var = (q * q) / 4
        if abs(mean - expect_mean) > 1e-6:
            ok = False
        if abs(var - expect_var) > 1e-4:
            ok = False
        if abs(var - 2 * expect_var) < 1e-4:
            ok = False  # fail-when-wrong: Var=q²/2 would match this
        rows[str(p)] = {
            "mean": float(mean.real),
            "var": var,
            "expect_mean": expect_mean,
            "expect_var": expect_var,
        }
    return {
        "proved": ok,
        "by_p": rows,
        "theorem": "Var(hatN)=q²/4 on Ω.  E[hatN]=q/2.  Not q²/2.",
    }


def prove_equicorrelation() -> dict:
    ok = True
    antipode_hit = 0
    rows = {}
    for p in (5, 7):
        F = _kernel_field(p)
        q = F["q"]
        q2 = float(q * q)
        N = load_N(p, F).astype(np.float64)
        xi = omega_set(F)[0]
        Fp = [a for a in range(1, q) if in_subfield(F, a)]
        H = np.stack([N @ e_vec(F, F["mul"][a][xi]) for a in Fp], axis=1)
        S = H.sum(axis=1)
        ES2 = float(np.mean(np.abs(S) ** 2))
        Q = live_Q_on_T(p)
        Qpp = float(np.mean([Q[r] / q2 for r in Q if merge_cls(F, r) == "++"]))
        pred = (p - 1) * q2 + (p - 1) * (p - 3) * q2 * Qpp / 16.0
        wrong = (p - 1) * q2 / 2 + (p - 1) * (p - 2) * q2 * Qpp / 16.0
        err = abs(ES2 - pred)
        if err > 1e-4 * max(1.0, ES2):
            ok = False
        if abs(ES2 - wrong) > 1.0:
            antipode_hit += 1
        rows[str(p)] = {
            "E_S2": ES2,
            "pred": pred,
            "wrong_no_antipode": wrong,
            "Qpp": Qpp,
            "err": err,
        }
    if antipode_hit == 0:
        ok = False
    return {
        "proved": ok,
        "antipode_hit": antipode_hit,
        "by_p": rows,
        "theorem": (
            "E[S²]=(p−1)q²+(p−1)(p−3)q² Q_{++}/16. "
            "Omitting antipodes (p−2 in place of p−3) fails."
        ),
    }


def prove_ns_lines() -> dict:
    ok = True
    chi_ok = True
    sq_spread_hit = 0
    rows = {}
    for p in (5, 7):
        F = _kernel_field(p)
        # Gauss sums along lines
        v = next(x for x in range(1, F["q"]) if F["chi"][x] == -1)
        off, onspan = [], []
        for w in range(F["q"]):
            sm = chi_line_sum(F, w, v)
            if in_span(F, w, v):
                onspan.append(sm)
            else:
                off.append(sm)
        if any(s != 1 for s in off):
            chi_ok = False
            ok = False
        if any(s != 1 - p for s in onspan):
            # on ⟨v⟩: ∑_{t≠0} χ(tv)=∑_{t≠0} χ(v)=-(p-1), χ(0)=0
            if any(s != -(p - 1) for s in onspan):
                chi_ok = False
                ok = False
        # intersections
        Y = np.load(YPATH[p])
        q = F["q"]
        D = np.sign(Y[Y[:, 0] > 0][:, 1 : 1 + q].astype(np.float64)) < 0
        ns_spreads = []
        sq_spreads = []
        target = (p - 1) / 2
        # vertical + slopes
        specs = [(p, "vert")] + [(1 + p * s, f"s{s}") for s in range(p)]
        for vec, _name in specs:
            chi = int(F["chi"][vec])
            if _name == "vert":
                lines = [[c + p * x1 for x1 in range(p)] for c in range(p)]
            else:
                s = (vec // p) % p
                lines = [
                    [x0 + p * ((s * x0 + c) % p) for x0 in range(p)] for c in range(p)
                ]
            ns = np.stack([D[:, ix].sum(axis=1) for ix in lines], axis=1)
            spr = float(ns.max() - ns.min())
            if chi == -1:
                ns_spreads.append(spr)
                if abs(ns.mean() - target) > 1e-9 or spr > 0:
                    ok = False
            else:
                sq_spreads.append(spr)
                if spr >= p - 1e-9:
                    sq_spread_hit += 1
        rows[str(p)] = {
            "chi_offspan": off[:3] if off else [],
            "ns_spreads": ns_spreads,
            "sq_spreads": sq_spreads,
            "target": target,
        }
    if not chi_ok or sq_spread_hit == 0:
        ok = False
    return {
        "proved": ok,
        "chi_line_sum_plus_one": chi_ok,
        "square_dir_varies": sq_spread_hit,
        "by_p": rows,
        "theorem": (
            "Nonsquare-direction lines meet D in (p−1)/2 points. "
            "∑_t χ(w+tv)=+1 off ⟨v⟩.  Square directions vary (spread p)."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "e1": bool(e1_closed_general()),
        "Q_tau_named_in_p": False,
        "lambda_ge_6": False,
        "note": (
            "ns-line theorem and E[S²]↔Q_{++} do not name E[(∑ n²)²] "
            "on a square parallel class, nor prove the 15.304 square-mass "
            "bound. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.305  Var(hatN)=q²/4; ns-lines meet D in (p−1)/2", flush=True)
    A = prove_variance()
    print(f"  A Var(hatN)=q²/4: {A['proved']}", flush=True)
    B = prove_equicorrelation()
    print(f"  B E[S²] vs Q++: {B['proved']} antipode={B['antipode_hit']}", flush=True)
    C = prove_ns_lines()
    print(
        f"  C ns-lines constant: {C['proved']} chi+1={C['chi_line_sum_plus_one']} "
        f"sq_varies={C['square_dir_varies']}",
        flush=True,
    )
    D = prove_open()
    print(f"  D floor: {D['proved']} phi_F={D['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.305",
        "title": "Var(hatN)=q²/4; ns-lines |D∩ℓ|=(p−1)/2; E[S²] names Q++ if E[(∑n²)²] named",
        "proved": {
            "hatN_variance": A["proved"],
            "Fp_orbit_ES2": B["proved"],
            "ns_lines_constant": C["proved"],
            "lambda_ge_6": False,
            "phi_F_ge_6_proved_general": D["phi_F_ge_6"],
            "Q_tau_named_in_p": False,
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
    path = ROOT / "evidence" / "e1_gmin_m4_prop15305.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
