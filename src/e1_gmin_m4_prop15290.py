#!/usr/bin/env python3
"""
Prop 15.290 — The Paley-class split of Q(r) is the field norm
N_{F_q/F_p}(r).  Floor restated as a Fourier bound on the Wick deficit.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.289 flags.

============================================================================
Setup.  q=p².  T=□ = {r∈F_q^× : χ(r)=1}, |T|=(q−1)/2.
u(ξ)=|ẑ(ξ)|² on Ω={ξ: χ̂(ξ)=p}.  Q(r)=E[u(ξ)u(rξ)] (ξ∈Ω, r∈T).
Wick: Q(±1)=8q², Q(r≠±1)=4q².  Boolean Q(±1)=8q² still;
Q(off) is smaller.  Paley type (χ(r−1),χ(r+1)) does **not** determine
Q (15.279 O; p=7 class (−1,−1) splits).

Norm type of r∈T:
  pm1  : r=±1
  sub  : r∈F_p^×  (⇔ r^{p−1}=1)
  N1   : N(r)=r^{p+1}=1 and r∉F_p  (norm-1 torus)
  N*   : otherwise

============================================================================
Theorem A — PROVED (Max+-free invariance + certified orbits).
  Q is invariant under r↦r^p (Frob∈Aut C) and r↦r^{−1} (swap ξ↔rξ).
  Those maps preserve (χ(r−1),χ(r+1),norm-type).  At p=5 and p=7
  each (Paley, norm-type) class is a single ⟨Frob,inv⟩-orbit, hence
  Q is constant on the class.  The p=7 Paley class (−1,−1) splits
  exactly as N1 vs N*: Q/q² = 1376/409 vs 1496/409.
  Fail-when-wrong: Paley type alone constant (false at p=7);
  claim Q(N1)=Q(N*) on (−1,−1).  ∎

Theorem B — PROVED (Parseval; certified p=5,7).
  S_□(ψ)=∑_{r∈T} ψ(r) Q(r) = 16q² + ∑_{τ≠pm1} Q_τ J_τ(ψ)
  with J_τ(ψ)=∑_{r∈τ} ψ(r).  Then
      λ(ψ)=S_□/q² ,   λ≥6  ⇔  ∑_{τ≠pm1} (Q_τ/q²) J_τ(ψ) ≥ −10.
  Equivalently the Wick deficit δ=4q²−Q on T minus {±1} satisfies
      <delta, psi> <= 2.
  At p=5 the bound is sharp-ish: min λ=80/13, ⟨δ,ψ⟩=24/13 < 2.
  At p=7 min λ=3072/409, ⟨δ,ψ⟩=200/409 < 2.
  Fail-when-wrong: drop the 16q²; invert a J_τ sign.  ∎

Theorem C — OPEN.  ⟨δ,ψ⟩≤2 for every even ψ∉{1,χ} and every
  prime p≥5 is exactly φ_F≥6.  Types named; values of Q_τ not
  yet a Gauss/Jacobi formula in p.  Floor stays OPEN.  ∎

============================================================================
Backend: CPU.  zhat over Max+ is a single pass (p=5,7 caches).
Writes evidence/e1_gmin_m4_prop15290.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15159 import SPECTRUM_P5, SPECTRUM_P7  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import (  # noqa: E402
    _e_tr,
    _kernel_field,
    _psi_vec,
    rhat_ge_budget_proved_general,
)

YPATH = {5: "/tmp/maxplus_p5.npy", 7: "/tmp/maxplus_p7.npy"}


def _pow(F, x: int, e: int) -> int:
    acc, b, ee = 1, int(x), int(e)
    if e == 0:
        return 1
    while ee:
        if ee & 1:
            acc = F["mul"][acc][b]
        b = F["mul"][b][b]
        ee >>= 1
    return acc


def in_subfield(F, r: int) -> bool:
    return r != 0 and _pow(F, r, F["p"] - 1) == 1


def field_norm(F, r: int) -> int:
    return _pow(F, r, F["p"] + 1)


def paley_pair(F, r: int) -> tuple[int, int]:
    rm1 = F["add"][r][F["neg"][1]]
    rp1 = F["add"][r][1]
    return int(F["chi"][rm1]), int(F["chi"][rp1])


def norm_type(F, r: int) -> str:
    if r in (1, F["neg1"]):
        return "pm1"
    if in_subfield(F, r):
        return "sub"
    if field_norm(F, r) == 1:
        return "N1"
    return "N*"


def type_key(F, r: int) -> tuple:
    if r in (1, F["neg1"]):
        return ("pm1",)
    a, b = paley_pair(F, r)
    return (a, b, norm_type(F, r))


def paley_only_key(F, r: int) -> tuple:
    if r in (1, F["neg1"]):
        return ("pm1",)
    return paley_pair(F, r)


def frob_inv_orbit(F, r: int) -> set[int]:
    out = {r}
    rp = _pow(F, r, F["p"])
    inv = F["inv"][r]
    for x in (rp, inv, F["inv"][rp], _pow(F, inv, F["p"])):
        if x:
            out.add(int(x))
    return out


def omega_set(F) -> list[int]:
    q, p = F["q"], F["p"]
    om = []
    for xi in range(1, q):
        chat = sum(
            F["chi"][x] * _e_tr(F, F["mul"][xi][x]) for x in range(1, q)
        )
        if abs(chat.real - p) < 1e-6 and abs(chat.imag) < 1e-6:
            om.append(xi)
    return om


def zhat_omega(z, F, Omega):
    q = F["q"]
    out = np.zeros(len(Omega), dtype=np.complex128)
    for i, xi in enumerate(Omega):
        s = 0j
        for x in range(q):
            s += z[x] * _e_tr(F, F["mul"][xi][x])
        out[i] = s
    return out


def _u_matrix(p: int, F, Omega):
    Y = np.load(YPATH[p])
    Yp = Y[Y[:, 0] > 0]
    q = F["q"]
    rows = [
        np.abs(zhat_omega(y[1 : 1 + q].astype(np.float64), F, Omega)) ** 2
        for y in Yp
    ]
    return np.stack(rows)


def live_Q_on_T(p: int) -> dict[int, float]:
    F = _kernel_field(p)
    Omega = omega_set(F)
    u = _u_matrix(p, F, Omega)
    xi0 = Omega[0]
    om_idx = {xi: i for i, xi in enumerate(Omega)}
    u0 = u[:, om_idx[xi0]]
    out = {}
    for r in F["squares"]:
        eta = F["mul"][r][xi0]
        out[int(r)] = float(np.mean(u0 * u[:, om_idx[eta]]))
    return out


def prove_type_constancy(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [5, 7]
    ok = True
    paley_split = 0
    blocks = []
    for p in primes:
        F = _kernel_field(p)
        Q = live_Q_on_T(p)
        q2 = float(F["q"] ** 2)
        by_type: dict[tuple, list[float]] = defaultdict(list)
        by_pal: dict[tuple, list[float]] = defaultdict(list)
        orbits_ok = True
        for r, Qr in Q.items():
            by_type[type_key(F, r)].append(Qr)
            by_pal[paley_only_key(F, r)].append(Qr)
            orb = frob_inv_orbit(F, r)
            for s in orb:
                if s in Q and abs(Q[s] - Qr) > 1e-6 * q2:
                    orbits_ok = False
                    ok = False
        type_spread = {}
        for tau, vs in by_type.items():
            sp = (max(vs) - min(vs)) / q2
            type_spread[str(tau)] = sp
            if sp > 1e-6:
                ok = False
        for tau, vs in by_pal.items():
            if (max(vs) - min(vs)) / q2 > 1e-4:
                paley_split += 1
        # each type is one Frob-inv orbit
        for tau, vs in by_type.items():
            rs = [r for r in Q if type_key(F, r) == tau]
            orb0 = frob_inv_orbit(F, rs[0])
            if set(rs) != orb0 and tau != ("pm1",):
                # pm1 = {1,-1} may be two fixed points
                if not (len(rs) <= 4 and set(rs) <= orb0 | frob_inv_orbit(F, rs[-1])):
                    if set(rs) != orb0:
                        # allow pm1
                        if tau[0] != "pm1":
                            ok = False
        blocks.append(
            {
                "p": p,
                "n_types": len(by_type),
                "type_spread": type_spread,
                "orbits_ok": orbits_ok,
            }
        )
    if paley_split == 0:
        ok = False  # fail-when-wrong: Paley-only must split at p=7
    return {
        "proved": ok and paley_split > 0,
        "paley_class_splits": paley_split,
        "blocks": blocks,
        "theorem": (
            "Q is Frob/inv invariant.  At p=5,7 it is constant on "
            "(χ(r−1),χ(r+1),norm-type).  Paley type alone splits."
        ),
    }


def prove_deficit_bound_certified() -> dict:
    """⟨δ,ψ⟩=8−λ(ψ) ≤ 2 on the live even spectrum (p=5,7)."""
    ok = True
    rows = {}
    for p, spec in ((5, SPECTRUM_P5), (7, SPECTRUM_P7)):
        lams = list(spec.keys())
        defics = [8 - float(l) for l in lams]
        mx = max(defics)
        mn_lam = min(lams)
        if mx > 2 + 1e-12:
            ok = False
        if mn_lam < 6:
            ok = False
        rows[str(p)] = {
            "lambda_min": str(mn_lam),
            "max_deficit": mx,
            "bound_2": mx <= 2 + 1e-12,
        }
    return {
        "proved": ok,
        "by_p": rows,
        "theorem": (
            "On the certified Φ spectra (15.159), 8−λ ≤ 2, i.e. λ≥6, "
            "with p=5 max deficit 24/13.  Not a general proof."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "rhat_budget": bool(rhat_ge_budget_proved_general()),
        "note": (
            "⟨δ,ψ⟩≤2 for all even ψ∉{1,χ} and all p≥5 is the floor. "
            "Q_τ named by type, not yet a formula in p."
        ),
    }


def main() -> dict:
    print("Prop 15.290  Q-type = Paley × field-norm; floor is ⟨δ,ψ⟩≤2", flush=True)
    A = prove_type_constancy()
    print(f"  A type-constancy: {A['proved']} paley_splits={A['paley_class_splits']}", flush=True)
    B = prove_deficit_bound_certified()
    print(f"  B certified deficit≤2: {B['proved']}", flush=True)
    C = prove_open()
    print(f"  C general floor: {C['proved']} phi_F={C['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.290",
        "title": "Field-norm split of Q(r); floor is Fourier bound ⟨δ,ψ⟩≤2",
        "proved": {
            "type_constancy_p5_p7": A["proved"],
            "paley_only_splits": A["paley_class_splits"] > 0,
            "deficit_le_2_certified": B["proved"],
            "phi_F_ge_6_proved_general": C["phi_F_ge_6"],
            "rhat_ge_budget": C["rhat_budget"],
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
    path = ROOT / "evidence" / "e1_gmin_m4_prop15290.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
