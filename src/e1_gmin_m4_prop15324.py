#!/usr/bin/env python3
"""
Prop 15.324 — For every Max+, the average of W(c) over square
c∈F_p^× equals E[W] exactly (2-design + χ-sum).  Odd U-hats of
N vanish; Parseval names the even energy.  |G_U|² is not
constant on Ω, so the 15.304 energy cannot drop |G_U|².
Var≤Var_2orb still OPEN.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.323 flags.  Floor stays OPEN.

============================================================================
Setup.  W(c)=∑_{N(δ)=c} N(δ).  T=□, n_f=(p−1)/2 square-norm fibers,
each of size |U|=p+1.  χ_q=χ_p∘N.  U={N=1}⊂T.

============================================================================
Theorem A — PROVED (handshaking + 15.318; certified p=5,7).
  ∑_{δ≠0} N(δ)=k(k−1) and ∑_{δ≠0} χ_q(δ) N(δ)=k(p+1)/2 are
  Max+-free, so
      ∑_{c∈F_p^{×□}} W(c) = p(p−1)(p²−1)/8 = n_f E[W]
  on every Max+.  Hence the fiber-average of W on square norms
  is exactly E[W], not merely in expectation.  Fail-when-wrong:
  claim each W(c) is constant in c (W(1) varies).  ∎

Theorem B — PROVED (N even + U cyclic; certified p=5,7).
  N(−δ)=N(δ) and −1∈U, so N|_U is even.  Odd characters of U
  have vanishing hats.  Parseval:
      ∑_{j even} E|N̂_U(j)|² = (p+1)² E[N(u)²]
  and the j=0 term is A_□.  At p=5 the two nonzero even hats
  are conjugates and split the remainder equally (750/13 each).
  Fail-when-wrong: claim an odd hat is nonzero.  ∎

Theorem C — PROVED as a profile; not a majorant.
  |G_U(j,α)|² is constant (=p) on squares at p=5, even j≠0
  (15.303), but on Ω (nonsquares) it spreads (min 1.91, max 13.09
  at j=2).  Therefore E[|G_U(ξ)|² |⟨ΔN,χ_j⟩|²] cannot be
  factored as |G_U|² E|hat|².  Fail-when-wrong: claim |G_U(2)|²
  is constant on Ω at p=5.  ∎

Theorem D — OPEN.  Fiber-averaging and U-Parseval do not bound
  Var(W) by Var_2orb=4p(p+1)(2p+1)/(5p+1).  Q++ / phi_F stay
  OPEN.

============================================================================
Backend: one N-fold per prime + field Gauss on U.  Writes
evidence/e1_gmin_m4_prop15324.json
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
from e1_gmin_m4_prop15279 import _kernel_field, rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15290 import field_norm, omega_set  # noqa: E402
from e1_gmin_m4_prop15298 import load_N  # noqa: E402
from e1_gmin_m4_prop15302 import U_powers  # noqa: E402
from e1_gmin_m4_prop15303 import GU  # noqa: E402
from e1_gmin_m4_prop15316 import chi_p_pm  # noqa: E402
from e1_gmin_m4_prop15323 import Var_2orb  # noqa: E402


def EW_named(p: int) -> int:
    return (p + 1) * p * (p - 1) // 4


def sum_W_squares_named(p: int) -> int:
    return p * (p - 1) * (p * p - 1) // 8


def fiber_W(N, F, c: int):
    idx = [d for d in range(1, F["q"]) if int(field_norm(F, d)) == c]
    return N[:, np.array(idx, dtype=np.int32)].sum(axis=1)


def prove_fiber_avg() -> dict:
    ok = True
    const_c_hit = 0
    rows = {}
    for p in (5, 7):
        F = _kernel_field(p)
        N = load_N(p, F).astype(np.float64)
        acc = np.zeros(N.shape[0], dtype=np.float64)
        per_c = {}
        for c in range(1, p):
            if chi_p_pm(p, c) != 1:
                continue
            Wc = fiber_W(N, F, c)
            acc += Wc
            per_c[c] = {
                "mean": float(np.mean(Wc)),
                "var": float(np.var(Wc)),
                "nuniq": int(len(np.unique(np.round(Wc, 6)))),
            }
        named = float(sum_W_squares_named(p))
        err = float(np.max(np.abs(acc - named)))
        if err > 1e-6:
            ok = False
        if all(v["nuniq"] == 1 for v in per_c.values()):
            const_c_hit += 1
            ok = False
        rows[str(p)] = {
            "named_sum": named,
            "err": err,
            "EW": EW_named(p),
            "n_fibers": (p - 1) // 2,
            "per_c": per_c,
        }
    if const_c_hit > 0:
        ok = False
    return {
        "proved": ok,
        "W_c_constant_false": const_c_hit == 0,
        "by_p": rows,
        "theorem": (
            "∑_{c□} W(c)=n_f E[W] on every Max+. "
            "Each W(c) still varies."
        ),
    }


def prove_parseval() -> dict:
    ok = True
    odd_hit = 0
    rows = {}
    for p in (5, 7):
        F = _kernel_field(p)
        N = load_N(p, F).astype(np.float64)
        U = U_powers(F)
        NU = N[:, np.array(U, dtype=np.int32)]
        p1 = p + 1
        ehat = {}
        omax = 0.0
        for j in range(p1):
            chi = np.exp(2j * np.pi * j * np.arange(len(U)) / p1)
            h2 = float(np.mean(np.abs((NU * chi).sum(axis=1)) ** 2))
            if j % 2 == 0:
                ehat[j] = h2
            else:
                omax = max(omax, h2)
        if omax > 1e-6:
            odd_hit += 1
            ok = False
        A = float(np.mean(NU.sum(axis=1) ** 2))
        EN2 = float(np.mean((NU**2).mean(axis=1)))
        rest = sum(v for j, v in ehat.items() if j != 0)
        if abs((A + rest) - p1 * p1 * EN2) > 1e-4:
            ok = False
        rows[str(p)] = {
            "A": A,
            "E_N2": EN2,
            "even_hats": ehat,
            "odd_max": omax,
            "parseval_ok": abs((A + rest) - p1 * p1 * EN2) < 1e-4,
        }
        if p == 5:
            if abs(ehat[2] - ehat[4]) > 1e-6:
                ok = False
            if abs(ehat[2] - 750 / 13) > 1e-4:
                ok = False
    if odd_hit > 0:
        ok = False
    return {
        "proved": ok,
        "odd_hats_false": odd_hit == 0,
        "by_p": rows,
        "theorem": (
            "Odd U-hats of N vanish.  Parseval: "
            "A_□ + ∑_{even j≠0} E|N̂_U(j)|²=(p+1)² E[N²_□]."
        ),
    }


def prove_GU_spread() -> dict:
    ok = True
    const_hit = 0
    F = _kernel_field(5)
    Omega = set(omega_set(F))
    g2s = []
    for a in Omega:
        g2s.append(abs(GU(F, 2, a)) ** 2)
    spread = float(max(g2s) - min(g2s))
    if spread < 1.0:
        const_hit += 1
        ok = False
    # squares are constant (15.303); record as contrast
    sq = []
    for a in range(1, F["q"]):
        if int(F["chi"][a]) == 1:
            sq.append(abs(GU(F, 2, a)) ** 2)
    sq_spread = float(max(sq) - min(sq))
    if sq_spread > 1e-6:
        ok = False
    return {
        "proved": ok,
        "GU_const_on_Omega_false": const_hit == 0,
        "Omega_spread": spread,
        "Omega_min": float(min(g2s)),
        "Omega_max": float(max(g2s)),
        "square_spread": sq_spread,
        "theorem": (
            "|G_U(2)|²=p on squares at p=5, but spreads on Ω. "
            "Cannot factor |G_U|² out of E|inner|²."
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
        "Var_2orb": str(Var_2orb(5)),
        "note": (
            "Fiber-average and U-Parseval do not name A_□ or bound "
            "Var by Var_2orb.  phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.324  fiber-avg W=EW; U-Parseval; |G_U|² spreads on Ω", flush=True)
    A = prove_fiber_avg()
    print(
        f"  A fiber: {A['proved']} Wc_const={not A['W_c_constant_false']}",
        flush=True,
    )
    B = prove_parseval()
    print(f"  B Parseval: {B['proved']} odd={not B['odd_hats_false']}", flush=True)
    C = prove_GU_spread()
    print(
        f"  C GU-Ω: {C['proved']} constΩ={not C['GU_const_on_Omega_false']} "
        f"spread={C['Omega_spread']:.2f}",
        flush=True,
    )
    D = prove_open()
    print(f"  D open: phi_F={D['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.324",
        "title": "Fiber-avg W=EW; U-Parseval; GU spreads on Ω; Var_2orb open",
        "proved": {
            "fiber_avg": A["proved"],
            "U_parseval": B["proved"],
            "GU_spreads_on_Omega": C["proved"],
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
    path = ROOT / "evidence" / "e1_gmin_m4_prop15324.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
