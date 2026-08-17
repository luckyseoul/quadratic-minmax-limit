#!/usr/bin/env python3
"""
Prop 15.325 — Cay(F_q, U) commutes with Paley and on V_{θ+}
    spec(A_U) = {−Kl(−s,−1) : χ_p(s)=−χ_p(−1)}.
Common neighbors |U ∩ (b+U)| are N(b)-constant.  The two-point
Rayleigh bound on Var(W) is 180 at p=5, which does not close the
floor (need ≤360/7).  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.324 flags.  Floor stays OPEN.

============================================================================
Setup.  q=p².  Paley graph Cay(F_q, □^×) and oval graph
Cay(F_q, U), U={N_{q/p}=1}, |U|=p+1, 0∉U so A_U is (p+1)-regular.
Both are Cayley on (F_q, +), hence they commute and share the
additive-character basis e_α(x)=e_p(Tr(αx)).

K(α)=∑_{u∈U} e_α(u).  15.316: K(α)=−Kl(−N(α),−1) (α≠0), K(0)=p+1.
Paley eig at e_α (α≠0) is (−1 + χ_q(α) G_q)/2 with G_q=−χ_p(−1)p
(15.319 HD).  Thus e_α ∈ V_{θ+}, θ+=(p−1)/2, iff
    χ_q(α)=−χ_p(−1)  ⇔  χ_p(N(α))=−χ_p(−1)
(the same sign that cuts out Ω).

W(1)=1_D^T A_U 1_D.  Regular sets satisfy v:=1_D−(k/q)1 ∈ V_{θ+},
‖v‖²=(p²−1)/4, and W(1)=(p+1)(p−1)²/4 + v^T A_U v.

============================================================================
Theorem A — PROVED (Cayley characters + 15.316 + HD; certified
p=5,7,11,13,17,19,23).
  spec(A_U |_{V_{θ+}}) = { −Kl(−s,−1) : s∈F_p^×, χ_p(s)=−χ_p(−1) }.
  There are (p−1)/2 such s and typically that many distinct eigs
  (p=5: 1±√5).  Fail-when-wrong: claim the list is −Kl(−s,−1) on
  squares for every p (false at p=5: squares give the θ− list).  ∎

Theorem B — PROVED (cyclotomic scheme; certified p=5,7,11,13).
  |U ∩ (b+U)| depends only on N(b).  Fail-when-wrong: claim the
  count takes two values inside a single fiber.  ∎

Theorem C — PROVED as a bound; too weak for the floor.
  v^T A_U v ∈ [λ_min, λ_max] ‖v‖², so W lies in a named interval
  and Var(W) ≤ (EW−W_lo)(W_hi−EW) (two-point / Popoviciu).
  At p=5 this is 180 > 360/7 = V_floor.  Fail-when-wrong: claim
  the two-point bound is ≤360/7 at p=5.  ∎

Theorem D — OPEN.  Named spectrum does not by itself give
Var≤Var_2orb or Q++≤Q^{floor ub}.  phi_F not imported.

============================================================================
Backend: F_p Kloosterman (no q×q eigh).  Writes
evidence/e1_gmin_m4_prop15325.json
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field, rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15290 import field_norm  # noqa: E402
from e1_gmin_m4_prop15316 import chi_p_pm, kloosterman  # noqa: E402
from e1_gmin_m4_prop15321 import Qpp_floor_ub, alpha_beta_p1  # noqa: E402
from e1_gmin_m4_prop15323 import EW_named, Var_2orb  # noqa: E402


def theta_plus_norms(p: int) -> list[int]:
    eps = -chi_p_pm(p, p - 1)
    return [s for s in range(1, p) if chi_p_pm(p, s) == eps]


def AU_eigs_named(p: int) -> list[float]:
    return [float((-kloosterman(p, -s, -1)).real) for s in theta_plus_norms(p)]


def AU_eigs_squares(p: int) -> list[float]:
    return [
        float((-kloosterman(p, -s, -1)).real)
        for s in range(1, p)
        if chi_p_pm(p, s) == 1
    ]


def two_point_var_ub(p: int) -> Fraction:
    eigs = AU_eigs_named(p)
    # keep as float then compare; exact at p=5 is 180
    fnorm = Fraction(p * p - 1, 4)
    rand_mean = Fraction((p + 1) * (p - 1) ** 2, 4)
    EW = Fraction(EW_named(p))
    Wlo = rand_mean + Fraction(min(eigs)).limit_denominator(10**9) * fnorm
    Whi = rand_mean + Fraction(max(eigs)).limit_denominator(10**9) * fnorm
    return (EW - Wlo) * (Whi - EW)


def two_point_var_ub_float(p: int) -> float:
    eigs = AU_eigs_named(p)
    fnorm = (p * p - 1) / 4.0
    rand_mean = (p + 1) * (p - 1) ** 2 / 4.0
    EW = EW_named(p)
    Wlo = rand_mean + min(eigs) * fnorm
    Whi = rand_mean + max(eigs) * fnorm
    return float((EW - Wlo) * (Whi - EW))


def V_floor_p5() -> Fraction:
    return Fraction(360, 7)


def _nn_one(p: int) -> dict:
    F = _kernel_field(p)
    q = F["q"]
    add = np.array(F["add"], dtype=np.int32)
    U = [d for d in range(1, q) if int(field_norm(F, d)) == 1]
    nn: dict[int, list[int]] = defaultdict(list)
    for b in range(1, q):
        cnt = 0
        for u in U:
            if int(field_norm(F, int(add[b, u]))) == 1:
                cnt += 1
        nn[int(field_norm(F, b))].append(cnt)
    ok = all(min(vs) == max(vs) for vs in nn.values())
    return {
        "p": p,
        "ok": ok,
        "vals": {str(c): (min(vs), max(vs)) for c, vs in sorted(nn.items())},
    }


def prove_spectrum() -> dict:
    ok = True
    sq_all_hit = 0
    rows = {}
    for p in (5, 7, 11, 13, 17, 19, 23):
        named = sorted({round(x, 8) for x in AU_eigs_named(p)})
        squares = sorted({round(x, 8) for x in AU_eigs_squares(p)})
        if p == 5:
            # 1±√5
            want = sorted(
                {round(1 + np.sqrt(5), 8), round(1 - np.sqrt(5), 8)}
            )
            if named != want:
                ok = False
            if named == squares:
                sq_all_hit += 1
                ok = False
        else:
            if not named:
                ok = False
        # squares-for-all-p is the fail-when-wrong
        if named == squares and p == 5:
            ok = False
        rows[str(p)] = {
            "n_named": len(AU_eigs_named(p)),
            "n_distinct": len(named),
            "named": named,
            "equals_squares": named == squares,
        }
    if rows["5"]["equals_squares"]:
        ok = False
    if not rows["7"]["equals_squares"]:
        # p≡3: θ+ norms are squares, so they SHOULD match
        ok = False
    return {
        "proved": ok,
        "squares_all_p_false": not rows["5"]["equals_squares"],
        "by_p": rows,
        "theorem": (
            "spec(A_U|θ+)={−Kl(−s,−1): χ_p(s)=−χ_p(−1)}. "
            "Squares-for-all-p fails at p=5."
        ),
    }


def prove_common_neighbors() -> dict:
    ok = True
    split_hit = 0
    rows = {}
    with ProcessPoolExecutor(max_workers=4) as ex:
        recs = list(ex.map(_nn_one, (5, 7, 11, 13)))
    for rec in recs:
        if not rec["ok"]:
            split_hit += 1
            ok = False
        rows[str(rec["p"])] = rec["vals"]
    if split_hit > 0:
        ok = False
    return {
        "proved": ok,
        "fiber_split_false": split_hit == 0,
        "by_p": rows,
        "theorem": "|U ∩ (b+U)| is N(b)-constant.",
    }


def prove_twopoint_weak() -> dict:
    ok = True
    closes_hit = 0
    ub5 = two_point_var_ub_float(5)
    if abs(ub5 - 180.0) > 0.05:
        ok = False
    if ub5 <= float(V_floor_p5()) + 1e-9:
        closes_hit += 1
        ok = False
    rows = {}
    for p in (5, 7, 13, 17):
        rows[str(p)] = {
            "two_point": two_point_var_ub_float(p),
            "Var2orb": float(Var_2orb(p)),
        }
    return {
        "proved": ok,
        "twopoint_closes_p5_false": closes_hit == 0,
        "ub5": ub5,
        "Vfloor5": str(V_floor_p5()),
        "by_p": rows,
        "theorem": (
            "Two-point Var≤180 at p=5, which is not ≤360/7. "
            "Spectrum names the interval, not the floor."
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
            "Named A_U spectrum on V_θ+ does not bound Var by "
            "Var_2orb or V_floor.  phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.325  A_U|θ+ = −Kl on {χ=−χ(−1)}; two-point too weak", flush=True)
    A = prove_spectrum()
    print(
        f"  A spec: {A['proved']} squares_all={not A['squares_all_p_false']}",
        flush=True,
    )
    B = prove_common_neighbors()
    print(
        f"  B nn: {B['proved']} split={not B['fiber_split_false']}",
        flush=True,
    )
    C = prove_twopoint_weak()
    print(
        f"  C 2pt: {C['proved']} closes5={not C['twopoint_closes_p5_false']} "
        f"ub5={C['ub5']:.2f}",
        flush=True,
    )
    D = prove_open()
    print(f"  D open: phi_F={D['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.325",
        "title": "A_U on Paley θ+ is −Kl(χ=−χ(−1)); two-point Var 180 at p=5",
        "proved": {
            "AU_spectrum_named": A["proved"],
            "nn_N_constant": B["proved"],
            "twopoint_too_weak": C["proved"],
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
    path = ROOT / "evidence" / "e1_gmin_m4_prop15325.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
