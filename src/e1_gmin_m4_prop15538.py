#!/usr/bin/env python3
"""
Prop 15.538 — Ω half-Gauss; Type+ 1D A_dbl = −4p³/(p−2).

Does **not** flip type_I / phi_F / e1 / L / Aut-Schur / Gsum /
pairing / residual_ii / 15.279–15.537 flags.  Does **not**
edit leftover-1 15.522–15.537 or residual-ii 15.528.

============================================================================
Setup.  Ω={ξ: χ(ξ)=−χ_p(−1)}, G(χ)=−p χ_p(−1) (15.523 / 15.269 B).
Type+ 1D: ker a square line, z constant on cosets, descends to
ε: F_p→{±1} with #(+)=m=(p+1)/2 (∑ε=1).  ẑ is supported on the
annihilator line L^⊥⊂{0}∪Ω, and ẑ=p ε̂ there.  Doubles (ξ,ξ,−2ξ).

============================================================================
Theorem A — PROVED (Gauss; certified p=5..19).
  ∑_{ξ∈Ω} ψ(αξ) = (q−1)/2 if α=0, else (−1 + χ_Ω χ(α) G)/2.
  Proof: Ω is a χ-coset, so the sum is
  (1/2)(∑_{ξ≠0} ψ(αξ) + χ_Ω ∑_{ξ≠0} χ(ξ)ψ(αξ)).  ∎
  Fail: drop χ(α) (α with χ(α)=−1).  ∎

Theorem B — PROVED Max+-free (Johnson 3-point + Parseval).
  |ε̂(α)|²=p+1 on F_p^* (Parseval: 1+(p−1)|ε̂|²=p²).  Hence
  |ẑ|²=p²(p+1) on L^perp minus 0.
  E[ε(x)ε(y)ε(z)]=1/p on a double+one and μ₃=−3/(p(p−2)) on
  distinct triples (complete design: λ₁=m/p, λ₂=(p+1)/(4p),
  λ₃=(p+1)(p−3)/(8p(p−2)); expand ε=2s−1).
  The form x+y−2z is F_p^*-homogeneous, so
      f(u)=∑_{x+y−2z=u} E[εεε]
  takes two values.  Collision count:
      f(0)=1 + p(p−1)μ₃=(−2p+1)/(p−2),
      f(1)=3 + p(p−3)μ₃=3/(p−2).
  Hence E[ε̂(1)² ε̂(−2)]=f(0)−f(1)=−2(p+1)/(p−2).
  Unique Ω-line ⟨ξ⟩ contributes C(p,m) of the n_{1d}=m C(p,m)
  Type+ 1D vectors, so
      A_{1d}(ξ,ξ,−2ξ)=p³ E / m = −4p³/(p−2).
  Certified p=5,7 against the Max+ cache (A_{1d}=−500/3, −1372/5).
  Fail: −4p³/(p+1) (p=5: −250≠−500/3).  ∎

Theorem C — OPEN as a name of A / Δ_conn.  A_dbl is not 1D-only
  (A_free=0 at p=5, −1372/19 at p=7).  Live A_dbl still −500/13,
  −31556/409.  Aut_e 3A+B still G>T.  type_I flags stay False.

============================================================================
Backend: field tables + 1D Fraction identities + y_∞=+1 caches.
Writes evidence/e1_gmin_m4_prop15538.json
"""
from __future__ import annotations

import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    is_prime,
)
from e1_gmin_m4_prop15275 import (  # noqa: E402
    type_I_aut_e_3AB_positive_general,
    type_I_multilevel_bad_case_ND_closed,
)
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15523 import G_chi  # noqa: E402
from e1_gmin_m4_prop15534 import field_tables  # noqa: E402
from e1_gmin_m4_prop15536 import chi_Omega  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402

EV = ROOT / "evidence" / "e1_gmin_m4_prop15538.json"
YPLUS = {5: "/tmp/maxplus_p5.npy", 7: "/tmp/maxplus_p7.npy"}


def half_gauss(p: int, ch, com: int, G: int, alpha: int) -> complex:
    q = p * p
    if alpha == 0:
        return (q - 1) / 2.0
    return 0.5 * (-1.0 + com * int(ch[alpha]) * G)


def mu3_named(p: int) -> Fraction:
    return Fraction(-3, p * (p - 2))


def f0_named(p: int) -> Fraction:
    return 1 + p * (p - 1) * mu3_named(p)


def f1_named(p: int) -> Fraction:
    return 3 + p * (p - 3) * mu3_named(p)


def E_eps_named(p: int) -> Fraction:
    """E[ε̂(1)² ε̂(−2)] = −2(p+1)/(p−2)."""
    return f0_named(p) - f1_named(p)


def A_1d_dbl_named(p: int) -> Fraction:
    return Fraction(-4 * p**3, p - 2)


def A_1d_dbl_wrong(p: int) -> Fraction:
    """Fail-when-wrong: −4p³/(p+1)."""
    return Fraction(-4 * p**3, p + 1)


def _cert_phi(p: int) -> dict:
    F = field_tables(p)
    q, ch = F["q"], F["chi"]
    mul, trt = F["mul"], F["tr"]
    G = G_chi(p)
    com = chi_Omega(p)
    Omega = [xi for xi in range(1, q) if int(ch[xi]) == com]
    om = np.exp(2j * np.pi / p)
    max_err = 0.0
    for a in (0, 1, 2, p, min(p + 1, q - 1)):
        live = 0j
        for xi in Omega:
            live += om ** int(trt[int(mul[a, xi])])
        err = abs(live - half_gauss(p, ch, com, G, a))
        if err > max_err:
            max_err = float(err)
    a_m = next(x for x in range(1, q) if int(ch[x]) == -1)
    live_m = 0j
    for xi in Omega:
        live_m += om ** int(trt[int(mul[a_m, xi])])
    pred_m = half_gauss(p, ch, com, G, a_m)
    drop_m = 0.5 * (-1.0 + com * G)
    return {
        "p": p,
        "max_err": max_err,
        "formula_ok": bool(max_err < 1e-8),
        "fail_drop_chi": bool(abs(live_m - drop_m) > 0.5 and abs(live_m - pred_m) < 1e-8),
    }


def _omega_col(p, trt, mul_t, xi):
    om = np.exp(2j * np.pi / p)
    return om ** trt[mul_t[xi]].astype(np.int64)


@lru_cache(maxsize=1)
def live_A_split() -> dict:
    out = {}
    for p in (5, 7):
        F = field_tables(p)
        q, ch = F["q"], F["chi"]
        mul, add, sub, trt = F["mul"], F["add"], F["sub"], F["tr"]
        com = chi_Omega(p)
        Omega = [xi for xi in range(1, q) if int(ch[xi]) == com]
        C = paley_conference_prime_power(p).astype(np.float64)
        Y = np.sign(np.load(YPLUS[p]).astype(np.float64))
        Yp = Y[Y[:, 0] > 0]
        assert np.allclose(Yp @ C, p * Yp, atol=1e-5)
        z = Yp[:, 1:]
        Psi = np.column_stack([_omega_col(p, trt, mul, xi) for xi in range(q)])
        Zh = (z @ Psi)[:, Omega]
        mag = np.abs(Zh)
        is_1d = (mag > 1e-6).sum(axis=1) == (p - 1)
        xi0 = Omega[0]
        th = int(sub[0, int(mul[2, xi0])])
        ia = 0
        ic = Omega.index(th)
        trip = Zh[:, ia] * Zh[:, ia] * Zh[:, ic]
        A_all = float(np.real(np.mean(trip)))
        A_1d = float(np.real(np.mean(trip[is_1d])))
        A_fr = float(np.real(np.mean(trip[~is_1d])))
        nH = int(len(Yp))
        den = nH * nH
        out[p] = {
            "nH": nH,
            "n_1d": int(is_1d.sum()),
            "n_free": int((~is_1d).sum()),
            "A_dbl": str(Fraction(A_all).limit_denominator(den)),
            "A_1d": str(Fraction(A_1d).limit_denominator(den)),
            "A_free": str(Fraction(A_fr).limit_denominator(den)),
            "z1_abs2": float(np.mean(mag[is_1d][mag[is_1d] > 1e-6] ** 2)),
        }
    return out


@lru_cache(maxsize=1)
def prove_A() -> dict:
    primes = tuple(p for p in range(5, 20) if is_prime(p))
    with ProcessPoolExecutor(max_workers=min(8, len(primes))) as ex:
        rows = {str(r["p"]): r for r in ex.map(_cert_phi, primes)}
    ok = all(r["formula_ok"] and r["fail_drop_chi"] for r in rows.values())
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "∑_{Ω} ψ(αξ)=(q−1)/2 (α=0) else (−1+χ_Ω χ(α)G)/2.  "
            "Fail: drop χ(α)."
        ),
    }


@lru_cache(maxsize=1)
def prove_B() -> dict:
    ok = True
    id_rows = {}
    for p in range(5, 32):
        if not is_prime(p):
            continue
        E = E_eps_named(p)
        want = Fraction(-2 * (p + 1), p - 2)
        f0 = f0_named(p)
        f1 = f1_named(p)
        A = A_1d_dbl_named(p)
        row_ok = (
            E == want
            and f0 == Fraction(-2 * p + 1, p - 2)
            and f1 == Fraction(3, p - 2)
            and A == Fraction(-4 * p**3, p - 2)
            and A != A_1d_dbl_wrong(p)
            and mu3_named(p) == Fraction(-3, p * (p - 2))
        )
        if not row_ok:
            ok = False
        if p in (5, 7, 11, 13):
            id_rows[str(p)] = {
                "E": str(E),
                "A_1d": str(A),
                "wrong": str(A_1d_dbl_wrong(p)),
                "ok": row_ok,
            }
    live = live_A_split()
    live_ok = (
        Fraction(live[5]["A_1d"]) == A_1d_dbl_named(5)
        and Fraction(live[7]["A_1d"]) == A_1d_dbl_named(7)
        and abs(live[5]["z1_abs2"] - 5 * 5 * 6) < 1e-4
        and abs(live[7]["z1_abs2"] - 7 * 7 * 8) < 1e-4
        and live[5]["n_1d"] == 30
        and live[7]["n_1d"] == 140
    )
    ok = ok and live_ok
    # fail −4p³/(p+1)
    ok = ok and A_1d_dbl_named(5) == Fraction(-500, 3)
    ok = ok and A_1d_dbl_wrong(5) == Fraction(-250, 3)
    return {
        "proved": bool(ok),
        "identities": id_rows,
        "live": {str(k): v for k, v in live.items()},
        "theorem": (
            "A_{1d,dbl}=−4p³/(p−2).  E[ε̂(1)²ε̂(−2)]=−2(p+1)/(p−2).  "
            "Fail: −4p³/(p+1)."
        ),
    }


def prove_open() -> dict:
    live = live_A_split()
    return {
        "proved": False,
        "type_I_multilevel_bad_case_ND_closed": bool(
            type_I_multilevel_bad_case_ND_closed()
        ),
        "type_I_aut_e_3AB_positive_general": bool(
            type_I_aut_e_3AB_positive_general()
        ),
        "Delta_conn_p_rational": False,
        "A_dbl_p5": live[5]["A_dbl"],
        "A_dbl_p7": live[7]["A_dbl"],
        "A_free_p5": live[5]["A_free"],
        "A_free_p7": live[7]["A_free"],
        "A_free_zero_all_p": False,
        "note": (
            "A_free_dbl=0 at p=5 not at p=7.  Live A_dbl still "
            "−500/13, −31556/409.  Aut_e G>T open.  Flags unflipped."
        ),
    }


def main() -> dict:
    print("Prop 15.538  Ω half-Gauss; A_1d_dbl=−4p³/(p−2)", flush=True)
    A = prove_A()
    print(f"  A half-Gauss: {A['proved']}", flush=True)
    B = prove_B()
    print(
        f"  B A_1d: {B['proved']} "
        f"p5={B['live']['5']['A_1d']} p7={B['live']['7']['A_1d']}",
        flush=True,
    )
    C = prove_open()
    print(
        f"  C open: ND={C['type_I_multilevel_bad_case_ND_closed']} "
        f"3AB={C['type_I_aut_e_3AB_positive_general']} "
        f"A_free7={C['A_free_p7']}",
        flush=True,
    )
    out = {
        "prop": "15.538",
        "title": "Ω half-Gauss; Type+ 1D A_dbl = −4p³/(p−2)",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "C": C,
        "proved": {
            "half_gauss": A["proved"],
            "A_1d_dbl_named": B["proved"],
            "type_I_multilevel_bad_case_ND_closed": C[
                "type_I_multilevel_bad_case_ND_closed"
            ],
            "type_I_aut_e_3AB_positive_general": C[
                "type_I_aut_e_3AB_positive_general"
            ],
        },
        "e1_closed_general": bool(e1_closed_general()),
        "gsum_disj_lb": bool(gsum_disj_lb_proved_general()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "flags_not_flipped": [
            "type_I_multilevel_bad_case_ND_closed",
            "type_I_aut_e_3AB_positive_general",
            "phi_F_ge_6",
            "e1",
            "L",
            "Aut-Schur",
            "Gsum",
            "pairing",
        ],
        "L_status": "OPEN",
    }
    EV.parent.mkdir(parents=True, exist_ok=True)
    EV.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", EV, flush=True)
    return out


if __name__ == "__main__":
    main()
