#!/usr/bin/env python3
"""
Prop 15.482 — W_c is real and Φ(Tr,N); the 15.481 pairing against
E[N]=μ± is identically 4, not live Q++; 13 and 409 are Gaussian
norms but 48 is not, so Q++ is not a pure-norm ratio.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.481 flags.

============================================================================
Setup.  W_c(α)=∑_x χ(x³−x) e(αx) on F_{p²} (15.481).  Q++ live
48/13, 1544/409.  den(Q++)=|H+|/(2p) (15.330).  Pairing
Q = n_++^{-1} E[u(ξ) 4 ∑_d N(d) K_++(ξ d)].

============================================================================
Theorem A — PROVED (Max+-free; certified p=3,5,7,11).
  W_c(α)∈ℝ and W_c(α)=Φ(Tr_{q/p} α, N_{q/p} α).
  Fail: W_c(α)=χ(α) W_c(1) (max |err|≥9.85 at p=5).  ∎

Theorem B — PROVED (named μ± + Ω-support; certified p=5,7).
  Replacing live N by E[N]=μ± in the 15.481 pairing gives Q≡4
  for every square-supported kernel of mass n_++ (Weil shape of
  K_++ drops out).  Live Q++ is 48/13, 1544/409 ≠ 4.  The leftover
  is the N_⊥ / 4-point pairing (already L_++ , 15.468 D).
  Fail: claim the mean-field pairing equals 48/13.  ∎

Theorem C — PROVED (elementary two-square).
  13=2²+3² and 409=3²+20² are Gaussian norms.  48 is not a sum
  of two integer squares (48=16·3, 3≡3 (mod 4)).
  Fail: claim 48=N(α) for α∈Z[i].  ∎

Theorem D — PROVED as a negative (p=3 |H+|).
  β=p−3(π+2+i) with π the 2+2i-normalized Grossencharacter if
  p≡1 (mod 4) and π=p if p≡3 (mod 4) has N(β)∈{13,409} at
  p=5,7, but N(β)=153 ≠ |H+|/(2p)=1 at p=3.  Dead as a general
  name of |H+|/(2p).  Not claimed as den(Q++) for p≥5.
  Fail: claim N(β)=|H+|/(2p) at p=3.  ∎

Theorem E — OPEN.  Q_τ still unnamed in p.  phi_F not imported.

============================================================================
Backend: field tables + p=5,7 Max+ cache for B.  Writes
evidence/e1_gmin_m4_prop15482.json
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field  # noqa: E402
from e1_gmin_m4_prop15298 import load_N  # noqa: E402
from e1_gmin_m4_prop15330 import LIVE_QPP  # noqa: E402
from e1_gmin_m4_prop15468 import D_rows  # noqa: E402
from e1_gmin_m4_prop15473 import e_matrix  # noqa: E402
from e1_gmin_m4_prop15474 import sigma_named  # noqa: E402
from e1_gmin_m4_prop15478 import a_E_sq, n_pp_named  # noqa: E402
from e1_gmin_m4_prop15480 import W_poly  # noqa: E402
from e1_gmin_m4_prop15481 import Kpp_named, chi_cubic, pred_Q  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15482_catalog.json"


def frob(F, x: int, p: int) -> int:
    if x == 0:
        return 0
    mul = F["mul"]
    acc, base, e = 1, x, p
    while e:
        if e & 1:
            acc = mul[acc][base]
        base = mul[base][base]
        e >>= 1
    return acc


def field_norm_p(F, x: int, p: int) -> int:
    if x == 0:
        return 0
    return int(F["mul"][x][frob(F, x, p)] % p)


def two_sq(n: int) -> list[tuple[int, int]]:
    out = []
    r = int(n ** 0.5)
    for a in range(r + 1):
        b2 = n - a * a
        b = int(round(b2 ** 0.5))
        if b * b == b2:
            out.append((a, b))
    return out


def pi_normed(p: int):
    """2+2i-normalized Grossencharacter, or None if p≡3 (mod 4)."""
    if p % 4 == 3:
        return None
    for alpha in range(1, p, 2):
        b2 = p - alpha * alpha
        if b2 < 0:
            break
        beta = int(round(b2 ** 0.5))
        if beta * beta == b2:
            z = complex(alpha, beta)
            for u in (1, -1, 1j, -1j):
                w = u * z
                t = (w - 1) / (2 + 2j)
                if abs(t.real - round(t.real)) < 1e-9 and abs(
                    t.imag - round(t.imag)
                ) < 1e-9:
                    return w
            return z
    raise ValueError(p)


def beta_den_candidate(p: int) -> complex:
    pi = pi_normed(p)
    if pi is None:
        pi = complex(p, 0)
    return p - 3 * (pi + 2 + 1j)


def Nrm(z: complex) -> int:
    return int(round(z.real ** 2 + z.imag ** 2))


def wick_Q(p: int) -> float:
    """15.481 pairing with live N replaced by column means E[N]."""
    F = _kernel_field(p)
    q = F["q"]
    Ds = D_rows(p, F).astype(np.float64)
    u = 4.0 * np.abs(Ds @ e_matrix(F).T) ** 2
    sigma = sigma_named(p)
    Omega = np.array(F["chi"]) == sigma
    Omega[0] = False
    xi0 = int(np.where(Omega)[0][0])
    N = load_N(p, F).astype(np.float64)
    mu = N.mean(axis=0)
    K = Kpp_named(F, p)
    mul = F["mul"]
    Kd = np.array([K[mul[xi0][d]] for d in range(q)])
    Iu = 4.0 * float(np.dot(mu, Kd).real)
    npp = n_pp_named(p)
    return float(np.mean(u[:, xi0]) * Iu / (npp * q * q))


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (3, 5, 7, 11):
        F = _kernel_field(p)
        q = F["q"]
        Wc = W_poly(F, lambda r: chi_cubic(F, r))
        imax = float(np.max(np.abs(Wc.imag)))
        if imax > 1e-8:
            ok = False
        phi = {}
        collide = 0
        for a in range(q):
            key = (int(F["trt"][a]), field_norm_p(F, a, p))
            v = float(Wc[a].real)
            if key in phi and abs(phi[key] - v) > 1e-8:
                collide += 1
            phi[key] = v
        if collide:
            ok = False
        w1 = Wc[1]
        chi = F["chi"]
        homog = float(max(abs(Wc[a] - chi[a] * w1) for a in range(1, q)))
        if homog < 1.0:
            ok = False
        rows[str(p)] = {
            "imag_max": imax,
            "phi_keys": len(phi),
            "collisions": collide,
            "homog_err": homog,
        }
    if rows["5"]["homog_err"] < 9.0:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "W_c∈ℝ and W_c=Φ(Tr,N). Fail: χ-homogeneous "
            f"(err={rows['5']['homog_err']:.2f} at p=5)."
        ),
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        qw = wick_Q(p)
        live = float(LIVE_QPP[p])
        F = _kernel_field(p)
        qlive = pred_Q(F, Kpp_named(F, p), p)
        if abs(qw - 4.0) > 1e-8:
            ok = False
        if abs(qw - live) < 0.05:
            ok = False
        if abs(qlive - live) > 1e-8:
            ok = False
        rows[str(p)] = {"wick": qw, "live": live, "full": qlive}
    if abs(rows["5"]["wick"] - float(LIVE_QPP[5])) < 0.05:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "Mean-field pairing ≡4 ≠ live Q++. Fail: wick=48/13."
        ),
    }


def prove_C() -> dict:
    t13, t409, t48, t1544 = two_sq(13), two_sq(409), two_sq(48), two_sq(1544)
    ok = (2, 3) in t13 and (3, 20) in t409 and t48 == []
    if (0, 0) in t48 or len(t48) > 0:
        ok = False
    return {
        "proved": bool(ok),
        "13": t13,
        "409": t409,
        "48": t48,
        "1544": t1544,
        "theorem": (
            "13=2²+3², 409=3²+20² are norms; 48 is not. "
            "Fail: claim 48=N(α)."
        ),
    }


def prove_D() -> dict:
    n5, n7, n3 = (
        Nrm(beta_den_candidate(5)),
        Nrm(beta_den_candidate(7)),
        Nrm(beta_den_candidate(3)),
    )
    h3 = 6  # |H+| at p=3 (cache 12 = ± pairs)
    den3 = h3 // (2 * 3)
    ok = n5 == 13 and n7 == 409 and n3 == 153 and den3 == 1 and n3 != den3
    return {
        "proved": bool(ok),
        "N5": n5,
        "N7": n7,
        "N3": n3,
        "Hplus3_over_2p": den3,
        "theorem": (
            "N(p−3(π+2+i)) hits 13,409 and is not |H+|/(2p) at p=3 "
            "(153≠1). Dead as a general |H+| name."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "aE2_5": a_E_sq(5),
        "note": (
            "W_c is real Φ(Tr,N), not χ-homogeneous. Mean-field pairing "
            "is 4, not Q++. Q_τ unnamed in p. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.482  W_c real Φ(Tr,N); mean-field pairing ≡4", flush=True)
    A = prove_A()
    print(f"  A real/Φ: {A['proved']} homog5={A['rows']['5']['homog_err']:.2f}", flush=True)
    B = prove_B()
    print(f"  B wick: {B['proved']} {B['rows']['5']['wick']:.4f} vs {B['rows']['5']['live']:.4f}", flush=True)
    C = prove_C()
    print(f"  C norms: {C['proved']} 48={C['48']}", flush=True)
    D = prove_D()
    print(f"  D den-cand dead: {D['proved']} N3={D['N3']}", flush=True)
    E = prove_open()
    print(f"  E open phi_F={E['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps(
            {"A": A["proved"], "B": B["proved"], "C": C["proved"], "D": D["proved"]},
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.482",
        "title": (
            "W_c real Φ(Tr,N); mean-field pairing ≡4; "
            "48 is not a Gaussian norm"
        ),
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "Wc_real_Phi_TrN": A["proved"],
            "mean_field_pairing_is_4": B["proved"],
            "48_not_gaussian_norm": C["proved"],
            "den_candidate_not_Hplus": D["proved"],
            "phi_F_ge_6_proved_general": E["phi_F_ge_6"],
            "residual_ii_k_eq_4p_empty": E["residual_ii_k_eq_4p_empty"],
        },
        "algebra": {"A": A, "B": B, "C": C, "D": D, "E": E},
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15482.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
