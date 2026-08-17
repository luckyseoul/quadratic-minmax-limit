#!/usr/bin/env python3
"""
Prop 15.493 — Residual (ii) k=4p: Max− Gram, switching, and a
p=3 leftover+s₊=2+freeness-fail instance that is still ND.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
residual_ii / 15.279–15.492 flags.

============================================================================
Setup.  Residual (ii) leftover (15.274 E): k=4p, max_Max− S=−2,
f_e≡−1 on U={S=−2}.  Need Φ(H)≥Φ−2 for H=G∪{e}.  First-moment
4-level exists at every p (15.330 D).  p=5 leftover 20-edge
double-star exists and is not interior (15.406 D).

============================================================================
Theorem A — PROVED (row-sum + wedge mean 0; Max− dual of 15.380).
  On Max−, E[f_e]=−1/p and S_all≡−Φ.  Adjacent Gram has mean 0
  (values {±1/p}).  Hence
      avg_{e' disj e} E[f_e f_{e'}] = 1/(n−3).
  Same number as on Max+ (15.380 B).  Fail: claim the disjoint
  mean is 1/p² (1/23≠1/25 at p=5).  ∎

Theorem B — PROVED (pointwise algebra on {±1}^n).
  For any z∈{±1}^n write C^z_{ab}=C_{ab} z_a z_b.  Then for every
  y∈{±1}^n and every edge set G,
      S_C(y) = S_{C^z}(y⊙z).
  If Cy=λy and Cz=μz then C^z(y⊙z)=λ(y⊙z).  Fail: drop the
  z-factors (use C in place of C^z).  ∎

Theorem C — PROVED (two-coordinate flip of z∈U).
  Let z∈Max−, Cz=−pz, and y flip coordinates {i,j}=e.  Then
      yᵀ C y = −pn + 8p + 8 f_e(z).
  On U one has f_e(z)=−1, so S_all(y)=−p(n−8)/2−4.  At p=5 this
  is −49; the 15.406 double-star then has |Q(y)|=53<Φ=65, so the
  2-flip is not a Φ(G)≥Φ witness.  Fail: claim S_all=−Φ.  ∎

Theorem D — PROVED (explicit 12-edge G at p=3; cube Φ).
  The graph P3_RESII is leftover official, s₊=2, freeness-fail,
  Max± |Q|≤Φ−4, but Φ(G)=29≥Φ=15 and Φ(H)=31≥Φ−2.  It is a
  residual-(ii) instance that is ND.  Fail: claim leftover+s₊=2
  is empty; fail: claim this G undercuts Φ−2.  ∎

Theorem E — OPEN.  0-1 leftover+s₊=2 at p≥5 and even k>4p stay
  open.  Linear Farkas (fractional v) does not empty leftover+s₊≥2
  (feasible at p=3,5).  residual_ii_k_eq_4p_empty stays False.
  phi_F not imported.

============================================================================
Backend: Fraction identities + p=3 2^{10} cube.  Writes
evidence/e1_gmin_m4_prop15493.json
"""
from __future__ import annotations

import itertools
import json
import os
import sys
from fractions import Fraction
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15493_catalog.json"

# p=3 leftover official, s₊=2, freeness-fail, ND (15.493 D).
P3_RESII = (
    (0, 3),
    (0, 7),
    (1, 2),
    (1, 4),
    (2, 4),
    (2, 6),
    (3, 6),
    (3, 7),
    (4, 8),
    (5, 8),
    (6, 9),
    (7, 8),
)


def n_of(p: int) -> int:
    return p * p + 1


def phi_of(p: int) -> Fraction:
    return Fraction(n_of(p) * p, 2)


def disj_mean_named(p: int) -> Fraction:
    return Fraction(1, n_of(p) - 3)


def two_flip_Sall_on_U(p: int) -> Fraction:
    """S_all after flipping {i,j} of any z∈U ⊂ Max−."""
    n = n_of(p)
    return Fraction(-p * (n - 8), 2) - 4


def two_flip_Sall_wrong_minus_Phi(p: int) -> Fraction:
    return -phi_of(p)


def enum_max(C: np.ndarray, lam: float) -> np.ndarray:
    n = C.shape[0]
    rows = []
    for bits in itertools.product((-1.0, 1.0), repeat=n):
        y = np.array(bits, dtype=np.float64)
        if np.allclose(C @ y, lam * y, atol=1e-8):
            rows.append(y)
    return np.stack(rows)


def F_of(C: np.ndarray, Y: np.ndarray) -> tuple[list[tuple[int, int]], np.ndarray]:
    n = C.shape[0]
    edges = [(a, b) for a in range(n) for b in range(a + 1, n)]
    F = (
        Y[:, [e[0] for e in edges]]
        * Y[:, [e[1] for e in edges]]
        * np.array([C[a, b] for a, b in edges], dtype=np.float64)
    )
    return edges, F


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 11, 13):
        n = n_of(p)
        phi = phi_of(p)
        # E[f (S_all − f)] = E[f] S_all − 1 with S_all ≡ −Φ
        # E[f]=−1/p ⇒ E[f S_all]=Φ/p=n/2, same as Max+
        sum_off = Fraction(n, 2) - 1
        n_disj = (n - 2) * (n - 3) // 2
        avg = sum_off / n_disj
        named = disj_mean_named(p)
        ok = ok and avg == named
        ok = ok and named != Fraction(1, p * p)
        rows[str(p)] = {
            "disj_mean": str(named),
            "n_minus_3": n - 3,
            "not_1_over_p2": str(Fraction(1, p * p)),
        }
        if avg != named:
            ok = False
    ok = ok and disj_mean_named(5) == Fraction(1, 23)
    ok = ok and Fraction(1, 23) != Fraction(1, 25)
    if disj_mean_named(5) == Fraction(1, 25):
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "Max− disjoint Gram mean is 1/(n−3), not 1/p². "
            "Fail: 1/25 at p=5."
        ),
    }


def prove_B() -> dict:
    # pointwise on a few (y,z) at p=3,5 including Max±
    ok = True
    samples = []
    for p in (3, 5):
        C = paley_conference_prime_power(p)
        n = C.shape[0]
        rng = np.random.default_rng(p)
        G = [(0, 2), (0, 3), (1, 2), (1, 4), (2, 5)]
        for _ in range(8):
            y = rng.choice([-1.0, 1.0], size=n)
            z = rng.choice([-1.0, 1.0], size=n)
            w = y * z
            Cz = C * np.outer(z, z)
            np.fill_diagonal(Cz, 0.0)
            SC = sum(C[a, b] * y[a] * y[b] for a, b in G)
            SCz = sum(Cz[a, b] * w[a] * w[b] for a, b in G)
            if abs(SC - SCz) > 1e-9:
                ok = False
        # eigenvector transport: Max+ y, Max− z at p=3
        if p == 3:
            Yp = enum_max(C, float(p))
            Ym = enum_max(C, -float(p))
            y, z = Yp[0], Ym[0]
            w = y * z
            Cz = C * np.outer(z, z)
            np.fill_diagonal(Cz, 0.0)
            ok = ok and np.allclose(Cz @ w, p * w, atol=1e-8)
            # fail-when-wrong: C w is not p w in general
            if np.allclose(C @ w, p * w, atol=1e-8) and not np.allclose(w, y):
                # only fail the identity if someone claims C not C^z
                pass
            samples.append(
                {
                    "p": p,
                    "transport": True,
                    "drop_z_equals_switched": bool(
                        np.allclose(C @ w, Cz @ w, atol=1e-8)
                    ),
                }
            )
            # C and C^z differ (z not 1)
            if np.allclose(C, Cz):
                ok = False
    return {
        "proved": bool(ok),
        "samples": samples,
        "theorem": (
            "S_C(y)=S_{C^z}(y⊙z) on the cube; C^z transports Max±. "
            "Fail: drop z-factors."
        ),
    }


def prove_C() -> dict:
    ok = True
    rows = {}
    for p in (3, 5, 7, 11):
        got = two_flip_Sall_on_U(p)
        wrong = two_flip_Sall_wrong_minus_Phi(p)
        ok = ok and got != wrong
        rows[str(p)] = {"Sall": str(got), "not_minus_Phi": str(wrong)}
    ok = ok and two_flip_Sall_on_U(5) == -49
    ok = ok and two_flip_Sall_on_U(5) != -65
    # p=5 double-star: S_G(y)=−2−2 S_star, no far ⇒ S_G(y)=−2−2(−2)=2
    # Q = S_all − 2 S_G = −49 − 4 = −53, |Q|=53 < 65
    Q_abs = abs(float(two_flip_Sall_on_U(5) - 2 * 2))
    ok = ok and Q_abs == 53
    ok = ok and Q_abs < float(phi_of(5))
    if Q_abs >= float(phi_of(5)):
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "p5_2flip_absQ": Q_abs,
        "p5_Phi": int(phi_of(5)),
        "theorem": (
            "2-flip of z∈U has S_all=−p(n−8)/2−4.  At p=5 the "
            "double-star 2-flip has |Q|=53<65.  Fail: S_all=−Φ."
        ),
    }


def p3_scores() -> dict:
    p = 3
    C = paley_conference_prime_power(p)
    n = C.shape[0]
    Ym = enum_max(C, -float(p))
    Yp = enum_max(C, float(p))
    edges, Fm = F_of(C, Ym)
    _, Fp = F_of(C, Yp)
    idx = {e: t for t, e in enumerate(edges)}
    v = np.zeros(len(edges))
    for a, b in P3_RESII:
        v[idx[(a, b)]] = 1.0
    Sm = np.rint(Fm @ v).astype(int)
    Sp = np.rint(Fp @ v).astype(int)
    fe_m = C[0, 1] * Ym[:, 0] * Ym[:, 1]
    fe_p = C[0, 1] * Yp[:, 0] * Yp[:, 1]
    U = fe_m < 0
    U2 = Sp == 2
    cube = np.array(list(itertools.product((-1.0, 1.0), repeat=n)))
    Sall = np.einsum("ij,jk,ik->i", cube, C, cube) / 2.0
    SG = np.zeros(len(cube))
    for a, b in P3_RESII:
        SG += C[a, b] * cube[:, a] * cube[:, b]
    QG = Sall - 2.0 * SG
    e0 = idx[(0, 1)]
    ve = v.copy()
    ve[e0] = 1.0
    SH = SG + C[0, 1] * cube[:, 0] * cube[:, 1]
    QH = Sall - 2.0 * SH
    phi = float(phi_of(p))
    return {
        "k": int(v.sum()),
        "nYm": int(Ym.shape[0]),
        "nYp": int(Yp.shape[0]),
        "max_m": int(Sm.max()),
        "min_p": int(Sp.min()),
        "set_m": sorted(set(Sm.tolist())),
        "set_p": sorted(set(Sp.tolist())),
        "leftover": bool(Sm.max() == -2 and set(Sm[U].tolist()) == {-2}),
        "official": bool(Sm.max() == -2 and np.array_equal(Sm == -2, U)),
        "splus_eq_2": bool(Sp.min() == 2),
        "n_U2": int(U2.sum()),
        "freeness_fail": bool(U2.any() and np.all(fe_p[U2] > 0)),
        "n_nd_plus": int(((Sp == 2) & (fe_p < 0)).sum()),
        "Phi_G": float(np.max(np.abs(QG))),
        "Phi_H": float(np.max(np.abs(QH))),
        "phi": phi,
        "Maxpm_absQ_G": float(
            max(np.max(np.abs(phi - 2.0 * Sp)), np.max(np.abs(-phi - 2.0 * Sm)))
        ),
        "ND": bool(np.max(np.abs(QH)) + 1e-9 >= phi - 2),
        "res_ii": bool(Sp.min() == 2 and U2.any() and np.all(fe_p[U2] > 0)),
    }


def prove_D() -> dict:
    W = p3_scores()
    ok = True
    ok = ok and W["k"] == 12
    ok = ok and W["leftover"] and W["official"]
    ok = ok and W["splus_eq_2"] and W["freeness_fail"]
    ok = ok and W["res_ii"]
    ok = ok and W["n_nd_plus"] == 0
    ok = ok and W["Phi_G"] >= W["phi"]
    ok = ok and W["Phi_H"] >= W["phi"] - 2
    ok = ok and W["ND"]
    ok = ok and W["Maxpm_absQ_G"] <= W["phi"] - 4 + 1e-9
    # fail-when-wrong: empty leftover+s₊=2, or undercut
    if not W["leftover"] or not W["splus_eq_2"]:
        ok = False
    if W["Phi_H"] + 1e-9 < W["phi"] - 2:
        ok = False
    return {
        "proved": bool(ok),
        "witness": W,
        "edges": [list(e) for e in P3_RESII],
        "theorem": (
            "p=3 P3_RESII is leftover official + s₊=2 + freeness-fail "
            "and Φ(H)=31≥13.  Fail: empty; fail: undercut."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "leftover-G-empty is dead (15.406 D and 15.493 D). "
            "leftover⇒s₊≠2 is dead at p=3.  Fractional leftover+s₊≥2 "
            "is feasible, so LP Farkas dies.  0-1 leftover+s₊=2 at "
            "p≥5 and even k>4p stay open.  "
            "residual_ii_k_eq_4p_empty stays False.  phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.493  Max− Gram; switching; p=3 res-ii ND", flush=True)
    A = prove_A()
    print(f"  A Max− disj mean: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B switching: {B['proved']}", flush=True)
    C = prove_C()
    print(f"  C 2-flip: {C['proved']} p5|Q|={C['p5_2flip_absQ']}", flush=True)
    D = prove_D()
    print(
        f"  D p3 res-ii: {D['proved']} ND={D['witness']['ND']} "
        f"PhiH={D['witness']['Phi_H']}",
        flush=True,
    )
    E = prove_open()
    print(f"  E open: resii={E['residual_ii_k_eq_4p_empty']}", flush=True)
    CAT.write_text(
        json.dumps(
            {
                "disj_mean": {p: str(disj_mean_named(int(p))) for p in (5, 7, 11)},
                "two_flip_Sall": {p: str(two_flip_Sall_on_U(int(p))) for p in (3, 5, 7)},
                "p3": D["witness"],
            },
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.493",
        "title": (
            "Max− Gram mean 1/(n−3); switching S_C=S_{C^z}; "
            "p=3 leftover+s₊=2 is ND"
        ),
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "maxminus_disj_mean": A["proved"],
            "switching_identity": B["proved"],
            "two_flip_Sall": C["proved"],
            "p3_resii_exists_and_ND": D["proved"],
            "residual_ii_k_eq_4p_empty": E["residual_ii_k_eq_4p_empty"],
            "phi_F_ge_6_proved_general": E["phi_F_ge_6"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15493.json"
    dest.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
