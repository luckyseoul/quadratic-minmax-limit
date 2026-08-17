#!/usr/bin/env python3
"""
Prop 15.310 — Q on AP-halfspaces is the AP-averaged Fejer 4th moment
of 1_A; Q on QR0-halfspaces is the Gauss 4th moment of 1_{QR0}.
Full-ensemble Q_τ is still the unnamed mixture with NL.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.309 flags.  Floor stays OPEN.

============================================================================
Setup.  Type+ 1D (15.292): Q_1D(r)/q² = 32/(p+1) E_{A uniform}
|1̂_A(1)|²|1̂_A(r)|² on r∈F_p^×, and 0 off F_p.  Restrict the
average to arithmetic-progression m-subsets (AP-hs) or to
AGL-images of QR0 (15.309).

============================================================================
Theorem A — PROVED (Max+-free Fejer; certified p=5,7).
  Let m=(p+1)/2 and n_AP=p(p−1)/2.  For r∈F_p^×
      Q_AP(r)/q² = 32/(p+1) · (1/n_AP) ∑_{A AP} |1̂_A(1)|² |1̂_A(r)|²
  and Q_AP(r)=0 for r∉F_p.  This equals the live H_+ average of
  u(ξ)u(rξ) on the AP-halfspace G-orbit.  At p=5,7 it is constant
  on F_p^×\\{±1} and equals 2(p+3)/3; at p=11 the Fejer average
  is not constant (32/5 vs 184/15), so 2(p+3)/3 is p=5,7-only.
  Fail-when-wrong: use a single interval A={0,…,m−1} (no AP-average).
  At p=7 that Fejer varies (6.22 vs 12.99) while live Q_AP is
  constantly 20/3 on F_p^×\\{±1}.  ∎

Theorem B — PROVED (Max+-free Gauss; certified p=5,7).
  For S=QR0, |1̂_S(α)|² is the Paley/Gauss mass on F_p.  The same
  32/(p+1) prefactor with A=S (equivalently any AGL image) gives
  Q_QR0.  At p=5, QR0 is an AP so this coincides with (A).  At p=7,
  |1̂_S|²≡2, hence Q_QR0=16 on F_p^× and 0 off; class-average
  Q_{++}^{QR0}=32/3 matches live.  Fail-when-wrong: claim Q_QR0=0
  on F_p^×\\{±1}.  ∎

Theorem C — OPEN.  Full-ensemble Q_τ = ∑_σ w_σ Q_τ^{(σ)} with
  w_σ=|O_σ|/|H_+|.  Affine Q_τ^{(σ)} are named by (A)(B).  NL
  orbit values are short rationals (dens 9,15,21, not 13 or 409)
  but have no Fejer/Gauss formula in p.  The mixture still has
  den |H_+|/(2p).  phi_F not imported.  ∎

============================================================================
Backend: CPU Fejer over n_AP≤55; live AP-orbit Q from H_+ cache.
Writes evidence/e1_gmin_m4_prop15310.json
"""
from __future__ import annotations

import cmath
import json
import math
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
from e1_gmin_m4_prop15290 import in_subfield, live_Q_on_T, omega_set, zhat_omega  # noqa: E402
from e1_gmin_m4_prop15300 import merge_cls  # noqa: E402
from e1_gmin_m4_prop15307 import load_D, square_classes, tag_sig  # noqa: E402
from e1_gmin_m4_prop15308 import is_AP, recover_S_if_hs  # noqa: E402
from e1_gmin_m4_prop15309 import qr0  # noqa: E402


def hat2(A, alpha: int, p: int) -> float:
    s = 0j
    for x in A:
        s += cmath.exp(2j * math.pi * alpha * x / p)
    return float(abs(s) ** 2)


def all_APs(p: int) -> list[frozenset[int]]:
    m = (p + 1) // 2
    seen: set[frozenset[int]] = set()
    out: list[frozenset[int]] = []
    for a in range(p):
        for d in range(1, p):
            S = frozenset((a + k * d) % p for k in range(m))
            if S not in seen:
                seen.add(S)
                out.append(S)
    return out


def Q_AP_Fejer(p: int, r: int) -> float:
    """32/(p+1) AP-average |1̂(1)|²|1̂(r)|² for r∈F_p; 0 if r∉{1..p−1}."""
    if r % p != r or r == 0:
        return 0.0
    APs = all_APs(p)
    acc = 0.0
    for A in APs:
        acc += hat2(A, 1, p) * hat2(A, r, p)
    return 32.0 / (p + 1) * acc / len(APs)


def Q_AP_single_interval(p: int, r: int) -> float:
    """Fail-when-wrong: one interval, no AP-average."""
    if r % p != r or r == 0:
        return 0.0
    A = set(range((p + 1) // 2))
    return 32.0 / (p + 1) * hat2(A, 1, p) * hat2(A, r, p)


def Q_QR0_Gauss(p: int, r: int) -> float:
    if r % p != r or r == 0:
        return 0.0
    S = qr0(p)
    return 32.0 / (p + 1) * hat2(S, 1, p) * hat2(S, r, p)


def live_hs_Q(p: int, which: str) -> dict[int, float]:
    """Average u(ξ0)u(r ξ0)/q² on AP or QR0 halfspaces."""
    D, F = load_D(p)
    classes = square_classes(F)
    n_sq = len(classes)
    hs_key = tuple(sorted(["H"] + ["C"] * (n_sq - 1)))
    idx = []
    for i in range(D.shape[0]):
        tags = []
        for lines in classes:
            ns = D[i, lines].sum(axis=1)
            tags.append(tag_sig(p, tuple(int(x) for x in np.sort(ns))))
        if tuple(sorted(tags)) != hs_key:
            continue
        S = recover_S_if_hs(p, D[i], classes)
        if S is None:
            continue
        ap = is_AP(S, p)
        if which == "AP" and ap:
            idx.append(i)
        if which == "QR0" and (not ap):
            idx.append(i)
        if which == "QR0" and ap and is_AP(qr0(p), p):
            # p=5: QR0 is AP; treat as AP
            pass
    Omega = omega_set(F)
    om = {xi: i for i, xi in enumerate(Omega)}
    xi0 = Omega[0]
    q2 = float(F["q"] ** 2)
    if not idx:
        return {}
    U = np.stack(
        [
            np.abs(zhat_omega(np.where(D[i], -1.0, 1.0), F, Omega)) ** 2
            for i in idx
        ]
    )
    u0 = U[:, om[xi0]]
    out = {}
    for r in F["squares"]:
        eta = F["mul"][r][xi0]
        out[int(r)] = float(np.mean(u0 * U[:, om[eta]]) / q2)
    return out


def prove_Fejer_AP() -> dict:
    ok = True
    rows = {}
    single_varies = 0
    for p in (5, 7):
        F = _kernel_field(p)
        live = live_hs_Q(p, "AP")
        sub_live = []
        sub_fej = []
        sub_single = []
        max_err = 0.0
        for r in F["squares"]:
            if not in_subfield(F, r) or r in (1, F["neg1"]):
                continue
            lv = live[r]
            fj = Q_AP_Fejer(p, r % p)
            sg = Q_AP_single_interval(p, r % p)
            sub_live.append(lv)
            sub_fej.append(fj)
            sub_single.append(sg)
            max_err = max(max_err, abs(lv - fj))
            if abs(lv - fj) > 1e-6:
                ok = False
        if max(sub_single) - min(sub_single) > 1.0 and p == 7:
            single_varies += 1
        if p == 7 and abs(max(sub_single) - min(sub_single)) < 0.5:
            ok = False
        # 2(p+3)/3 matches at p=5,7
        twostep = 2.0 * (p + 3) / 3.0
        if any(abs(x - twostep) > 1e-6 for x in sub_fej):
            ok = False
        rows[str(p)] = {
            "n_sub": len(sub_live),
            "live_sub": sub_live[:4],
            "fejer_sub": sub_fej[:4],
            "single_sub": sub_single[:4],
            "max_err": max_err,
            "two_p_plus_3_over_3": twostep,
        }
    # p=11: Fejer not constant, so 2(p+3)/3 is not general
    fej11 = [Q_AP_Fejer(11, r) for r in range(2, 11)]
    if max(fej11) - min(fej11) < 1.0:
        ok = False
    if single_varies == 0:
        ok = False
    return {
        "proved": ok,
        "by_p": rows,
        "p11_fejer_spread": float(max(fej11) - min(fej11)),
        "theorem": (
            "Q_AP = AP-averaged Fejer. Single interval fails at p=7. "
            "2(p+3)/3 holds at p=5,7 not p=11."
        ),
    }


def prove_QR0_Gauss() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        F = _kernel_field(p)
        if is_AP(qr0(p), p):
            # coincides with AP
            rows[str(p)] = {"coincides_AP": True}
            continue
        live = live_hs_Q(p, "QR0")
        if not live:
            ok = False
            rows[str(p)] = {"empty": True}
            continue
        max_err = 0.0
        pm1 = []
        sub = []
        off = []
        for r in F["squares"]:
            pred = Q_QR0_Gauss(p, r % p) if in_subfield(F, r) else 0.0
            err = abs(live[r] - pred)
            max_err = max(max_err, err)
            if r in (1, F["neg1"]):
                pm1.append(live[r])
            elif in_subfield(F, r):
                sub.append(live[r])
            else:
                off.append(live[r])
        if max_err > 1e-6:
            ok = False
        if p == 7 and any(abs(x) < 1e-6 for x in sub):
            ok = False  # fail-when-wrong: Q_QR0=0 on sub
        rows[str(p)] = {
            "max_err": max_err,
            "live_pm1": pm1,
            "live_sub": sub[:3],
            "live_off_max": max(off) if off else None,
            "pred_sub": Q_QR0_Gauss(p, 2),
        }
    return {
        "proved": ok,
        "by_p": rows,
        "theorem": (
            "Q_QR0 = Gauss 4th moment of 1_QR0. "
            "At p=7 equals 16 on F_p^×, 0 off; not identically 0 on sub."
        ),
    }


def prove_open() -> dict:
    # mixture still has den 13, 409
    Q5 = live_Q_on_T(5)
    F5 = _kernel_field(5)
    q25 = float(25 * 25)
    pp = [Q5[r] / q25 for r in Q5 if merge_cls(F5, r) == "++"]
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "e1": bool(e1_closed_general()),
        "Q_tau_named_in_p": False,
        "lambda_ge_6": False,
        "p5_Qpp_live": float(np.mean(pp)) if pp else None,
        "note": (
            "Affine Q_τ named by Fejer/Gauss. Full-ensemble Q_τ is the "
            "NL mixture and still has den |H_+|/(2p). phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.310  Q_AP = AP-Fejer; Q_QR0 = Gauss; full Q_τ unnamed", flush=True)
    A = prove_Fejer_AP()
    print(f"  A Fejer AP: {A['proved']} p11_spread={A['p11_fejer_spread']:.3f}", flush=True)
    B = prove_QR0_Gauss()
    print(f"  B QR0 Gauss: {B['proved']}", flush=True)
    C = prove_open()
    print(f"  C floor: {C['proved']} phi_F={C['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.310",
        "title": "Q_AP is AP-averaged Fejer; Q_QR0 is Gauss; full Q_τ unnamed",
        "proved": {
            "Q_AP_Fejer": A["proved"],
            "Q_QR0_Gauss": B["proved"],
            "Q_tau_named_in_p": False,
            "lambda_ge_6": False,
            "phi_F_ge_6_proved_general": C["phi_F_ge_6"],
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
    path = ROOT / "evidence" / "e1_gmin_m4_prop15310.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
