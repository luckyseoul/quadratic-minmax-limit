#!/usr/bin/env python3
"""
Prop 15.353 — Ensemble Q_{++} is a mixture of 2 (p=5) or 4 (p=7)
per-y values.  The hs+triangle 2-type formula is exact at p=5 and
fails at p=7 (1796/409 ≠ 1544/409).  Q_τ still unnamed in p.
phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.352 flags.

============================================================================
Setup.  Q_y^{++} = Ω-average of u(ξ)u(rξ)/q² on ++.  Ensemble is
(1/|H_+|) ∑_y Q_y^{++}.  M=∑ n_ℓ² on a square class; 15.306 C
converts E[M²] to Q_{++}.

============================================================================
Theorem A — PROVED (cache fold; certified p=5).
  Q_y^{++} ∈ {16/9, 64/15} with weights n_1d=30 and p²(p−1)=100.
  The mix is 48/13.  Fail: claim a single value, or 2-type E[M²]
  with those weights fails to recover 48/13.  ∎

Theorem B — PROVED (cache fold; certified p=7).
  Q_y^{++} ∈ {40/9 (2730), 8/3 (1764), 32/9 (1176), 32/3 (56)}.
  The mix is 1544/409.  32/3 is the QR0 halfspace (15.310).
  Fail: claim the hs+triangle 2-type value 1796/409.  ∎

Theorem C — OPEN.  No Gauss/Jacobi formula in p hits both 48/13
  and 1544/409.  phi_F not imported.

============================================================================
Backend: ProcessPool zhat over H_+.  Writes
evidence/e1_gmin_m4_prop15353.json
"""
from __future__ import annotations

import json
import sys
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15290 import omega_set, zhat_omega  # noqa: E402
from e1_gmin_m4_prop15292 import n_1d  # noqa: E402
from e1_gmin_m4_prop15300 import merge_cls  # noqa: E402
from e1_gmin_m4_prop15307 import M_named, load_D  # noqa: E402
from e1_gmin_m4_prop15330 import LIVE_QPP  # noqa: E402


HPLUS = {5: 130, 7: 5726}
W = 86


def _qpp_shard(spec):
    p, lo, hi = spec
    D, F = load_D(p)
    Omega = omega_set(F)
    om = {xi: i for i, xi in enumerate(Omega)}
    q = F["q"]
    q2 = float(q * q)
    rpp = [r for r in F["squares"] if merge_cls(F, r) == "++"]
    out = []
    for i in range(lo, hi):
        z = np.where(D[i], -1.0, 1.0)
        u = np.abs(zhat_omega(z, F, Omega)) ** 2
        sm = 0.0
        for r in rpp:
            for xi in Omega:
                sm += u[om[xi]] * u[om[F["mul"][r][xi]]]
        qpp = sm / (len(Omega) * len(rpp) * q2)
        out.append(str(Fraction(qpp).limit_denominator(p ** 4)))
    return out


def qpp_hist(p: int) -> dict[str, int]:
    D, _ = load_D(p)
    H = int(D.shape[0])
    sw = max(1, (H + W - 1) // W)
    jobs = [(p, lo, min(H, lo + sw)) for lo in range(0, H, sw)]
    vals = []
    with ProcessPoolExecutor(max_workers=min(W, len(jobs))) as ex:
        for part in ex.map(_qpp_shard, jobs):
            vals.extend(part)
    return dict(Counter(vals))


def mix_of(hist: dict[str, int]) -> Fraction:
    n = sum(hist.values())
    return sum(Fraction(k) * v for k, v in hist.items()) / n


def twotype_Q(p: int) -> Fraction:
    n1 = n_1d(p)
    N = HPLUS[p]
    n_sq = (p + 1) // 2
    named = M_named(p)
    MH, MC, MR = named["H"], named["C"], named["R"]
    hs_m2 = (MH ** 2 + (n_sq - 1) * MC ** 2) / n_sq
    tri_m2 = (3 * MR ** 2 + (n_sq - 3) * MC ** 2) / n_sq
    EM2 = (n1 * Fraction(hs_m2) + (N - n1) * Fraction(tri_m2)) / N
    k = p * (p - 1) // 2
    q = p * p
    EM = Fraction(p * (q - 1), 4)
    ES2 = p * p * EM2 - 2 * p * k * k * EM + k ** 4
    return 16 * (ES2 / ((p - 1) * q * q) - 1) / (p - 3)


def prove_A() -> dict:
    hist = qpp_hist(5)
    mix = mix_of(hist)
    ok = (
        set(hist) == {"16/9", "64/15"}
        and hist["16/9"] == n_1d(5)
        and hist["64/15"] == 5 * 5 * 4
        and mix == LIVE_QPP[5]
        and twotype_Q(5) == LIVE_QPP[5]
    )
    if len(hist) == 1:
        ok = False
    return {
        "proved": bool(ok),
        "hist": hist,
        "mix": str(mix),
        "twotype": str(twotype_Q(5)),
        "theorem": "p=5: Q_y^{++}∈{16/9,64/15} mix to 48/13.",
    }


def prove_B() -> dict:
    hist = qpp_hist(7)
    mix = mix_of(hist)
    two = twotype_Q(7)
    ok = (
        set(hist) == {"40/9", "8/3", "32/9", "32/3"}
        and hist["40/9"] == 2730
        and hist["8/3"] == 1764
        and hist["32/9"] == 1176
        and hist["32/3"] == 56
        and mix == LIVE_QPP[7]
        and two == Fraction(1796, 409)
        and two != LIVE_QPP[7]
    )
    if two == LIVE_QPP[7]:
        ok = False
    return {
        "proved": bool(ok),
        "hist": hist,
        "mix": str(mix),
        "twotype": str(two),
        "theorem": (
            "p=7: four Q_y^{++} values mix to 1544/409.  "
            "hs+triangle 2-type is 1796/409, not the ensemble."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "Q_tau_named_in_p": False,
        "note": (
            "Mixture of named per-y values recovers the ensemble at p=5,7 "
            "but is not a Gauss/Jacobi formula in p.  phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.353  ensemble Q++ is a 2- or 4-level mix", flush=True)
    A = prove_A()
    print(f"  A p5: {A['proved']} {A['hist']} mix={A['mix']}", flush=True)
    B = prove_B()
    print(f"  B p7: {B['proved']} {B['hist']} mix={B['mix']} 2type={B['twotype']}", flush=True)
    C = prove_open()
    print(f"  C open: phi_F={C['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.353",
        "title": "Ensemble Q++ is a 2-/4-level mix; 2-type dies at p=7",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "p5_two_level_mix": A["proved"],
            "p7_four_level_mix_not_2type": B["proved"],
            "Q_tau_named_in_p": False,
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15353.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
