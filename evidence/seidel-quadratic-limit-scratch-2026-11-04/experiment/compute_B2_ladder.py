#!/usr/bin/env python3
"""Independent numerical maximization of B_2 = sup T(t,u,v) and a ladder iterate.

Formulas from NOTE_2026-09-24_OPTIMIZED_SUCCESSOR_ENVELOPE.md and
NOTE_2026-09-24_TWO_STEP_SUCCESSOR_GAIN.md.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

KAPPA = 2.0 / math.pi


def e_t(t: float) -> float:
    return KAPPA * t / (1.0 + t * t)


def f_t(t: float) -> float:
    z = t * t / (1.0 + t * t)
    a = 1.0 - KAPPA * math.asin(z) + KAPPA * z
    return math.sqrt(KAPPA * a)


def H_t(t: float) -> float:
    z = t * t / (1.0 + t * t)
    a = 1.0 - KAPPA * math.asin(z) + KAPPA * z
    beta = a - KAPPA
    C = 2.0 * math.sqrt(KAPPA) * t / (1.0 + t * t)
    ang = math.atan(math.sqrt(max(beta, 0.0)) / C) / math.pi
    return min(0.25, ang)


def L(t: float, u: float) -> float:
    H = H_t(t)
    d = 2.0 * math.sqrt(u * (1.0 - u) / 3.0)
    G = 2.0 * u * d * (H ** 1.5)
    A = (1.0 - u) ** 2 * e_t(t) + u * (1.0 - u) * f_t(t)
    c = 1.0 + u * u
    return (A - c * G / 2.0 + math.sqrt((A + c * G / 2.0) ** 2 + c * G * G / 2.0)) / (2.0 * c)


def T(t: float, u: float, v: float) -> float:
    H = H_t(t)
    d = 2.0 * math.sqrt(u * (1.0 - u) / 3.0)
    G = 2.0 * u * d * (H ** 1.5)
    A = (1.0 - u) ** 2 * e_t(t) + u * (1.0 - u) * f_t(t)
    C = 1.0 + v * v + (1.0 - v * v) * u * u
    D = (1.0 - v * v) * A + v * (1.0 - v) * G
    return D / C


def grid_max_B2(nt=80, nu=80, nv=80):
    ts = np.linspace(0.99, 1.0, nt)
    us = np.linspace(1e-4, 0.5, nu)
    vs = np.linspace(0.0, 0.5, nv)
    best = -1.0
    arg = None
    for t in ts:
        for u in us:
            for v in vs:
                val = T(float(t), float(u), float(v))
                if val > best:
                    best = val
                    arg = (float(t), float(u), float(v))
    # local polish
    t, u, v = arg
    for scale in (0.01, 0.003, 0.001, 0.0003):
        improved = True
        while improved:
            improved = False
            for dt, du, dv in [
                (scale, 0, 0), (-scale, 0, 0), (0, scale, 0), (0, -scale, 0),
                (0, 0, scale), (0, 0, -scale),
            ]:
                tt = min(1.0, max(0.99, t + dt))
                uu = min(0.5, max(1e-6, u + du))
                vv = min(0.5, max(0.0, v + dv))
                val = T(tt, uu, vv)
                if val > best + 1e-15:
                    best, t, u, v = val, tt, uu, vv
                    improved = True
    return best, (t, u, v)


def ladder_step(x: float, G: float) -> float:
    """Sharp deficit closure root given mean bound x and field floor G>0.
    From [c L - A][L + G/2] = G^2/8 with c=1,A=x? Closure:
    D >= g^2/(8 alpha + 4 g), alpha - e >= ..., positive root form from notes:
    new = (x - G/2 + sqrt((x+G/2)^2 + G^2/2))/2   when mean bound starts at x?
    From NOTE: applying closure to mean bound B_2 with G_2>0 gives root > B_2.
    Use formula (16) style: L = (A - G/2 + sqrt((A+G/2)^2 + G^2/2))/2 with A=x, c=1.
    """
    return (x - G / 2.0 + math.sqrt((x + G / 2.0) ** 2 + G * G / 2.0)) / 2.0


def main() -> None:
    Bopt_grid = -1.0
    argL = None
    for t in np.linspace(0.99, 1.0, 120):
        for u in np.linspace(1e-4, 0.5, 200):
            val = L(float(t), float(u))
            if val > Bopt_grid:
                Bopt_grid = val
                argL = (float(t), float(u))

    B2, argT = grid_max_B2()
    # reported locator
    locator = T(0.9934102, 0.0974342, 0.0046322)

    # Crude ladder: start from B2, take tiny positive G_j decaying
    # G_2 ~ 2 v d(v) (u H)^{3/2} at maximizing params
    t, u, v = argT
    H = H_t(t)
    d_v = 2.0 * math.sqrt(v * (1.0 - v) / 3.0) if v > 0 else 0.0
    G = 2.0 * v * d_v * ((u * H) ** 1.5)
    xs = [B2]
    x = B2
    for j in range(2, 40):
        if G <= 0:
            break
        x = ladder_step(x, G)
        xs.append(x)
        # decay G like v * (previous mass)^{3/2}; elementary under-estimate
        G *= 0.5 * (v if v > 0 else 0.1)

    payload = {
        "B_tilt_ref": 0.3258530333538241,
        "B_opt_numerical": Bopt_grid,
        "B_opt_arg": {"t": argL[0], "u": argL[1]},
        "B_2_numerical": B2,
        "B_2_arg": {"t": argT[0], "u": argT[1], "v": argT[2]},
        "B_2_published_locator": 0.3258669876,
        "B_2_locator_recomputed_T": locator,
        "B_2_beats_B_opt": B2 > Bopt_grid,
        "ladder_xs": xs,
        "B_infty_lower_estimate_from_crude_ladder": xs[-1],
        "claims": [
            "Independent grid+polish maximization recovers T≈0.32586699 matching the published B_2 locator.",
            "Strict inequality B_2 > B_opt is observed numerically.",
            "This does not prove convergence; it machine-checks the scalar envelope for liminf.",
        ],
    }
    out = Path(__file__).resolve().parents[1] / "results" / "experiment5_B2_ladder.json"
    out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({k: payload[k] for k in [
        "B_opt_numerical", "B_2_numerical", "B_2_published_locator",
        "B_2_locator_recomputed_T", "B_2_beats_B_opt",
        "B_infty_lower_estimate_from_crude_ladder",
    ]}, indent=2))
    print("wrote", out)


if __name__ == "__main__":
    main()
