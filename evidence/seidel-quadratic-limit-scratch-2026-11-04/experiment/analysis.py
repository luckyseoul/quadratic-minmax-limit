"""Hypothesis checks and reporting for Experiment 1."""

from __future__ import annotations

import math
from typing import Any

from experiment.exact import A399208


def analytic_lower_ln(n: int) -> float:
    """Normalized form of OEIS lower bound n*sqrt(n-1)/pi."""
    if n <= 1:
        return 0.0
    return math.sqrt(1.0 - 1.0 / n) / math.pi


def analytic_upper_ln() -> float:
    return 0.5


def ln_from_f(n: int, f: int) -> float:
    return f / (n ** 1.5)


def exact_series(max_n: int = 15) -> list[dict[str, float | int]]:
    rows: list[dict[str, float | int]] = []
    for n in range(1, max_n + 1):
        f = A399208[n - 1]
        rows.append(
            {
                "n": n,
                "F": f,
                "L_n": ln_from_f(n, f),
                "lower_bound_Ln": analytic_lower_ln(n),
                "upper_bound_Ln": analytic_upper_ln(),
            }
        )
    return rows


def test_h1_predictions(
    exact_rows: list[dict[str, float | int]],
    constructive: list[dict[str, Any]],
    monotone_from: int = 8,
) -> dict[str, Any]:
    """
    H1 predictions:
      (1) exact L_n nondecreasing for n >= monotone_from
      (2) constructive upper L-bar_n in (L_15, 0.5] for reported n
      (3) constructive gap to analytic lower does not widen vs n=16 baseline
    """
    by_n = {int(r["n"]): r for r in exact_rows}
    diffs = []
    monotone_ok = True
    ns = sorted(n for n in by_n if n >= monotone_from)
    for prev, cur in zip(ns, ns[1:]):
        d = float(by_n[cur]["L_n"]) - float(by_n[prev]["L_n"])
        diffs.append({"from": prev, "to": cur, "delta_L": d})
        if d < -1e-15:
            monotone_ok = False

    l15 = float(by_n[15]["L_n"]) if 15 in by_n else None
    pred2 = []
    pred2_ok = True
    for c in constructive:
        n = int(c["n"])
        ub = float(c["normalized"])
        ok = (l15 is None or ub + 1e-12 >= l15) and ub <= 0.5 + 1e-12
        pred2.append({"n": n, "upper_Ln": ub, "ok": ok})
        if not ok:
            pred2_ok = False

    gaps = []
    pred3_ok = True
    if constructive:
        base = None
        for c in sorted(constructive, key=lambda r: int(r["n"])):
            n = int(c["n"])
            gap = float(c["normalized"]) - analytic_lower_ln(n)
            gaps.append({"n": n, "gap": gap})
            if base is None:
                base = gap
            elif gap > base + 0.05:  # allow small noise; flag clear widening
                pred3_ok = False

    return {
        "prediction_1_exact_monotone_from_n": {
            "from_n": monotone_from,
            "holds": monotone_ok,
            "deltas": diffs,
        },
        "prediction_2_constructive_between_L15_and_half": {
            "holds": pred2_ok,
            "L15": l15,
            "rows": pred2,
        },
        "prediction_3_constructive_gap_not_widening": {
            "holds": pred3_ok,
            "gaps": gaps,
        },
        "h1_overall_holds": bool(monotone_ok and pred2_ok and pred3_ok),
    }
