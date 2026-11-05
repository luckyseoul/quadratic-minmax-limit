"""Run Experiment 1: trend of L_n and constructive upper envelopes."""

from __future__ import annotations

import argparse
import csv
import json
import time
from pathlib import Path

from experiment.analysis import (
    analytic_lower_ln,
    analytic_upper_ln,
    exact_series,
    test_h1_predictions,
)
from experiment.constructive import conference_like_upper, local_search_upper
from experiment.exact import A399208, verify_against_oeis


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"


def parse_ns(s: str) -> list[int]:
    return [int(x.strip()) for x in s.split(",") if x.strip()]


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--max-exact",
        type=int,
        default=7,
        help="Exhaustively verify OEIS values up to this n (default 7; 8 is slow).",
    )
    p.add_argument(
        "--constructive-ns",
        type=str,
        default="16,24,32,48,64",
        help="Comma-separated n for local-search upper bounds.",
    )
    p.add_argument("--seed", type=int, default=0)
    p.add_argument(
        "--n-outer",
        type=int,
        default=120,
        help="Simulated-annealing outer steps per constructive n.",
    )
    p.add_argument(
        "--skip-constructive",
        action="store_true",
        help="Only run exact verification + OEIS trend analysis.",
    )
    args = p.parse_args()

    RESULTS.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    print("== Exact verification vs OEIS A399208 ==")
    verify_ns = list(range(1, max(1, args.max_exact) + 1))
    verified = verify_against_oeis(verify_ns)
    for n, row in verified.items():
        status = "OK" if row["match"] else "MISMATCH"
        print(f"  n={n:2d}  computed={row['computed']:3d}  oeis={row['oeis']:3d}  [{status}]")
    if not all(r["match"] for r in verified.values()):
        raise SystemExit("Exact verification failed; aborting.")

    print("\n== Exact normalized series L_n = F(n)/n^{3/2} (A399208) ==")
    exact_rows = exact_series(15)
    for r in exact_rows:
        print(
            f"  n={int(r['n']):2d}  F={int(r['F']):3d}  L_n={r['L_n']:.6f}  "
            f"lower={r['lower_bound_Ln']:.6f}  upper={r['upper_bound_Ln']:.6f}"
        )

    constructive_rows: list[dict] = []
    if not args.skip_constructive:
        print("\n== Constructive upper bounds ==")
        for n in parse_ns(args.constructive_ns):
            # Prefer Paley when available; always also run local search and take min.
            candidates = []
            paley = conference_like_upper(n, seed=args.seed)
            if paley is not None:
                candidates.append(paley)
                print(
                    f"  n={n:3d}  paley_upper={paley.upper_f:4d}  "
                    f"Lbar={paley.normalized:.6f}"
                )
            # Scale search effort mildly with n
            outer = args.n_outer if n <= 32 else max(40, args.n_outer // 2)
            restarts = 16 if n <= 32 else 8
            sa = local_search_upper(
                n, seed=args.seed + n, n_outer=outer, n_restarts=restarts
            )
            candidates.append(sa)
            print(
                f"  n={n:3d}  sa_upper={sa.upper_f:4d}  Lbar={sa.normalized:.6f}  "
                f"evals={sa.evals}"
            )
            best = min(candidates, key=lambda c: c.upper_f)
            constructive_rows.append(
                {
                    "n": best.n,
                    "upper_F": best.upper_f,
                    "normalized": best.normalized,
                    "method": best.method,
                    "seed": best.seed,
                    "evals": best.evals,
                    "all_candidates": [
                        {
                            "method": c.method,
                            "upper_F": c.upper_f,
                            "normalized": c.normalized,
                        }
                        for c in candidates
                    ],
                }
            )

    print("\n== Hypothesis H1 checks ==")
    h1 = test_h1_predictions(exact_rows, constructive_rows)
    for key in (
        "prediction_1_exact_monotone_from_n",
        "prediction_2_constructive_between_L15_and_half",
        "prediction_3_constructive_gap_not_widening",
    ):
        print(f"  {key}: holds={h1[key]['holds']}")
    print(f"  H1 overall: holds={h1['h1_overall_holds']}")

    payload = {
        "hypothesis": (
            "Along exact F(n) for n<=15, L_n is eventually monotone nondecreasing "
            "for n>=8, and constructive upper envelopes for n<=64 stay in a "
            "band between L_15 and 1/2 without clear gap widening."
        ),
        "open_problem": "existence of lim L_n (Ivanisvili, MathOverflow 413935)",
        "oeis": "A399208",
        "exact_verification": {str(k): v for k, v in verified.items()},
        "exact_series": exact_rows,
        "constructive_upper_bounds": constructive_rows,
        "h1": h1,
        "A399208": A399208,
        "elapsed_sec": round(time.time() - t0, 3),
        "limitations": [
            "Does not prove or disprove existence of the limit.",
            "Constructive values are approximate upper bounds (inner max over x is heuristic).",
            "n<=64 is far from a conclusive asymptotic regime.",
        ],
    }

    json_path = RESULTS / "experiment1.json"
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    csv_path = RESULTS / "experiment1_series.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "n",
                "F_exact",
                "F_upper",
                "L_n",
                "L_upper",
                "analytic_lower",
                "analytic_upper",
                "source",
            ]
        )
        for r in exact_rows:
            w.writerow(
                [
                    r["n"],
                    r["F"],
                    "",
                    f"{r['L_n']:.10f}",
                    "",
                    f"{r['lower_bound_Ln']:.10f}",
                    f"{r['upper_bound_Ln']:.10f}",
                    "A399208",
                ]
            )
        for c in constructive_rows:
            n = c["n"]
            w.writerow(
                [
                    n,
                    "",
                    c["upper_F"],
                    "",
                    f"{c['normalized']:.10f}",
                    f"{analytic_lower_ln(n):.10f}",
                    f"{analytic_upper_ln():.10f}",
                    c["method"],
                ]
            )

    print(f"\nWrote {json_path}")
    print(f"Wrote {csv_path}")
    print(f"Elapsed {payload['elapsed_sec']}s")


if __name__ == "__main__":
    main()
