#!/usr/bin/env python3
"""Count scaled dual vectors through a bound without storing them."""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from minmax_quadratic import paley_conference_prime_power  # noqa: E402


def paley_conference(p: int) -> np.ndarray:
    return np.rint(paley_conference_prime_power(p)).astype(np.int64)


def gp_matrix(matrix: np.ndarray) -> str:
    rows = [",".join(str(int(value)) for value in row) for row in matrix]
    return "[" + ";".join(rows) + "]"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p", type=int, required=True)
    parser.add_argument("--scaled-bound", type=int, required=True)
    args = parser.parse_args()
    p = args.p
    C = paley_conference(p)
    A = C - p * np.eye(len(C), dtype=np.int64)
    program = f"""
default(parisize,2G);
A={gp_matrix(A)};
B=matkerint(A);
G=B~*B;
Q={2*p}*G^-1;
if(denominator(Q)!=1,error("nonintegral scaled dual Gram"));
gettime();
DM=qfminim(Q,{args.scaled_bound},0);
print("P=",{p});
print("BOUND=",{args.scaled_bound});
print("SIGNED_COUNT=",DM[1]);
print("MAXNORM=",DM[2]);
print("ELAPSED_MS=",gettime());
quit;
"""
    proc = subprocess.run(
        ["gp", "-fq", "-s", "2G"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
    )
    print(proc.stdout, end="")
    if proc.stderr.strip():
        print(proc.stderr, end="")


if __name__ == "__main__":
    main()
