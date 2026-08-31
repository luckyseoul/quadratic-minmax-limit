#!/usr/bin/env python3
"""Compute exact ordinary Paley-dual theta counts without ``mffromqf``.

PARI's ``mffromqf`` always enumerates through the full Sturm bound of the
ambient modular-form space.  The p=11 R1 cone only needs a shorter prefix.
For an even Gram matrix Q, ``qfrep(Q, bound, 1)`` returns half the number of
vectors with Q-norm 2, 4, ..., 2*bound, so doubling and prepending the zero
vector gives the required exact theta coefficients directly.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess

import numpy as np

from r1_multigaussian_window import gp_matrix, paley_conference


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p", type=int, required=True)
    parser.add_argument("--bound", type=int, required=True)
    parser.add_argument("--stack", default="8G")
    parser.add_argument("--remote")
    parser.add_argument("--gp", default="gp")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    p = args.p
    conference = paley_conference(p)
    difference = conference - p * np.eye(len(conference), dtype=np.int64)
    program = f"""
default(parisize,{args.stack});
A={gp_matrix(difference)};
B=matkerint(A);
G=B~*B;
Q={4*p}*G^-1;
if(denominator(Q)!=1,error("nonintegral even Gram"));
if(!qfiseven(Q),error("Gram is not even"));
print("P={p} RANK=",matsize(Q)[1]," BOUND={args.bound}");
gettime(); V=qfrep(Q,{args.bound},1); print("QFREP_MS=",gettime());
print("THETA_COEFS=",concat([1],2*V));
quit;
"""
    command = [args.gp, "-fq", "-s", args.stack]
    if args.remote:
        command = [
            "ssh",
            args.remote,
            *command,
        ]
    result = subprocess.run(
        command,
        input=program,
        text=True,
        capture_output=True,
        check=False,
    )
    if args.output is not None:
        args.output.write_text(result.stdout)
    print(result.stdout, end="")
    if result.stderr:
        print(result.stderr, end="")
    raise SystemExit(result.returncode)


if __name__ == "__main__":
    main()
