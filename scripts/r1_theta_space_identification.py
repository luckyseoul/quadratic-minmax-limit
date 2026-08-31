#!/usr/bin/env python3
"""Ask PARI to identify the exact modular space of the ordinary dual theta."""
from __future__ import annotations

import argparse
import subprocess

import numpy as np

from r1_multigaussian_window import gp_matrix, paley_conference


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p", type=int, required=True)
    parser.add_argument("--remote")
    parser.add_argument("--coefficients", type=int, default=8)
    parser.add_argument("--stack", default="2G")
    parser.add_argument("--gp", default="gp")
    parser.add_argument("--threads", type=int)
    parser.add_argument(
        "--half-cusp-terms",
        type=int,
        help="also print this many exact coefficients at the cusp 1/2",
    )
    parser.add_argument(
        "--standard-cusps",
        action="store_true",
        help="also print the geometrically relevant initial blocks at 0, 1/4, and 1/p",
    )
    args = parser.parse_args()
    p = args.p
    C = paley_conference(p)
    A = C - p * np.eye(len(C), dtype=np.int64)
    thread_setup = (
        f"default(nbthreads,{args.threads});" if args.threads is not None else ""
    )
    half_cusp = ""
    if args.half_cusp_terms is not None:
        if args.half_cusp_terms < 1:
            raise ValueError("--half-cusp-terms must be positive")
        last = args.half_cusp_terms - 1
        half_cusp = f"""
params=0;
C=mfslashexpansion(T[1],T[2],[1,0;2,1],{last},1,&params);
print("CUSP_HALF_PARAMS=",params);
print("CUSP_HALF_COEFS=",C);
"""
    standard_cusps = ""
    if args.standard_cusps:
        zero_last = (p + 1) // 2 - 1
        quarter_last = 2 * (p + 1) - 1
        p_last = (p - 1) // 2 - 1
        standard_cusps = f"""
params0=0; C0=mfslashexpansion(T[1],T[2],[0,-1;1,0],{zero_last},1,&params0);
print("CUSP_ZERO_PARAMS=",params0); print("CUSP_ZERO_COEFS=",C0);
params4=0; C4=mfslashexpansion(T[1],T[2],[1,0;4,1],{quarter_last},1,&params4);
print("CUSP_QUARTER_PARAMS=",params4); print("CUSP_QUARTER_COEFS=",C4);
paramsp=0; Cp=mfslashexpansion(T[1],T[2],[1,0;{p},1],{p_last},1,&paramsp);
print("CUSP_P_PARAMS=",paramsp); print("CUSP_P_COEFS=",Cp);
"""
    program = f"""
default(parisize,{args.stack});
{thread_setup}
A={gp_matrix(A)};
B=matkerint(A);
G=B~*B;
Q={4*p}*G^-1;
if(denominator(Q)!=1,error("nonintegral even Gram"));
if(sum(i=1,matsize(Q)[1],Q[i,i]%2)!=0,error("Gram is not even"));
print("P={p} METHOD=ordinary");
gettime(); T=mffromqf(Q); print("MFFROMQF_MS=",gettime());
print("RANK=",matsize(Q)[1]);
print("DET=",matdet(Q));
print("MF_PARAMS=",mfparams(T[1]));
print("FORM_PARAMS=",mfparams(T[2]));
print("STURM=",mfsturm(T[1]));
print("COEFS=",mfcoefs(T[2],{args.coefficients}));
{half_cusp}
{standard_cusps}
quit;
"""
    command = [args.gp, "-fq", "-s", args.stack]
    if args.remote:
        command = [
            "ssh",
            "-F",
            "/dev/null",
            "-o",
            "StrictHostKeyChecking=no",
            "-o",
            "UserKnownHostsFile=/dev/null",
            args.remote,
            *command,
        ]
    proc = subprocess.run(
        command,
        input=program,
        text=True,
        check=False,
    )
    raise SystemExit(proc.returncode)


if __name__ == "__main__":
    main()
