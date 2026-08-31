#!/usr/bin/env python3
"""Compute one exact Paley-dual harmonic theta series with PARI/GP.

This is an R1 reconnaissance tool.  It uses the projected tensor of one
square circle, so the polynomial is in the actual admissible R1 tensor
space rather than in the full ambient harmonic-polynomial space.
"""
from __future__ import annotations

import argparse
import subprocess

import numpy as np

from r1_multigaussian_window import gp_matrix, paley_conference

from e1_gmin_m4_prop15634 import explicit_square_circles


def gp_vector(values: np.ndarray) -> str:
    return "[" + ",".join(str(int(x)) for x in values) + "]~"


def square_circle_word(p: int) -> np.ndarray:
    _blocks, words = explicit_square_circles(p)
    return words[0]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p", type=int, required=True)
    parser.add_argument("--remote")
    parser.add_argument("--coefficients", type=int, default=40)
    args = parser.parse_args()
    p = args.p
    C = paley_conference(p)
    w = square_circle_word(p)
    n = len(C)
    d = n // 2
    variables = ",".join(f"x{i}" for i in range(1, d + 1))
    program = f"""
default(parisize,20G);
C={gp_matrix(C)};
w={gp_vector(w)};
n={n}; p={p}; d={d};
Pm=(p*matid(n)+C)/(2*p);
vv=vector(n,i,w[i]^2)~;
aa=(4*p^2*vv-2*sum(i=1,n,vv[i])*vector(n,i,1)~)/(p^2-1);
W=w*w~;
for(i=1,n,W-=aa[i]*Pm[,i]*Pm[,i]~);
if(W!=Pm*W*Pm,error("projection failure"));
if(vector(n,i,W[i,i])!=vector(n),error("diagonal failure"));
A=C-p*matid(n);
B=matkerint(A);
G=B~*B;
D=B*G^-1;
Q=4*p*G^-1;
if(denominator(Q)!=1,error("nonintegral Q"));
M=D~*W*D;
N=D~*W*W*D;
R=G^-1;
F=trace(W*W);
xx=[{variables}];
qm=xx*M*xx~;
qn=xx*N*xx~;
qr=xx*R*xx~;
H=qm^2-4*qr*qn/(d+4)+2*F*qr^2/((d+2)*(d+4));
print("P=",p," RANK=",d," FROB2=",F);
gettime(); T=mffromqf(Q,H); print("MFFROMQF_MS=",gettime());
print("MF_PARAMS=",mfparams(T[1]));
print("FORM_PARAMS=",mfparams(T[2]));
print("COEFS=",mfcoefs(T[2],{args.coefficients}));
quit;
"""
    command = ["gp", "-fq", "-s", "20G"]
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
    result = subprocess.run(command, input=program, text=True, check=False)
    raise SystemExit(result.returncode)


if __name__ == "__main__":
    main()
