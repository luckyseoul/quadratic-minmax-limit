#!/usr/bin/env python3
"""Compute the multiplicity-weighted trace of every R1 shell operator.

For an orthonormal basis ``W_a`` of

    Z = {W : PWP=W, diag(W)=0},

this constructs the harmonic polynomial ``sum_a H_{W_a}`` without building
that basis.  Its coefficient at a dual shell is the trace of the harmonic
shell operator.  Adding the universal radial correction gives the trace of
the positive raw-quartic operator, which couples all PSL channels.

The default calculation uses a symmetry reduction.  Coordinate transitivity
makes the theta series of all coordinate-zonal harmonic quartics equal, and
the trace polynomial is their orbit sum.  PARI therefore sees one compact
zonal quartic instead of an expanded sum of ``p^2+1`` fourth powers.  The
older direct polynomial remains available as an independent cross-check.
"""
from __future__ import annotations

import argparse
import subprocess

import numpy as np

from r1_multigaussian_window import gp_matrix, paley_conference


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p", type=int, required=True)
    parser.add_argument("--remote")
    parser.add_argument("--coefficients", type=int, default=40)
    parser.add_argument("--stack", default="12G")
    parser.add_argument(
        "--gp",
        default="gp",
        help="GP executable on the local or --remote host",
    )
    parser.add_argument(
        "--threads",
        type=int,
        help="PARI worker-thread count (defaults to the host build setting)",
    )
    parser.add_argument(
        "--method",
        choices=("zonal", "direct"),
        default="zonal",
        help="use the transitive one-coordinate reduction or the full orbit sum",
    )
    args = parser.parse_args()

    p = args.p
    C = paley_conference(p)
    n = p * p + 1
    d = n // 2
    zdim = n * (n - 6) // 8
    variables = ",".join(f"x{i}" for i in range(1, d + 1))
    if args.method == "zonal":
        # For ell_i(u)=u_i on range(P), ||ell_i||^2=P_ii=1/2.  Its harmonic
        # quartic is ell_i^4-3*r*ell_i^2/(d+4)+3*r^2/[4(d+2)(d+4)].
        # Summing over the n transitive coordinates and simplifying the
        # radial coefficient gives H_tr=-4p^2/(p^2-1)*sum_i Z_i.
        polynomial = f"""
ell=u[1];
Hone=ell^4-3*qr*ell^2/(d+4)+3*qr^2/(4*(d+2)*(d+4));
tfac=-4*p^2*n/(p^2-1);
H=Hone;
"""
        coefficient_expression = "tfac*mfcoefs(T[2],%d)" % args.coefficients
    else:
        polynomial = """
sum4=sum(i=1,n,u[i]^4);
proj4=qr^2-4*p^2*sum4/(p^2-1)+2*qr^2/(p^2-1);
H=proj4-4*zdim*qr^2/(d*(d+4))+2*zdim*qr^2/((d+2)*(d+4));
"""
        coefficient_expression = "mfcoefs(T[2],%d)" % args.coefficients

    thread_setup = (
        f"default(nbthreads,{args.threads});" if args.threads is not None else ""
    )
    program = f"""
default(parisize,{args.stack});
{thread_setup}
C={gp_matrix(C)};
n={n}; p={p}; d={d}; zdim={zdim};
A=C-p*matid(n); B=matkerint(A); G=B~*B; D=B*G^-1;
R=G^-1; Q=4*p*R;
if(denominator(Q)!=1,error("nonintegral Q"));
xx=[{variables}]; u=D*xx~; qr=xx*R*xx~;
{polynomial}
print("P=",p," RANK=",d," ZDIM=",zdim," METHOD={args.method}");
gettime(); T=mffromqf(Q,H); print("MFFROMQF_MS=",gettime());
print("MF_PARAMS=",mfparams(T[1]));
print("FORM_PARAMS=",mfparams(T[2]));
print("TRACE_HARMONIC_COEFS=",{coefficient_expression});
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
    result = subprocess.run(command, input=program, text=True, check=False)
    raise SystemExit(result.returncode)


if __name__ == "__main__":
    main()
