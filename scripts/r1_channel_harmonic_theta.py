#!/usr/bin/env python3
"""Resolve Paley harmonic theta forms by square-circle PSL channel.

At p=5 the three eigenspaces of the square-circle tensor operator are
exactly the three irreducible constituents of the admissible tensor space.
This script constructs one rational tensor in each eigenspace and asks PARI
for its normalized harmonic theta series.  It is a structural reconnaissance
tool; larger p split some of these three spaces further.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15634 import explicit_square_circles
from r1_multigaussian_window import gp_matrix, paley_conference


def gp_vector(values: np.ndarray) -> str:
    return "[" + ",".join(str(int(value)) for value in values) + "]~"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p", type=int, default=5)
    parser.add_argument("--coefficients", type=int, default=30)
    parser.add_argument("--remote")
    parser.add_argument("--stack", default="4G")
    args = parser.parse_args()

    p = args.p
    C = paley_conference(p)
    _blocks, words = explicit_square_circles(p)
    n = p * p + 1
    d = n // 2
    variables = ",".join(f"x{i}" for i in range(1, d + 1))
    seed = np.zeros(n, dtype=np.int64)
    seed[:5] = (1, 2, -1, 3, -2)

    program = f"""
default(parisize,{args.stack});
C={gp_matrix(C)};
Words={gp_matrix(words)};
z={gp_vector(seed)};
n={n}; p={p}; d={d}; b=matsize(Words)[1];
Pm=(p*matid(n)+C)/(2*p);
projZ(v)={{my(vv,aa,W); vv=vector(n,i,v[i]^2)~;
  aa=(4*p^2*vv-2*sum(i=1,n,vv[i])*vector(n,i,1)~)/(p^2-1);
  W=v*v~; for(i=1,n,W-=aa[i]*Pm[,i]*Pm[,i]~); W;}};
Blist=vector(b,s,projZ(Words[s,]~));
op(A)={{my(R=matrix(n,n)); for(s=1,b,R+=trace(Blist[s]*A)*Blist[s]); R;}};
ll=p^3*(p-1); lh=p^3*(p+1);
split_low(A)=(lh*op(A)-op(op(A)))/(ll*(lh-ll));
split_high(A)=(op(op(A))-ll*op(A))/(lh*(lh-ll));
B0=Blist[1];
Wl=split_low(B0); Wh=split_high(B0);
A0=projZ(Pm*z); Wk=A0-split_low(A0)-split_high(A0);
Ws=[Wk,Wl,Wh]; labels=["circle-kernel","circle-low","circle-high"];
for(j=1,3,{{
  W=Ws[j];
  if(W!=Pm*W*Pm || vector(n,i,W[i,i])!=vector(n),error("bad W",j));
  ev=if(j==1,0,if(j==2,ll,lh));
  if(op(W)!=ev*W,error("bad channel",j));
  F=trace(W*W); if(F==0,error("zero channel",j));
  A=C-p*matid(n); B=matkerint(A); G=B~*B; D=B*G^-1; Q=4*p*G^-1;
  M=D~*W*D; N=D~*W*W*D; R=G^-1;
  xx=[{variables}]; qm=xx*M*xx~; qn=xx*N*xx~; qr=xx*R*xx~;
  H=qm^2-4*qr*qn/(d+4)+2*F*qr^2/((d+2)*(d+4));
  T=mffromqf(Q,H); coeff=mfcoefs(T[2],{args.coefficients})/F;
  print("CHANNEL=",labels[j]," F=",F," PARAMS=",mfparams(T[2]));
  print("NORMALIZED_COEFS=",coeff);
  sh=mfshimura(T[1],T[2]/F,p);
  print("SHIMURA_PARAMS=",mfparams(sh[2]));
  print("SHIMURA_COEFS=",mfcoefs(sh[2],20));
  print("SHIMURA_NEW=",mftonew(sh[1],sh[2]));
}});
quit;
"""
    command = ["gp", "-fq", "-s", args.stack]
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
