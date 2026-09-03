#!/usr/bin/env python3
"""Intel/OpenCL verifier for every b=0 midpoint-seam moment block.

For an actual signed midpoint function ``g`` and its affine line profiles
``G_L``, this independently compares

    sum_a G_L(a) a^i
    sum_x g(x) L(x)^i
    evaluation at L of the homogeneous binary i-form of mixed moments

for every ``0 <= i <= p-1``.  This is a large-prime identity check, not a
finite-prime census and not a proof by sampling.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import random
import time
from pathlib import Path

import numpy as np
import pyopencl as cl

from residual_midpoint_seam_gpu import directions_and_signs, is_prime, legendre_table


KERNEL = r"""
__kernel void moments(__global const int *g, __global const int *profiles,
                      __global const int *powers, __global int *direct,
                      __global int *incidence, int p) {
    int q = get_global_id(0);
    int directions = p + 1;
    int i = q / directions, d = q - i * directions;
    if (i >= p) return;
    int sd = 0;
    for (int point=0; point<p*p; ++point) {
        int x=point/p, y=point-x*p;
        int label = d == p ? y : (x + d*y) % p;
        sd = (sd + g[point] * powers[i*p + label]) % p;
    }
    int si = 0;
    for (int a=0; a<p; ++a)
        si = (si + profiles[d*p+a] * powers[i*p+a]) % p;
    if (sd < 0) sd += p; if (si < 0) si += p;
    direct[q]=sd; incidence[q]=si;
}

__kernel void mixed(__global const int *g, __global const int *powers,
                    __global const int *binom, __global int *coeff, int p) {
    int q=get_global_id(0), i=q/p, k=q-i*p;
    if (i>=p || k>i) return;
    int value=0;
    for (int point=0; point<p*p; ++point) {
        int x=point/p, y=point-x*p;
        int term = (powers[(i-k)*p+x] * powers[k*p+y]) % p;
        value = (value + g[point] * term) % p;
    }
    value = (value * binom[i*p+k]) % p;
    if (value<0) value+=p;
    coeff[q]=value;
}
"""


def make_signed_edge_fields(
    p: int, edge_count: int, seed: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray, str]:
    _directions, _eta, nonsquare = directions_and_signs(p)
    legendre = legendre_table(p)
    rng = random.Random(seed)
    edges: set[tuple[int, int]] = set()
    while len(edges) < edge_count:
        u, v = rng.sample(range(p*p), 2)
        edges.add((min(u, v), max(u, v)))
    array = np.asarray(sorted(edges), dtype=np.int32)
    ux, uy = array[:, 0]//p, array[:, 0]%p
    vx, vy = array[:, 1]//p, array[:, 1]%p
    dx, dy = (vx-ux)%p, (vy-uy)%p
    tau = legendre[(dx*dx-nonsquare*dy*dy)%p]
    inv2=pow(2,-1,p)
    midpoint=(((ux+vx)*inv2)%p)*p+((uy+vy)*inv2)%p
    g=np.zeros(p*p,dtype=np.int32)
    np.add.at(g,midpoint,tau)
    signed_degree=np.zeros(p*p,dtype=np.int32)
    np.add.at(signed_degree,array[:,0],tau)
    np.add.at(signed_degree,array[:,1],tau)
    kappa=signed_degree-2*g
    return g, signed_degree, kappa, hashlib.sha256(array.tobytes()).hexdigest()


def line_profiles(g: np.ndarray, p: int) -> np.ndarray:
    profiles=np.zeros((p+1,p),dtype=np.int32)
    x=np.repeat(np.arange(p,dtype=np.int32),p)
    y=np.tile(np.arange(p,dtype=np.int32),p)
    for d in range(p):
        np.add.at(profiles[d],(x+d*y)%p,g)
    np.add.at(profiles[p],y,g)
    return profiles


def main() -> None:
    parser=argparse.ArgumentParser()
    parser.add_argument("--p",type=int,default=131)
    parser.add_argument("--edges",type=int,default=12289)
    parser.add_argument("--seed",type=int,default=15765)
    parser.add_argument("--mode",choices=("midpoint","curvature"),default="midpoint")
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    p=args.p
    if not is_prime(p) or p<29:
        raise ValueError("need a held-out prime p>=29")
    g,signed_degree,kappa,edge_hash=make_signed_edge_fields(p,args.edges,args.seed)
    field=g if args.mode=="midpoint" else kappa
    profiles=line_profiles(field,p)
    powers=np.empty((p,p),dtype=np.int32)
    powers[0]=1
    for i in range(1,p): powers[i]=(powers[i-1]*np.arange(p,dtype=np.int32))%p
    binom=np.zeros((p,p),dtype=np.int32)
    for i in range(p):
        for k in range(i+1): binom[i,k]=math.comb(i,k)%p

    platforms=cl.get_platforms()
    devices=[d for plat in platforms for d in plat.get_devices() if d.type & cl.device_type.GPU]
    if not devices: raise RuntimeError("no OpenCL GPU")
    device=devices[0]
    context=cl.Context([device]); queue=cl.CommandQueue(context)
    mf=cl.mem_flags
    g_b=cl.Buffer(context,mf.READ_ONLY|mf.COPY_HOST_PTR,hostbuf=field)
    profiles_b=cl.Buffer(context,mf.READ_ONLY|mf.COPY_HOST_PTR,hostbuf=profiles)
    powers_b=cl.Buffer(context,mf.READ_ONLY|mf.COPY_HOST_PTR,hostbuf=powers)
    binom_b=cl.Buffer(context,mf.READ_ONLY|mf.COPY_HOST_PTR,hostbuf=binom)
    direct=np.empty(p*(p+1),dtype=np.int32); incidence=np.empty_like(direct)
    coeff=np.zeros(p*p,dtype=np.int32)
    direct_b=cl.Buffer(context,mf.WRITE_ONLY,direct.nbytes)
    incidence_b=cl.Buffer(context,mf.WRITE_ONLY,incidence.nbytes)
    coeff_b=cl.Buffer(context,mf.READ_WRITE|mf.COPY_HOST_PTR,hostbuf=coeff)
    program=cl.Program(context,KERNEL).build()
    started=time.time()
    program.moments(queue,(direct.size,),None,g_b,profiles_b,powers_b,direct_b,incidence_b,np.int32(p))
    program.mixed(queue,(coeff.size,),None,g_b,powers_b,binom_b,coeff_b,np.int32(p))
    cl.enqueue_copy(queue,direct,direct_b); cl.enqueue_copy(queue,incidence,incidence_b)
    cl.enqueue_copy(queue,coeff,coeff_b); queue.finish()
    elapsed=time.time()-started
    direct=direct.reshape(p,p+1); incidence=incidence.reshape(p,p+1); coeff=coeff.reshape(p,p)
    evaluated=np.zeros_like(direct)
    for i in range(p):
        for lam in range(p):
            value=0
            for k in range(i,-1,-1): value=(value*lam+int(coeff[i,k]))%p
            evaluated[i,lam]=value
        evaluated[i,p]=coeff[i,i]
    errors_a=(direct-incidence)%p; errors_b=(direct-evaluated)%p
    profile_energy=int(np.dot(profiles.astype(np.int64).ravel(),profiles.astype(np.int64).ravel()))
    field_energy=int(np.dot(field.astype(np.int64),field.astype(np.int64)))
    boundary_size=int(np.count_nonzero(signed_degree%2))
    result={
        "classification":"exact all-degree large-prime verification; symbolic identity candidate, not proof by sampling",
        "host":platform.node(),"architecture":platform.machine(),
        "gpu_backend":"OpenCL","gpu_device":device.name.strip(),
        "p":p,"random_simple_edges":args.edges,"edge_list_sha256":edge_hash,"seed":args.seed,
        "mode":args.mode,
        "degree_blocks_checked":p,"projective_directions":p+1,
        "incidence_vs_direct_nonzero_errors":int(np.count_nonzero(errors_a)),
        "homogeneous_form_vs_direct_nonzero_errors":int(np.count_nonzero(errors_b)),
        "extra_b0_p_torsion_constraints_beyond_equal_totals":p*(p-1)//2,
        "full_b0_cokernel_including_equal_totals":p*(p+1)//2,
        "identity":"sum_a R_L(f)(a)a^i = sum_x f(x)L(x)^i = F_i(L)",
        "field":"g" if args.mode=="midpoint" else "kappa=d_z-2*g",
        "field_total":int(field.sum()),
        "field_energy":field_energy,
        "profile_energy":profile_energy,
        "point_Radon_parseval_exact_for_zero_sum":bool(
            int(field.sum())!=0 or profile_energy==p*field_energy
        ),
        "graph_boundary_size":boundary_size,
        "curvature_parity_equals_graph_boundary":bool(
            args.mode!="curvature" or np.array_equal(kappa%2,signed_degree%2)
        ),
        "curvature_energy_dominates_p_boundary":bool(
            args.mode!="curvature" or profile_energy>=p*boundary_size
        ),
        "all_degree_seam_validated":bool(not np.any(errors_a) and not np.any(errors_b)),
        "elapsed_seconds":elapsed,
        "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps(result,indent=2,sort_keys=True),flush=True)


if __name__=="__main__": main()
