#!/usr/bin/env python3
"""Exact integer midpoint-seam model for one Proposition 15.758 template.

Eliminate the hidden diagonal profiles by an integer midpoint function g:

    D_L(a) = eta_L * sum_{x:Lx=a} g(x) - A_L(a) >= 0.

An actual h-edge graph also forces sum(g)=T and ||g||_1<=h.  This model
minimizes ||g||_1.  Feasibility is only the uncollapsed midpoint seam; it is
not a full simple-edge lift and is never reported as residual closure.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import time
from pathlib import Path

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import csr_matrix, hstack, vstack

from residual_midpoint_seam_gpu import compact_template


def point_line_matrix(p: int) -> csr_matrix:
    rows: list[int] = []
    cols: list[int] = []
    for direction in range(p):
        for label in range(p):
            row = direction * p + label
            for y in range(p):
                x = (label - direction * y) % p
                rows.append(row); cols.append(x * p + y)
    for label in range(p):
        row = p * p + label
        for x in range(p):
            rows.append(row); cols.append(x * p + label)
    return csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(p*(p+1), p*p))


def main() -> None:
    parser=argparse.ArgumentParser()
    parser.add_argument("--p",type=int,required=True)
    parser.add_argument("--branch",choices=("B","C"),required=True)
    parser.add_argument("--endpoint",choices=("lower","upper"),default="lower")
    parser.add_argument("--template-seed",type=int,required=True)
    parser.add_argument("--time-limit",type=float,default=120.0)
    parser.add_argument("--solver",choices=("highs","cpsat"),default="cpsat")
    parser.add_argument("--workers",type=int,default=16)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    template=compact_template(args.p,args.branch,args.endpoint,args.template_seed)
    p=args.p; n=p*p; h=int(template["edge_count"]); T=int(template["signed_total"])
    eta=np.asarray(template["eta"],dtype=np.int64)
    anti=np.asarray(template["anti"],dtype=np.int64)
    R=point_line_matrix(p)
    # eta*R*(gplus-gminus) >= anti.
    signed_R=R.multiply(np.repeat(eta,p)[:,None])
    inequality=hstack([-signed_R,signed_R],format="csr")
    total=csr_matrix(np.r_[np.ones(n),-np.ones(n)][None,:])
    l1=csr_matrix(np.ones((1,2*n)))
    matrix=vstack([inequality,total,l1],format="csr")
    lower=np.r_[np.full(p*(p+1),-np.inf),T,-np.inf]
    upper=np.r_[-anti.ravel(),T,h]
    started=time.time()
    if args.solver == "highs":
        result=milp(
            np.ones(2*n),integrality=np.ones(2*n),bounds=Bounds(0,np.inf),
            constraints=LinearConstraint(matrix,lower,upper),
            options={"time_limit":args.time_limit,"mip_rel_gap":0.0,"presolve":True},
        )
        result_x = result.x
        result_status = int(result.status)
        result_message = str(result.message)
        result_fun = result.fun
        result_bound = getattr(result,"mip_dual_bound",None)
        result_gap = getattr(result,"mip_gap",None)
    else:
        from ortools.sat.python import cp_model
        model=cp_model.CpModel()
        g_vars=[model.new_int_var(-h,h,f"g_{i}") for i in range(n)]
        abs_vars=[model.new_int_var(0,h,f"a_{i}") for i in range(n)]
        for gv,av in zip(g_vars,abs_vars): model.add_abs_equality(av,gv)
        model.add(sum(g_vars)==T); model.add(sum(abs_vars)<=h)
        for direction in range(p+1):
            sign=int(eta[direction])
            for label in range(p):
                if direction==p:
                    indices=[x*p+label for x in range(p)]
                else:
                    indices=[((label-direction*y)%p)*p+y for y in range(p)]
                model.add(sign*sum(g_vars[index] for index in indices)>=int(anti[direction,label]))
        model.minimize(sum(abs_vars))
        solver=cp_model.CpSolver(); solver.parameters.max_time_in_seconds=args.time_limit
        solver.parameters.num_search_workers=args.workers
        solver.parameters.random_seed=args.template_seed
        status=solver.solve(model)
        feasible=status in (cp_model.FEASIBLE,cp_model.OPTIMAL)
        result_x=None
        if feasible:
            g_found=np.array([solver.value(var) for var in g_vars],dtype=np.int64)
            result_x=np.r_[np.maximum(g_found,0),np.maximum(-g_found,0)]
        result_status=int(status); result_message=solver.status_name(status)
        result_fun=float(solver.objective_value) if feasible else None
        result_bound=float(solver.best_objective_bound) if feasible else None
        result_gap=(None if not feasible or result_fun==0 else (result_fun-result_bound)/abs(result_fun))
    elapsed=time.time()-started
    payload: dict[str,object]={
        "classification":"exact integer midpoint-seam optimization for one randomized template; not full graph lift",
        "host":platform.node(),"architecture":platform.machine(),"p":p,
        "branch":args.branch,"endpoint":args.endpoint,"template_seed":args.template_seed,
        "template_sha256":template["template_sha256"],"edge_count_h":h,"signed_total_T":T,
        "solver":args.solver,"status":result_status,
        "message":result_message,"elapsed_seconds":elapsed,
        "mip_dual_bound":None if result_bound is None else float(result_bound),
        "mip_gap":None if result_gap is None else float(result_gap),
    }
    if result_x is not None:
        plus=np.rint(result_x[:n]).astype(np.int64); minus=np.rint(result_x[n:]).astype(np.int64)
        g=plus-minus; line=np.asarray(R@g).ravel(); D=(np.repeat(eta,p)*line-anti.ravel()).reshape(p+1,p)
        l1_value=int(np.abs(g).sum())
        exact=bool(int(g.sum())==T and int(D.min())>=0 and l1_value<=h)
        g_hash=hashlib.sha256(g.astype("<i8").tobytes()).hexdigest()
        payload.update({
            "incumbent_found":True,"minimum_l1_reported":float(result_fun),
            "exact_replayed_l1":l1_value,"l1_room_h_minus_l1":h-l1_value,
            "midpoint_function_min":int(g.min()),"midpoint_function_max":int(g.max()),
            "midpoint_function_support":int(np.count_nonzero(g)),"midpoint_function_sha256":g_hash,
            "minimum_hidden_diagonal_count":int(D.min()),"maximum_hidden_diagonal_count":int(D.max()),
            "hidden_diagonal_sums_match_parallel_counts":bool(np.array_equal(D.sum(axis=1),template["parallel"])),
            "necessary_midpoint_cancellation_parity":int((h-l1_value)%2),
            "integer_midpoint_seam_feasible":exact,
            "full_simple_edge_lift_proved":False,
            "g_values":g.tolist() if exact else None,
        })
    else:
        payload.update({"incumbent_found":False,"integer_midpoint_seam_feasible":False,"full_simple_edge_lift_proved":False})
    payload["script_sha256"]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    args.output.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps(payload,indent=2,sort_keys=True),flush=True)


if __name__=="__main__": main()
