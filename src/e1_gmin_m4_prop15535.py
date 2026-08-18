#!/usr/bin/env python3
"""
Prop 15.535 — Twisted involutions σ(x)=a x̄+b, a=±v^{1-p},
N(a)=1, b+a x̄(b)=0, freely pair T_ns-orbits at p=7 (all 14)
and at no p=5 (0/10; fix∈{2,30}, inv∈{2,6}).  Fail: a=+
is free at p=5.  Inversion+χ preserves H_+ but does not
send T_ns-orbits to T_ns-orbits.  D=n_orb/2 is a free
⟨T_ns,σ⟩ count only for p≡3.  Q_τ unnamed.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.534 flags.  Does **not** overwrite 15.522–15.534.
"""
from __future__ import annotations

import json
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import _kernel_field
from e1_gmin_m4_prop15526 import load_Hplus
from e1_gmin_m4_prop15532 import DIRS

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / "evidence" / "e1_gmin_m4_prop15535.json"

W = 48
NS = {5: DIRS[5]["ns"], 7: DIRS[7]["ns"]}


def _key(y: np.ndarray) -> bytes:
    return np.ascontiguousarray(y.astype(np.int8)).tobytes()


def _pow_el(F, x: int, e: int) -> int:
    if x == 0:
        return 0
    r, b = 1, int(x)
    ee = int(e) % (F["q"] - 1)
    if ee < 0:
        ee += F["q"] - 1
    while ee:
        if ee & 1:
            r = F["mul"][r][b]
        b = F["mul"][b][b]
        ee >>= 1
    return r


def _frob(F, x: int) -> int:
    return _pow_el(F, x, F["p"])


def _norm(F, x: int) -> int:
    return F["mul"][x][_frob(F, x)]


def a_pm(F, v: int) -> tuple[int, int]:
    """a_+=v^{1-p}=v/v^p, a_-=-a_+."""
    vp = _frob(F, v)
    a_plus = F["mul"][v][_pow_el(F, vp, F["q"] - 2)]
    a_minus = F["mul"][F["neg"][1]][a_plus]
    return a_plus, a_minus


def b_solutions(F, a: int) -> list[int]:
    out = []
    for b in range(F["q"]):
        if F["add"][b][F["mul"][a][_frob(F, b)]] == 0:
            out.append(b)
    return out


def _apply_sigma(y: np.ndarray, F, a: int, b: int) -> np.ndarray:
    q = F["q"]
    z = np.empty_like(y)
    z[0] = 1.0
    for x in range(q):
        dest = F["add"][F["mul"][a][_frob(F, x)]][b]
        z[1 + dest] = y[1 + x]
    return z


def _sigma_job(spec) -> dict:
    p, lab, b = spec
    F = _kernel_field(p)
    H = load_Hplus(p)
    va, vb = NS[p]
    v = va + vb * p
    ap, am = a_pm(F, v)
    a = ap if lab == "+" else am
    table = {_key(H[i]): i for i in range(len(H))}
    perm_T = [F["add"][x][v] for x in range(F["q"])]

    def push_T(y):
        z = np.empty_like(y)
        z[0] = 1.0
        for x in range(F["q"]):
            z[1 + perm_T[x]] = y[1 + x]
        return z

    oid = -np.ones(len(H), dtype=np.int32)
    n_orb = 0
    for i in range(len(H)):
        if oid[i] >= 0:
            continue
        j = i
        while oid[j] < 0:
            oid[j] = n_orb
            j = table[_key(push_T(H[j]))]
        n_orb += 1
    img = []
    nfix = 0
    for i in range(len(H)):
        z = _apply_sigma(H[i], F, a, b)
        j = table[_key(z)]
        img.append(j)
        if j == i:
            nfix += 1
    inv = set()
    for i in range(len(H)):
        if oid[i] == oid[img[i]]:
            inv.add(int(oid[i]))
    return {
        "p": p,
        "lab": lab,
        "a": a,
        "b": b,
        "N_a": _norm(F, a),
        "n_orb": n_orb,
        "nfix": nfix,
        "n_inv": len(inv),
        "free": nfix == 0 and len(inv) == 0,
    }


def _jobs_for(p: int) -> list[tuple]:
    F = _kernel_field(p)
    va, vb = NS[p]
    v = va + vb * p
    ap, am = a_pm(F, v)
    out = []
    for lab, a in (("+", ap), ("-", am)):
        for b in b_solutions(F, a):
            out.append((p, lab, b))
    return out


def prove_A() -> dict:
    """Field identities: N(a)=1, p solutions b, involution."""
    ok = True
    rec = {}
    for p in (5, 7, 11, 13):
        F = _kernel_field(p)
        va, vb = (1, 2) if p == 5 else (1, 1)
        if p >= 11:
            # any ns: packed 1+s p with χ=-1
            v = None
            for s in range(p):
                cand = 1 + s * p
                if int(F["chi"][cand]) == -1:
                    v = cand
                    break
        else:
            v = va + vb * p
        ap, am = a_pm(F, v)
        bs_p = b_solutions(F, ap)
        bs_m = b_solutions(F, am)
        n_ok = (
            _norm(F, ap) == 1
            and _norm(F, am) == 1
            and len(bs_p) == p
            and len(bs_m) == p
        )
        # involution on a sample of F_q
        a, b0 = ap, bs_p[0]

        def fn(x):
            return F["add"][F["mul"][a][_frob(F, x)]][b0]

        id_ok = all(fn(fn(x)) == x for x in range(F["q"]))
        if not (n_ok and id_ok):
            ok = False
        rec[str(p)] = {
            "N_plus": _norm(F, ap),
            "N_minus": _norm(F, am),
            "n_b_plus": len(bs_p),
            "n_b_minus": len(bs_m),
            "involution": id_ok,
        }
    if rec["5"]["N_plus"] != 1:
        ok = False
    return {
        "proved": bool(ok),
        "by_p": rec,
        "theorem": (
            "a=±v^{1-p} has N(a)=1; #{b: b+a b̄=0}=p; σ is an involution. "
            "Fail N(a)≠1."
        ),
    }


def prove_B() -> dict:
    jobs = _jobs_for(5) + _jobs_for(7)
    with ProcessPoolExecutor(max_workers=min(W, len(jobs))) as ex:
        rows = list(ex.map(_sigma_job, jobs))
    byp = {5: [], 7: []}
    for r in rows:
        byp[r["p"]].append(r)
    ok = True
    rec = {}
    for p in (5, 7):
        rs = byp[p]
        n_free = sum(1 for r in rs if r["free"])
        fixes = sorted({r["nfix"] for r in rs})
        invs = sorted({r["n_inv"] for r in rs})
        if p == 5:
            row_ok = n_free == 0 and fixes == [2, 30] and invs == [2, 6]
        else:
            row_ok = n_free == len(rs) == 14 and fixes == [0] and invs == [0]
        if not row_ok:
            ok = False
        rec[str(p)] = {
            "n_maps": len(rs),
            "n_free": n_free,
            "fix_set": fixes,
            "inv_set": invs,
            "ok": row_ok,
        }
    # fail-when-wrong: a=+ free at p=5
    plus5 = [r for r in byp[5] if r["lab"] == "+"]
    if any(r["free"] for r in plus5):
        ok = False
    return {
        "proved": bool(ok),
        "by_p": rec,
        "theorem": (
            "All 14 twisted σ are free at p=7; none of 10 at p=5 "
            "(fix∈{2,30}, inv∈{2,6}). Fail a=+ free at p=5."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "D_named_as_orbits": False,
        "free_pairing_p5": False,
        "Q_tau_named_in_p": False,
        "phi_F_imported": False,
        "note": (
            "⟨T_ns,σ⟩ is free only for p≡3 in the named twisted family. "
            "p=5 has no free Aut_∞-normalizer of T_ns. D and Q_τ unnamed. "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("15.535 twisted σ free on T_ns-orbits at p=7, not p=5", flush=True)
    A, B, C = prove_A(), prove_B(), prove_open()
    out = {
        "prop": "15.535",
        "title": "Twisted σ=a x̄+b freely pairs T_ns at p=7, not p=5",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "C": C,
        "e1_closed_general": bool(e1_closed_general()),
        "gsum_disj_lb": bool(gsum_disj_lb_proved_general()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "flags_not_flipped": [
            "phi_F_ge_6",
            "e1",
            "L",
            "Aut-Schur",
            "Gsum",
            "pairing",
        ],
        "L_status": "OPEN",
    }
    EV.parent.mkdir(parents=True, exist_ok=True)
    EV.write_text(json.dumps(out, indent=2) + "\n")
    print(f"A={A['proved']} B={B['proved']} phi_F={out['phi_F_ge_6']}", flush=True)
    return out


if __name__ == "__main__":
    main()
