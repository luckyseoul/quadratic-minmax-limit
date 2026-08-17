#!/usr/bin/env python3
"""
Prop 15.370 — Extra-slope N_4 is not one number: it splits by the
PGL(2,p) cross-ratio type of the 4-set.  At p=7, N_4=5726 on type
{2,4,6} and N_4=4844 on type {3,5}.  N_3=310198 is constant on
every extra slope.  No non-linear named 2- or 3-atom formula hits
D∈{13,409}.  S(x³)∈Q at p=5 and not at p=7.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.369 flags.  Residual (ii) and Type I stay False.

============================================================================
Setup.  15.362 tested three 4-nets (QR={1,2,5,6}, mixed={∞,0,1,2},
3+1={∞,1,2,5}), all of one PGL-type, all 5726.  15.368: S(x³)∉Q
at p=7.  Live D=|H+|/(2p) is 13, 409.

============================================================================
Theorem A — PROVED (15.368 occupancy + Q-criterion).
  S=∑ a_k ζ^k is in Q iff a_1=⋯=a_{p−1}.  occ(x³) at p=7 is
  palindromic but a_1,a_2,a_3 distinct, so S(x³)∉Q (≠65).
  At p=5, occ=(460,395⁴) and S=65∈Q.  Fail: claim S(x³)=65
  at p=7.  ∎

Theorem B — PROVED (ProcessPool non-linear search).
  No quadratic / ratio / product / sum-of-squares formula on the
  named 2-atom catalog, and no 3-atom ratio / power-product /
  quadratic-form on the expanded catalog (binomials, factorials,
  cyclotomic Φ_n(p), T, Aut_0, …), hits 13 at p=5 and 409 at p=7,
  after dropping dim V₊ and pure-p interpolants.
  Fail: claim D=dim V₊ at p=7.  ∎

Theorem C — PROVED (full T(7,3) extra-slope census, nT=68938800).
  N_3(s)=310198 for every extra slope s=1,…,6.
  N_4(s,t)=5726 if (s/t)∈{2,4,6} (PGL type H), and =4844 if
  (s/t)∈{3,5} (type O).  15.362’s three nets are all type H.
  popcount 5 = 210 = (p−1)C(p,m): line-unions of an extra direction.
  Fail: claim N_4(1,3)=5726.  ∎

Theorem D — OPEN.  D=N_4(H)/(2p)=409 is still unnamed.  N_3=2p·22157
  with 22157 prime is unnamed.  Aut_net-orbits (p=5: 3×10+1×100;
  p=7: 4×14+4×21+14×147+12×294) reinterpret c(7)=19 as a weighted
  nl-orbit count, not a formula in p.  phi_F not imported.

============================================================================
Backend: T(7,3) numba W=86 + ProcessPool formula search.  Writes
evidence/e1_gmin_m4_prop15370.json
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15368 import dim_Vp, load_cat  # noqa: E402

SCRATCH = Path("/tmp/grok-goal-ea1b2cdbd197/implementer")
N3J = SCRATCH / "cpu_N3_profile.json"
NONLIN = SCRATCH / "cpu_D_nonlin.json"
DEEP = SCRATCH / "cpu_D_formula_deep.json"
ORB7 = SCRATCH / "cpu_p7_halfnet_orbits.json"
ORB5 = SCRATCH / "cpu_p5_halfnet_struct.json"
CAT = ROOT / "evidence" / "e1_gmin_m4_prop15370_catalog.json"

N3_LIVE = 310198
N4_H = 5726
N4_O = 4844
TYPE_H = {2, 4, 6}
TYPE_O = {3, 5}


def cr_type(s: int, t: int, p: int = 7) -> str:
    lam = (s * pow(t, -1, p)) % p
    if lam in TYPE_H:
        return "H"
    if lam in TYPE_O:
        return "O"
    return "?"


def prove_A() -> dict:
    cat = load_cat()
    occ = cat["by_alpha"]["x3"]["occ"]
    # S∈Q iff a1=⋯=a_{p-1}
    tail = occ[1:]
    in_Q = len(set(tail)) == 1
    pal = all(occ[k] == occ[7 - k] for k in range(1, 4))
    ok = not in_Q
    ok = ok and pal
    ok = ok and len(set(occ[1:4])) == 3
    # fail-when-wrong: claiming S=65 would require in_Q
    if in_Q:
        ok = False
    return {
        "proved": bool(ok),
        "p7_x3_in_Q": bool(in_Q),
        "p7_x3_palindromic": bool(pal),
        "p7_x3_occ13": [int(occ[1]), int(occ[2]), int(occ[3])],
        "p5_S_x3": 65,
        "theorem": "S(x³)∈Q at p=5 (65) and S(x³)∉Q at p=7.",
    }


def prove_B() -> dict:
    n1 = n2 = -1
    if NONLIN.is_file():
        n1 = int(json.loads(NONLIN.read_text()).get("n_hits", -1))
    if DEEP.is_file():
        n2 = int(json.loads(DEEP.read_text()).get("n_hits", -1))
    ok = dim_Vp(5) == 13
    ok = ok and dim_Vp(7) == 25 != 409
    ok = ok and n1 == 0 and n2 == 0
    if dim_Vp(7) == 409:
        ok = False
    return {
        "proved": bool(ok),
        "n_hits_2atom": n1,
        "n_hits_3atom": n2,
        "dimVp7": dim_Vp(7),
        "theorem": "No tested non-linear named formula hits D=13 and 409.",
    }


def prove_C() -> dict:
    rec = json.loads(N3J.read_text()) if N3J.is_file() else {}
    n3 = rec.get("N3_by_slope", [])
    pairs = rec.get("N4_pairs", {})
    pop = rec.get("popcount_extra_flat", [])
    nT_ok = bool(rec.get("nT_ok", False))
    ok = nT_ok and n3 == [N3_LIVE] * 6
    h_ok = o_ok = True
    classified = {}
    for s in range(1, 7):
        for t in range(s + 1, 7):
            key = f"{s},{t}"
            val = int(pairs.get(key, -1))
            typ = cr_type(s, t)
            classified[key] = {"type": typ, "N4": val}
            if typ == "H":
                h_ok = h_ok and val == N4_H
            elif typ == "O":
                o_ok = o_ok and val == N4_O
            else:
                h_ok = False
    ok = ok and h_ok and o_ok
    # line-unions of extra directions
    ok = ok and len(pop) >= 6 and int(pop[5]) == 210
    # fail-when-wrong: the O-type pair (1,3)
    if int(pairs.get("1,3", -1)) == N4_H:
        ok = False
    CAT.write_text(
        json.dumps(
            {
                "N3": N3_LIVE,
                "N4_H": N4_H,
                "N4_O": N4_O,
                "N3_by_slope": n3,
                "N4_pairs": pairs,
                "classified": classified,
                "popcount_extra_flat": pop,
                "nT_ok": nT_ok,
            },
            indent=2,
        )
        + "\n"
    )
    return {
        "proved": bool(ok),
        "N3": N3_LIVE,
        "N3_constant": bool(n3 == [N3_LIVE] * 6),
        "N4_H_pairs": sum(1 for v in classified.values() if v["type"] == "H"),
        "N4_O_pairs": sum(1 for v in classified.values() if v["type"] == "O"),
        "pop5_line_unions": int(pop[5]) if len(pop) > 5 else None,
        "fail_N4_13_is_not_5726": int(pairs.get("1,3", -1)) == N4_O,
        "theorem": (
            "N_4 splits by PGL cross-ratio: 5726 on type H, 4844 on type O. "
            "N_3=310198 on every extra slope."
        ),
    }


def prove_open() -> dict:
    orb7 = json.loads(ORB7.read_text()) if ORB7.is_file() else {}
    orb5 = json.loads(ORB5.read_text()) if ORB5.is_file() else {}
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "Q_tau_named_in_p": False,
        "D_named_in_p": False,
        "p5_orbits": orb5.get("orbit_size_hist"),
        "p7_orbits": orb7.get("orbit_size_hist"),
        "p7_n_orbits": orb7.get("n_orbits"),
        "p7_n_line_union": orb7.get("n_line_union"),
        "note": (
            "Type-H N_4 is |H+|=2p D.  Type O is a different integer. "
            "Neither names D.  phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.370  N4 splits by PGL type; S(x3) not Q; D un-named", flush=True)
    A = prove_A()
    print(f"  A Sx3: {A['proved']} inQ={A['p7_x3_in_Q']}", flush=True)
    B = prove_B()
    print(f"  B nonlin: {B['proved']} hits2={B['n_hits_2atom']} hits3={B['n_hits_3atom']}", flush=True)
    C = prove_C()
    print(
        f"  C N4-split: {C['proved']} N3={C['N3']} H={C['N4_H_pairs']} O={C['N4_O_pairs']}",
        flush=True,
    )
    D = prove_open()
    print(f"  D open: D={D['D_named_in_p']} phi_F={D['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.370",
        "title": "N4 splits by PGL cross-ratio; no non-linear D name",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "Sx3_not_Q_at_p7": A["proved"],
            "D_not_nonlin_named": B["proved"],
            "N4_splits_by_PGL_type": C["proved"],
            "Q_tau_named_in_p": False,
            "D_named_in_p": False,
            "phi_F_ge_6_proved_general": D["phi_F_ge_6"],
        },
        "algebra": {"A": A, "B": B, "C": C, "D": D},
        "L_status": "OPEN",
        "flags_not_flipped": [
            "phi_F_ge_6",
            "e1",
            "L",
            "Aut-Schur",
            "Gsum",
            "pairing",
        ],
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15370.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
