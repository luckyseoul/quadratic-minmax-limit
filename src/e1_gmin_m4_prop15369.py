#!/usr/bin/env python3
"""
Prop 15.369 — The only named α with S(α) given by one formula at both
p=5 and p=7 is affine: S=T.  Deg-2 constancy S(x²)=S(x²+x) holds at
p=5 (both 40) and dies at p=7 (occ0 179539290 ≠ 230044675).  D is
not a 3-term linear form in the expanded named catalog.  phi_F not
imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.368 flags.  Residual (ii) and Type I stay False.

============================================================================
Setup.  15.363: affine α ⇒ ⟨α,Diag−m⟩≡0 on T, so S=T.  15.368:
S(x²)∉Q at p=7.  Live D=|H+|/(2p) is 13, 409.  dim V₊ dies.

============================================================================
Theorem A — PROVED (15.363 algebra, all p; certified p=5 census).
  For α(c)=uc+v, ⟨α,Diag−m⟩=0 on every A∈T, hence S(α)=T(p,m).
  At p=5, occ(const)=occ(x)=(2040,0,0,0,0).  Fail: claim S(x)=40.  ∎

Theorem B — PROVED (p=5 census + p=7 T(7,3) occupancy).
  S(x²)=S(x²+x)=40 at p=5 (identical occ).  At p=7,
  occ(x²)[0]=179539290 ≠ occ(x²+x)[0]=230044675, so deg-2
  constancy dies.  Fail: claim S(x²)=S(x²+x) at p=7.  ∎

Theorem C — PROVED (ProcessPool 3-term search).
  No aX+bY+cZ with |a,b,c|≤2 on the named list
  {n_1d,n_pp,A,k,G,dim V₊,q,p,U,|T|,m,(p+1)/2,n_Hoff,Φ_3,Φ_4,Φ_6,Φ_12}
  equals 13 at p=5 and 409 at p=7.  dim V₊=25≠409.  Fail: claim
  D=dim V₊ at p=7.  ∎

Theorem D — OPEN.  Affine S=T does not name D.  Deg-2 is not a
  constant-S type in p.  phi_F not imported.

============================================================================
Backend: p=5 T(5,2) + 15.368 catalog.  Writes
evidence/e1_gmin_m4_prop15369.json
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
from e1_gmin_m4_prop15363 import gen_T, occupancy, poly_alpha, wrap_diag_sums  # noqa: E402
from e1_gmin_m4_prop15368 import dim_Vp  # noqa: E402

CAT7 = ROOT / "evidence" / "e1_gmin_m4_prop15368_catalog.json"
DDEEP = Path("/tmp/grok-goal-ea1b2cdbd197/implementer/cpu_D_deep.json")
CAT = ROOT / "evidence" / "e1_gmin_m4_prop15369_catalog.json"


def p5_occ(coeffs) -> tuple:
    p, m = 5, 2
    mats = gen_T(p, m)
    diags = [wrap_diag_sums(g, 1) for g in mats]
    return occupancy(diags, poly_alpha(coeffs, p), m), len(mats)


def prove_A() -> dict:
    occ_x, nT = p5_occ((0, 1))
    occ_1, _ = p5_occ((1,))
    ok = occ_x == (nT,) + (0,) * 4
    ok = ok and occ_1 == occ_x
    if occ_x[0] == 40:
        ok = False
    return {
        "proved": bool(ok),
        "p5_S_x_is_T": occ_x[0] == nT and set(occ_x[1:]) == {0},
        "nT": nT,
        "theorem": "Affine S=T at p=5; algebra of 15.363 gives it for all p.",
    }


def prove_B() -> dict:
    occ_x2, _ = p5_occ((0, 0, 1))
    occ_x2x, _ = p5_occ((0, 1, 1))
    cat7 = json.loads(CAT7.read_text())
    o2 = cat7["by_alpha"]["x2"]["occ0"]
    o2x = cat7["by_alpha"]["x2+x"]["occ0"]
    ok = occ_x2 == occ_x2x
    ok = ok and occ_x2[0] - occ_x2[1] == 40
    ok = ok and o2 == 179539290
    ok = ok and o2x == 230044675
    ok = ok and o2 != o2x
    if o2 == o2x:
        ok = False
    return {
        "proved": bool(ok),
        "p5_x2_eq_x2x": occ_x2 == occ_x2x,
        "p5_S": occ_x2[0] - occ_x2[1],
        "p7_x2_occ0": o2,
        "p7_x2x_occ0": o2x,
        "theorem": "S(x²)=S(x²+x)=40 at p=5; occupancies split at p=7.",
    }


def prove_C() -> dict:
    rec = {}
    if DDEEP.is_file():
        rec = json.loads(DDEEP.read_text())
    n_hits = int(rec.get("n_hits", -1))
    ok = dim_Vp(5) == 13
    ok = ok and dim_Vp(7) == 25 != 409
    ok = ok and n_hits == 0
    if dim_Vp(7) == 409:
        ok = False
    # persist a small catalog
    CAT.write_text(
        json.dumps(
            {
                "n_hits": n_hits,
                "dimVp": {5: dim_Vp(5), 7: dim_Vp(7)},
                "source": "cpu_D_deep.json",
            },
            indent=2,
        )
        + "\n"
    )
    return {
        "proved": bool(ok),
        "n_hits": n_hits,
        "dimVp7": dim_Vp(7),
        "theorem": "No 3-term named linear form hits D=13 and 409; dim V+ dies.",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "Q_tau_named_in_p": False,
        "D_named_in_p": False,
        "note": (
            "Affine S=T is the only tested type with one formula at "
            "both primes.  It does not name D.  phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.369  affine S=T only; deg-2 constancy dies", flush=True)
    A = prove_A()
    print(f"  A affine: {A['proved']} S_x=T {A['p5_S_x_is_T']}", flush=True)
    B = prove_B()
    print(
        f"  B deg2: {B['proved']} p5eq={B['p5_x2_eq_x2x']} p7 {B['p7_x2_occ0']}≠{B['p7_x2x_occ0']}",
        flush=True,
    )
    C = prove_C()
    print(f"  C D-lin: {C['proved']} hits={C['n_hits']} dim7={C['dimVp7']}", flush=True)
    D = prove_open()
    print(f"  D open: D={D['D_named_in_p']} phi_F={D['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.369",
        "title": "Only affine S=T is constant at p=5 and p=7; deg-2 dies",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "affine_S_is_T": A["proved"],
            "deg2_constancy_dies_p7": B["proved"],
            "D_not_3term_named": C["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15369.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
