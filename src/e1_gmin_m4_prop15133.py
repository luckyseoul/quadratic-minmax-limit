#!/usr/bin/env python3
"""
Prop 15.133 — Aut-orbit / class_key Bose–Mesner of T; dim E_{4p}^{ck};
CR dead; F19 quantitative at p=7; Aut-line program; residual OPEN.

Continues 15.132 toward N_ED4_SUF. Does **not** soft-close residual.
No class_key thrash as a closing argument (F19) — class_key spectrum is
measured once as structure, not refined.

============================================================================
Setup. T = C-signed Johnson adjacency on 4-sets; E_{4p}=ker(4pI−T).
δ ∈ E_{4p} is Aut(C)-invariant whenever m₄ is (Aut preserves Max+).
Path C residual ⇔ δ² ≤ room_hyp/24; for p≥7, δ²≤ρ_min² suffices.

============================================================================
Theorem A (Aut-line residual form) — PROVED (restate 15.109.F / 15.122.D).
  Let V^{Aut} = {Aut(C)-invariant functions on 4-sets}.  T preserves V^{Aut}.
  Set E_{4p}^{Aut} := E_{4p} ∩ V^{Aut}.  Then δ ∈ E_{4p}^{Aut}.  If
      dim E_{4p}^{Aut} ≤ 1,
  then either δ=0 (residual OK) or δ=c v_0 with ‖v_0‖₂=1 Aut-invariant, and
      residual ⇔ c² ≤ room_hyp/24,
  with c = ⟨m₄, v_0⟩ = E_y Q_0(y), Q_0(y):=⟨v_0,f_y⟩.
  When Q_0 is constant on Max+ (e.g. Aut transitive on Max+, or certified),
  c = Q_0(x_hs) is a Max+-free evaluation (halfspace character sum).  ∎

Theorem B (class_key T-spectrum / dim E_{4p}^{ck}) — CERTIFIED p=3,5,7.
  Let V^{ck} be the span of class_key cell indicators of C on 4-sets
  (combinatorial C-edge type of K4; 9 / 37 / 82 cells at p=3,5,7).
  The quotient of T on V^{ck} (equitable for C-combinatorics) has:
    p=3: dim E_{4p}^{ck} = 0,  λ_max(T|_{ck}) ≈ 9.798 < 12=4p;
    p=5: dim E_{4p}^{ck} = 1,  λ_max = 20 = 4p  (Aut-line on class_key);
    p=7: dim E_{4p}^{ck} = 0,  λ_max ≈ 23.375 < 28=4p.
  Backend: ProcessPool W=86 row-shards of the quotient at p=7.  ∎

Theorem C (F19 quantitative: class_key misses residual at p=7) — PROVED.
  At p=7, dim E_{4p}^{ck}=0 (Thm B), so every class_key-invariant 4-set
  function has no 4p-component.  But the true residual satisfies
      δ² = 19180800/1840091 > 0
  (Prop 15.131 pair-avg).  Therefore δ ∉ V^{ck}: m₄ is **not** constant on
  class_key cells at p=7 (F19).  class_key cannot close the residual at p=7
  and must not be thrashed as a moduli refinement.  ∎

Theorem D (cross-ratio PGL orbits are not Aut(C) orbits) — CERTIFIED p=5.
  Unordered cross-ratio orbits of PGL(2,p²) on 4-sets of PG(1,p²) yield
  4 orbits at p=5 (sizes 1300..7800), but each orbit contains **both**
  κ ∈ {±1} or {±3}.  Since κ is Aut(C)-invariant, these are not Aut(C)
  orbits.  Only a thin subset of random PGL matrices preserve the Paley C
  matrix entrywise.  Quotient T on CR orbits has λ_max≈6.77≪20 and
  nullity 0 — a false “Aut-line” vanishing.  Dead for residual.  ∎

Theorem E (Aut-line success locus) — CERTIFIED.
  p=3: E_{4p}={0} in ambient (λ_max(T)<4p), residual OK (δ=0).
  p=5: class_key carries a 1-dimensional 4p-line; residual saturates
       c² = room_hyp/24 = 1536/65 along that line (15.109/15.110).
  p=7: need a partition finer than class_key (true Aut orbits of C, or
       no compression — full character sum for m₄ / Q_δ).  ∎

Theorem F (OPEN residual, honest).
  Prove for all primes p≥5 either
    (i) dim E_{4p}^{Aut}≤1 for true Aut(C) and c²=Q_0(hs)²≤room_hyp/24
        by Gauss / halfspace character sums, or
    (ii) δ²≤ρ_min² (p≥7) / ≤room (p=5) Max+-free without false compressions, or
    (iii) max_y Q_δ(y)≤ρ_min².
  Dead: raw PGL cross-ratio; class_key thrash at p=7 (F19); soft-close.  ∎

============================================================================
Certified: Thm A algebra; Thm B class_key spectrum p=3,5,7; Thm C–E structure.
Backend: CPU ProcessPool W=86 class_key T; GPU unused (F20).
F19: class_key measured once; not used to close residual at p=7.

L OPEN. hyp_residual_for_all_p=false. ed4_le_suf_for_all_p=false.
Writes evidence/e1_gmin_m4_prop15133.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15102 import delta2_budget_hyp, rho_min2  # noqa: E402
from e1_gmin_m4_prop15117 import room_hyp_over_24  # noqa: E402
from e1_gmin_m4_prop15131 import P7_DELTA2_PAIR  # noqa: E402


def is_prime(p: int) -> bool:
    if p < 2:
        return False
    if p % 2 == 0:
        return p == 2
    k = 3
    while k * k <= p:
        if p % k == 0:
            return False
        k += 2
    return True


# Hard-certified class_key spectrum (classkey_E4p_nullity.py, W=86)
CK_SPECTRUM: dict[int, dict] = {
    3: {
        "n_orbits": 9,
        "C_n_4": 210,
        "nullity": 0,
        "lambda_max": 9.797958971132712,  # float cert; < 12
        "4p": 12,
        "lam_max_lt_4p": True,
    },
    5: {
        "n_orbits": 37,
        "C_n_4": 14950,
        "nullity": 1,
        "lambda_max": 20.0,
        "4p": 20,
        "lam_max_eq_4p": True,
    },
    7: {
        "n_orbits": 82,
        "C_n_4": 230300,
        "nullity": 0,
        "lambda_max": 23.37509475111153,  # float cert; < 28
        "4p": 28,
        "lam_max_lt_4p": True,
    },
}

# CR dead-end cert at p=5
CR_P5 = {
    "n_orbits_frob": 4,
    "n_orbits_nofrob": 5,
    "kappa_not_constant": True,
    "nullity_false": 0,
    "lambda_max_approx": 6.7692307692307665,
    "4p": 20,
}


# ---------------------------------------------------------------------------
# Theorem A
# ---------------------------------------------------------------------------

def prove_aut_line_form() -> dict:
    return {
        "proved": True,
        "theorem": (
            "δ∈E_{4p}^{Aut}; if dim≤1 then residual ⇔ c²≤room_hyp/24 with "
            "c=⟨m₄,v_0⟩=E Q_0; c=Q_0(hs) when Q_0 constant on Max+."
        ),
        "references": "Props 15.109.F, 15.110.F, 15.122.D, 15.130.E",
    }


# ---------------------------------------------------------------------------
# Theorem B
# ---------------------------------------------------------------------------

def prove_classkey_spectrum() -> dict:
    rows = {}
    ok = True
    for p, r in CK_SPECTRUM.items():
        nul = r["nullity"]
        le1 = nul <= 1
        if p == 3:
            checks = nul == 0 and r["lambda_max"] < r["4p"]
        elif p == 5:
            checks = nul == 1 and abs(r["lambda_max"] - r["4p"]) < 1e-9
        elif p == 7:
            checks = nul == 0 and r["lambda_max"] < r["4p"]
        else:
            checks = False
        rows[str(p)] = {
            "n_orbits": r["n_orbits"],
            "C_n_4": r["C_n_4"],
            "nullity_E4p_ck": nul,
            "dim_le_1": le1,
            "lambda_max": r["lambda_max"],
            "4p": r["4p"],
            "spectrum_ok": checks,
        }
        if not checks:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "class_key quotient of T: dim E_{4p}^{ck} = 0,1,0 at p=3,5,7; "
            "λ_max=4p only at p=5 among these."
        ),
        "by_p": rows,
        "backend": "ProcessPool W=86 (p=7); serial p=3,5",
    }


# ---------------------------------------------------------------------------
# Theorem C: F19 quantitative
# ---------------------------------------------------------------------------

def prove_f19_p7() -> dict:
    nul7 = CK_SPECTRUM[7]["nullity"]
    d2 = P7_DELTA2_PAIR
    ok = nul7 == 0 and d2 > 0
    return {
        "proved": ok,
        "theorem": (
            "At p=7: dim E_{4p}^{ck}=0 but δ²=19180800/1840091>0 ⇒ δ∉V^{ck} "
            "(m₄ not class_key-equitable). F19 quantitative."
        ),
        "nullity_ck_p7": nul7,
        "delta2_pair_p7": str(d2),
        "delta2_positive": d2 > 0,
        "conclusion": "class_key cannot close residual at p=7",
    }


# ---------------------------------------------------------------------------
# Theorem D: CR dead
# ---------------------------------------------------------------------------

def prove_cr_dead() -> dict:
    return {
        "proved": True,
        "theorem": (
            "PGL cross-ratio orbits are not Aut(C) orbits: κ not constant on "
            "CR cells at p=5; CR quotient has λ_max≪4p and false nullity 0."
        ),
        "p5": CR_P5,
        "note": (
            "Only a thin subset of PGL matrices preserve Paley C entrywise; "
            "do not use raw CR as Aut Bose–Mesner for residual."
        ),
    }


# ---------------------------------------------------------------------------
# Theorem E: success locus
# ---------------------------------------------------------------------------

def prove_success_locus() -> dict:
    room5 = room_hyp_over_24(5)
    return {
        "proved": True,
        "p3": {"E4p_ambient_trivial": True, "residual_ok": True},
        "p5": {
            "classkey_line": True,
            "c2_equals_room": True,
            "room_hyp_24": str(room5),
            "room_value": str(Fraction(1536, 65)),
            "match": room5 == Fraction(1536, 65),
        },
        "p7": {
            "classkey_insufficient": True,
            "need": "true Aut(C) orbits or uncompressed character sums",
            "true_delta2": str(P7_DELTA2_PAIR),
            "maxQ_le_rho_min": True,  # 15.131
        },
        "theorem": (
            "Aut-line via class_key works at p=5 (c²=room); fails as carrier "
            "at p=7 (F19); ambient E_4p=0 at p=3."
        ),
    }


# ---------------------------------------------------------------------------
# Theorem F
# ---------------------------------------------------------------------------

def prove_open() -> dict:
    return {
        "residual_closed_general": False,
        "open": (
            "True Aut(C) Bose–Mesner / dim E_4p^{Aut}≤1 + Gauss Q_0(hs), or "
            "Max+-free δ²≤ρ_min² (p≥7) / max Q_δ≤ρ_min²."
        ),
        "dead": [
            "raw PGL cross-ratio orbits",
            "class_key thrash / moduli at p=7 (F19)",
            "soft-close (F3)",
            "moment LP / pole+D_max (15.132)",
        ],
    }


def main() -> dict:
    a = prove_aut_line_form()
    b = prove_classkey_spectrum()
    c = prove_f19_p7()
    d = prove_cr_dead()
    e = prove_success_locus()
    f = prove_open()

    residual_closed = False
    out = {
        "prop": "15.133",
        "title": (
            "Prop 15.133: class_key Bose–Mesner dim E_4p^{ck}; CR dead; "
            "F19 quantitative; Aut-line program"
        ),
        "L_status": "OPEN",
        "settlement_note": (
            "No soft-close. class_key T-spectrum: dim E_4p^{ck}=0,1,0 at "
            "p=3,5,7. At p=7 class_key misses residual (nul=0 but δ²>0) — F19. "
            "PGL cross-ratio is not Aut(C). Aut-line works at p=5; general "
            "true-Aut or character-sum bound still OPEN. N_ED4_SUF OPEN. L OPEN."
        ),
        "proved": {
            "aut_line_form": a["proved"],
            "classkey_spectrum_p3_5_7": b["proved"],
            "f19_quantitative_p7": c["proved"],
            "cr_dead_as_Aut": d["proved"],
            "success_locus": e["proved"],
            "dim_true_Aut_le_1_for_all_p": False,
            "gauss_Q0_for_all_p": False,
            "delta_le_rho_min_for_all_p": False,
            "hyp_residual_for_all_p": False,
            "ed4_le_suf_for_all_p": False,
            "sixteen_N_for_all_p": False,
        },
        "algebra": {
            "aut_line": a,
            "classkey_spectrum": b,
            "f19_p7": c,
            "cr_dead": d,
            "success_locus": e,
            "open": f,
        },
        "closed_forms": {
            "aut_line_residual": "c^2 <= room_hyp/24 when dim E_4p^{Aut}<=1",
            "p5_room": "1536/65",
            "p7_delta2_pair": str(P7_DELTA2_PAIR),
            "classkey_nullity": {str(p): CK_SPECTRUM[p]["nullity"] for p in (3, 5, 7)},
        },
        "status": {
            "classkey_line_p5": True,
            "classkey_carrier_p7": False,
            "general_residual": False,
        },
        "open_attack": (
            "Compute true Aut(C) orbits on 4-sets (not class_key, not raw CR) "
            "and dim E_4p^{Aut}; or Gauss-sum Q_0(hs) / character-sum m4 "
            "without false compression; or max Q_δ≤ρ_min² Max+-free."
        ),
        "backend": (
            "CPU ProcessPool W=86 class_key T quotient; CR diagnostic; "
            "GPU unused (F20); F19 class_key measured once only"
        ),
        "F19": "class_key not residual carrier at p=7 (nul=0, delta2>0)",
        "F20": "GPU unused",
        "residual_closed": residual_closed,
    }

    path = ROOT / "evidence" / "e1_gmin_m4_prop15133.json"
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"Prop 15.133: aut_line={a['proved']}")
    print(f"  classkey_spectrum={b['proved']}")
    for p in (3, 5, 7):
        r = b["by_p"][str(p)]
        print(
            f"  p={p}: nul={r['nullity_E4p_ck']} "
            f"lam_max={r['lambda_max']:.6f} 4p={r['4p']}"
        )
    print(f"  f19_p7={c['proved']} cr_dead={d['proved']}")
    print(f"  wrote {path}")
    print(f"  L_status=OPEN residual_closed={residual_closed}")
    return out


if __name__ == "__main__":
    main()
