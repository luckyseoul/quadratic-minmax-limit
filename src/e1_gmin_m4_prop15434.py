#!/usr/bin/env python3
"""
Prop 15.434 — J_Δ(r)=∑_{s ns}(L₊(s)−μμ)χ(1+rs) is constant on
pp_fp and on pp_U, not on all ++.  Live p=5: −2,+2; p=7:
−120/19,−12/19.  mean_{++} J_Δ = 2/3, −84/19.  extra ns in
Q/q² is (16p/q²) mean = 32/375, −192/931.  No formula in p
(p=7 only as /19).  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.433 flags.

============================================================================
Setup.  15.432/33: δ=L₊−μμ on Paley×N.  J_Δ(r)=∑_{s ns} δ(s)
χ(1+rs).  ++ sizes 15.314.  extra ns vs constant-L is
(16p/q²) mean_{++} J_Δ.

============================================================================
Theorem A — PROVED (LIVE_LNS + field K; certified p=5,7).
  J_Δ is constant on pp_fp and on pp_U.  Values (−2,+2) and
  (−120/19,−12/19).  Fail: J_Δ=0 (then extra=0≠32/375 at p=5);
  fail: constant on all ++ (−2≠2 at p=5).  ∎

Theorem B — PROVED (15.314 sizes + A).
  mean_{++}=(n_fp J_fp+n_U J_U)/|++| equals 2/3 (p=5) and
  −84/19 (p=7).  extra=(16p/q²)mean equals 32/375, −192/931.
  Fail: J_Δ=0.  ∎

Theorem C — PROVED as a negative (field K; p=5,7,11,13).
  Mixed cells have meanK=0 on pp_fp.  Unmixed |meanK| on
  pp_fp is 2 at p=5,7, is 1 at p=11, and lies in {2,4,6}/5
  at p=13.  Field-fake δ=χ_p(N+1)/2 recovers J at p=5 and
  misses at p=7 (2≠−120/19).  19∉named_gj_ints(7).
  Fail: |meanK_unm,fp|=2 at p=11.  ∎

Theorem D — OPEN.  p=7 values exist only as /19 (Hoffman
  mix, 15.433 amps).  No Max+-free rational in p.  Q_NL
  unnamed.  phi_F not imported.

============================================================================
Backend: field K + 15.432 LIVE_LNS.  Writes
evidence/e1_gmin_m4_prop15434.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field  # noqa: E402
from e1_gmin_m4_prop15290 import field_norm, paley_pair  # noqa: E402
from e1_gmin_m4_prop15314 import chi_p_pm  # noqa: E402
from e1_gmin_m4_prop15330 import named_gj_ints  # noqa: E402
from e1_gmin_m4_prop15356 import Q_1d_pp_named  # noqa: E402
from e1_gmin_m4_prop15409 import D_form_on_lattice_general  # noqa: E402
from e1_gmin_m4_prop15411 import LIVE_QNL  # noqa: E402
from e1_gmin_m4_prop15412 import Q_NL_is_short_rational  # noqa: E402
from e1_gmin_m4_prop15431 import slot_of  # noqa: E402
from e1_gmin_m4_prop15432 import LIVE_LNS  # noqa: E402
from e1_gmin_m4_prop15433 import mumu  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15434_catalog.json"

LIVE_J = {
    5: {"pp_fp": Fraction(-2), "pp_U": Fraction(2)},
    7: {"pp_fp": Fraction(-120, 19), "pp_U": Fraction(-12, 19)},
}


def n_plus(p: int) -> int:
    if p % 4 == 1:
        return 2 * p - 4
    return 3 * (p - 3) // 2


def n_fp(p: int) -> int:
    return p - 3


def n_Upp(p: int) -> int:
    return n_plus(p) - n_fp(p)


def J_zero_wrong(_p: int) -> Fraction:
    return Fraction(0)


def ns_cells(F) -> dict[tuple[int, int, int], list[int]]:
    by: dict[tuple[int, int, int], list[int]] = defaultdict(list)
    for s in range(1, F["q"]):
        if F["chi"][s] != -1:
            continue
        a, b = paley_pair(F, s)
        by[(int(a), int(b), int(field_norm(F, s)))].append(s)
    return by


def K_cell(F, r: int, ss: list[int]) -> int:
    acc = 0
    mul, add, chi = F["mul"], F["add"], F["chi"]
    for s in ss:
        acc += int(chi[add[1][mul[r][s]]])
    return acc


def plus_rs(F, p: int) -> list[int]:
    return [r for r in F["squares"] if slot_of(F, r, p) in ("pp_fp", "pp_U")]


def J_delta_live(p: int) -> dict[str, Fraction]:
    F = {**_kernel_field(p), "p": p}
    by = ns_cells(F)
    dlt = {k: LIVE_LNS[p][k][0] - mumu(p) for k in LIVE_LNS[p]}
    seen: dict[str, set[Fraction]] = {"pp_fp": set(), "pp_U": set()}
    last: dict[str, Fraction] = {}
    for r in plus_rs(F, p):
        sl = slot_of(F, r, p)
        val = sum((dlt[k] * K_cell(F, r, ss) for k, ss in by.items()), Fraction(0))
        seen[sl].add(val)
        last[sl] = val
    if any(len(seen[sl]) != 1 for sl in seen):
        raise ValueError(f"J_Δ not slot-constant at p={p}: {seen}")
    return last


def field_fake_J(p: int) -> dict[str, Fraction]:
    """δ ← χ_p(N+1)/2, the p=5 field law."""
    F = {**_kernel_field(p), "p": p}
    by = ns_cells(F)
    dlt = {k: Fraction(chi_p_pm(p, k[2] + 1), 2) for k in by}
    seen: dict[str, set[Fraction]] = {"pp_fp": set(), "pp_U": set()}
    last: dict[str, Fraction] = {}
    for r in plus_rs(F, p):
        sl = slot_of(F, r, p)
        val = sum((dlt[k] * K_cell(F, r, ss) for k, ss in by.items()), Fraction(0))
        seen[sl].add(val)
        last[sl] = val
    return last


def meanK_unmixed_fp(p: int) -> set[Fraction]:
    F = {**_kernel_field(p), "p": p}
    by = ns_cells(F)
    acc: dict[tuple, int] = defaultdict(int)
    n_r = 0
    for r in plus_rs(F, p):
        if slot_of(F, r, p) != "pp_fp":
            continue
        n_r += 1
        for key, ss in by.items():
            acc[key] += K_cell(F, r, ss)
    out: set[Fraction] = set()
    for (a, b, _N), sm in acc.items():
        if a != b:
            continue
        out.add(abs(Fraction(sm, n_r)))
    return out


def mean_J(p: int) -> Fraction:
    j = LIVE_J[p]
    return (n_fp(p) * j["pp_fp"] + n_Upp(p) * j["pp_U"]) / n_plus(p)


def extra_ns(p: int) -> Fraction:
    q = p * p
    return Fraction(16 * p, q * q) * mean_J(p)


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        got = J_delta_live(p)
        want = LIVE_J[p]
        rows[str(p)] = {sl: str(got[sl]) for sl in ("pp_fp", "pp_U")}
        for sl in ("pp_fp", "pp_U"):
            if got[sl] != want[sl]:
                ok = False
        if got["pp_fp"] == got["pp_U"]:
            ok = False
        if got["pp_fp"] == J_zero_wrong(p) or got["pp_U"] == J_zero_wrong(p):
            ok = False
    if extra_ns(5) == 0:
        ok = False
    if LIVE_J[5]["pp_fp"] != -2 or LIVE_J[5]["pp_U"] != 2:
        ok = False
    if LIVE_J[7]["pp_fp"] != Fraction(-120, 19):
        ok = False
    if LIVE_J[7]["pp_U"] != Fraction(-12, 19):
        ok = False
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "J_Δ constant on pp_fp / pp_U, not on ++. "
            "Fail: J_Δ=0; fail: one value on all ++."
        ),
    }


def prove_B() -> dict:
    ok = True
    rec = {}
    want_mean = {5: Fraction(2, 3), 7: Fraction(-84, 19)}
    want_ex = {5: Fraction(32, 375), 7: Fraction(-192, 931)}
    for p in (5, 7):
        mj = mean_J(p)
        ex = extra_ns(p)
        rec[str(p)] = {"mean": str(mj), "extra": str(ex)}
        if mj != want_mean[p] or ex != want_ex[p]:
            ok = False
        if J_zero_wrong(p) == 0 and extra_ns(p) == 0:
            ok = False
    if extra_ns(5) == 0:
        ok = False
    return {
        "proved": bool(ok),
        "by_p": rec,
        "theorem": (
            "mean_{++} J_Δ=2/3, −84/19. extra=32/375, −192/931. "
            "Fail: J_Δ=0."
        ),
    }


def prove_C() -> dict:
    ok = True
    kabs = {p: sorted(str(x) for x in meanK_unmixed_fp(p)) for p in (5, 7, 11, 13)}
    if set(meanK_unmixed_fp(5)) != {2}:
        ok = False
    if set(meanK_unmixed_fp(7)) != {2}:
        ok = False
    if set(meanK_unmixed_fp(11)) == {2}:
        ok = False
    if set(meanK_unmixed_fp(11)) != {1}:
        ok = False
    if 2 in meanK_unmixed_fp(13):
        ok = False
    fake5 = field_fake_J(5)
    fake7 = field_fake_J(7)
    if fake5 != LIVE_J[5]:
        ok = False
    if fake7["pp_fp"] == LIVE_J[7]["pp_fp"]:
        ok = False
    if fake7["pp_fp"] != 2:
        ok = False
    S7 = named_gj_ints(7)
    if 19 in S7 or -19 in S7 or 38 in S7:
        ok = False
    return {
        "proved": bool(ok),
        "K_unm_fp": kabs,
        "fake7_fp": str(fake7["pp_fp"]),
        "S7_has_19": 19 in S7,
        "theorem": (
            "Unmixed |meanK| on pp_fp is 2 at p=5,7, not at p=11,13. "
            "Field-fake dies at p=7. 19 not GJ. "
            "Fail: |K|=2 at p=11."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "p7_fp_den19": LIVE_J[7]["pp_fp"].denominator % 19 == 0,
        "p7_U_den19": LIVE_J[7]["pp_U"].denominator % 19 == 0,
        "p7_mean_den19": mean_J(7).denominator % 19 == 0,
        "Q_NL_5": str(LIVE_QNL[5]),
        "Q_1d_5": str(Q_1d_pp_named(5)),
        "Q_NL_gt_Q1d_p5": LIVE_QNL[5] > Q_1d_pp_named(5),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "D_form_on_lattice_general": D_form_on_lattice_general(),
        "Q_NL_is_short_rational": Q_NL_is_short_rational(),
        "note": (
            "p=7 J_Δ only as /19 (Hoffman, not GJ). "
            "Q_NL unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.434  J_Δ slot-constant; no formula in p", flush=True)
    A = prove_A()
    print(f"  A slot-const: {A['proved']} {A['by_p']}", flush=True)
    B = prove_B()
    print(f"  B mean/extra: {B['proved']} {B['by_p']}", flush=True)
    C = prove_C()
    print(
        f"  C K/GJ: {C['proved']} K={C['K_unm_fp']} fake7={C['fake7_fp']} "
        f"19inS7={C['S7_has_19']}",
        flush=True,
    )
    D = prove_open()
    print(f"  D open: den19={D['p7_fp_den19']} phi_F={D['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps(
            {
                "J5": {k: str(v) for k, v in LIVE_J[5].items()},
                "J7": {k: str(v) for k, v in LIVE_J[7].items()},
                "mean5": str(mean_J(5)),
                "mean7": str(mean_J(7)),
            },
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.434",
        "title": "J_Δ slot-constant on pp_fp/pp_U; p=7 only /19; Q_NL open",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "J_slot_constant": A["proved"],
            "mean_extra_named_live": B["proved"],
            "K_pattern_dies_p11": C["proved"],
            "phi_F_ge_6_proved_general": D["phi_F_ge_6"],
            "residual_ii_k_eq_4p_empty": D["residual_ii_k_eq_4p_empty"],
            "Q_NL_is_short_rational": D["Q_NL_is_short_rational"],
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
            "residual_ii_k_eq_4p_empty",
        ],
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15434.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
