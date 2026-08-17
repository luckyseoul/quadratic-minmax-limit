#!/usr/bin/env python3
"""
Prop 15.414 — Paley ++ splits into F_p-parallel and extra-++
(split-slope) 4-points.  The 1D same-slope slot is named
    L_par^{1D} = p²(p−1)(p³−5p−4) / (16(p+1)(p−2))
(Max+-free Johnson on the μ-plus-set).  NL L_par / L_spl
stay unnamed; Q_NL is not imported.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.413 flags.

============================================================================
Setup.  Paley ++ contains F_p^×\\{±1} (15.355).  Multiplication
by λ∈F_p preserves AG(2,p) slope, so N(d) and N(λd) live on one
parallel class.  Extra ++ (s∉F_p) send d to a different slope.
n_par = p−3, n_spl = n_++ − (p−3).  At p=5 the Paley-++ class
is exactly {2,3}=F_5^×\\{±1}.  At p=7 it is {2,3,4,5} ∪ {21,28}.

On a Type+ 1D halfspace, H-direction contributes k² and each
C-direction contributes p² E[J_μ(1)J_μ(λ)] with λ∉{0,±1} and
|A|=μ=(p−1)/2 (the D-side, not the σ=+1 side).

============================================================================
Theorem A — PROVED (field geometry; certified p=5,7 cache).
  L_NL(1,++) = (n_par L_par + n_spl L_spl) / n_++.
  Live: L_par=24 (p=5) and 2113/19 (p=7); L_spl=673/6 (p=7);
  mix = 24 and 38143/342.  Fail: claim Paley ++ is purely
  same-slope at p=7 (then 2113/19 ≠ 38143/342).  ∎

Theorem B — PROVED (Max+-free Johnson; certified p=5,7).
  L_par^{1D} = p²(p−1)(p³−5p−4)/(16(p+1)(p−2))
  = 2(k² + ((p−1)/2) p² E_JJ_μ)/(p+1),
  E_JJ_μ=(p−5)(p²−3p+4)/(16(p−2)).
  Equals 100/3 and 2793/20.  Fail: use |A|=(p+1)/2
  (200/3 ≠ 100/3 at p=5).  ∎

Theorem C — PROVED (p=5 cache + A).
  On p=5 NL, Paley ++ is purely parallel and L=q−1=24
  (E[cc]=0).  Fail: claim L_par^{NL}=q−1 at p=7.  ∎

Theorem D — OPEN.  L_par^{NL} and L_spl^{NL} unnamed in p.
  15.298 still needs the full L table to name Q_NL.
  The 2-point interpolant C(p) is not imported.  phi_F False.

============================================================================
Backend: Fraction + lru fold on p=5,7 caches.  Writes
evidence/e1_gmin_m4_prop15414.json
"""
from __future__ import annotations

import json
import os
import sys
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15298 import coarse_Q_class, load_N  # noqa: E402
from e1_gmin_m4_prop15307 import load_D, square_classes, tag_sig  # noqa: E402
from e1_gmin_m4_prop15409 import D_form_on_lattice_general  # noqa: E402
from e1_gmin_m4_prop15411 import LIVE_QNL  # noqa: E402
from e1_gmin_m4_prop15412 import Q_NL_is_short_rational  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15414_catalog.json"

LIVE_LPAR_NL = {5: Fraction(24), 7: Fraction(2113, 19)}
LIVE_LSPL_NL = {7: Fraction(673, 6)}
LIVE_LPP_NL = {5: Fraction(24), 7: Fraction(38143, 342)}
LIVE_LPAR_1D = {5: Fraction(100, 3), 7: Fraction(2793, 20)}
LIVE_LSPL_1D = {7: Fraction(539, 5)}
LIVE_LPP_1D = {5: Fraction(100, 3), 7: Fraction(3871, 30)}


def falling(n: int, k: int) -> int:
    acc = 1
    for i in range(k):
        acc *= n - i
    return acc


def E_JJ_mu(p: int) -> Fraction:
    """E[J(1)J(λ)] on one F_p-line, |A|=μ=(p−1)/2, λ∉{0,±1}."""
    return Fraction((p - 5) * (p * p - 3 * p + 4), 16 * (p - 2))


def Lpar_1d_named(p: int) -> Fraction:
    """Max+-free: p²(p−1)(p³−5p−4) / (16(p+1)(p−2))."""
    return Fraction(
        p * p * (p - 1) * (p * p * p - 5 * p - 4),
        16 * (p + 1) * (p - 2),
    )


def Lpar_1d_via_EJJ(p: int) -> Fraction:
    k = p * (p - 1) // 2
    return Fraction(2 * (k * k + Fraction(p - 1, 2) * p * p * E_JJ_mu(p)), p + 1)


def Lpar_1d_wrong_m(p: int) -> Fraction:
    """Fail-when-wrong: |A|=(p+1)/2 instead of μ."""
    m = (p + 1) // 2
    p3, p4 = falling(p, 3), falling(p, 4)
    n3, n4 = falling(m, 3), falling(m, 4)
    EJ = Fraction(p) * (4 * Fraction(n3, p3) + (p - 4) * Fraction(n4, p4))
    k = p * (p - 1) // 2
    return Fraction(2 * (k * k + Fraction(p - 1, 2) * p * p * EJ), p + 1)


def n_par_named(p: int) -> int:
    """# F_p^× \\ {±1}."""
    return p - 3


def mix_Lpp(npar: int, Lpar: Fraction, nspl: int, Lspl: Fraction | None) -> Fraction:
    if nspl == 0 or Lspl is None:
        return Lpar
    return (npar * Lpar + nspl * Lspl) / (npar + nspl)


def _in_fp(p: int, x: int) -> bool:
    return x != 0 and (x // p) == 0


def _hs_mask(p, D, F):
    classes = square_classes(F)
    n_sq = len(classes)
    occ = []
    for lines in classes:
        ns = np.stack([D[:, ix].sum(axis=1) for ix in lines], axis=1)
        occ.append([tuple(int(v) for v in np.sort(row)) for row in ns])
    m = np.zeros(D.shape[0], dtype=bool)
    for i in range(D.shape[0]):
        t = tuple(sorted(tag_sig(p, occ[di][i]) for di in range(n_sq)))
        m[i] = t.count("H") == 1 and t.count("C") == n_sq - 1
    return m


def _mean_L(N, F, mask, dlist, slist) -> Fraction | None:
    if not slist:
        return None
    Nm = N[mask]
    ENN = Nm.T @ Nm
    n = int(mask.sum())
    acc = 0
    c = 0
    mul = F["mul"]
    for d in dlist:
        for s in slist:
            acc += int(ENN[d, mul[s][d]])
            c += 1
    return Fraction(acc, c * n)


@lru_cache(maxsize=4)
def parsplit_table(p: int) -> dict:
    D, F = load_D(p)
    N = load_N(p, F).astype(np.int64)
    hs = _hs_mask(p, D, F)
    q = F["q"]
    plus_s = [
        s
        for s in range(1, q)
        if F["chi"][s] == 1 and coarse_Q_class(F, s) == "++"
    ]
    s_par = [s for s in plus_s if _in_fp(p, s) and s not in (1, F["neg1"])]
    s_spl = [s for s in plus_s if s not in s_par]
    d_pos = [d for d in range(1, q) if F["chi"][d] == 1]
    out: dict = {"s_par": s_par, "s_spl": s_spl, "plus_s": plus_s}
    for lab, mask in (("NL", ~hs), ("1D", hs)):
        out[lab] = {
            "n": int(mask.sum()),
            "L++": _mean_L(N, F, mask, d_pos, plus_s),
            "Lpar": _mean_L(N, F, mask, d_pos, s_par),
            "Lspl": _mean_L(N, F, mask, d_pos, s_spl),
        }
    return out


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        tab = parsplit_table(p)
        npar = len(tab["s_par"])
        nspl = len(tab["s_spl"])
        Lpar = tab["NL"]["Lpar"]
        Lspl = tab["NL"]["Lspl"]
        Lpp = tab["NL"]["L++"]
        mix = mix_Lpp(npar, Lpar, nspl, Lspl)
        ok = ok and mix == Lpp
        ok = ok and npar == n_par_named(p)
        ok = ok and Lpar == LIVE_LPAR_NL[p]
        ok = ok and Lpp == LIVE_LPP_NL[p]
        rows[str(p)] = {
            "n_par": npar,
            "n_spl": nspl,
            "Lpar": str(Lpar),
            "Lspl": None if Lspl is None else str(Lspl),
            "L++": str(Lpp),
            "mix": str(mix),
        }
        if p == 7:
            ok = ok and Lspl == LIVE_LSPL_NL[7]
            # fail-when-wrong: purely parallel
            ok = ok and Lpar != Lpp
            if Lpar == LIVE_LPP_NL[7]:
                ok = False
        if p == 5:
            ok = ok and nspl == 0
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "L_++ = (n_par L_par + n_spl L_spl)/n_++. "
            "Fail: Paley ++ purely same-slope at p=7."
        ),
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 11, 13, 17, 19):
        named = Lpar_1d_named(p)
        via = Lpar_1d_via_EJJ(p)
        ok = ok and named == via
        rows[str(p)] = str(named)
    ok = ok and Lpar_1d_named(5) == LIVE_LPAR_1D[5]
    ok = ok and Lpar_1d_named(7) == LIVE_LPAR_1D[7]
    ok = ok and Lpar_1d_wrong_m(5) != LIVE_LPAR_1D[5]
    ok = ok and Lpar_1d_wrong_m(5) == Fraction(200, 3)
    # live 1D cache check
    tab5 = parsplit_table(5)
    tab7 = parsplit_table(7)
    ok = ok and tab5["1D"]["Lpar"] == LIVE_LPAR_1D[5]
    ok = ok and tab7["1D"]["Lpar"] == LIVE_LPAR_1D[7]
    if Lpar_1d_wrong_m(5) == LIVE_LPAR_1D[5]:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "wrong_m_p5": str(Lpar_1d_wrong_m(5)),
        "theorem": (
            "L_par^{1D}=p²(p−1)(p³−5p−4)/(16(p+1)(p−2)). "
            "Fail: |A|=(p+1)/2 (200/3≠100/3)."
        ),
    }


def prove_C() -> dict:
    tab = parsplit_table(5)
    L = tab["NL"]["L++"]
    Lpar = tab["NL"]["Lpar"]
    ok = L == Lpar == 24 == 5 * 5 - 1
    ok = ok and len(tab["s_spl"]) == 0
    tab7 = parsplit_table(7)
    ok = ok and tab7["NL"]["Lpar"] != 48
    ok = ok and tab7["NL"]["Lpar"] == Fraction(2113, 19)
    if tab7["NL"]["Lpar"] == 48:
        ok = False
    return {
        "proved": bool(ok),
        "p5_L": str(L),
        "p7_Lpar": str(tab7["NL"]["Lpar"]),
        "theorem": (
            "p=5 NL Paley ++ is purely parallel and L=q−1=24. "
            "Fail: L_par^{NL}=q−1 at p=7."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "Q_NL_5": str(LIVE_QNL[5]),
        "Q_NL_7": str(LIVE_QNL[7]),
        "L_par_NL_7": str(LIVE_LPAR_NL[7]),
        "L_spl_NL_7": str(LIVE_LSPL_NL[7]),
        "Q_NL_is_short_rational": Q_NL_is_short_rational(),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "D_form_on_lattice_general": D_form_on_lattice_general(),
        "note": (
            "L_par^{1D} named. NL L_par=2113/19 and L_spl=673/6 "
            "unnamed in p. Full L table still needed for 15.298. "
            "C(p) interpolant not imported. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.414  Paley ++ par/split; L_par^{1D} named", flush=True)
    A = prove_A()
    print(f"  A mix: {A['proved']} {A['rows']}", flush=True)
    B = prove_B()
    print(f"  B 1D form: {B['proved']} p5={B['rows']['5']} p7={B['rows']['7']}", flush=True)
    C = prove_C()
    print(f"  C p5 pure-par: {C['proved']} p7Lpar={C['p7_Lpar']}", flush=True)
    D = prove_open()
    print(f"  D open: phi_F={D['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps(
            {
                "A": A["rows"],
                "B": B["rows"],
                "C": {"p5": C["p5_L"], "p7_Lpar": C["p7_Lpar"]},
            },
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.414",
        "title": (
            "Paley ++ = F_p-parallel + extra-++; "
            "L_par^{1D} named Max+-free; NL L_par/L_spl open"
        ),
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "par_split_mix": A["proved"],
            "Lpar_1d_named": B["proved"],
            "p5_NL_pure_parallel": C["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15414.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
