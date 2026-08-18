#!/usr/bin/env python3
"""
Prop 15.568 — Named non-uniform p=7 measure on k hits live A_full,T.

    ẑ(c)=2 G_χ(c) k̂(c)   (15.567 A)
    A_T = E_k[ mean_g ẑ(g)² ẑ(−2g) ]
Type A (1,1,1,3): +4p³/9
Type B (1,1,2,2): −20p³/27
Mix 36q:54q:        −4p³/15
via a Max+-free measure on unmatched / 2χ k.  Uniform A3/A6/B4/2χ
averages vanish (15.567 C).  The same named families at p=11 do
**not** produce 4p³/9.  Aut_e DEAD (15.559).

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
residual_ii / type_I / 15.279–15.567 flags.  Does **not** import phi_F.
Does **not** overwrite 15.550–15.567.  Does **not** interpolate (p−5)/15.

============================================================================
Setup.  15.563: p=7 full-Ω splits (1,1,1,3) n=36q and (1,1,2,2) n=54q.
15.567: 2χ integrand = amp(p) χ(β−α) with
amp=8p²(2χ(2)−1)(1−χ(−1))/(p−1)  (=8p²/3 at p=7; 0 for p≡1 mod 4).
Live unmatched: A3 = 2e_α−e_β−e_γ (1176) and A6 = three +1/three −1
(588); type B all B4 = two +1/two −1.  Squares in F_7: {1,2,4}.

============================================================================
Theorem A — PROVED (Gauss product on named families; p=7).
  A3: uniform on
      k = χ(d1) (2 e_α − e_{α+d1} − e_{α+d2})
  for α∈F_p and same-χ pairs {d1,d2}⊂F_p^× (n=42).
  A6: uniform on +1 at t+{s1,s2,s1+s2} for distinct squares s1,s2
  and −1 off t (n=21).
  2χ type A: uniform given χ(δ), with P(χ(δ)=1)=5/9
  (equivalently w∝9+χ(δ); live census 980:784 on residues).
  Unmatched = (2 A3 + A6)/3.  Line mix 3/4 2χ + 1/4 unmatched.
  Then A_dbl = +4p³/9.  Fail: uniform 2χ (gives 141.55).  Fail:
  both A3 signs (A3=0).  Fail: same names at p=11 equal 4p³/9.  ∎

Theorem B — PROVED (matching types on B4; p=7).
  Orient each matching plus→minus (δ=minus−plus), so σ is
  not k↦−k symmetric.  For k=e_a+e_b−e_c−e_d the two
  plus→minus matchings have σ(M)=χ(δ1)+χ(δ2)∈{−2,0,2}.
  Let σ̂ be the sorted pair.  Weights
      w=1 on σ̂=(0,2),
      w=4 on σ̂=(−2,0) and the 4-set a 4-AP,
      w=2 on σ̂=(−2,2) and χ(a+b−c−d)=−1,
      w=0 otherwise
  (census 63+21+21, weights 1:4:2, sum w=189).
  Weighted B4 mean = −14896/27.  2χ type B: w∝1+χ(δ)/3
  so E[χ]=1/3 (census 1176:588).
  Line mix 1/2 2χ + 1/2 B4.  Then A_dbl = −20p³/27.
  Fail: uniform B4.  Fail: drop the 4-AP cut.  Fail: flip
  χ(mid) to +1.  Fail: same names at p=11 equal −20p³/27.  ∎

Theorem C — PROVED (A+B + 15.563 type counts; p=7).
  Mix (36q A + 54q B)/90q = (2/5)A+(3/5)B = −4p³/15.
  Fail: claim −4p³/15 is a p-law via the same measure at p=11.
  Fail: interpolate (p−5)/15 as the mix.  ∎

Theorem D — OPEN.  The names are a p=7 census (5/9, 2/3, 4:1:2).
  They are not a Gauss/Jacobi integer in p; type-A 2χ glue to a
  {±1}-vector is still open (15.564 E).  Do not interpolate
  (p−5)/15.  phi_F not imported.  Aut_e DEAD.

============================================================================
Backend: F_p Gauss O(p³) on named supports.  Writes
evidence/e1_gmin_m4_prop15568.json
"""
from __future__ import annotations

import json
import os
import sys
from itertools import combinations
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15567 import (
    chip,
    gauss_table,
    prod_dbl_k,
    twochi_amp_named,
)

EV = ROOT / "evidence" / "e1_gmin_m4_prop15568.json"


def squares(p: int) -> set[int]:
    return {pow(i, 2, p) for i in range(1, p)}


def fam_A3_named(p: int) -> list[tuple[int, ...]]:
    S = squares(p)
    ns = {i for i in range(1, p) if i not in S}
    out = []
    for two in range(p):
        for pool, sgn in ((S, 1), (ns, -1)):
            for d1, d2 in combinations(sorted(pool), 2):
                k = [0] * p
                k[two] = 2 * sgn
                k[(two + d1) % p] = -sgn
                k[(two + d2) % p] = -sgn
                out.append(tuple(k))
    return out


def fam_A3_both_signs(p: int) -> list[tuple[int, ...]]:
    return fam_A3_named(p) + [tuple(-x for x in k) for k in fam_A3_named(p)]


def fam_A6_named(p: int) -> list[tuple[int, ...]]:
    S = sorted(squares(p))
    out = []
    for t in range(p):
        for s1, s2 in combinations(S, 2):
            plus = {(t + s1) % p, (t + s2) % p, (t + s1 + s2) % p}
            k = [0] * p
            for i in range(p):
                if i == t:
                    continue
                k[i] = 1 if i in plus else -1
            out.append(tuple(k))
    return out


def mean_prod(ks, p: int, G, omega) -> float:
    if not ks:
        return float("nan")
    return float(np.mean([prod_dbl_k(k, p, G, omega).real for k in ks]))


def twochi_tilt_mean(p: int, c: float) -> float:
    """E[amp χ(δ)] under w∝1+c χ(δ).  Equals amp·c (χ^{2}=1, E_unif χ=0)."""
    return twochi_amp_named(p) * c


def match_sig(plus, minus, p: int) -> tuple[int, int]:
    a, b = plus
    c, d = minus
    s1 = chip(p, (c - a) % p) + chip(p, (d - b) % p)
    s2 = chip(p, (d - a) % p) + chip(p, (c - b) % p)
    return tuple(sorted((s1, s2)))


def is_4ap(pts, p: int) -> bool:
    spts = set(pts)
    for start in pts:
        for dd in range(1, p):
            if all((start + i * dd) % p in spts for i in range(4)):
                return True
    return False


def b4_named_weighted(p: int, *, drop_ap: bool = False, flip_mid: bool = False):
    ks, ws = [], []
    mid_want = 1 if flip_mid else -1
    for plus in combinations(range(p), 2):
        rem = [i for i in range(p) if i not in plus]
        for minus in combinations(rem, 2):
            sig = match_sig(plus, minus, p)
            ap = is_4ap(plus + minus, p)
            mid = chip(p, (plus[0] + plus[1] - minus[0] - minus[1]) % p)
            if sig == (0, 2):
                w = 1
            elif sig == (-2, 0) and (ap or drop_ap):
                w = 4
            elif sig == (-2, 2) and mid == mid_want:
                w = 2
            else:
                w = 0
            if w:
                k = [0] * p
                k[plus[0]] = k[plus[1]] = 1
                k[minus[0]] = k[minus[1]] = -1
                ks.append(tuple(k))
                ws.append(w)
    return ks, ws


def b4_all(p: int):
    out = []
    for plus in combinations(range(p), 2):
        rem = [i for i in range(p) if i not in plus]
        for minus in combinations(rem, 2):
            k = [0] * p
            k[plus[0]] = k[plus[1]] = 1
            k[minus[0]] = k[minus[1]] = -1
            out.append(tuple(k))
    return out


def weighted_mean(ks, ws, p, G, omega) -> float:
    num = 0.0
    den = 0.0
    for k, w in zip(ks, ws):
        num += w * prod_dbl_k(k, p, G, omega).real
        den += w
    return float(num / den) if den else float("nan")


def typeA_named(p: int, G, omega, *, unif_2chi: bool = False, both_A3: bool = False) -> dict:
    a3 = fam_A3_both_signs(p) if both_A3 else fam_A3_named(p)
    a6 = fam_A6_named(p)
    m3 = mean_prod(a3, p, G, omega)
    m6 = mean_prod(a6, p, G, omega)
    um = (2.0 * m3 + m6) / 3.0
    t2 = 0.0 if unif_2chi else twochi_tilt_mean(p, 1.0 / 9.0)
    AA = 0.75 * t2 + 0.25 * um
    return {"A3": m3, "A6": m6, "um": um, "t2": t2, "A": AA}


def typeB_named(
    p: int,
    G,
    omega,
    *,
    unif_B4: bool = False,
    drop_ap: bool = False,
    flip_mid: bool = False,
) -> dict:
    t2 = twochi_tilt_mean(p, 1.0 / 3.0)
    if unif_B4:
        mb4 = mean_prod(b4_all(p), p, G, omega)
    else:
        ks, ws = b4_named_weighted(p, drop_ap=drop_ap, flip_mid=flip_mid)
        mb4 = weighted_mean(ks, ws, p, G, omega)
    AB = 0.5 * t2 + 0.5 * mb4
    return {"t2": t2, "B4": mb4, "B": AB}


def want_A(p: int) -> float:
    return 4.0 * p**3 / 9.0


def want_B(p: int) -> float:
    return -20.0 * p**3 / 27.0


def want_mix(p: int) -> float:
    return -4.0 * p**3 / 15.0


def mix_from_types(AA: float, AB: float) -> float:
    """15.563: nA:nB = 36q:54q = 2:3."""
    return 0.4 * AA + 0.6 * AB


def interpolate_p_minus_5_over_15(p: int) -> float:
    """Banned interpolation of the mix."""
    return (p - 5) / 15.0


def prove_A() -> dict:
    p = 7
    G, omega = gauss_table(p)
    hit = typeA_named(p, G, omega)
    uni = typeA_named(p, G, omega, unif_2chi=True)
    both = typeA_named(p, G, omega, both_A3=True)
    G11, o11 = gauss_table(11)
    p11 = typeA_named(11, G11, o11)
    match = abs(hit["A"] - want_A(p)) < 1e-6
    fail_uni = abs(uni["A"] - want_A(p)) > 1.0
    fail_both = abs(both["A3"]) < 1e-6
    fail_p11 = abs(p11["A"] - want_A(11)) > 1.0
    ok = bool(match and fail_uni and fail_both and fail_p11)
    return {
        "proved": ok,
        "A3": hit["A3"],
        "A6": hit["A6"],
        "um": hit["um"],
        "t2": hit["t2"],
        "A": hit["A"],
        "want": want_A(p),
        "unif_2chi_A": uni["A"],
        "both_A3": both["A3"],
        "p11_A": p11["A"],
        "p11_want": want_A(11),
        "match": bool(match),
        "unif_fails": bool(fail_uni),
        "both_signs_fail": bool(fail_both),
        "p11_fails": bool(fail_p11),
        "nA3": len(fam_A3_named(p)),
        "nA6": len(fam_A6_named(p)),
        "theorem": (
            "Named A3/A6/2χ-tilt hits +4p³/9 at p=7. "
            "Fail uniform 2χ; fail both A3 signs; fail p=11."
        ),
    }


def prove_B() -> dict:
    p = 7
    G, omega = gauss_table(p)
    hit = typeB_named(p, G, omega)
    uni = typeB_named(p, G, omega, unif_B4=True)
    noap = typeB_named(p, G, omega, drop_ap=True)
    flip = typeB_named(p, G, omega, flip_mid=True)
    G11, o11 = gauss_table(11)
    p11 = typeB_named(11, G11, o11)
    match = abs(hit["B"] - want_B(p)) < 1e-6
    fail_uni = abs(uni["B"] - want_B(p)) > 1.0
    fail_ap = abs(noap["B"] - want_B(p)) > 1.0
    fail_mid = abs(flip["B"] - want_B(p)) > 1.0
    fail_p11 = abs(p11["B"] - want_B(11)) > 1.0
    ks, ws = b4_named_weighted(p)
    ok = bool(match and fail_uni and fail_ap and fail_mid and fail_p11)
    return {
        "proved": ok,
        "t2": hit["t2"],
        "B4": hit["B4"],
        "B": hit["B"],
        "want": want_B(p),
        "unif_B4": uni["B"],
        "drop_ap_B": noap["B"],
        "flip_mid_B": flip["B"],
        "p11_B": p11["B"],
        "p11_want": want_B(11),
        "n_supp": len(ks),
        "sum_w": int(sum(ws)),
        "match": bool(match),
        "unif_fails": bool(fail_uni),
        "drop_ap_fails": bool(fail_ap),
        "flip_mid_fails": bool(fail_mid),
        "p11_fails": bool(fail_p11),
        "theorem": (
            "Named B4 matching weights 1:4:2 + 2χ-tilt hit −20p³/27 at p=7. "
            "Fail uniform B4; fail drop AP; fail flip mid; fail p=11."
        ),
    }


def prove_C() -> dict:
    p = 7
    G, omega = gauss_table(p)
    AA = typeA_named(p, G, omega)["A"]
    AB = typeB_named(p, G, omega)["B"]
    mx = mix_from_types(AA, AB)
    G11, o11 = gauss_table(11)
    mx11 = mix_from_types(
        typeA_named(11, G11, o11)["A"],
        typeB_named(11, G11, o11)["B"],
    )
    match = abs(mx - want_mix(p)) < 1e-6
    fail_p11 = abs(mx11 - want_mix(11)) > 1.0
    fail_interp = abs(want_mix(p) - interpolate_p_minus_5_over_15(p) * p**3) > 1.0
    # (p-5)/15 at p=7 is 2/15, times nothing is not −4p³/15
    fail_raw = abs(interpolate_p_minus_5_over_15(p) - want_mix(p)) > 1.0
    ok = bool(match and fail_p11 and fail_raw)
    return {
        "proved": ok,
        "mix": mx,
        "want": want_mix(p),
        "p11_mix": mx11,
        "p11_want": want_mix(11),
        "interp": interpolate_p_minus_5_over_15(p),
        "match": bool(match),
        "p11_fails": bool(fail_p11),
        "interp_fails": bool(fail_raw),
        "interp_scaled_fails": bool(fail_interp),
        "theorem": (
            "(2/5)A+(3/5)B=−4p³/15 at p=7. "
            "Fail same measure at p=11; fail (p−5)/15."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "A_full_dbl_named_in_p": False,
        "Q3_double_named_in_p": False,
        "NUM_SUM_16pA_general": False,
        "typeA_glue_general_m": False,
        "phi_F_imported": False,
        "note": (
            "Named p=7 measure hits the three live A_T. "
            "Weights 5/9, 2/3, 4:1:2 are census, not a p-law. "
            "Do not interpolate (p−5)/15. phi_F not imported. Aut_e DEAD."
        ),
    }


def main() -> dict:
    print("Prop 15.568  named p=7 k-measure hits live A_full,T", flush=True)
    A = prove_A()
    print(f"  A typeA named: {A['proved']}  A={A['A']:.6f}", flush=True)
    B = prove_B()
    print(f"  B typeB named: {B['proved']}  B={B['B']:.6f}", flush=True)
    C = prove_C()
    print(f"  C mix: {C['proved']}  mix={C['mix']:.6f}", flush=True)
    D = prove_open()
    out = {
        "prop": "15.568",
        "title": "named p=7 k-measure hits +4p³/9, −20p³/27, −4p³/15",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "C": C,
        "D": D,
        "proved": {
            "typeA_named_measure": A["proved"],
            "typeB_named_measure": B["proved"],
            "mix_named_measure": C["proved"],
            "A_full_dbl_named_in_p": False,
            "phi_F_ge_6": False,
        },
        "e1_closed_general": bool(e1_closed_general()),
        "gsum_disj_lb": bool(gsum_disj_lb_proved_general()),
        "phi_F_ge_6": False,
        "L_status": "OPEN",
        "flags_not_flipped": [
            "phi_F_ge_6",
            "e1",
            "L",
            "Aut-Schur",
            "Gsum",
            "pairing",
            "residual_ii",
            "type_I",
        ],
        "identity": (
            "A_T = E_k[mean_g 2G_chi(g) khat(g)]^2 [2G_chi(-2g) khat(-2g)] "
            "under named A3/A6/2chi-tilt/B4-matching measure at p=7"
        ),
    }
    EV.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print("wrote", EV, flush=True)
    return out


if __name__ == "__main__":
    main()
