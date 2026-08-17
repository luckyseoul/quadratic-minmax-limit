#!/usr/bin/env python3
"""
Prop 15.448 — Every Type+ 1D slot of L is named in p.  The
split is (χ(s), χ(s−1), χ(s+1), 1_{s∈F_p}), not Paley-Tor/N*.
Ensemble L_τ / Q_τ still carry |H+|.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.447 flags.

============================================================================
Setup.  Type+ 1D: y=σ∘L, ker L a square F_p-line, |σ^{-1}(+1)|=(p+1)/2
(15.279 R, 15.292).  Affine N(d)=k if L(d)=0 else p J_μ(L(d)), with
|B|=μ=(p−1)/2 the minus set and J the Johnson 1-point on F_p
(15.414–15.415).  L_1D(s)=E[N(d)N(sd)] over χ(d)=+1 and Type+ lifts.

Johnson (named):
    EJ   = (p−3)/4,
    EJ2  = (p−3)(p²−4p+7)/(16(p−2))   = E[J(1)²]=E[J(1)J(−1)],
    EJJ  = (p−5)(p²−3p+4)/(16(p−2))   = E[J(1)J(λ)], λ∉{0,±1}.

χ|_{F_p^×}≡+1 because (q−1)/2=(p²−1)/2 is a multiple of p−1.
A Type+ line is therefore all F_q-squares.

============================================================================
Theorem A — PROVED (field character; all odd p).
  χ|_{F_p^×}≡+1.  Every Type+ kernel line is contained in the
  F_q-squares, so #{d∈ker^× : χ(d)=+1}=p−1.  Fail: claim
  χ|_{F_p}=χ_p (then χ(3)=−1 at p=5).  ∎

Theorem B — PROVED (A + line geometry; all odd p).
  Average over Type+ forms and square d, n_sq=(q−1)/2.  Fractions:
      n00 = 2/(p+1) on s∈F_p, else 0
      n01 = 2/(p+1) on s∉F_p, else 0
      n10 = 2/(p+1) on squares s∉F_p, else 0
      npm = 2/(p+1)·(1_{χ(s−1)=+1}+1_{χ(s+1)=+1}) on s∉F_p
          = (p−1)/(p+1) on s=±1, else 0 on F_p
      ng  = 1 − n00 − n01 − n10 − npm.
  Fail: put n10=n01 on nonsquares (kills ns_mm at p=5).  ∎

Theorem C — PROVED (B + Johnson; certified p=5,7 lifts).
  L_1D is the Johnson combination of those weights, hence
      L_pm1     = p²(p−1)(p²−5)/(16(p−2)),
      L_Sub     = L_par^{1D}   (15.414),
      L_sq(1,1) = L_spl^{1D}   (15.415; now all odd p, empty at p=5),
      L_sq mix  = p²(p³−4p²+p+10)/(16(p−2)),
      L_sq(−,−) = p²(p−3)(p²−p−4)/(16(p−2)),
      L_ns(1,1) = p²(p−1)(p²−5p+8)/(16(p−2)),
      L_ns mix  = p²(p−1)(p−3)/16   (=μ₊μ₋),
      L_ns(−,−) = p²(p−1)²(p−4)/(16(p−2)).
  Fail: invert p−2 on L_pm1 (125/4≠125/3 at p=5);
  fail: claim L_spl^{1D}(5)=25 is a live 1D value.  ∎

Theorem D — OPEN.  Ensemble L_τ = w L_1D + (1−w) L_NL still
  has w=n_1d/|H+| and L_NL unnamed.  Q_τ unnamed.  phi_F False.

============================================================================
Backend: Fraction + Type+ lifts at p=5,7.  Writes
evidence/e1_gmin_m4_prop15448.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import defaultdict
from fractions import Fraction
from itertools import combinations
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import (  # noqa: E402
    _kernel_field,
    _linear_forms,
    is_type_plus_form,
    lift_sigma_vector,
)
from e1_gmin_m4_prop15292 import n_1d  # noqa: E402
from e1_gmin_m4_prop15298 import mu_minus, mu_plus  # noqa: E402
from e1_gmin_m4_prop15414 import E_JJ_mu, Lpar_1d_named  # noqa: E402
from e1_gmin_m4_prop15415 import Lspl_1d_named  # noqa: E402
from e1_gmin_m4_prop15445 import pow_f  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15448_catalog.json"
CERT_LIFTS = (5, 7)
CERT_WT = (5, 7, 11, 13)


def EJ_mu(p: int) -> Fraction:
    return Fraction(p - 3, 4)


def EJ2_mu(p: int) -> Fraction:
    return Fraction((p - 3) * (p * p - 4 * p + 7), 16 * (p - 2))


def Lpm1_1d_named(p: int) -> Fraction:
    return Fraction(p * p * (p - 1) * (p * p - 5), 16 * (p - 2))


def Lpm1_1d_inverted(p: int) -> Fraction:
    """Fail-when-wrong: p−1 in the denominator."""
    return Fraction(p * p * (p - 1) * (p * p - 5), 16 * (p - 1))


def Lsq_11_named(p: int) -> Fraction:
    """Same polynomial as 15.415 L_spl, now for every odd p."""
    return Lspl_1d_named(p)


def Lsq_mix_named(p: int) -> Fraction:
    return Fraction(p * p * (p**3 - 4 * p * p + p + 10), 16 * (p - 2))


def Lsq_mm_named(p: int) -> Fraction:
    return Fraction(p * p * (p - 3) * (p * p - p - 4), 16 * (p - 2))


def Lns_11_named(p: int) -> Fraction:
    return Fraction(p * p * (p - 1) * (p * p - 5 * p + 8), 16 * (p - 2))


def Lns_mix_named(p: int) -> Fraction:
    return Fraction(p * p * (p - 1) * (p - 3), 16)


def Lns_mm_named(p: int) -> Fraction:
    return Fraction(p * p * (p - 1) * (p - 1) * (p - 4), 16 * (p - 2))


def L1d_from_weights(
    p: int, n00: Fraction, n01: Fraction, n10: Fraction, npm: Fraction, ng: Fraction
) -> Fraction:
    k = p * (p - 1) // 2
    return (
        n00 * k * k
        + (n01 + n10) * k * p * EJ_mu(p)
        + npm * p * p * EJ2_mu(p)
        + ng * p * p * E_JJ_mu(p)
    )


def weights_named(p: int, slot: str) -> dict[str, Fraction]:
    d = Fraction(p + 1)
    two = Fraction(2, d)
    table = {
        "pm1": (two, Fraction(0), Fraction(0), Fraction(p - 1, d), Fraction(0)),
        "Sub": (two, Fraction(0), Fraction(0), Fraction(0), Fraction(p - 1, d)),
        "sq_11": (
            Fraction(0),
            two,
            two,
            Fraction(4, d),
            Fraction(p - 7, d),
        ),
        "sq_mix": (
            Fraction(0),
            two,
            two,
            two,
            Fraction(p - 5, d),
        ),
        "sq_mm": (Fraction(0), two, two, Fraction(0), Fraction(p - 3, d)),
        "ns_11": (
            Fraction(0),
            two,
            Fraction(0),
            Fraction(4, d),
            Fraction(p - 5, d),
        ),
        "ns_mix": (Fraction(0), two, Fraction(0), two, Fraction(p - 3, d)),
        "ns_mm": (Fraction(0), two, Fraction(0), Fraction(0), Fraction(p - 1, d)),
    }
    n00, n01, n10, npm, ng = table[slot]
    return {"n00": n00, "n01": n01, "n10": n10, "npm": npm, "ng": ng}


def L1d_named(p: int, slot: str) -> Fraction:
    return {
        "pm1": Lpm1_1d_named,
        "Sub": Lpar_1d_named,
        "sq_11": Lsq_11_named,
        "sq_mix": Lsq_mix_named,
        "sq_mm": Lsq_mm_named,
        "ns_11": Lns_11_named,
        "ns_mix": Lns_mix_named,
        "ns_mm": Lns_mm_named,
    }[slot](p)


def slot_of(F, p: int, s: int) -> str:
    if s in (1, F["neg1"]):
        return "pm1"
    infp = pow_f(F, s, p) == s
    if infp:
        return "Sub"
    sm1 = F["add"][s][F["neg"][1]]
    sp1 = F["add"][s][1]
    cm = int(F["chi"][sm1])
    cp = int(F["chi"][sp1])
    ch = int(F["chi"][s])
    pair = (cm, cp)
    if ch == 1:
        if pair == (1, 1):
            return "sq_11"
        if pair == (-1, -1):
            return "sq_mm"
        return "sq_mix"
    if pair == (1, 1):
        return "ns_11"
    if pair == (-1, -1):
        return "ns_mm"
    return "ns_mix"


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in CERT_WT + (17, 19):
        F = _kernel_field(p)
        fp = [int(F["chi"][i]) for i in range(1, p)]
        if any(c != 1 for c in fp):
            ok = False
        n_plus = 0
        n_chk = 0
        for ab in _linear_forms(p):
            if not is_type_plus_form(p, *ab):
                continue
            n_plus += 1
            nchi = 0
            nker = 0
            for x in range(1, F["q"]):
                a, b = x % p, x // p
                if (ab[0] * a + ab[1] * b) % p == 0:
                    nker += 1
                    if int(F["chi"][x]) == 1:
                        nchi += 1
            if nker != p - 1 or nchi != p - 1:
                ok = False
            n_chk += 1
        if n_plus != (p + 1) // 2 or n_chk != n_plus:
            ok = False
        rows[str(p)] = {"chi_Fp": fp[:4], "n_plus": n_plus}
    # fail-when-wrong: χ_p(3)=−1 at p=5, but F_q-χ(3)=+1
    F5 = _kernel_field(5)
    if int(F5["chi"][3]) != 1:
        ok = False
    if pow(3, (5 - 1) // 2, 5) != 4:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "χ|_{F_p^×}≡+1 and every Type+ line is all squares. "
            "Fail: χ|_{F_p}=χ_p (χ_p(3)=−1 at p=5)."
        ),
    }


def _count_weights(p: int) -> dict[str, dict[str, Fraction]]:
    F = _kernel_field(p)
    forms = [ab for ab in _linear_forms(p) if is_type_plus_form(p, *ab)]
    raw: dict[str, dict[str, int]] = defaultdict(
        lambda: {"n00": 0, "n01": 0, "n10": 0, "npm": 0, "ng": 0, "n": 0}
    )
    for s in range(1, F["q"]):
        sl = slot_of(F, p, s)
        for alpha, beta in forms:

            def lin(x: int, al: int = alpha, be: int = beta) -> int:
                return (al * (x % p) + be * (x // p)) % p

            for d in range(1, F["q"]):
                if F["chi"][d] != 1:
                    continue
                raw[sl]["n"] += 1
                ad = lin(d)
                bd = lin(F["mul"][s][d])
                if ad == 0 and bd == 0:
                    raw[sl]["n00"] += 1
                elif ad == 0:
                    raw[sl]["n01"] += 1
                elif bd == 0:
                    raw[sl]["n10"] += 1
                else:
                    rho = (bd * pow(ad, p - 2, p)) % p
                    if rho in (1, p - 1):
                        raw[sl]["npm"] += 1
                    else:
                        raw[sl]["ng"] += 1
    out: dict[str, dict[str, Fraction]] = {}
    for sl, c in raw.items():
        n = c["n"]
        out[sl] = {k: Fraction(c[k], n) for k in ("n00", "n01", "n10", "npm", "ng")}
    return out


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in CERT_WT:
        got = _count_weights(p)
        for sl, wt in got.items():
            want = weights_named(p, sl)
            if sl == "sq_11" and p == 5:
                ok = False if sl in got else ok
                continue
            if wt != want:
                ok = False
        if p == 5 and "sq_11" in got:
            ok = False
        rows[str(p)] = {sl: {k: str(v) for k, v in wt.items()} for sl, wt in got.items()}
        # fail: force n10=n01 on ns_mm
        bad = dict(weights_named(p, "ns_mm"))
        bad["n10"] = bad["n01"]
        live = L1d_from_weights(p, **weights_named(p, "ns_mm"))
        forced = L1d_from_weights(p, **bad)
        if p == 5 and forced == live:
            ok = False
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "n00,n01,n10,npm named by χ(s), χ(s±1), 1_{s∈F_p}. "
            "Fail: n10=n01 on nonsquares."
        ),
    }


def _L_from_Y(Y: np.ndarray, F) -> dict[int, Fraction]:
    Ds = Y[:, 1 : 1 + F["q"]] < 0
    add = np.array(F["add"], dtype=np.int32)
    N = np.stack(
        [(Ds & Ds[:, add[:, d]]).sum(axis=1) for d in range(F["q"])], axis=1
    )
    ENN = (N.astype(np.float64).T @ N) / N.shape[0]
    out: dict[int, Fraction] = {}
    for s in range(1, F["q"]):
        acc = 0.0
        n = 0
        for d in range(1, F["q"]):
            if F["chi"][d] != 1:
                continue
            acc += ENN[d, F["mul"][s][d]]
            n += 1
        out[s] = Fraction(acc / n).limit_denominator(10**7)
    return out


def _build_1d(p: int):
    F = _kernel_field(p)
    m = (p + 1) // 2
    forms = [ab for ab in _linear_forms(p) if is_type_plus_form(p, *ab)]
    plus = list(combinations(range(p), m))
    rows = [
        lift_sigma_vector(p, a, b, set(A)) for a, b in forms for A in plus
    ]
    Y = np.array(rows, dtype=np.float64)
    if len(Y) != n_1d(p):
        raise RuntimeError(f"n1d {len(Y)} != {n_1d(p)}")
    return F, Y


def prove_C() -> dict:
    ok = True
    # algebra: Johnson combo equals the closed form
    alg = {}
    for p in (5, 7, 11, 13, 17, 19, 23):
        row = {}
        for sl in ("pm1", "Sub", "sq_11", "sq_mix", "sq_mm", "ns_11", "ns_mix", "ns_mm"):
            wt = weights_named(p, sl)
            via = L1d_from_weights(p, **wt)
            named = L1d_named(p, sl)
            row[sl] = str(named)
            if via != named:
                ok = False
        if L1d_named(p, "Sub") != Lpar_1d_named(p):
            ok = False
        if L1d_named(p, "sq_11") != Lspl_1d_named(p):
            ok = False
        if Lns_mix_named(p) != mu_plus(p) * mu_minus(p):
            ok = False
        alg[str(p)] = row
    if Lpm1_1d_inverted(5) == Lpm1_1d_named(5):
        ok = False
    if Lpm1_1d_inverted(5) != Fraction(125, 4):
        ok = False
    if Lsq_11_named(5) != 25:
        ok = False
    # live Type+ lifts
    live = {}
    for p in CERT_LIFTS:
        F, Y = _build_1d(p)
        L1 = _L_from_Y(Y, F)
        by = defaultdict(set)
        for s, val in L1.items():
            by[slot_of(F, p, s)].add(val)
        live[str(p)] = {sl: [str(v) for v in sorted(vs)] for sl, vs in by.items()}
        for sl, vs in by.items():
            if sl == "sq_11" and p == 5:
                ok = False
                continue
            if vs != {L1d_named(p, sl)}:
                ok = False
        if p == 5:
            if 25 in {v for vs in by.values() for v in vs}:
                ok = False
            if "sq_11" in by:
                ok = False
    return {
        "proved": bool(ok),
        "algebra": alg,
        "live": live,
        "wrong_pm1_p5": str(Lpm1_1d_inverted(5)),
        "empty_sq11_p5": str(Lsq_11_named(5)),
        "theorem": (
            "Eight L_1D slots named. Fail: invert p−2 on L_pm1 "
            f"({Lpm1_1d_inverted(5)}≠{Lpm1_1d_named(5)}); "
            "fail: 25 is live at p=5."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "All Type+ L_1D slots named. Ensemble L_τ still has "
            "w=n_1d/|H+| and unnamed L_NL. Q_τ unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.448  every Type+ L_1D slot named in p", flush=True)
    A = prove_A()
    print(f"  A χ on F_p: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B weights: {B['proved']}", flush=True)
    C = prove_C()
    print(
        f"  C named: {C['proved']} pm1_p5={C['algebra']['5']['pm1']} "
        f"fail={C['wrong_pm1_p5']}",
        flush=True,
    )
    D = prove_open()
    print(f"  D open phi_F={D['phi_F_ge_6']} e1={D['e1']}", flush=True)
    CAT.write_text(
        json.dumps(
            {
                "A": A["proved"],
                "B": B["proved"],
                "C": C["proved"],
                "pm1_5": C["algebra"]["5"]["pm1"],
                "wrong": C["wrong_pm1_p5"],
            },
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.448",
        "title": (
            "Every Type+ L_1D slot is named in p by Johnson weights "
            "on (χ(s), χ(s−1), χ(s+1), 1_{s∈F_p})"
        ),
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "chi_Fp_trivial": A["proved"],
            "weights_named": B["proved"],
            "L1d_slots_named": C["proved"],
            "phi_F_ge_6_proved_general": D["phi_F_ge_6"],
            "residual_ii_k_eq_4p_empty": D["residual_ii_k_eq_4p_empty"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15448.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
