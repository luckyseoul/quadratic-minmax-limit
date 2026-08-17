#!/usr/bin/env python3
"""
Prop 15.495 — L_χ on the nonlinear Max+ family is not a
named function of p, so 15.298 does not name Q_NL.

At p=5 (NL = Hoffman) χ=+1 square L is n or q−1.  At p=7
the same slots have dens {6,19,38,57,114}, not 409, and
lie outside the Gauss/Jacobi rational closure (except the
already-named ns mean μ₊μ₋).  Live L_NL recovers Q_NL via
15.298; the p=5 names do not.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum /
pairing / 15.279–15.494 flags.  Does **not** overwrite
15.493/494.

============================================================================
Setup.  L_χ(s)=E[N(d)N(s d)] on y_∞=+1 Max+ with D not an
F_p-line (NL = Max+ \\ Type+ 1D).  n=p²+1, q=p².
15.298: Q is the half-Gauss image of the L table of the
ensemble.  Q_NL := Q_{++}^{NL}/q² (15.356).  Live
Q_NL ∈ {64/15, 632/171}.  named_gj_ints as in 15.330.

============================================================================
Theorem A — PROVED (NL cache fold; certified p=5,7).
  χ=+1 square L_NL at p=5:
      pm1 = N1 = n = 26,   sub = N* = q−1 = 24.
  At p=7 the same slots are
      pm1 4435/38,  sub 2113/19,  ++N1 673/6,
      --N1 6131/57,  N* same 12349/114,  N* cross 6269/57
  with dens {6,19,38,57,114} not 409.
  Fail: L_NL(pm1)=n at p=7 (4435/38 ≠ 50);
  fail: L_NL(sub)=q−1 at p=7 (2113/19 ≠ 48);
  fail: those dens divide 409.  ∎

Theorem B — PROVED (named_gj_ints(7) + field catalog).
  Except ns-mean = μ₊μ₋, no p=7 NL χ=+1 slot is in the
  rational closure of named_gj_ints(7) ∪ {n, q−1, μ±, L_par(R)}.
  19,38,57,114,2113,673,4435 ∉ named_gj_ints(7).
  Fail: 19∈S7; 2113/19 = L_par(R)=133.  ∎

Theorem C — PROVED (15.298 on the NL ensemble).
  Live L_NL recovers Q_NL (64/15, 632/171).
  Replacing square L by the p=5 names n / q−1 and ns by μ₊μ₋
  misses both (1568/375 ≠ 64/15; negative ≠ 632/171).
  Constant-L μ₊μ₋ recovers 4(q−1)/q, not Q_NL.
  Fail: claim the p=5 L-law names Q_NL via 15.298.  ∎

Theorem D — OPEN.  L_NL unnamed in p ⇒ 15.298 does not
  import Q_NL.  Do not interpolate {26, 4435/38}.
  Mixture Q++ = w Q_1d + (1−w) Q_NL still carries unnamed
  N+ unless Q_1d=Q_NL (false: 16/9≠64/15, 104/15≠632/171).
  phi_F stays False.

============================================================================
Backend: integer ENN on /tmp/maxplus_p{5,7}.npy; 15.298
predict_Q.  Writes evidence/e1_gmin_m4_prop15495.json
"""
from __future__ import annotations

import json
import os
import sys
from fractions import Fraction
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field  # noqa: E402
from e1_gmin_m4_prop15292 import n_1d  # noqa: E402
from e1_gmin_m4_prop15294 import _fp_dirs, _is_1d  # noqa: E402
from e1_gmin_m4_prop15298 import (  # noqa: E402
    coarse_Q_class,
    load_N,
    mu_minus,
    mu_plus,
    ns_key,
    predict_Q,
)
from e1_gmin_m4_prop15330 import named_gj_ints  # noqa: E402
from e1_gmin_m4_prop15356 import Q_1d_pp_named  # noqa: E402
from e1_gmin_m4_prop15409 import D_form_on_lattice_general  # noqa: E402
from e1_gmin_m4_prop15411 import LIVE_QNL  # noqa: E402
from e1_gmin_m4_prop15412 import Q_NL_is_short_rational  # noqa: E402
from e1_gmin_m4_prop15430 import Lpar_R_named  # noqa: E402

EV = ROOT / "evidence" / "e1_gmin_m4_prop15495.json"
YPATH = {5: "/tmp/maxplus_p5.npy", 7: "/tmp/maxplus_p7.npy"}

LIVE_LNL_PLUS = {
    5: {
        ("pm1",): Fraction(26),
        (1, 1, "sub"): Fraction(24),
        (-1, -1, "N*"): Fraction(24),
        (-1, 1, "N1"): Fraction(26),
        (1, -1, "N1"): Fraction(26),
    },
    7: {
        ("pm1",): Fraction(4435, 38),
        (1, 1, "sub"): Fraction(2113, 19),
        (1, 1, "N1"): Fraction(673, 6),
        (-1, -1, "N1"): Fraction(6131, 57),
        (-1, -1, "N*"): Fraction(12349, 114),
        (-1, 1, "N*"): Fraction(6269, 57),
        (1, -1, "N*"): Fraction(6269, 57),
    },
}

P5_SQ_DENS = {1}
P7_SQ_DENS = {6, 19, 38, 57, 114}
P7_DENS = {2, 6, 19, 38, 57, 114}


def p5_square_name(tk, chi: int, p: int) -> Fraction:
    """p=5 law: n on pm1/N1, q−1 on sub/N*; χ=−1 is /4."""
    n = Fraction(p * p + 1)
    qm1 = Fraction(p * p - 1)
    if tk == ("pm1",) or (isinstance(tk, tuple) and "N1" in str(tk)):
        base = n
    else:
        base = qm1
    return base if chi == 1 else base / 4


def field_catalog(p: int) -> dict[str, Fraction]:
    q = p * p
    return {
        "n": Fraction(q + 1),
        "q-1": Fraction(q - 1),
        "mu+mu-": mu_plus(p) * mu_minus(p),
        "mu+2": mu_plus(p) ** 2,
        "Lpar_R": Lpar_R_named(p),
    }


def gj_closure(p: int) -> set[Fraction]:
    S = named_gj_ints(p)
    rats: set[Fraction] = set()
    base = [a for a in S if a]
    for a in base:
        rats.add(Fraction(a))
        for b in base:
            rats.add(Fraction(a, b))
    for v in field_catalog(p).values():
        rats.add(v)
        for a in base:
            rats.add(v * a)
            rats.add(v / a)
    return rats


def fold_NL(p: int) -> dict:
    F = _kernel_field(p)
    N = load_N(p, F)
    Y = np.load(YPATH[p])
    Yp = Y[Y[:, 0] > 0]
    z = np.sign(Yp[:, 1 : 1 + F["q"]].astype(np.float64))
    dirs = _fp_dirs(p)
    one = np.array([_is_1d(zi, dirs, F["add"], p) for zi in z])
    mask = ~one
    Nm = N[mask].astype(np.int64)
    G = Nm.T @ Nm
    n = int(Nm.shape[0])
    q = F["q"]
    L: dict = {}
    for d in range(1, q):
        for s in range(1, q):
            key = (int(F["chi"][d]), ns_key(F, s))
            val = Fraction(int(G[d, F["mul"][s][d]]), n)
            prev = L.get(key)
            if prev is None:
                L[key] = val
            elif prev != val:
                raise RuntimeError(f"L not constant on {key} at p={p}")
    return {
        "F": {**F, "p": p},
        "L": L,
        "EN": N[mask].astype(np.float64).mean(axis=0),
        "n_NL": n,
        "n_1d": int(one.sum()),
    }


def L_plus_squares(L: dict) -> dict:
    out = {}
    for (chi, tk), v in L.items():
        if chi != 1:
            continue
        if isinstance(tk, tuple) and tk and tk[0] == "ns":
            continue
        out[tk] = v
    return out


def apply_p5_law(L: dict, p: int) -> dict:
    mumu = mu_plus(p) * mu_minus(p)
    out = {}
    for (chi, tk), _v in L.items():
        if isinstance(tk, tuple) and tk and tk[0] == "ns":
            out[(chi, tk)] = mumu
        else:
            out[(chi, tk)] = p5_square_name(tk, chi, p)
    return out


def Qpp_from_pred(pred: dict, F) -> Fraction:
    q2 = float(F["q"] ** 2)
    vals = []
    for r, Qr in pred.items():
        if r in (1, F["neg1"]):
            continue
        c = coarse_Q_class(F, r)
        if c == "++" or (F["p"] == 5 and str(c).startswith("mixed")):
            vals.append(Qr / q2)
    return Fraction(float(np.mean(vals))).limit_denominator(4000)


def prove_A(folds: dict | None = None) -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        pack = (folds or {}).get(p) or fold_NL(p)
        plus = L_plus_squares(pack["L"])
        live = LIVE_LNL_PLUS[p]
        if set(plus) != set(live):
            ok = False
        for tk, val in live.items():
            if plus.get(tk) != val:
                ok = False
        dens = {v.denominator for v in plus.values()}
        rows[p] = {str(k): str(v) for k, v in plus.items()}
        rows[p]["dens"] = sorted(dens)
        rows[p]["n_NL"] = pack["n_NL"]
        rows[p]["n_1d"] = pack["n_1d"]
        if pack["n_1d"] != n_1d(p):
            ok = False
        if p == 5:
            if plus[("pm1",)] != Fraction(p * p + 1):
                ok = False
            if plus[(1, 1, "sub")] != Fraction(p * p - 1):
                ok = False
            if dens != P5_SQ_DENS:
                ok = False
        if p == 7:
            if plus[("pm1",)] == Fraction(p * p + 1):
                ok = False
            if plus[(1, 1, "sub")] == Fraction(p * p - 1):
                ok = False
            if dens != P7_SQ_DENS:
                ok = False
            if 409 in dens:
                ok = False
            for d in (19, 38, 57, 114):
                if 409 % d == 0 and d != 1:
                    ok = False
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "p=5 NL L_+=n or q−1; p=7 dens {6,19,38,57,114} not 409. "
            "Fail: L(pm1)=n or L(sub)=q−1 at p=7."
        ),
    }


def prove_B() -> dict:
    S7 = named_gj_ints(7)
    clos = gj_closure(7)
    ok = True
    hits = {}
    for tk, val in LIVE_LNL_PLUS[7].items():
        hit = val in clos
        hits[str(tk)] = hit
        if hit:
            ok = False
    mumu = mu_plus(7) * mu_minus(7)
    if mumu not in clos:
        ok = False
    if Lpar_R_named(7) != 133:
        ok = False
    if LIVE_LNL_PLUS[7][(1, 1, "sub")] == Lpar_R_named(7):
        ok = False
    for bad in (19, 38, 57, 114, 2113, 673, 4435, 409):
        if bad in S7:
            ok = False
    return {
        "proved": bool(ok),
        "slot_in_gj": hits,
        "S7_has_19": 19 in S7,
        "Lpar_R_7": str(Lpar_R_named(7)),
        "theorem": (
            "p=7 NL square L not in GJ closure. "
            "Fail: 19∈S7 or 2113/19=L_par(R)."
        ),
    }


def prove_C(folds: dict | None = None) -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        pack = (folds or {}).get(p) or fold_NL(p)
        F = pack["F"]
        live = {k: float(v) for k, v in pack["L"].items()}
        pred = predict_Q(p, F, live, pack["EN"], naive_ns=False)
        Qlive = Qpp_from_pred(pred, F)
        law = {k: float(v) for k, v in apply_p5_law(pack["L"], p).items()}
        Qlaw = Qpp_from_pred(predict_Q(p, F, law, pack["EN"], naive_ns=False), F)
        mm = {k: float(mu_plus(p) * mu_minus(p)) for k in pack["L"]}
        Qmm = Qpp_from_pred(predict_Q(p, F, mm, pack["EN"], naive_ns=False), F)
        target = LIVE_QNL[p]
        if abs(float(Qlive) - float(target)) > 1e-8:
            ok = False
        if abs(float(Qlaw) - float(target)) < 1e-6:
            ok = False
        if abs(float(Qmm) - float(target)) < 1e-6:
            ok = False
        rows[p] = {
            "live": str(Qlive),
            "target": str(target),
            "p5law": str(Qlaw),
            "mumu": str(Qmm),
        }
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "15.298(live L_NL)=Q_NL. p=5 L-names miss. "
            "Fail: p=5 law recovers 64/15 or 632/171."
        ),
    }


def prove_open() -> dict:
    q1d5, q1d7 = Q_1d_pp_named(5), Q_1d_pp_named(7)
    return {
        "proved": False,
        "Q_tau_named_in_p": False,
        "phi_F_imported": False,
        "Q_1d_eq_Q_NL": bool(
            q1d5 == LIVE_QNL[5] or q1d7 == LIVE_QNL[7]
        ),
        "Q_1d_pp": {5: str(q1d5), 7: str(q1d7)},
        "note": (
            "L_NL unnamed in p; 15.298 cannot import Q_NL. "
            "Q_1d≠Q_NL so the mixture still has unnamed N+. "
            "⟨δ,ψ⟩≤2 stays OPEN."
        ),
    }


def main() -> dict:
    print("Prop 15.495  L_NL not named in p; 15.298 does not name Q_NL", flush=True)
    folds = {p: fold_NL(p) for p in (5, 7)}
    A = prove_A(folds)
    print(f"  A L table: {A['proved']} {A['by_p']}", flush=True)
    B = prove_B()
    print(f"  B GJ miss: {B['proved']} S7_has_19={B['S7_has_19']}", flush=True)
    C = prove_C(folds)
    print(f"  C 15.298: {C['proved']} {C['by_p']}", flush=True)
    D = prove_open()
    print(f"  D floor: {D['proved']} phi_F imported={D['phi_F_imported']}", flush=True)
    out = {
        "prop": "15.495",
        "title": "L_χ^{NL} is not named in p; 15.298 does not name Q_NL",
        "proved": {
            "L_NL_p5_is_n_or_qm1": A["proved"],
            "L_NL_p7_not_gj": B["proved"],
            "gauss_image_needs_unnamed_L": C["proved"],
            "phi_F_ge_6_proved_general": bool(phi_F_ge_6_proved_general()),
            "Q_tau_named_in_p": False,
        },
        "algebra": {"A": A, "B": B, "C": C, "D": D},
        "L_status": "OPEN",
        "e1_closed_general": bool(e1_closed_general()),
        "gsum_disj_lb": bool(gsum_disj_lb_proved_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "D_form_on_lattice_general": bool(D_form_on_lattice_general()),
        "Q_NL_is_short_rational": bool(Q_NL_is_short_rational()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "flags_not_flipped": [
            "phi_F_ge_6",
            "e1",
            "L",
            "Aut-Schur",
            "Gsum",
            "pairing",
        ],
    }
    EV.parent.mkdir(parents=True, exist_ok=True)
    EV.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print("wrote", EV, flush=True)
    return out


if __name__ == "__main__":
    main()
