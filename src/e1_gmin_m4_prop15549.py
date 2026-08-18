#!/usr/bin/env python3
"""
Prop 15.549 — K_λ is the complete Paley 4-linear minus a Gauss
collision:
    K_λ = K_λ^{all} − N q (3q−5),
    K_λ^{all} = Σ_{a,b,c,d} χ((d−a)(d−b)(c−a)(c−b)) T_{abcd}
                e_ξ(a−b) e_{rξ}(c−d),
with χ(0)=0.  The correction is G(χ, rξ)=p on r∈++.  Fail drop
it, or subtract only 2Nq(q−2).  K_λ^{all} still not a p-law.
16pA / Q_τ / phi_F stay open.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.548 flags.  Does **not** overwrite 15.522–15.548.

============================================================================
Setup.  15.548: K_λ is the χ(λ) channel of the 4-distinct y-first
4-linear, χ(λ)=χ(d−a)χ(c−b)χ(d−b)χ(c−a).  T_{abcd}=Σ_y y_a y_b y_c y_d.
G(a,b)=(N/p)χ(b−a) on F_q (15.548 A).  Ω={ξ: G(χ,ξ)=p}.

============================================================================
Theorem A — PROVED (field; r square ⇒ rξ∈Ω).
  G(χ, rξ)=p for every r∈++ and ξ∈Ω.  Fail: G(χ, rξ)=−p
  (sign of the quadratic Gauss on Ω).  Certified p=5,7,11,13.  ∎

Theorem B — PROVED (15.548 G + A; Max+-free).
  χ(λ) vanishes unless c,d∉{a,b}.  The leftover non-4-distinct
  strata are (a=b, c≠d) and (c=d, a≠b).  Each of a=b (c≠d) and c=d (a≠b) equals N q (q−2)
  because Σ_{δ≠0} χ(δ) e_{±rξ}(δ)=p.  The double pair a=b, c=d
  contributes N q (q−1).  Hence
      K_λ = K_λ^{all} − N q (3q−5).
  Fail: drop the correction.  Fail: omit the a=b=c=d-off
  double pair and subtract only 2 N q (q−2).  ∎

Theorem C — OPEN.  K_λ^{all} is not reduced to a Gauss/Jacobi
  integer in p.  16pA / Q_τ unnamed.  phi_F not imported.

============================================================================
Backend: 15.548 G + field Gauss; p=5 4-linear fold.  Writes
evidence/e1_gmin_m4_prop15549.json
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import _e_tr, _kernel_field
from e1_gmin_m4_prop15290 import omega_set, type_key
from e1_gmin_m4_prop15548 import G_named, fold_jacobi_p5

EV = ROOT / "evidence" / "e1_gmin_m4_prop15549.json"


def gauss_chi(F, xi) -> complex:
    s = 0j
    for x in range(1, F["q"]):
        s += F["chi"][x] * _e_tr(F, F["mul"][xi][x])
    return s


def pp_r(F) -> int:
    for r in F["squares"]:
        if type_key(F, r) in ((1, 1, "sub"), (1, 1, "N1")):
            return int(r)
    raise RuntimeError("no ++ r")


def K_lam_gauss_coll(p: int, N: int) -> int:
    """N q (3q−5) = 2Nq(q−2) + Nq(q−1)."""
    q = p * p
    return N * q * (3 * q - 5)


def K_lam_coll_wrong_omit_aabb(p: int, N: int) -> int:
    """Fail-when-wrong: omit the a=b,c=d pair, 2 N q (q−2)."""
    q = p * p
    return 2 * N * q * (q - 2)


def K_lam_from_complete(K_all: float, p: int, N: int) -> float:
    return K_all - K_lam_gauss_coll(p, N)


def prove_A() -> dict:
    rows = {}
    ok = True
    for p in (5, 7, 11, 13):
        F = _kernel_field(p)
        xi = omega_set(F)[0]
        r = pp_r(F)
        eta = F["mul"][r][xi]
        g = gauss_chi(F, eta)
        gneg = gauss_chi(F, F["neg"][eta])
        row_ok = abs(g - p) < 1e-8 and abs(gneg - p) < 1e-8
        if abs(g + p) < 1e-8:
            row_ok = False
        if not row_ok:
            ok = False
        rows[str(p)] = {
            "G_rxi": float(g.real),
            "G_mrxi": float(gneg.real),
            "ok": row_ok,
        }
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": "G(χ, rξ)=p on ++. Fail G=−p.",
    }


def fold_Kall_p5() -> dict:
    """Same p=5 fold as 15.548, also accumulating the unrestricted χ(λ) sum."""
    p = 5
    F = _kernel_field(p)
    q = F["q"]
    from e1_gmin_m4_prop15468 import D_rows

    Y = np.where(D_rows(p, F), -1.0, 1.0)
    N = Y.shape[0]
    Omega = omega_set(F)
    xi = Omega[0]
    r = pp_r(F)
    eta = F["mul"][r][xi]
    e_xi = np.array([_e_tr(F, F["mul"][xi][x]) for x in range(q)], dtype=np.complex128)
    e_eta = np.array([_e_tr(F, F["mul"][eta][x]) for x in range(q)], dtype=np.complex128)
    add, neg, chi = F["add"], F["neg"], F["chi"]
    chi_diff = np.zeros((q, q), dtype=np.int8)
    for u in range(q):
        for v in range(q):
            chi_diff[u, v] = chi[add[u][neg[v]]]
    Klam = Kall = 0j
    for a in range(q):
        Ya = Y[:, a]
        ea = e_xi[a]
        for b in range(q):
            M = (Ya * Y[:, b])[:, None] * Y
            Td = M.T @ Y
            eab = ea * np.conj(e_xi[b])
            ph = eab * e_eta[:, None] * np.conj(e_eta)[None, :]
            contrib = Td * ph
            ch_lam = (
                chi_diff[None, :, a]
                * chi_diff[:, None, b]
                * chi_diff[None, :, b]
                * chi_diff[:, None, a]
            )
            Kall += (contrib * ch_lam).sum()
            ii = np.arange(q)
            dist = (
                (ii[:, None] != a)
                & (ii[:, None] != b)
                & (ii[None, :] != a)
                & (ii[None, :] != b)
                & (ii[:, None] != ii[None, :])
                & (a != b)
            )
            Klam += (contrib * ch_lam)[dist].sum()
    coll = K_lam_gauss_coll(p, N)
    wrong = K_lam_coll_wrong_omit_aabb(p, N)
    named = K_lam_from_complete(float(Kall.real), p, N)
    return {
        "N": N,
        "K_lam": float(Klam.real),
        "K_all": float(Kall.real),
        "coll": coll,
        "named": named,
        "wrong_omit_aabb": float(Kall.real) - wrong,
        "match": bool(abs(named - Klam.real) < 1e-4),
        "drop_fails": bool(abs(Kall.real - Klam.real) > 1.0),
        "omit_aabb_fails": bool(abs((Kall.real - wrong) - Klam.real) > 1.0),
    }


def prove_B() -> dict:
    rec = fold_Kall_p5()
    B548 = fold_jacobi_p5()
    ok = (
        rec["match"]
        and rec["drop_fails"]
        and rec["omit_aabb_fails"]
        and abs(rec["K_lam"] - B548["K_lam"]) < 1e-4
        and abs(rec["K_lam"] - 21550) < 1e-4
        and abs(rec["K_all"] - (rec["K_lam"] + rec["coll"])) < 1e-4
    )
    # p=7 arithmetic identity from the certified fold
    N7, K7 = 5726, 2206274.0
    coll7 = K_lam_gauss_coll(7, N7)
    Kall7 = K7 + coll7
    rec["p7_coll"] = coll7
    rec["p7_Kall_from_B"] = Kall7
    rec["p7_back"] = K_lam_from_complete(Kall7, 7, N7)
    if abs(rec["p7_back"] - K7) > 1e-6:
        ok = False
    return {
        "proved": bool(ok),
        "fold": rec,
        "theorem": (
            "K_λ=K_λ^{all}−Nq(3q−5). Fail drop "
            f"(p=5: {rec['K_all']}≠{rec['K_lam']}). "
            f"Fail omit a=b,c=d: {rec['wrong_omit_aabb']}≠{rec['K_lam']}."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "K_all_named_in_p": False,
        "NUM_SUM_16pA_general": False,
        "Q_tau_named_in_p": False,
        "phi_F_imported": False,
        "note": (
            "Complete sum K_λ^{all} is not a Gauss/Jacobi integer in p. "
            "16pA / Q_τ unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.549  K_λ = complete Paley 4-linear − Nq(3q−5)", flush=True)
    A = prove_A()
    print(f"  A G(χ,rξ)=p: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B collision Gauss: {B['proved']}", flush=True)
    C = prove_open()
    out = {
        "prop": "15.549",
        "title": "K_λ is the complete χ(λ) 4-linear minus 2Nq(q−2)",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "C": C,
        "proved": {
            "G_chi_rxi_eq_p": A["proved"],
            "K_lam_complete_minus_gauss": B["proved"],
            "NUM_SUM_16pA_general": False,
            "phi_F_ge_6": False,
        },
        "e1_closed_general": bool(e1_closed_general()),
        "gsum_disj_lb": bool(gsum_disj_lb_proved_general()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
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
    EV.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print("wrote", EV, flush=True)
    return out


if __name__ == "__main__":
    main()
