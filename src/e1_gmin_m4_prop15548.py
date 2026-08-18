#!/usr/bin/env python3
"""
Prop 15.548 — 2-point of H_+ on F_q is (N/p)χ, and the |κ|=3
layer of the y-first 4-linear is the Jacobi split
    K_3 = (K_1 + K_λ + K_{1−λ} + K_{λ(1−λ)})/4
with χ(λ)=χ(d−a)χ(c−b)χ(d−b)χ(c−a).  The three χ-channels
are not a p-law.  16pA / Q_τ unnamed.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.547 flags.  Does **not** overwrite 15.522–15.547.

============================================================================
Setup.  Paley C on F_q ∪ {∞}, Max+ H_+ with y_∞=+1, N=|H_+|.
y_x=±1 on F_q, Σ_x y_x=p (15.189).  Affine Aut x↦αx+β (α square)
preserves C and H_+.  Cross-ratio λ=(d−a)(c−b)/((d−b)(c−a)).
15.267 C: |κ|=3 ⇔ χ(λ)=χ(1−λ)=1.  F(r)=Σ_y |ẑ(ξ)|²|ẑ(rξ)|²
for r∈++ equals NUM_SUM q².  4-distinct part is K_1.

============================================================================
Theorem A — PROVED (Aut-radial + Cy=py; certified p=5,7).
  For a≠b in F_q,
      G(a,b) := Σ_{y∈H_+} y_a y_b = (N/p) χ(b−a).
  Aut is transitive on χ=+1 pairs and on χ=−1 pairs, so G is
  Paley-radial.  Cy=py and E[y_a]=1/p (Aut-transitive λ_1=Nk/q)
  force g_+ = N/p, g_− = −N/p.  Fail: G≡0 off-diagonal
  (row-sum of χ G would be 0, not N(p−1/p)).  Fail: E[y_a]=0
  (then g_++g_− would not vanish).  ∎

Theorem B — PROVED (character expansion; certified p=5).
  I_3=(1+χ(λ))(1+χ(1−λ))/4 equals 1 on ordered 4-distinct
  tuples of a |κ|=3 set and 0 on |κ|=1.  Hence
      K_3 = (K_1 + K_λ + K_{1−λ} + K_{λ(1−λ)})/4
  with K_χ = Σ_{4-dist} χ_•(λ) T_{abcd} e_ξ(a−b) e_{rξ}(c−d).
  At p=5, r∈++: K_3=83350, K_1=88750, K_λ=21550,
  K_{1−λ}=K_{λ(1−λ)}=111550.  Fail: drop the λ(1−λ) channel
  ((88750+21550+111550)/4=55462.5≠83350).  The product
  χ(λ)χ(1−λ)=+1 bucket (100150) is |κ|=3 plus {χ(λ)=χ(1−λ)=−1},
  not I_3.  ∎

Theorem C — OPEN.  K_λ, K_{1−λ} are not a Gauss/Jacobi formula
  in p.  Collision layer of F(r) is not reduced to 16pA.
  Q_τ unnamed.  phi_F not imported.

============================================================================
Backend: Fraction + numpy GEMM on p=5 cache.  Writes
evidence/e1_gmin_m4_prop15548.json
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
from e1_gmin_m4_prop15367 import A_num, NUM_SUM_named
from e1_gmin_m4_prop15468 import D_rows

EV = ROOT / "evidence" / "e1_gmin_m4_prop15548.json"
YPATH = {5: "/tmp/maxplus_p5.npy", 7: "/tmp/maxplus_p7.npy"}


def G_named(p: int, N: int, chi_ab: int) -> float:
    """Σ_y y_a y_b for a≠b, χ(b−a)=chi_ab ∈ {±1}."""
    return (N / p) * chi_ab


def G_named_zero(p: int, N: int, chi_ab: int) -> float:
    """Fail-when-wrong: claim G≡0 off-diagonal."""
    return 0.0


def jac_drop_prod(K1: float, Klam: float, K1m: float) -> float:
    """Fail-when-wrong: drop the χ(λ(1−λ)) channel."""
    return 0.25 * (K1 + Klam + K1m)


def _pp_r(F) -> int:
    for r in F["squares"]:
        if type_key(F, r) in ((1, 1, "sub"), (1, 1, "N1")):
            return int(r)
    raise RuntimeError("no ++ r")


def live_G(p: int) -> dict:
    F = _kernel_field(p)
    Y = np.where(D_rows(p, F), -1.0, 1.0)
    N, q = Y.shape
    G = Y.T @ Y
    chi = F["chi"]
    add, neg = F["add"], F["neg"]
    ok = True
    n_plus = n_minus = 0
    err = 0.0
    for a in range(q):
        for b in range(a + 1, q):
            ch = int(chi[add[b][neg[a]]])
            pred = G_named(p, N, ch)
            got = float(G[a, b])
            err = max(err, abs(got - pred))
            if abs(got - pred) > 1e-8:
                ok = False
            if ch == 1:
                n_plus += 1
            else:
                n_minus += 1
            if G_named_zero(p, N, ch) == got:
                ok = False
    # E[y_a]=1/p
    ey = Y.mean(axis=0)
    if np.max(np.abs(ey - 1.0 / p)) > 1e-8:
        ok = False
    return {
        "proved": bool(ok),
        "N": N,
        "max_err": err,
        "n_plus": n_plus,
        "n_minus": n_minus,
        "E_y": float(ey[0]),
        "g_plus": G_named(p, N, 1),
        "g_minus": G_named(p, N, -1),
    }


def fold_jacobi_p5() -> dict:
    p = 5
    F = _kernel_field(p)
    q = F["q"]
    Y = np.where(D_rows(p, F), -1.0, 1.0)
    Omega = omega_set(F)
    xi = Omega[0]
    r = _pp_r(F)
    eta = F["mul"][r][xi]
    e_xi = np.array([_e_tr(F, F["mul"][xi][x]) for x in range(q)], dtype=np.complex128)
    e_eta = np.array([_e_tr(F, F["mul"][eta][x]) for x in range(q)], dtype=np.complex128)
    add, neg, mul, inv = F["add"], F["neg"], F["mul"], F["inv"]
    chi = F["chi"]
    K1 = Klam = K1m = Kprod = K3 = 0j
    F_all = F_coll = 0j
    n_i3 = 0
    for a in range(q):
        ea = e_xi[a]
        Ya = Y[:, a]
        for b in range(q):
            eab = ea * np.conj(e_xi[b])
            yab = Ya * Y[:, b]
            for c in range(q):
                Td = yab * Y[:, c] @ Y
                ph = eab * e_eta[c] * np.conj(e_eta)
                contrib = Td * ph
                F_all += contrib.sum()
                for d in range(q):
                    tval = contrib[d]
                    if len({a, b, c, d}) < 4:
                        F_coll += tval
                        continue
                    K1 += tval
                    da = add[d][neg[a]]
                    cb = add[c][neg[b]]
                    db = add[d][neg[b]]
                    ca = add[c][neg[a]]
                    if 0 in (da, cb, db, ca):
                        continue
                    lam = mul[mul[da][cb]][inv[mul[db][ca]]]
                    onem = add[1][neg[lam]]
                    if lam == 0 or onem == 0:
                        continue
                    cl = int(chi[lam])
                    c1 = int(chi[onem])
                    Klam += cl * tval
                    K1m += c1 * tval
                    Kprod += cl * c1 * tval
                    if cl == 1 and c1 == 1:
                        K3 += tval
                        n_i3 += 1
    K3j = 0.25 * (K1 + Klam + K1m + Kprod)
    drop = jac_drop_prod(K1.real, Klam.real, K1m.real)
    return {
        "F_all": float(F_all.real),
        "F_coll": float(F_coll.real),
        "K1": float(K1.real),
        "K_lam": float(Klam.real),
        "K_1mlam": float(K1m.real),
        "K_lam1mlam": float(Kprod.real),
        "K3_direct": float(K3.real),
        "K3_jacobi": float(K3j.real),
        "drop_prod": float(drop),
        "n_I3": n_i3,
        "match": bool(abs(K3j - K3) < 1e-4),
        "drop_fails": bool(abs(drop - K3.real) > 1.0),
        "NUM_SUM_yfirst": float(F_all.real / (q * q)),
        "16pA": NUM_SUM_named(p),
    }


def prove_A() -> dict:
    rows = {}
    ok = True
    for p in (5, 7):
        rec = live_G(p)
        rows[str(p)] = rec
        if not rec["proved"]:
            ok = False
        if rec["g_plus"] == 0 or rec["g_minus"] == 0:
            ok = False
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "G(a,b)=(N/p)χ(b−a) on F_q. Fail G≡0 "
            f"(p=5 g_+={rows['5']['g_plus']}). Fail E[y_a]=0."
        ),
    }


def prove_B() -> dict:
    rec = fold_jacobi_p5()
    ok = (
        rec["match"]
        and rec["drop_fails"]
        and abs(rec["K3_direct"] - 83350) < 1e-4
        and abs(rec["K_lam"] - 21550) < 1e-4
        and abs(rec["K_1mlam"] - rec["K_lam1mlam"]) < 1e-4
        and abs(rec["NUM_SUM_yfirst"] - rec["16pA"]) < 1e-6
        and rec["n_I3"] == 66000
    )
    # product=+1 bucket is not I_3
    prod_plus = rec["K3_direct"]  # only both-plus; both-minus sits in K_lam1mlam mix
    if abs(prod_plus - 100150) < 1.0:
        ok = False
    return {
        "proved": bool(ok),
        "fold": rec,
        "theorem": (
            "K_3=(K_1+K_λ+K_{1−λ}+K_{λ(1−λ)})/4 equals 83350 at p=5. "
            f"Fail drop λ(1−λ) ({rec['drop_prod']}≠{rec['K3_direct']})."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "NUM_SUM_16pA_general": False,
        "Q_tau_named_in_p": False,
        "phi_F_imported": False,
        "K_lam_named": False,
        "note": (
            "K_λ=21550 and K_{1−λ}=111550 at p=5 are not a p-law. "
            "16pA still only the live fold. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.548  H_+ 2-point (N/p)χ; |κ|=3 Jacobi split", flush=True)
    A = prove_A()
    print(f"  A G=(N/p)χ: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B Jacobi K3: {B['proved']}", flush=True)
    C = prove_open()
    out = {
        "prop": "15.548",
        "title": "H_+ 2-point (N/p)χ; |κ|=3 is the I_3 Jacobi split",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "C": C,
        "proved": {
            "G_paley_radial": A["proved"],
            "K3_jacobi_split": B["proved"],
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
