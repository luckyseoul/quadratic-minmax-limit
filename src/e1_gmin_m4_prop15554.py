#!/usr/bin/env python3
"""
Prop 15.554 — ŷ k-point on Ω is supported on Σ α_i=0, and
K_λ^{all}=Term0+Cross+ΩΩ.  Q3 is constant on 3-line triples
and on each 1-line class (double vs generic).  Q4≠Wick.
Ω-bulk values still not a p-law.  16pA / Q_τ / phi_F open.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
residual_ii / type_I / 15.279–15.553 flags.  Does **not** import phi_F.

============================================================================
Setup.  Translations x↦x+b preserve H_+.  ŷ^{(b)}(α)=e_α(b)ŷ(α).
15.550 W=q^{-1}[pJ+Σ_Ω ŷ J(α−η)].  15.553 names Term0 (ŷ(0) piece).
15.551: Ω-triple is 1-line iff r=η/ξ∈F_p, else 3-line.  A double
is a 1-line triple with a repeated frequency: r∈{1,−2,−1/2}.

============================================================================
Theorem A — PROVED (translation Aut of H_+; certified p=5,7).
  For any k≥1 and α_1,…,α_k∈F_q,
      Σ_{y∈H_+} ∏_i ŷ(α_i) = 0 unless Σ_i α_i = 0.
  Fail: Q3 on a nont-triple equals Q3 on a triple.  ∎

Theorem B — PROVED (A + cache; p=5,7).
  On Ω, Q2(α,−α)=2Nq and Q2=0 off antipodes.
  Q3 lives only on additive triples.  It is constant on:
    3-line (p=5: −25000; p=7: −672280),
    1-line doubles r∈{1,−2,−1/2} (p=5: −5000; p=7: −441784),
    generic 1-line (p=7: −1248520; empty at p=5).
  Fail: Q3≡0 on triples.  Fail: 3-line = double.  ∎

Theorem C — OPEN.  K_all=Term0+Cross+ΩΩ
  (p=5: 650+7040+241360=249050; p=7: 578326+1040032+40429424).
  Fail: Wick-ΩΩ as the bulk (p=5 WickOO=31839600≠248400).
  Q3 constants and Q4 are not a p-law.  16pA / Q_τ unnamed.

============================================================================
Backend: p=5,7 Max+ caches + field tables.  Writes
evidence/e1_gmin_m4_prop15554.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import defaultdict
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15279 import _kernel_field
from e1_gmin_m4_prop15550 import _e_tab, omega_set_local, pp_r

EV = ROOT / "evidence" / "e1_gmin_m4_prop15554.json"
YPATH = {5: "/tmp/maxplus_p5.npy", 7: "/tmp/maxplus_p7.npy"}

LIVE_Q3 = {
    5: {"double": -5000.0, "generic": None, "n3": -25000.0},
    7: {"double": -441784.0, "generic": -1248520.0, "n3": -672280.0},
}
LIVE_SPLIT = {
    5: {"Term0": 650.0, "Cross": 7040.0, "OO": 241360.0, "Kall": 249050.0},
    7: {"Term0": 578326.0, "Cross": 1040032.0, "OO": 40429424.0, "Kall": 42047782.0},
}


def load_Y(p: int, q: int) -> np.ndarray:
    raw = np.load(YPATH[p])
    block = raw[raw[:, 0] > 0][:, 1 : 1 + q].astype(np.float64)
    return np.where(np.sign(block) < 0, -1.0, 1.0)


def triple_kind(p: int, r: int) -> str:
    """double / generic / n3 / degenerate."""
    if r < p:
        if (1 + r) % p == 0:
            return "degenerate"
        half = pow(2, p - 2, p)  # 1/2 in F_p
        neg2 = (-2) % p
        neghalf = (-half) % p
        if r % p in (1, neg2, neghalf):
            return "double"
        return "generic"
    return "n3"


def yhat_omega(p: int):
    F = _kernel_field(p)
    q = F["q"]
    Y = load_Y(p, q)
    Omega = omega_set_local(F)
    mul = F["mul"]
    e_tab = np.stack(
        [np.array([_e_tab(F)[mul[om][t]] for t in range(q)], dtype=np.complex128) for om in Omega]
    )
    yhat = Y @ e_tab.T
    return F, Y, Omega, yhat


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        F, Y, Omega, yhat = yhat_omega(p)
        add, neg = F["add"], F["neg"]
        idx = {om: i for i, om in enumerate(Omega)}
        nO = len(Omega)
        # k=1: Σ ŷ(α) = 0 on Ω
        m1 = max(abs(float(np.sum(yhat[:, i]).real)) for i in range(nO))
        # k=2 off hyperplane
        m2_off = 0.0
        m2_on = 0.0
        want2 = 2 * Y.shape[0] * F["q"]
        for i, a in enumerate(Omega):
            for j, b in enumerate(Omega):
                val = abs(float(np.sum(yhat[:, i] * yhat[:, j]).real))
                if add[a][b] == 0:
                    m2_on = max(m2_on, abs(val - want2))
                else:
                    m2_off = max(m2_off, val)
        # k=3: nont vs triple
        q3_nont = 0.0
        q3_trip = 0.0
        n_nont = n_trip = 0
        for i, a in enumerate(Omega):
            for j, b in enumerate(Omega):
                c = neg[add[a][b]]
                s = yhat[:, i] * yhat[:, j]
                if c in idx:
                    k = idx[c]
                    val = abs(float(np.sum(s * yhat[:, k]).real))
                    q3_trip = max(q3_trip, val)
                    n_trip += 1
                # nont: pick a frequency not completing the sum
                k2 = (i + j + 1) % nO
                if Omega[k2] != c:
                    val = abs(float(np.sum(s * yhat[:, k2]).real))
                    q3_nont = max(q3_nont, val)
                    n_nont += 1
        row_ok = (
            m1 < 1e-6
            and m2_off < 1e-6
            and m2_on < 1e-4
            and q3_nont < 1e-4
            and q3_trip > 1.0
        )
        # fail-when-wrong: nont == trip
        fail_ok = abs(q3_nont - q3_trip) > 1.0
        if not (row_ok and fail_ok):
            ok = False
        rows[str(p)] = {
            "nOmega": nO,
            "m1_max": m1,
            "m2_off_max": m2_off,
            "m2_on_err": m2_on,
            "q3_nont_max": q3_nont,
            "q3_trip_max": q3_trip,
            "n_trip": n_trip,
            "ok": bool(row_ok and fail_ok),
        }
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": "Σ ∏ŷ(α_i)=0 unless Σα_i=0. Fail Q3(nont)=Q3(triple).",
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        F, Y, Omega, yhat = yhat_omega(p)
        add, neg, mul, inv = F["add"], F["neg"], F["mul"], F["inv"]
        idx = {om: i for i, om in enumerate(Omega)}
        N, q = Y.shape[0], F["q"]
        want2 = 2 * N * q
        buckets: dict[str, list[float]] = defaultdict(list)
        for i, a in enumerate(Omega):
            for j, b in enumerate(Omega):
                c = neg[add[a][b]]
                if c not in idx:
                    continue
                k = idx[c]
                r = mul[b][inv[a]]
                kind = triple_kind(p, r)
                if kind == "degenerate":
                    continue
                q3 = float(np.sum(yhat[:, i] * yhat[:, j] * yhat[:, k]).real)
                buckets[kind].append(q3)
        rec = {}
        live = LIVE_Q3[p]
        for kind, want in (("double", live["double"]), ("n3", live["n3"]), ("generic", live["generic"])):
            vals = buckets.get(kind, [])
            if want is None:
                rec[kind] = {"n": len(vals), "ok": len(vals) == 0}
                if vals:
                    ok = False
                continue
            if not vals:
                rec[kind] = {"n": 0, "ok": False}
                ok = False
                continue
            arr = np.array(vals)
            kind_ok = (
                abs(arr.max() - arr.min()) < 1e-4
                and abs(arr.mean() - want) < 1e-4
            )
            rec[kind] = {
                "n": int(len(vals)),
                "mean": float(arr.mean()),
                "min": float(arr.min()),
                "max": float(arr.max()),
                "want": want,
                "ok": bool(kind_ok),
            }
            if not kind_ok:
                ok = False
        # fail 3-line = double
        if rec.get("n3", {}).get("ok") and rec.get("double", {}).get("ok"):
            if abs(rec["n3"]["mean"] - rec["double"]["mean"]) < 1.0:
                ok = False
        # Q2
        q2_ok = True
        for i, a in enumerate(Omega):
            got = float(np.sum(yhat[:, i] * yhat[:, idx[neg[a]]]).real)
            if abs(got - want2) > 1e-4:
                q2_ok = False
        rec["Q2"] = {"want": float(want2), "ok": bool(q2_ok)}
        if not q2_ok:
            ok = False
        rec["ok"] = bool(
            rec["Q2"]["ok"]
            and rec.get("double", {}).get("ok", False)
            and rec.get("n3", {}).get("ok", False)
            and rec.get("generic", {}).get("ok", False)
        )
        if not rec["ok"]:
            ok = False
        rows[str(p)] = rec
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "Q2=2Nq on antipodes. Q3 constant on 3-line / double / "
            "generic 1-line. Fail Q3≡0; fail 3-line=double."
        ),
    }


def split_Kall(p: int) -> dict:
    F, Y, Omega, yhat = yhat_omega(p)
    q = F["q"]
    N = Y.shape[0]
    xi = int(Omega[0])
    eta = int(F["mul"][pp_r(F)][xi])
    add = np.array(F["add"], dtype=np.int32)
    neg = np.array(F["neg"], dtype=np.int32)
    mul = np.array(F["mul"], dtype=np.int32)
    chi = np.array(F["chi"], dtype=np.float64)
    e1 = _e_tab(F)
    e_all = e1[mul]
    ta = add[:, neg]
    prod = mul[ta[:, :, None], ta[:, None, :]]
    Q = chi[prod].transpose(1, 2, 0)
    J = np.einsum("abt,gt->gab", Q, e_all)
    eta_m = int(neg[eta])
    Jpos = np.stack([J[int(add[eta, neg[int(om)]])] for om in Omega], axis=0)
    Jneg = np.stack([J[int(add[eta_m, neg[int(om)]])] for om in Omega], axis=0)
    Wpos = np.tensordot(yhat, Jpos, axes=(1, 0))
    Wneg = np.tensordot(yhat, Jneg, axes=(1, 0))
    invq2 = 1.0 / (q * q)
    phase = e_all[xi][:, None] * np.conj(e_all[xi])[None, :]
    YY = Y[:, :, None] * Y[:, None, :]
    cross = (p * invq2) * (J[eta][None] * Wneg + J[eta_m][None] * Wpos)
    oo = invq2 * Wpos * Wneg
    Cross = np.einsum("nab,nab,ab->", YY, cross, phase)
    OO = np.einsum("nab,nab,ab->", YY, oo, phase)
    G = Y.T @ Y
    term0 = (p * p * invq2) * np.einsum("ab,ab,ab->", G, phase, J[eta] * J[eta_m])
    wick_oo = 0j
    want2 = 2 * N * q
    for i, om in enumerate(Omega):
        Jm = J[int(add[eta_m, int(om)])]
        wick_oo += np.einsum("ab,ab,ab->", G, phase, Jpos[i] * Jm)
    wick_oo *= invq2 * want2
    return {
        "N": int(N),
        "Term0": float(np.real(term0)),
        "Cross": float(np.real(Cross)),
        "OO": float(np.real(OO)),
        "Kall": float(np.real(term0 + Cross + OO)),
        "WickOO": float(np.real(wick_oo)),
        "bulk": float(np.real(Cross + OO)),
    }


def prove_open() -> dict:
    splits = {}
    ok_split = True
    for p, want in LIVE_SPLIT.items():
        rec = split_Kall(p)
        match = (
            abs(rec["Term0"] - want["Term0"]) < 1e-2
            and abs(rec["Cross"] - want["Cross"]) < 1e-2
            and abs(rec["OO"] - want["OO"]) < 1e-2
            and abs(rec["Kall"] - want["Kall"]) < 1e-2
        )
        wick_fails = abs(rec["WickOO"] - rec["bulk"]) > 1.0
        rec["match"] = bool(match)
        rec["wick_fails"] = bool(wick_fails)
        if not (match and wick_fails):
            ok_split = False
        splits[str(p)] = rec
    return {
        "proved": False,
        "split_certified": bool(ok_split),
        "K_all_named_in_p": False,
        "Q3_named_in_p": False,
        "Q4_is_Wick": False,
        "NUM_SUM_16pA_general": False,
        "Q_tau_named_in_p": False,
        "phi_F_imported": False,
        "splits": splits,
        "note": (
            "K_all=Term0+Cross+ΩΩ certified p=5,7. Q3 constants live "
            "but not a p-law. Q4≠Wick. 16pA / Q_τ unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.554  ŷ k-point support; Q3 by line type; Ω-bulk split", flush=True)
    A = prove_A()
    print(f"  A support: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B Q3 classes: {B['proved']}", flush=True)
    C = prove_open()
    print(f"  C split certified: {C['split_certified']}", flush=True)
    out = {
        "prop": "15.554",
        "title": "ŷ k-point supported on Σα=0; Q3 constant on line types; K_all=Term0+Cross+ΩΩ",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "C": C,
        "proved": {
            "ypoint_support": A["proved"],
            "Q3_constant_on_types": B["proved"],
            "K_all_split": C["split_certified"],
            "K_all_named_in_p": False,
            "NUM_SUM_16pA_general": False,
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
    }
    EV.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print("wrote", EV, flush=True)
    return out


if __name__ == "__main__":
    main()
