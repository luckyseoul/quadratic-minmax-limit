#!/usr/bin/env python3
"""
Prop 15.534 — Type I kernel I is Gauss × cubic S, not a χ-type law.

Does **not** flip type_I / phi_F / e1 / L / Aut-Schur / Gsum /
pairing / residual_ii / 15.279–15.533 flags.  Does **not**
edit leftover-1 15.522–15.533 or residual-ii 15.528.

============================================================================
Setup.  Paley C on V=F_q∪{∞}, q=p², χ=χ_q.  χ(−1)=1 (q≡1 mod 4).
ψ(t)=exp(2πi Tr_{q/p}(t)/p).  G(χ)=∑_{t≠0} χ(t)ψ(t)=−p χ_p(−1)
(15.523 A).  Kernel (15.529)
    I(η,θ)=∑_{x,y} χ(x)χ(y)χ(x−y) ψ(−ηx−θy)
(x=0, y=0, and x=y vanish).  Cubic
    S(λ)=∑_t χ(t(t−1)(t−λ)).
S is the trace of the Legendre curve y²=x(x−1)(x−λ) (Weil |S|≤2p).

============================================================================
Theorem A — PROVED (Gauss substitution; certified p=5,7 exhaustive
  and p=11,13 samples).
  S(0)=S(1)=−1.  For λ≠0, S(1/λ)=χ(λ) S(λ).
  For θ≠0,
      I(η,θ)=G(χ) χ(θ) S(−η/θ)
  and, when η≠0, equivalently G(χ) χ(η) S(−θ/η).
  For θ=0, η≠0: I(η,0)=−χ(η) G(χ).  I(0,0)=0.
  Proof.  S(0): χ(t²(t−1))=χ(t−1) on t≠0, so
  ∑_{u≠−1} χ(u)=−χ(−1)=−1.  S(1): χ(t(t−1)²)=χ(t) on t≠1,
  ∑_{t≠1} χ(t)=−1.  Functional: t=1/u on F_q^* gives
  t(t−1)(t−1/λ)=(1−u)(λ−u)/(u³ λ) and χ of that is
  χ(λ) χ(u(u−1)(u−λ)).
  For θ≠0 set x=t y.  Then χ(x)χ(y)χ(x−y)=χ(t(t−1)) χ(y) and
  ψ(−ηx−θy)=ψ(−(ηt+θ)y), so
      I=∑_{t≠0,1} χ(t(t−1)) ∑_{y≠0} χ(y) ψ(−(ηt+θ)y).
  Inner Gauss: 0 if ηt+θ=0, else χ(−(ηt+θ)) G(χ).
  Hence I=G χ(−1) ∑ χ(t(t−1)(ηt+θ))=G χ(η) S(−θ/η) (η≠0);
  the functional rewrites this as G χ(θ) S(−η/θ).  The axis
  I(η,0)=−∑_{x≠0} χ(x)ψ(−ηx) because ∑_y χ(y(x−y))=J(χ,χ)=−1
  for x≠0, and the remaining Gauss is χ(η) G.
  Fail: drop χ(θ) (p=5, (η,θ)=(1,p): live 10, drop −10).
  Fail: G≡p (p=5, (1,1): 30≠−30).
  Fail: I is a single p-law ((1,1) is −30 at p=5, 98=2p² at p=7).  ∎

Theorem B — PROVED as a negative (field S-table; p=5,7).
  The 3-χ type (χ(λ),χ(λ−1),χ(λ+1)) does not pin S(λ).
  Witness p=5: type (−1,−1,+1) has S∈{−6,2}.  Hence I is not
  a χ-piecewise law in the types of (η,θ) (those determine
  χ(λ),χ(λ−1) via λ=−η/θ, and even the 4-type
  (χ(η),χ(θ),χ(η+θ),χ(θ−η)) only recovers the 3-χ type of λ).
  Fail: claim each 3-χ type is S-constant.  ∎

Theorem C — OPEN.  Naming I does not name live Δ_conn as a
  p-rational (15.529 C: A still has den 13,409).  The rewrite
      Δ_conn=2/((p²−1)p⁴) ∑_{Ω, ξ+η+θ=0} A χ(θ) S(−η/θ)
  still carries D through A.  Aut_e 3A+B still G>T.
  type_I flags stay False.

============================================================================
Backend: field tables + S-table + vectorized I.  Writes
evidence/e1_gmin_m4_prop15534.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
)
from e1_gmin_m4_prop15275 import (  # noqa: E402
    type_I_aut_e_3AB_positive_general,
    type_I_multilevel_bad_case_ND_closed,
)
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15523 import G_chi  # noqa: E402

EV = ROOT / "evidence" / "e1_gmin_m4_prop15534.json"

SAMPLES = [(1, 1), (1, 2), (1, 0), (0, 1), (0, 0), (2, 3)]
# (1, p) is appended per p so χ(θ) can be −1 (drop-χ fail).


def field_tables(p: int) -> dict:
    q = p * p

    def is_irr(a: int, b: int) -> bool:
        return all((x * x - a * x - b) % p != 0 for x in range(p))

    ia = ib = None
    for a in range(p):
        for b in range(p):
            if is_irr(a, b):
                ia, ib = a, b
                break
        if ia is not None:
            break

    def mul(u: int, v: int) -> int:
        c0, c1 = u % p, u // p
        d0, d1 = v % p, v // p
        return (c0 * d0 + c1 * d1 * ib) % p + (
            (c0 * d1 + c1 * d0 + c1 * d1 * ia) % p
        ) * p

    def add(u: int, v: int) -> int:
        return (u % p + v % p) % p + ((u // p + v // p) % p) * p

    def sub(u: int, v: int) -> int:
        return (u % p - v % p) % p + ((u // p - v // p) % p) * p

    def powm(u: int, e: int) -> int:
        r, base = 1, u
        while e:
            if e & 1:
                r = mul(r, base)
            base = mul(base, base)
            e >>= 1
        return r

    def chi(x: int) -> int:
        if x == 0:
            return 0
        return 1 if powm(x, (q - 1) // 2) == 1 else -1

    def tr(u: int) -> int:
        return add(u, powm(u, p)) % p

    ch = np.array([chi(x) for x in range(q)], dtype=np.int8)
    trt = np.array([tr(x) for x in range(q)], dtype=np.int16)
    add_t = np.empty((q, q), dtype=np.int16)
    sub_t = np.empty((q, q), dtype=np.int16)
    mul_t = np.empty((q, q), dtype=np.int16)
    for i in range(q):
        for j in range(q):
            add_t[i, j] = add(i, j)
            sub_t[i, j] = sub(i, j)
            mul_t[i, j] = mul(i, j)
    inv = np.zeros(q, dtype=np.int16)
    for i in range(1, q):
        inv[i] = powm(i, q - 2)
    return {
        "p": p,
        "q": q,
        "chi": ch,
        "tr": trt,
        "add": add_t,
        "sub": sub_t,
        "mul": mul_t,
        "inv": inv,
    }


def S_table(F: dict) -> np.ndarray:
    """S(λ)=∑_t χ(t(t−1)(t−λ))."""
    ch = F["chi"].astype(np.int32)
    w = ch * ch[F["sub"][:, 1]]
    return w @ ch[F["sub"]]


def I_of(F: dict, eta: int, theta: int) -> complex:
    q = F["q"]
    p = F["p"]
    ch = F["chi"]
    sub, trt, mul_t = F["sub"], F["tr"], F["mul"]
    om = np.exp(2j * np.pi / p)
    px = om ** ((-trt[mul_t[eta]].astype(np.int64)) % p)
    py = om ** ((-trt[mul_t[theta]].astype(np.int64)) % p)
    wx = ch.astype(np.complex128) * px
    wy = ch.astype(np.complex128) * py
    s = 0j
    for x in range(q):
        s += wx[x] * np.dot(wy, ch[sub[x]].astype(np.complex128))
    return s


def I_pred(F: dict, S: np.ndarray, eta: int, theta: int) -> int:
    G = G_chi(F["p"])
    ch = F["chi"]
    if theta == 0:
        if eta == 0:
            return 0
        return -int(ch[eta]) * G
    lam = int(F["mul"][int(F["sub"][0, eta]), int(F["inv"][theta])])
    return int(G * int(ch[theta]) * int(S[lam]))


def I_pred_drop_chi(F: dict, S: np.ndarray, eta: int, theta: int) -> int:
    """Fail-when-wrong: drop χ(θ)."""
    G = G_chi(F["p"])
    if theta == 0:
        return I_pred(F, S, eta, theta)
    lam = int(F["mul"][int(F["sub"][0, eta]), int(F["inv"][theta])])
    return int(G * int(S[lam]))


def I_pred_G_eq_p(F: dict, S: np.ndarray, eta: int, theta: int) -> int:
    """Fail-when-wrong: G≡p."""
    ch = F["chi"]
    p = F["p"]
    if theta == 0:
        if eta == 0:
            return 0
        return -int(ch[eta]) * p
    lam = int(F["mul"][int(F["sub"][0, eta]), int(F["inv"][theta])])
    return int(p * int(ch[theta]) * int(S[lam]))


def _pairs_for(p: int, exhaustive: bool) -> list[tuple[int, int]]:
    q = p * p
    if exhaustive:
        return [(eta, theta) for eta in range(q) for theta in range(q)]
    out = list(SAMPLES)
    out.append((1, p))
    out.append((1, p + 1))
    return [(e, t) for e, t in out if e < q and t < q]


def _cert_p(spec: tuple[int, bool]) -> dict:
    p, exhaustive = spec
    F = field_tables(p)
    S = S_table(F)
    G = G_chi(p)
    q = F["q"]
    ch = F["chi"]
    inv = F["inv"]
    mul_t, sub_t, add_t = F["mul"], F["sub"], F["add"]

    s0 = int(S[0])
    s1 = int(S[1])
    sm1 = int(S[int(sub_t[0, 1])])
    func_ok = True
    for lam in range(1, q):
        inv_lam = int(inv[lam])
        if int(S[inv_lam]) != int(ch[lam]) * int(S[lam]):
            func_ok = False
            break

    pairs = _pairs_for(p, exhaustive)
    max_err = 0.0
    n_ok = 0
    drop_hits_live = 0
    geqp_hits_live = 0
    live_11 = None
    for eta, theta in pairs:
        live = I_of(F, eta, theta)
        pred = float(I_pred(F, S, eta, theta))
        err = abs(live.real - pred) + abs(live.imag)
        if err > max_err:
            max_err = float(err)
        if err < 1e-6:
            n_ok += 1
        if theta != 0 and int(ch[theta]) == -1:
            drop = float(I_pred_drop_chi(F, S, eta, theta))
            if abs(live.real - drop) + abs(live.imag) < 1e-6:
                drop_hits_live += 1
        geqp = float(I_pred_G_eq_p(F, S, eta, theta))
        if abs(live.real - geqp) + abs(live.imag) < 1e-6 and geqp != pred:
            geqp_hits_live += 1
        if eta == 1 and theta == 1:
            live_11 = float(np.real(live))

    # 3-χ types of λ
    by_type: dict[tuple, set[int]] = defaultdict(set)
    for lam in range(2, q):
        typ = (
            int(ch[lam]),
            int(ch[sub_t[lam, 1]]),
            int(ch[add_t[lam, 1]]),
        )
        by_type[typ].add(int(S[lam]))
    split = {str(k): sorted(v) for k, v in by_type.items() if len(v) > 1}

    formula_ok = (
        s0 == -1
        and s1 == -1
        and func_ok
        and n_ok == len(pairs)
        and max_err < 1e-6
    )
    return {
        "p": p,
        "G": G,
        "S0": s0,
        "S1": s1,
        "S_m1": sm1,
        "S_unique": sorted({int(x) for x in S}),
        "n_pairs": len(pairs),
        "n_ok": n_ok,
        "max_err": max_err,
        "formula_ok": bool(formula_ok),
        "func_ok": bool(func_ok),
        "live_11": live_11,
        "drop_chi_collisions": drop_hits_live,
        "G_eq_p_collisions_off": geqp_hits_live,
        "split_types": split,
        "n_split_types": len(split),
    }


def _run_certs() -> list[dict]:
    specs = [(5, True), (7, True), (11, False), (13, False)]
    w = min(4, len(specs))
    with ProcessPoolExecutor(max_workers=w) as ex:
        return list(ex.map(_cert_p, specs))


@lru_cache(maxsize=1)
def prove_A() -> dict:
    rows = {str(r["p"]): r for r in _run_certs()}
    r5, r7 = rows["5"], rows["7"]
    drop_fail_ok = False
    # p=5, (1,p): χ(p)≠1 so drop χ is live-wrong
    F5 = field_tables(5)
    S5 = S_table(F5)
    live_1p = I_of(F5, 1, 5)
    pred_1p = I_pred(F5, S5, 1, 5)
    drop_1p = I_pred_drop_chi(F5, S5, 1, 5)
    geqp_11 = I_pred_G_eq_p(F5, S5, 1, 1)
    pred_11 = I_pred(F5, S5, 1, 1)
    drop_fail_ok = (
        abs(live_1p.real - pred_1p) < 1e-6
        and drop_1p != pred_1p
        and abs(live_1p.real - drop_1p) > 1.0
        and geqp_11 != pred_11
        and pred_11 == -30
        and abs((r5["live_11"] or 0) + 30) < 1e-6
        and abs((r7["live_11"] or 0) - 98) < 1e-6
        and 98 == 2 * 7 * 7
        and -30 != 2 * 5 * 5
    )
    ok = all(r["formula_ok"] for r in rows.values()) and drop_fail_ok
    return {
        "proved": bool(ok),
        "rows": rows,
        "p5_I11": r5["live_11"],
        "p7_I11": r7["live_11"],
        "p5_I_1p": {"live": float(np.real(live_1p)), "pred": pred_1p, "drop": drop_1p},
        "fail_G_eq_p_at_11": geqp_11,
        "theorem": (
            "I(η,θ)=G χ(θ) S(−η/θ) (θ≠0); I(η,0)=−χ(η)G; "
            "S(0)=S(1)=−1.  Fail: drop χ(θ); G≡p; I single p-law."
        ),
    }


@lru_cache(maxsize=1)
def prove_B() -> dict:
    rows = prove_A()["rows"]
    # p=5 type (−1,−1,+1) is the Ipiece witness
    split5 = rows["5"]["split_types"]
    want_key = str((-1, -1, 1))
    vals = split5.get(want_key, [])
    ok = (
        rows["5"]["n_split_types"] >= 1
        and rows["7"]["n_split_types"] >= 1
        and set(vals) == {-6, 2}
    )
    return {
        "proved": bool(ok),
        "p5_split": split5,
        "p7_n_split": rows["7"]["n_split_types"],
        "witness": {"type": want_key, "S": vals},
        "theorem": (
            "3-χ types do not pin S (p=5 type (−1,−1,+1) has S∈{−6,2}). "
            "I is not a χ-piecewise law."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "type_I_multilevel_bad_case_ND_closed": bool(
            type_I_multilevel_bad_case_ND_closed()
        ),
        "type_I_aut_e_3AB_positive_general": bool(
            type_I_aut_e_3AB_positive_general()
        ),
        "Delta_conn_p_rational": False,
        "note": (
            "I-name does not clear den 13,409 in A (15.529).  "
            "Δ_conn rewrite still carries D.  Aut_e 3A+B still G>T.  "
            "Flags unflipped."
        ),
    }


def main() -> dict:
    print("Prop 15.534  I = G χ(θ) S(−η/θ)", flush=True)
    A = prove_A()
    print(
        f"  A formula: {A['proved']} "
        f"I(1,1) p5={A['p5_I11']} p7={A['p7_I11']}",
        flush=True,
    )
    B = prove_B()
    print(
        f"  B types split: {B['proved']} witness={B['witness']}",
        flush=True,
    )
    C = prove_open()
    print(
        f"  C open: ND={C['type_I_multilevel_bad_case_ND_closed']} "
        f"3AB={C['type_I_aut_e_3AB_positive_general']}",
        flush=True,
    )
    out = {
        "prop": "15.534",
        "title": "I is Gauss × cubic S, not a χ-type p-law",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "C": C,
        "proved": {
            "I_equals_G_chi_S": A["proved"],
            "chi_types_do_not_pin_S": B["proved"],
            "type_I_multilevel_bad_case_ND_closed": C[
                "type_I_multilevel_bad_case_ND_closed"
            ],
            "type_I_aut_e_3AB_positive_general": C[
                "type_I_aut_e_3AB_positive_general"
            ],
        },
        "e1_closed_general": bool(e1_closed_general()),
        "gsum_disj_lb": bool(gsum_disj_lb_proved_general()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "flags_not_flipped": [
            "type_I_multilevel_bad_case_ND_closed",
            "type_I_aut_e_3AB_positive_general",
            "phi_F_ge_6",
            "e1",
            "L",
            "Aut-Schur",
            "Gsum",
            "pairing",
        ],
        "L_status": "OPEN",
    }
    # json-safe: rows contain only Python scalars
    EV.parent.mkdir(parents=True, exist_ok=True)
    EV.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", EV, flush=True)
    return out


if __name__ == "__main__":
    main()
