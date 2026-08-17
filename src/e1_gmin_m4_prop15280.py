#!/usr/bin/env python3
"""
Prop 15.280 — Weil expansion of N̂ on DΔT; circle |N̂_T|²; Hoffman
support.  Aut_∞ mixture / floor still OPEN (Cross not evaluated;
N+(p) unnamed).  Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur /
Gsum / pairing / 15.279 inequality.

============================================================================
Setup.  Even multiplicative ψ_k of F_q^×, q=p², ψ(0)=0.
N̂_S(ψ)=∑_{x,y∈S} ψ(y−x)=∑_{δ≠0} N_S(δ) ψ(δ).
Cross_ψ(A,B)=∑_{x∈A,y∈B} ψ(y−x).  Regular sets: |D|=p(p−1)/2.
NL Type+ are D_A Δ T, T a Hoffman evec-support of C'=D_{y_A} C D_{y_A}
(15.138–15.141).  Floor: E|N̂|² ≥ 3q(q−1)/16 on even Aut_L-orbits.

============================================================================
Theorem A (DΔT expansion) — PROVED algebra, any D,T.
  I=D∩T.  For even ψ (so Cross(A,B)=Cross(B,A)),
      N̂_{DΔT} = N̂_D + N̂_T + 2 Cross(D,T)
               − 4 Cross(D,I) − 4 Cross(T,I) + 4 N̂_I.
  Equivalently, with K=D\\I and J=T\\I,
      N̂_{DΔT} = N̂_K + N̂_J + 2 Cross(K,J).
  Naive N̂_D+N̂_T±2 Cross holds iff I=∅ or T⊂D (or D⊂T).
  Fail-when-wrong: invert the Cross sign.  ∎

Theorem B (norm-circle N̂_T) — PROVED.
  H=ker(N:F_q^×→F_p^×), |H|=p+1.  Every S_{t,c} is a translate of
  a scale of H.  N̂ is translation-invariant; N̂_{λH}=ψ(λ) N̂_H.
  N̂_H=(∑_{x∈H} ψ(x)) (∑_{z∈H} ψ(z−1)).
  ∑_H ψ = (p+1) if (p+1)|k, else 0.
  Ratio count: r(c)=#{(z,w)∈(H\\{1})²:(z−1)/(w−1)=c} equals
      r(1)=p,  r(c)=0 for c∈F_p\\{1},  r(c)=1 for c∉F_p.
  Proof of the F_p-vanishing: (z−1)=c(w−1) with c∈F_p\\{0,1} and
  N(z)=N(w)=1 forces w+w^{-1}=2, hence w=1, excluded.
  Hence for even ψ∉{1,χ} with (p+1)|k (so (p-1)∤k, else ψ∈{1,χ}),
      |∑_{z∈H} ψ(z−1)|² = p + ∑_{c∉F_p} ψ(c) = p
  because ∑_{F_p^×} ψ=0.  Therefore
      |N̂_T|² = 0           if (p+1)∤k,
      |N̂_T|² = p(p+1)²     if (p+1)|k.
  Fail-when-wrong: invert the (p+1)|k slot; put r=1 on F_p\\{1}.  ∎

Theorem C (Hoffman constraints) — PROVED structure.
  T a coclique of C'=D_{y_A} C D_{y_A} (no ∞ in T) ⇒
  χ(b−a)=y_A(a) y_A(b) for a≠b in T.
  |DΔT|=|D| ⇒ |I|=|T|/2.  Then N_I, N_J are supported on squares
  and C_{I,J} on nonsquares.
  Fail-when-wrong: non-Hoffman T leaks N_I onto nonsquares.
  Norm-circle Hoffman T of C_hs is empty at p=13 (15.138).  ∎

Theorem D (OPEN mixture / floor).
  Cross(D,T) is not a linear function of (N̂_D, N̂_T, N̂_I) for
  general p (lstsq fails at p=7,11).  N_I=(N̂_T+N̂_{χψ})/4 is not
  Hoffman-automatic (fails on some p=7 types and on non-Hoffman T).
  N+(p) unnamed (6,130,5726 at p=3,5,7; 409 prime; type-list
  ∑|G|/stab is 15.144-dead for p≥11).  1D-only fails the floor at
  p=7 k=8 (60.38<441).  Per-orbit positivity is false (NL orbits
  with N̂=0).  phi_F_ge_6 stays False.  ∎

============================================================================
Backend: CPU field arithmetic.  GPU unused.
Writes evidence/e1_gmin_m4_prop15280.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15138 import YSTAR_CERT, construct_ystar  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field, _psi_vec  # noqa: E402


def even_ks(p: int) -> list[int]:
    q = p * p
    return [k for k in range(2, q - 1, 2) if k != (q - 1) // 2]


def _add_neg(F: dict):
    q = F["q"]
    add = np.array(F["add"], dtype=np.int32)
    neg = np.array(F["neg"], dtype=np.int32)
    return add, neg, q


def N_vec(d: np.ndarray, add: np.ndarray) -> np.ndarray:
    q = d.shape[0]
    return np.array(
        [int((d * d[add[:, delta]]).sum()) for delta in range(q)], dtype=np.int32
    )


def nhat(d: np.ndarray, add: np.ndarray, psi: np.ndarray) -> complex:
    N = N_vec(d, add)
    return complex(np.dot(N[1:].astype(np.complex128), psi[1:]))


def cross(a: np.ndarray, b: np.ndarray, psi: np.ndarray, add: np.ndarray, neg: np.ndarray) -> complex:
    ia = np.flatnonzero(a)
    ib = np.flatnonzero(b)
    if ia.size == 0 or ib.size == 0:
        return 0j
    return complex(psi[add[np.ix_(ib, neg[ia])]].sum())


def expansion_rhs(d, t, add, neg, psi, sign: float = 1.0) -> complex:
    I = d & t
    return (
        nhat(d, add, psi)
        + nhat(t, add, psi)
        + (2.0 * sign) * cross(d, t, psi, add, neg)
        - 4.0 * cross(d, I, psi, add, neg)
        - 4.0 * cross(t, I, psi, add, neg)
        + 4.0 * nhat(I, add, psi)
    )


def unit_circle(F: dict) -> list[int]:
    q = F["q"]
    mul = F["mul"]
    # N(u)=u^{p+1}=1 ⇔ u ≠ 0 and the field norm is 1
    # Use chi of the F_p-valued norm via u * u^p
    p = F["p"]

    def frob(u: int) -> int:
        # u^p via repeated square
        r, b, e = 1, u, p
        while e:
            if e & 1:
                r = mul[r][b]
            b = mul[b][b]
            e >>= 1
        return r

    pts = []
    for u in range(1, q):
        if mul[u][frob(u)] == 1:
            pts.append(u)
    return pts


def ratio_r(F: dict, H: list[int]) -> np.ndarray:
    q = F["q"]
    add, neg = np.array(F["add"], dtype=np.int32), np.array(F["neg"], dtype=np.int32)
    mul = F["mul"]
    inv = F["inv"]
    r = np.zeros(q, dtype=np.int32)
    one = 1
    for w in H:
        if w == one:
            continue
        wm1 = add[w, neg[one]]
        iw = inv[int(wm1)]
        for z in H:
            if z == one:
                continue
            zm1 = add[z, neg[one]]
            r[mul[int(zm1)][iw]] += 1
    return r


def prove_expansion(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [5, 7]
    rows = {}
    ok = True
    inv_hits = 0
    nchk = 0
    for p in primes:
        F = _kernel_field(p)
        add, neg, q = _add_neg(F)
        ks = even_ks(p)
        rng = np.random.default_rng(p)
        id_ok = inv_ok = 0
        n = 0
        for _ in range(4):
            d = (rng.random(q) < 0.4).astype(np.int8)
            t = (rng.random(q) < 0.15).astype(np.int8)
            dlt = np.bitwise_xor(d, t)
            for k in ks[:: max(1, len(ks) // 4)]:
                psi = np.array(_psi_vec(F, k), dtype=np.complex128)
                lhs = nhat(dlt, add, psi)
                rhs = expansion_rhs(d, t, add, neg, psi, 1.0)
                rhs_inv = expansion_rhs(d, t, add, neg, psi, -1.0)
                n += 1
                if abs(lhs - rhs) < 1e-6:
                    id_ok += 1
                if abs(lhs - rhs_inv) < 1e-6:
                    inv_ok += 1
        rows[str(p)] = {"ok": id_ok, "n": n, "inv_ok": inv_ok}
        if id_ok != n:
            ok = False
        inv_hits += inv_ok
        nchk += n
    return {
        "proved": ok,
        "theorem": (
            "N̂_{DΔT}=N̂_D+N̂_T+2 Cross(D,T)−4 Cross(D,I)−4 Cross(T,I)+4 N̂_I "
            "(even ψ, any D,T).  Invert-Cross is not an identity."
        ),
        "by_p": rows,
        "invert_cross_is_identity": inv_hits == nchk and nchk > 0,
        "invert_hits": inv_hits,
        "n_checks": nchk,
    }


def prove_circle_nhat(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [5, 7, 11, 13]
    rows = {}
    ok = True
    for p in primes:
        F = _kernel_field(p)
        add, neg, q = _add_neg(F)
        H = unit_circle(F)
        if len(H) != p + 1:
            ok = False
            rows[str(p)] = {"|H|": len(H), "bad": True}
            continue
        r = ratio_r(F, H)
        r_ok = int(r[1]) == p and int(r[0]) == 0
        # F_p is 0..p-1 in the a+b·p encoding
        for c in range(1, q):
            in_fp = c < p
            want = (p if c == 1 else 0) if in_fp else 1
            if int(r[c]) != want:
                r_ok = False
        tH = np.zeros(q, dtype=np.int8)
        tH[H] = 1
        van_ok = slot_ok = 0
        van_n = slot_n = 0
        inv_slot_fail = 0
        for k in even_ks(p):
            psi = np.array(_psi_vec(F, k), dtype=np.complex128)
            val = abs(nhat(tH, add, psi)) ** 2
            if k % (p + 1) == 0:
                slot_n += 1
                if abs(val - p * (p + 1) ** 2) < 1e-4:
                    slot_ok += 1
                if abs(val) < 1e-4:
                    inv_slot_fail += 1  # inverted vanishing would want 0 here
            else:
                van_n += 1
                if val < 1e-6:
                    van_ok += 1
        # inverted slot: claiming vanish ON (p+1)|k fails
        rows[str(p)] = {
            "|H|": len(H),
            "r_formula": r_ok,
            "vanish_off": f"{van_ok}/{van_n}",
            "slot_p(p+1)^2": f"{slot_ok}/{slot_n}",
            "pred": p * (p + 1) ** 2,
            "inverted_slot_vanishes": inv_slot_fail == slot_n and slot_n > 0,
        }
        if not (r_ok and van_ok == van_n and slot_ok == slot_n):
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "On every norm circle T, |N̂_T|²=0 unless (p+1)|k, then p(p+1)². "
            "r(c)=p 1_{c=1}+1_{c∉F_p}; |∑_H ψ(z−1)|²=p."
        ),
        "by_p": rows,
    }


def prove_hoffman_support() -> dict:
    """χ-support + |I|=|T|/2 on ystar (p=5,7); p=13 empty; leak on a non-H set."""
    from e1_gmin_m4_prop15139 import _field_ops, seidel_switch
    from minmax_quadratic import halfspace_boolean_vector, paley_conference_prime_power

    rows = {}
    ok = True
    for p in (5, 7):
        rec = construct_ystar(p)
        F = _kernel_field(p)
        add, neg, q = _add_neg(F)
        chi = np.array(F["chi"], dtype=np.int8)
        hs = halfspace_boolean_vector(p)
        d = (np.rint(hs[1:]) < 0).astype(np.int8)
        T = np.zeros(q, dtype=np.int8)
        for s in rec["S"]:
            T[s - 1] = 1
        I = d & T
        J = T & (1 - d)
        bal = 2 * int(I.sum()) == int(T.sum())
        NI, NJ = N_vec(I, add), N_vec(J, add)
        ia = np.flatnonzero(I)
        leak_I = leak_J = leak_C = 0
        for delta in range(1, q):
            cij = int(J[add[ia, delta]].sum()) if ia.size else 0
            if chi[delta] == -1:
                leak_I += int(NI[delta])
                leak_J += int(NJ[delta])
            else:
                leak_C += cij
        pair = True
        yf = np.rint(hs[1:]).astype(np.int8)
        idx = np.flatnonzero(T)
        for i, a in enumerate(idx):
            for b in idx[i + 1 :]:
                diff = add[int(b), int(neg[int(a)])]
                if int(chi[diff]) != int(yf[int(a)]) * int(yf[int(b)]):
                    pair = False
        rows[str(p)] = {
            "found": rec["found"],
            "balanced": bal,
            "pair_chi": pair,
            "leak_I_ns": leak_I,
            "leak_J_ns": leak_J,
            "leak_C_sq": leak_C,
        }
        if not (rec["found"] and bal and pair and leak_I == 0 and leak_J == 0 and leak_C == 0):
            ok = False

    # non-Hoffman: first norm circle of C_hs that fails coclique, p=5
    p = 5
    F = _kernel_field(p)
    add, neg, q = _add_neg(F)
    chi = np.array(F["chi"], dtype=np.int8)
    mul, addf, negf, frob, qq = _field_ops(p)
    C = paley_conference_prime_power(p).astype(np.float64)
    hs = halfspace_boolean_vector(p)
    Cp = seidel_switch(C, hs)
    d = (np.rint(hs[1:]) < 0).astype(np.int8)
    leaked = False
    n_non = 0
    for t in range(q):
        for c in range(1, p):
            pts = [u for u in range(q) if mul(addf(u, negf(t)), frob(addf(u, negf(t)))) == c]
            if len(pts) != p + 1:
                continue
            S = [1 + u for u in pts]
            if all(Cp[S[i], S[j]] > 0 for i in range(p + 1) for j in range(i + 1, p + 1)):
                continue
            n_non += 1
            T = np.zeros(q, dtype=np.int8)
            T[pts] = 1
            I = d & T
            NI = N_vec(I, add)
            if any(int(NI[delta]) and chi[delta] == -1 for delta in range(1, q)):
                leaked = True
                break
        if leaked:
            break
    p13_empty = YSTAR_CERT[13].get("found") is False
    return {
        "proved": ok and leaked and p13_empty,
        "theorem": (
            "Hoffman coclique ⇒ χ(b−a)=y_A(a)y_A(b), |I|=|T|/2, "
            "N_I,N_J on □, C_IJ on ⊠.  Non-Hoffman leaks.  p=13 nc empty."
        ),
        "ystar": rows,
        "nonhoffman_leaks_p5": leaked,
        "n_nonhoffman_scanned": n_non,
        "p13_norm_circle_empty": p13_empty,
    }


def prove_open() -> dict:
    return {
        "floor_closed": False,
        "cross_closed": False,
        "Nplus_named": False,
        "per_orbit_positive": False,
        "one_d_only_meets_floor": False,
        "open": (
            "Evaluate Cross_ψ(D,T) Max+-free on every Hoffman type; name N+(p); "
            "write the Aut_∞ mixture ≥ 3q(q−1)/16.  1D-only fails at p=7 k=8."
        ),
        "dead": [
            "1D-only mixture",
            "per-orbit positivity",
            "norm-circle T as the only T (empty at p=13)",
            "15.144 type-list for N+ at p≥11",
            "census deg≤3 through {3,5,7} as a formula for N+ (409 prime)",
        ],
    }


def main() -> dict:
    a = prove_expansion()
    b = prove_circle_nhat()
    c = prove_hoffman_support()
    d = prove_open()
    out = {
        "prop": "15.280",
        "title": (
            "Prop 15.280: N̂ on DΔT Weil expansion; circle |N̂_T|²=p(p+1)²; "
            "Hoffman support; mixture/floor OPEN"
        ),
        "L_status": "OPEN",
        "settlement_note": (
            "Named Max+-free identities A–C.  Cross and N+(p) unnamed; "
            "Aut_∞ mixture does not meet 3q(q−1)/16 without them.  "
            "No flag flip."
        ),
        "proved": {
            "d_symdiff_expansion": a["proved"],
            "circle_Nhat_eval": b["proved"],
            "hoffman_support": c["proved"],
            "invert_cross_is_identity": a["invert_cross_is_identity"],
            "cross_closed": False,
            "Nplus_named": False,
            "mixture_meets_floor": False,
            "phi_F_ge_6_proved_general": False,
            "inequality_proved": False,
        },
        "algebra": {"expansion": a, "circle": b, "hoffman": c, "open": d},
        "closed_forms": {
            "expansion": "N̂_{DΔT}=N̂_D+N̂_T+2X_DT-4X_DI-4X_TI+4 N̂_I",
            "circle": "|N̂_T|²=p(p+1)² if (p+1)|k else 0",
            "r": "r(1)=p, r|F_p\\{1}=0, r|F_q\\F_p=1",
        },
        "fail_when_wrong": [
            "invert Cross sign (not an identity)",
            "invert (p+1)|k vanishing of N̂_T",
            "r=1 on F_p\\{1} (kills |S|²=p)",
            "Hoffman support on non-Hoffman T",
            "only-norm-circle T at p=13 (empty)",
        ],
        "status": {
            "identities_A_B_C": bool(a["proved"] and b["proved"] and c["proved"]),
            "floor": False,
        },
        "backend": "CPU field + character sums; GPU unused",
        "residual_closed": False,
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15280.json"
    path.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"Prop 15.280 expansion={a['proved']} circle={b['proved']} hoffman={c['proved']}")
    print(f"  invert_cross_identity={a['invert_cross_is_identity']} floor OPEN")
    print(f"  wrote {path}")
    return out


if __name__ == "__main__":
    main()
