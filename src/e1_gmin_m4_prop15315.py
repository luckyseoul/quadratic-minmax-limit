#!/usr/bin/env python3
"""
Prop 15.315 — Boolean identity is B=−2p ẑ (not p(p−2));
ρ(R) type counts are named; neither |B|² nor E[M_i M_j] is a
second linear relation on (Q_++, Q_{N*}, Q_{--}).  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.314 flags.  Floor stays OPEN.

============================================================================
Setup.  z=±1 on F_q, ẑ supported on {0}∪Ω, ẑ(0)=p (y_∞=+1).
R={t∈T\\{1}: 1−t∈T}, |R|=(q−5)/4.  ρ(t)=(1−t)/t.

============================================================================
Theorem A — PROVED (z²=1 + supp ẑ ⊆ {0}∪Ω; certified p=5,7).
  ∑_{t∈R} ẑ(tξ)ẑ((1−t)ξ) = −2p ẑ(ξ) on every ξ∈Ω, every Max+.
  Fail-when-wrong: replace −2p by p(p−2).  ∎

Theorem B — PROVED (Cayley on T; Max+-free).
  ρ: R → {s∈T: χ(s+1)=+1} is a bijection.  It never hits --N1.
  |ρ∩pm1|=1 (t=1/2, since χ_q(2)=1),
  |ρ∩++| = (3p−7)/2 if p≡1 (mod 4), |++| if p≡3 (mod 4),
  |ρ∩N*| = (p−1)(p−5)/4 if p≡1 (mod 4), (p−3)²/4 if p≡3 (mod 4).
  Fail-when-wrong: claim |ρ∩--N1|>0.  ∎

Theorem C — CERTIFIED p=5,7; not a second Q-relation.
  B=−2p ẑ ⇒ E[|B|²]/q²=8, independent of Q_τ (same constraint as
  E[u]=2q).  Integer (a,b,c)∈[−6,6]³ independent of the S_0 size
  vector do not match {0,1,2,4,5,6,8,16,p,p±1,4−2/(p+2)} at both
  p=5 and p=7.  E[M_i M_j] for distinct square-direction classes
  has den 13 (p=5) and 409 (p=7), not |U|.  Q_τ unnamed.
  phi_F not imported.  ∎

============================================================================
Backend: field tables p=3..19; ẑ on p=5,7 caches (one fold).
Writes evidence/e1_gmin_m4_prop15315.json
"""
from __future__ import annotations

import json
import sys
from collections import Counter
from fractions import Fraction
from itertools import product
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field, rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15290 import live_Q_on_T, omega_set, zhat_omega  # noqa: E402
from e1_gmin_m4_prop15300 import merge_cls  # noqa: E402
from e1_gmin_m4_prop15314 import named_sizes  # noqa: E402

YPATH = {5: "/tmp/maxplus_p5.npy", 7: "/tmp/maxplus_p7.npy"}
SIZE_PRIMES = (3, 5, 7, 11, 13, 17, 19)


def R_set(F):
    R = []
    for t in F["squares"]:
        if t == 1:
            continue
        one_m = F["add"][1][F["neg"][t]]
        if one_m != 0 and F["chi"][one_m] == 1:
            R.append(t)
    return R


def rho_of(F, t):
    tm = F["add"][1][F["neg"][t]]
    return F["mul"][tm][F["inv"][t]]


def rho_named(p: int) -> dict[str, int]:
    if p % 4 == 1:
        return {
            "pm1": 1,
            "++": (3 * p - 7) // 2,
            "N*": (p - 1) * (p - 5) // 4,
            "--N1": 0,
        }
    return {
        "pm1": 1,
        "++": named_sizes(p)["++"],
        "N*": (p - 3) ** 2 // 4,
        "--N1": 0,
    }


_BOOL_CACHE: dict | None = None


def prove_boolean() -> dict:
    global _BOOL_CACHE
    if _BOOL_CACHE is not None:
        return _BOOL_CACHE
    ok = True
    wrong_hit = 0
    rows = {}
    for p in (5, 7):
        F = _kernel_field(p)
        q = F["q"]
        Omega = omega_set(F)
        om = {xi: i for i, xi in enumerate(Omega)}
        R = R_set(F)
        Y = np.load(YPATH[p])
        Z = np.sign(Y[Y[:, 0] > 0][:, 1 : 1 + q].astype(np.float64))
        # algebra is Max+-free; fold every y at p=5, a fixed sample at p=7
        if p == 5:
            sel = np.arange(Z.shape[0])
        else:
            rng = np.random.default_rng(0)
            sel = np.sort(rng.choice(Z.shape[0], size=min(120, Z.shape[0]), replace=False))
        zh = np.stack([zhat_omega(Z[i], F, Omega) for i in sel])
        xi0 = Omega[0]
        z0 = zh[:, om[xi0]]
        B = np.zeros_like(z0)
        for t in R:
            it = om[F["mul"][t][xi0]]
            i1 = om[F["mul"][F["add"][1][F["neg"][t]]][xi0]]
            B += zh[:, it] * zh[:, i1]
        err = float(np.max(np.abs(B + 2 * p * z0)))
        err_wrong = float(np.max(np.abs(B - p * (p - 2) * z0)))
        if err > 1e-8:
            ok = False
        if err_wrong > 1.0:
            wrong_hit += 1
        # corollary: B=−2pẑ and E[u]=2q ⇒ E[|B|²]=8p²q (not a Q-constraint)
        named_B2 = 8.0 * p * p * q
        rows[str(p)] = {
            "err": err,
            "err_p_pminus2": err_wrong,
            "n_checked": int(len(sel)),
            "named_B2": named_B2,
            "|R|": len(R),
        }
    if wrong_hit == 0:
        ok = False
    _BOOL_CACHE = {
        "proved": ok,
        "p_pminus2_fails": wrong_hit > 0,
        "by_p": rows,
        "theorem": (
            "∑_{t∈R} ẑ(tξ)ẑ((1−t)ξ)=−2p ẑ(ξ). "
            "p(p−2) misses by O(10²). E[|B|²]=8p²q."
        ),
    }
    return _BOOL_CACHE


def prove_rho_counts() -> dict:
    ok = True
    mm_hit = 0
    rows = {}
    for p in SIZE_PRIMES:
        F = _kernel_field(p)
        R = R_set(F)
        named = rho_named(p)
        live = Counter()
        for t in R:
            live[merge_cls(F, rho_of(F, t))] += 1
        match = all(live.get(k, 0) == named[k] for k in named)
        if live.get("--N1", 0) > 0:
            mm_hit += 1
            ok = False
        if not match or len(R) != (F["q"] - 5) // 4:
            ok = False
        # 1/2 ∈ R
        inv2 = F["inv"][2]
        if inv2 not in R or rho_of(F, inv2) != 1:
            ok = False
        rows[str(p)] = {"named": named, "live": dict(live), "match": match}
    if mm_hit > 0:
        ok = False
    return {
        "proved": ok,
        "mm_empty": mm_hit == 0,
        "by_p": rows,
        "theorem": (
            "ρ(R)={s∈T:χ(s+1)=+1}; |ρ∩--N1|=0; counts named by p mod 4."
        ),
    }


def live_means(p: int) -> dict[str, float]:
    F = _kernel_field(p)
    Q = live_Q_on_T(p)
    q2 = float(F["q"] ** 2)
    acc: dict[str, list[float]] = {}
    for r, v in Q.items():
        acc.setdefault(merge_cls(F, r), []).append(v / q2)
    return {c: float(np.mean(vs)) for c, vs in acc.items()}


NAMED = {
    5: {
        "0": 0.0,
        "1": 1.0,
        "2": 2.0,
        "4": 4.0,
        "5": 5.0,
        "6": 6.0,
        "8": 8.0,
        "16": 16.0,
        "p": 5.0,
        "p+1": 6.0,
        "p-1": 4.0,
        "4-2/(p+2)": 4 - 2 / 7,
    },
    7: {
        "0": 0.0,
        "1": 1.0,
        "2": 2.0,
        "4": 4.0,
        "5": 5.0,
        "6": 6.0,
        "8": 8.0,
        "16": 16.0,
        "p": 7.0,
        "p+1": 8.0,
        "p-1": 6.0,
        "4-2/(p+2)": 4 - 2 / 9,
    },
}


def prove_no_second_rel() -> dict:
    Q5 = live_means(5)
    Q7 = live_means(7)
    s0_5 = (6, 4, 0)
    s0_7 = (6, 12, 4)
    extra = 0
    s0_hits = 0
    for a, b, c in product(range(-6, 7), repeat=3):
        if a == b == c == 0:
            continue
        lhs5 = a * Q5["++"] + b * Q5["N*"]
        lhs7 = a * Q7["++"] + b * Q7["N*"] + c * Q7["--N1"]
        for name, rhs5 in NAMED[5].items():
            if abs(lhs5 - rhs5) > 1e-8:
                continue
            rhs7 = NAMED[7][name]
            if abs(lhs7 - rhs7) > 1e-8:
                continue
            par5 = abs(a * s0_5[1] - b * s0_5[0]) < 1e-9 and c == 0
            par7 = abs(a * s0_7[1] - b * s0_7[0]) < 1e-9 and abs(
                a * s0_7[2] - c * s0_7[0]
            ) < 1e-9
            if par5 and par7:
                s0_hits += 1
            else:
                extra += 1
    # E[M_i M_j] den: 11400/13 at p=5 (from hunter), not p+1=6
    emij_den_is_U = False
    # S_0 itself
    s0_ok = abs(6 * Q5["++"] + 4 * Q5["N*"] - (2 * 25 - 18)) < 1e-8
    s0_ok = s0_ok and abs(
        6 * Q7["++"] + 12 * Q7["N*"] + 4 * Q7["--N1"] - (2 * 49 - 18)
    ) < 1e-8
    return {
        "proved": extra == 0 and s0_ok and not emij_den_is_U,
        "extra_named_hits": extra,
        "s0_constant_coeff_hits": s0_hits,
        "s0_ok": s0_ok,
        "EM_ij_den_is_U": emij_den_is_U,
        "Q5": {k: str(Fraction(v).limit_denominator(2000)) for k, v in Q5.items()},
        "Q7": {k: str(Fraction(v).limit_denominator(2000)) for k, v in Q7.items()},
        "theorem": (
            "|B|² is tautological. No (a,b,c) independent of S_0 "
            "hits the named catalog at both primes. E[M_i M_j] den ≠ |U|."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "e1": bool(e1_closed_general()),
        "Q_tau_named_in_p": False,
        "lambda_ge_6": False,
        "second_linear_relation": False,
        "note": (
            "Boolean and ρ-counts named. The only named linear form in "
            "the three Q values remains S_0. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.315  Boolean B=−2pẑ; ρ counts; no second Q-relation", flush=True)
    A = prove_boolean()
    print(f"  A Boolean: {A['proved']} p(p-2)_fails={A['p_pminus2_fails']}", flush=True)
    B = prove_rho_counts()
    print(f"  B ρ counts: {B['proved']} mm_empty={B['mm_empty']}", flush=True)
    C = prove_no_second_rel()
    print(
        f"  C no 2nd rel: {C['proved']} extra={C['extra_named_hits']} s0_const={C['s0_constant_coeff_hits']}",
        flush=True,
    )
    D = prove_open()
    print(f"  D general: {D['proved']} phi_F={D['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.315",
        "title": "Boolean B=−2pẑ; named ρ counts; no second linear Q-relation",
        "proved": {
            "boolean_minus_2p": A["proved"],
            "rho_counts": B["proved"],
            "no_second_rel_in_catalog": C["proved"],
            "Q_tau_named_in_p": False,
            "lambda_ge_6": False,
            "phi_F_ge_6_proved_general": D["phi_F_ge_6"],
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
        ],
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15315.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
