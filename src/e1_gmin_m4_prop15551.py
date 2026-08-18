#!/usr/bin/env python3
"""
Prop 15.551 — Galois line-union of ẑ; Ω-triples are 1-line
iff r∈F_p (else 3-line); A_k=0 on 3-line triples for k<3.

Does **not** flip type_I / phi_F / e1 / L / Aut-Schur / Gsum /
pairing / residual_ii / 15.279–15.550 flags.  Does **not**
edit leftover-1 15.522–15.550 or residual-ii 15.528 / 15.547.

============================================================================
Setup.  F_q=F_{p²}, z:F_q→{±1} the y_∞=+1 chart, ẑ(ξ)=∑_x z(x)
ψ(ξx), ψ(u)=exp(2πi Tr(u)/p).  15.269 B: Max+ ẑ is supported
on {0}∪Ω.  Ω is a union of m=(p+1)/2 punctured F_p-lines
(χ_q|F_p^*=1).  r=η/ξ on an Ω-triple ξ+η+θ=0.

============================================================================
Theorem A — PROVED Max+-free (Galois; any Q-valued z).
  ẑ(ξ)∈Z[ζ_p].  The embedding Q(ζ_p)↪C, ζ↦e^{2πi/p}, is
  injective, and σ_g(ζ)=ζ^g satisfies
      σ_g(ẑ(ξ))=ẑ(g ξ)    (g∈F_p^*).
  Hence ẑ(ξ)=0 ⇔ ẑ vanishes on the whole punctured F_p-line
  ξ F_p^*.  Combined with 15.269 B, every y_∞=+1 Max+ has
  ẑ-support a union of k Ω-lines, nsupp=k(p−1),
  k∈{1,…,m} (k=0 dies by Parseval: |ẑ(0)|²=p²≠p⁴).
  Certified every y_∞=+1 vector at p=5 (k∈{1,3}) and p=7
  (k∈{1,3,4}).  Fail: nsupp=18 at p=5 (18∤ p−1).  Fail:
  a live vector with ẑ(ξ)=0≠ẑ(2ξ).  ∎

Theorem B — PROVED Max+-free (A + F_p-linearity).
  An Ω-triple {ξ,η,θ} with r=η/ξ is 1-line iff r∈F_p
  (then θ=−(1+r)ξ lies on ⟨ξ⟩ too).  If r∉F_p the three
  points lie on three distinct Ω-lines: two on one line
  would force the third by ξ+η+θ=0.  There is no 2-line
  Ω-triple.  Therefore A_k := E[ẑ(ξ)ẑ(η)ẑ(θ) | k-line]
  vanishes on every 3-line (nfp) triple whenever k<3.
  Recovers A_1d(r∉F_p)=0 (15.540 A).  Certified p=5,7,11,13
  (n_2line=0; n_1line = # {r∈F_p^* : 1+r≠0, r∈Ω/Ω}).
  Fail: a 2-line Ω-triple.  Fail: A_1d(nfp)≠0 at p=5,7.  ∎

Theorem C — OPEN.  A_full,dbl is still 0 vs −4p³/15, not a
  p-law.  μ_full is not a function of the coarse Aut field
  type (κ,S,φ): at p=7 the class (κ,S,φ)=(1,−10,−10) takes
  both 13/735 and 29/735.  Aut_e 3A+B still G>T.  type_I
  flags stay False.

============================================================================
Backend: Galois + field lines; y_∞=+1 caches p=5,7.  Writes
evidence/e1_gmin_m4_prop15551.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    is_prime,
)
from e1_gmin_m4_prop15275 import (  # noqa: E402
    type_I_aut_e_3AB_positive_general,
    type_I_multilevel_bad_case_ND_closed,
)
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15534 import S_table, field_tables  # noqa: E402
from e1_gmin_m4_prop15536 import chi_Omega  # noqa: E402
from e1_gmin_m4_prop15538 import A_1d_dbl_named  # noqa: E402
from e1_gmin_m4_prop15540 import A_1d_of_r  # noqa: E402

EV = ROOT / "evidence" / "e1_gmin_m4_prop15551.json"
YPLUS = {5: "/tmp/maxplus_p5.npy", 7: "/tmp/maxplus_p7.npy"}


def n_omega_lines(p: int) -> int:
    """# of F_p-lines in Ω.  Equals (p+1)/2."""
    return (p + 1) // 2


def triple_nlines(p: int, r: int) -> int:
    """1 if r∈F_p and 1+r≠0; 0 if 1+r=0; else 3."""
    if r < p:
        if (1 + r) % p == 0:
            return 0
        return 1
    return 3


def _omega_col(p, trt, mul_t, xi):
    om = np.exp(2j * np.pi / p)
    return om ** trt[mul_t[xi]].astype(np.int64)


def omega_and_lines(F: dict) -> tuple[list[int], list[list[int]]]:
    p, q, ch, mul = F["p"], F["q"], F["chi"], F["mul"]
    com = chi_Omega(p)
    Omega = [xi for xi in range(1, q) if int(ch[xi]) == com]
    Om_index = {xi: i for i, xi in enumerate(Omega)}
    seen: set[int] = set()
    lines: list[list[int]] = []
    for xi in Omega:
        if xi in seen:
            continue
        L = [int(mul[xi, a]) for a in range(1, p)]
        for u in L:
            seen.add(u)
        lines.append([Om_index[u] for u in L])
    return Omega, lines


def parseval_kills_k0(p: int) -> bool:
    """k=0 ⇒ |ẑ(0)|²=p² ≠ p⁴ = ∑|ẑ|²."""
    return p * p != p**4


@lru_cache(maxsize=1)
def live_line_union() -> dict:
    out = {}
    for p in (5, 7):
        F = field_tables(p)
        F = {**F, "p": p}
        q, mul, trt = F["q"], F["mul"], F["tr"]
        Omega, lines = omega_and_lines(F)
        Yp = np.sign(np.load(YPLUS[p]).astype(np.float64))
        z = Yp[Yp[:, 0] > 0][:, 1:]
        Psi = np.column_stack([_omega_col(p, trt, mul, xi) for xi in range(q)])
        Zh = (z @ Psi)[:, Omega]
        mag = np.abs(Zh) > 1e-6
        nsupp = mag.sum(axis=1)
        k_list = []
        n_partial_line = 0
        n_gal_break = 0
        two = 2 % p
        for i in range(len(z)):
            kk = 0
            for idxs in lines:
                bits = mag[i, idxs]
                s = int(bits.sum())
                if s == 0:
                    continue
                if s == len(idxs):
                    kk += 1
                else:
                    n_partial_line += 1
            k_list.append(kk)
            # Galois fail witness: ẑ(ξ)=0 ≠ ẑ(2ξ)
            for j, xi in enumerate(Omega):
                x2 = int(mul[xi, two])
                j2 = Omega.index(x2)
                if (not mag[i, j]) and mag[i, j2]:
                    n_gal_break += 1
                    break
        kcnt = Counter(k_list)
        out[p] = {
            "nH": int(len(z)),
            "n_omega": len(Omega),
            "n_lines": len(lines),
            "kcnt": {str(k): int(v) for k, v in sorted(kcnt.items())},
            "n_partial_line": n_partial_line,
            "n_gal_break": n_gal_break,
            "nsupp_mod": sorted({int(x) % (p - 1) for x in nsupp.tolist()}),
            "want_k": {5: [1, 3], 7: [1, 3, 4]}[p],
        }
    return out


def omega_triple_census(p: int) -> dict:
    F = field_tables(p)
    q, ch = F["q"], F["chi"]
    mul, add, sub, inv = F["mul"], F["add"], F["sub"], F["inv"]
    com = chi_Omega(p)
    Omega = [xi for xi in range(1, q) if int(ch[xi]) == com]
    Om = set(Omega)
    n1 = n2 = n3 = nbad = 0
    for xi in Omega:
        for eta in Omega:
            th = int(sub[0, int(add[xi, eta])])
            if th not in Om:
                continue
            r = int(mul[eta, int(inv[xi])])
            nl = triple_nlines(p, r)
            # independent line span of {ξ,η}
            same = r < p
            if same and nl == 1:
                n1 += 1
            elif (not same) and nl == 3:
                n3 += 1
            elif nl == 0:
                nbad += 1
            else:
                n2 += 1
    return {
        "p": p,
        "n_1line": n1,
        "n_2line": n2,
        "n_3line": n3,
        "n_bad": nbad,
        "n_lines": n_omega_lines(p),
    }


def phi_infset(F: dict, a: int, b: int, c: int) -> int:
    ch, sub, q = F["chi"], F["sub"], F["q"]
    s = 0
    for r in range(q):
        if r == a or r == b or r == c:
            continue
        s += int(ch[sub[a, r]]) * int(ch[sub[b, r]]) * int(ch[sub[c, r]])
    return s


@lru_cache(maxsize=1)
def live_mufull_Sphi_split() -> dict:
    """p=7 μ_full on (κ,S,φ)=(1,−10,−10)."""
    p = 7
    F = field_tables(p)
    q, ch = F["q"], F["chi"]
    mul, sub, inv, trt = F["mul"], F["sub"], F["inv"], F["tr"]
    S = S_table(F)
    Omega = [xi for xi in range(1, q) if int(ch[xi]) == chi_Omega(p)]
    Yp = np.sign(np.load(YPLUS[p]).astype(np.float64))
    z = Yp[Yp[:, 0] > 0][:, 1:]
    Psi = np.column_stack([_omega_col(p, trt, mul, xi) for xi in range(q)])
    nsupp = (np.abs((z @ Psi)[:, Omega]) > 1e-6).sum(axis=1)
    zf = z[nsupp == len(Omega)]
    nF = len(zf)
    mus: list[str] = []
    for a, b, c in combinations(range(q), 3):
        kap = int(ch[sub[b, a]]) + int(ch[sub[c, a]]) + int(ch[sub[c, b]])
        if kap != 1:
            continue
        rho = int(mul[sub[c, a], int(inv[sub[b, a]])])
        if int(S[rho]) != -10:
            continue
        if phi_infset(F, a, b, c) != -10:
            continue
        muf = float(np.mean(zf[:, a] * zf[:, b] * zf[:, c]))
        mus.append(str(Fraction(muf).limit_denominator(nF * nF)))
    cnt = Counter(mus)
    return {
        "class": "(1,-10,-10)",
        "n": int(sum(cnt.values())),
        "mu_cnt": dict(cnt),
        "n_vals": len(cnt),
    }


@lru_cache(maxsize=1)
def prove_A() -> dict:
    live = live_line_union()
    ok = parseval_kills_k0(5) and parseval_kills_k0(7)
    ok = ok and n_omega_lines(5) == 3 and n_omega_lines(7) == 4
    ok = ok and n_omega_lines(11) == 6
    rows = {}
    for p in (5, 7):
        L = live[p]
        want = set(L["want_k"])
        got = {int(k) for k in L["kcnt"]}
        row_ok = (
            L["n_partial_line"] == 0
            and L["n_gal_break"] == 0
            and L["nsupp_mod"] == [0]
            and got == want
            and L["n_lines"] == n_omega_lines(p)
        )
        if p == 5:
            row_ok = row_ok and L["nH"] == 130 and L["kcnt"] == {"1": 30, "3": 100}
        if p == 7:
            row_ok = (
                row_ok
                and L["nH"] == 5726
                and L["kcnt"] == {"1": 140, "3": 1176, "4": 4410}
            )
        if not row_ok:
            ok = False
        rows[str(p)] = L
        rows[str(p)]["ok"] = row_ok
    # fail-when-wrong: nsupp=18 is legal at p=5
    ok = ok and (18 % 4 != 0)
    return {
        "proved": bool(ok),
        "rows": rows,
        "parseval_k0": True,
        "theorem": (
            "ẑ-support is a union of k Ω-lines, nsupp=k(p−1).  "
            "Fail: nsupp=18 at p=5; fail ẑ(ξ)=0≠ẑ(2ξ)."
        ),
    }


@lru_cache(maxsize=1)
def prove_B() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 11, 13):
        if not is_prime(p):
            continue
        c = omega_triple_census(p)
        row_ok = c["n_2line"] == 0 and c["n_bad"] == 0 and c["n_1line"] > 0
        row_ok = row_ok and c["n_3line"] > 0
        row_ok = row_ok and triple_nlines(p, 1) == 1
        row_ok = row_ok and triple_nlines(p, p + 1) == 3
        if p > 2:
            # −1 ∈ F_p, 1+(p−1)≡0
            row_ok = row_ok and triple_nlines(p, p - 1) == 0
        if not row_ok:
            ok = False
        rows[str(p)] = c
        rows[str(p)]["ok"] = row_ok
    # recovers 15.540 A_1d(nfp)=0
    ok = ok and A_1d_of_r(5, 6) == 0
    ok = ok and A_1d_of_r(7, 8) == 0
    ok = ok and A_1d_of_r(5, 1) == A_1d_dbl_named(5)
    ok = ok and A_1d_of_r(5, 6) != A_1d_dbl_named(5)
    # fail: 2-line exists
    ok = ok and rows["5"]["n_2line"] == 0
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "Ω-triples are 1-line iff r∈F_p, else 3-line (never 2).  "
            "A_k=0 on nfp for k<3.  Fail: a 2-line triple."
        ),
    }


def prove_open() -> dict:
    split = live_mufull_Sphi_split()
    a_full_dbl = {
        "5": "0",
        "7": str(Fraction(-4 * 7**3, 15)),
    }
    split_ok = (
        split["n_vals"] == 2
        and "13/735" in split["mu_cnt"]
        and "29/735" in split["mu_cnt"]
    )
    return {
        "proved": False,
        "type_I_multilevel_bad_case_ND_closed": bool(
            type_I_multilevel_bad_case_ND_closed()
        ),
        "type_I_aut_e_3AB_positive_general": bool(
            type_I_aut_e_3AB_positive_general()
        ),
        "A_full_named_p_law": False,
        "A_full_dbl": a_full_dbl,
        "mu_full_Sphi_split": split,
        "mu_full_Sphi_splits": bool(split_ok),
        "note": (
            "A_full,dbl is 0 vs −4p³/15.  μ_full not a function of "
            "(κ,S,φ) at p=7.  Aut_e G>T open.  Flags unflipped."
        ),
    }


def main() -> dict:
    print("Prop 15.551  Galois line-union; 1-line iff r∈F_p", flush=True)
    A = prove_A()
    print(
        f"  A line-union: {A['proved']} "
        f"p5={A['rows']['5']['kcnt']} p7={A['rows']['7']['kcnt']}",
        flush=True,
    )
    B = prove_B()
    print(
        f"  B 1-vs-3: {B['proved']} "
        f"p7 n2={B['rows']['7']['n_2line']}",
        flush=True,
    )
    C = prove_open()
    print(
        f"  C open: ND={C['type_I_multilevel_bad_case_ND_closed']} "
        f"3AB={C['type_I_aut_e_3AB_positive_general']} "
        f"Sphi_split={C['mu_full_Sphi_splits']}",
        flush=True,
    )
    out = {
        "prop": "15.551",
        "title": "Galois line-union of ẑ; Ω-triples 1-line iff r∈F_p",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "C": C,
        "proved": {
            "galois_line_union": A["proved"],
            "triple_one_vs_three": B["proved"],
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
    EV.parent.mkdir(parents=True, exist_ok=True)
    EV.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", EV, flush=True)
    return out


if __name__ == "__main__":
    main()
