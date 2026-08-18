#!/usr/bin/env python3
"""
Prop 15.536 — On Ω-triples, I = p S(r/(1+r)); doubles are CM.

Does **not** flip type_I / phi_F / e1 / L / Aut-Schur / Gsum /
pairing / residual_ii / 15.279–15.535 flags.  Does **not**
edit leftover-1 15.522–15.535 or residual-ii 15.528.

============================================================================
Setup.  15.534: I(η,θ)=G(χ) χ(θ) S(−η/θ) for θ≠0, with
S(λ)=∑_t χ(t(t−1)(t−λ)) and G(χ)=−p χ_p(−1).  Ω={ξ≠0: χ̂(ξ)=p}
={ξ: χ(ξ)=−χ_p(−1)} (15.269 B, 15.523 A).  On ξ+η+θ=0 in Ω
write r=η/ξ (15.269 H).  Then θ=−ξ(1+r) and −η/θ=r/(1+r).

============================================================================
Theorem A — PROVED (15.269 B + 15.534 A; certified p=5,7,11,13).
  χ_Ω=−χ_p(−1), so G χ(θ)=p for every Ω-triple.  Hence
      I(η,θ)=p S(r/(1+r)),     r=η/ξ.
  Also χ(r/(1+r))=χ(r/(1+r)−1)=1.  Fail: drop χ_Ω and write
  I=G S (p=5 doubles: −30≠30).  Fail: I=2p² (p=5: 30≠50).  ∎

Theorem B — PROVED (CM j=1728 + field S; certified p=5..19).
  Doubles are (ξ,ξ,−2ξ), r∈{1,−2,−1/2}, λ∈{1/2,2,−1}.
  These three λ have j(E_λ)=1728, so S(1/2)=S(2)=S(−1).
  Thus I_dbl=p S(−1).  For p≡3 (mod 4), a_p(y²=x³−x)=0 so
  a_{p²}=−2p and S(−1)=2p, hence I_dbl=2p².  For p≡1 (mod 4),
  S(−1)=2p−a_p² with p=a²+b², a odd, a_p²=4a²; I_dbl≠2p²
  (p=5: S=6, I=30).  Fail: I_dbl=2p² at all p≥5.  ∎

Theorem C — OPEN as a name of Δ_conn.  Distinct Ω-triples have
  generic S not a single p-law (p=5: −10; p=7: −2; p=11 splits
  {−10,6,22}).  Even the 3-χ type (1,1,1) splits at p=7
  (S∈{−2,14}).  The contraction still carries den 13,409
  through A.  Aut_e 3A+B still G>T.  type_I flags stay False.

============================================================================
Backend: 15.534 field tables + S-table.  Writes
evidence/e1_gmin_m4_prop15536.json
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
from e1_gmin_m4_prop15523 import G_chi, chi_p_minus_1  # noqa: E402
from e1_gmin_m4_prop15534 import S_table, field_tables  # noqa: E402

EV = ROOT / "evidence" / "e1_gmin_m4_prop15536.json"


def chi_Omega(p: int) -> int:
    """χ on Ω equals −χ_p(−1)."""
    return -chi_p_minus_1(p)


def ap2_cm(p: int) -> int:
    """a_p(y²=x³−x)²: 0 if p≡3 (mod 4), else 4a² for p=a²+b², a odd."""
    if p % 4 == 3:
        return 0
    for a in range(1, p):
        b2 = p - a * a
        if b2 <= 0:
            break
        b = int(b2**0.5)
        if b * b != b2:
            continue
        odd = a if (a % 2 == 1) else b
        return 4 * odd * odd
    raise ValueError(f"no two-square split of {p}")


def S_minus1_named(p: int) -> int:
    """S(−1)=2p−a_p² on y²=x³−x."""
    return 2 * p - ap2_cm(p)


def _cert_p(p: int) -> dict:
    F = field_tables(p)
    S = S_table(F)
    q, ch = F["q"], F["chi"]
    mul, add, sub, inv = F["mul"], F["add"], F["sub"], F["inv"]
    G = G_chi(p)
    com = chi_Omega(p)
    Omega = [xi for xi in range(1, q) if int(ch[xi]) == com]
    Om = set(Omega)
    two = 2
    half = int(inv[two])
    m1 = int(sub[0, 1])
    S_half, S_2, S_m1 = int(S[half]), int(S[two]), int(S[m1])

    n_tri = n_ok = n_dbl = 0
    fail_GS = 0
    I_dbl: set[int] = set()
    I_dist: set[int] = set()
    S_dist: set[int] = set()
    type111: set[int] = set()
    Gchi_ok = True
    for xi in Omega:
        for eta in Omega:
            th = int(sub[0, int(add[xi, eta])])
            if th not in Om:
                continue
            n_tri += 1
            r = int(mul[eta, int(inv[xi])])
            opr = int(add[1, r])
            lam = int(mul[r, int(inv[opr])])
            lam2 = int(mul[int(sub[0, eta]), int(inv[th])])
            I534 = int(G * int(ch[th]) * int(S[lam2]))
            IpS = int(p * int(S[lam]))
            I_GS = int(G * int(S[lam]))
            if int(G * int(ch[th])) != p:
                Gchi_ok = False
            if I534 == IpS and lam == lam2:
                n_ok += 1
            if I_GS == I534:
                fail_GS += 1
            dbl = len({xi, eta, th}) < 3
            if dbl:
                n_dbl += 1
                I_dbl.add(I534)
            else:
                I_dist.add(I534)
                S_dist.add(int(S[lam]))
            if (
                int(ch[lam]) == 1
                and int(ch[int(sub[lam, 1])]) == 1
                and int(ch[int(add[lam, 1])]) == 1
            ):
                type111.add(int(S[lam]))

    named = S_minus1_named(p)
    return {
        "p": p,
        "G": G,
        "chi_Omega": com,
        "n_tri": n_tri,
        "n_dbl": n_dbl,
        "n_ok": n_ok,
        "I_eq_pS": n_ok == n_tri and Gchi_ok,
        "n_GS_collisions": fail_GS,
        "S_half": S_half,
        "S_2": S_2,
        "S_m1": S_m1,
        "S_cm_equal": S_half == S_2 == S_m1,
        "S_m1_named": named,
        "S_m1_hits_named": S_m1 == named,
        "I_dbl": sorted(I_dbl),
        "I_dist": sorted(I_dist),
        "S_dist": sorted(S_dist),
        "I_dbl_is_p_Sm1": I_dbl == {p * S_m1},
        "I_dbl_is_2p2": I_dbl == {2 * p * p},
        "type111_S": sorted(type111),
    }


def _run_certs(primes: tuple[int, ...]) -> list[dict]:
    w = min(8, len(primes))
    with ProcessPoolExecutor(max_workers=w) as ex:
        return list(ex.map(_cert_p, primes))


@lru_cache(maxsize=1)
def prove_A() -> dict:
    rows = {str(r["p"]): r for r in _run_certs((5, 7, 11, 13))}
    r5 = rows["5"]
    ok = all(r["I_eq_pS"] for r in rows.values())
    # fail G S: p=5 doubles I=30, G S(1/2)=−30
    fail_GS = r5["G"] * r5["S_half"] == -30 and r5["I_dbl"] == [30]
    fail_2p2 = 30 != 2 * 5 * 5 and r5["I_dbl"] == [30]
    ok = ok and fail_GS and fail_2p2 and r5["n_GS_collisions"] == 0
    return {
        "proved": bool(ok),
        "rows": rows,
        "p5_I_dbl": r5["I_dbl"],
        "p5_GS": r5["G"] * r5["S_half"],
        "theorem": (
            "On Ω-triples I=p S(r/(1+r)).  "
            "Fail: I=G S (p=5 dbl −30≠30); I=2p² (30≠50)."
        ),
    }


@lru_cache(maxsize=1)
def prove_B() -> dict:
    primes = tuple(p for p in range(5, 20) if is_prime(p))
    extra = {str(r["p"]): r for r in _run_certs(primes) if r["p"] not in (5, 7, 11, 13)}
    rows = dict(prove_A()["rows"])
    rows.update(extra)
    ok = True
    rec = {}
    for p in primes:
        r = rows[str(p)]
        want_2p2 = p % 4 == 3
        row_ok = (
            r["S_cm_equal"]
            and r["S_m1_hits_named"]
            and r["I_dbl_is_p_Sm1"]
            and r["I_dbl_is_2p2"] is want_2p2
            and r["S_m1"] == S_minus1_named(p)
        )
        if not row_ok:
            ok = False
        rec[str(p)] = {
            "S_m1": r["S_m1"],
            "named": S_minus1_named(p),
            "ap2": ap2_cm(p),
            "I_dbl": r["I_dbl"],
            "is_2p2": r["I_dbl_is_2p2"],
            "ok": row_ok,
        }
    # fail-when-wrong: 2p² at p=5
    ok = ok and rec["5"]["I_dbl"] == [30] and not rec["5"]["is_2p2"]
    return {
        "proved": bool(ok),
        "by_p": rec,
        "theorem": (
            "I_dbl=p S(−1)=p(2p−a_p²).  "
            "Equals 2p² iff p≡3 (mod 4).  Fail: 2p² at p=5."
        ),
    }


def prove_open() -> dict:
    rows = prove_A()["rows"]
    split11 = len(rows["11"]["S_dist"]) > 1
    split7_111 = set(rows["7"]["type111_S"]) == {-2, 14}
    return {
        "proved": False,
        "type_I_multilevel_bad_case_ND_closed": bool(
            type_I_multilevel_bad_case_ND_closed()
        ),
        "type_I_aut_e_3AB_positive_general": bool(
            type_I_aut_e_3AB_positive_general()
        ),
        "Delta_conn_p_rational": False,
        "S_dist_p5": rows["5"]["S_dist"],
        "S_dist_p7": rows["7"]["S_dist"],
        "S_dist_p11": rows["11"]["S_dist"],
        "type111_p7": rows["7"]["type111_S"],
        "generic_S_is_p_law": False,
        "split_p11": bool(split11),
        "split_type111_p7": bool(split7_111),
        "note": (
            "Generic S on distinct Ω-triples is not a p-law "
            "(−10, −2, split at p=11).  Δ_conn still den 13,409.  "
            "Aut_e 3A+B still G>T.  Flags unflipped."
        ),
    }


def main() -> dict:
    print("Prop 15.536  I_Ω = p S(r/(1+r)); doubles CM", flush=True)
    A = prove_A()
    print(
        f"  A I=pS: {A['proved']} p5_Idbl={A['p5_I_dbl']} GS={A['p5_GS']}",
        flush=True,
    )
    B = prove_B()
    print(
        f"  B CM doubles: {B['proved']} "
        f"p5={B['by_p']['5']['I_dbl']} p7={B['by_p']['7']['I_dbl']}",
        flush=True,
    )
    C = prove_open()
    print(
        f"  C open: ND={C['type_I_multilevel_bad_case_ND_closed']} "
        f"3AB={C['type_I_aut_e_3AB_positive_general']} "
        f"S_dist11={C['S_dist_p11']}",
        flush=True,
    )
    out = {
        "prop": "15.536",
        "title": "On Ω-triples I=p S(r/(1+r)); doubles are CM S(−1)",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "C": C,
        "proved": {
            "I_Omega_eq_p_S": A["proved"],
            "I_dbl_CM": B["proved"],
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
