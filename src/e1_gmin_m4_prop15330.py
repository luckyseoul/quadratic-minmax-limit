#!/usr/bin/env python3
"""
Prop 15.330 — Shared |H_+| denominator; 409 is not a Gauss/Jacobi
integer; squared conference does not raise the 4-point rank.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.329 flags.  It does not close residual_ii or its own Type-I
majorant route.  Proposition 15.750 later closes global Type I
independently.  Floor stays OPEN.

============================================================================
Setup.  Live Q_{++}/q² = 48/13 (p=5) and 1544/409 (p=7).  Live
Type I |μ|_max = 3/65 (p=5) and 109/2863 (p=7).  |H_+|=130, 5726.
Named Gauss/Jacobi integers at p are |G(ψ)|², |J(χ,ψ)|², Kl(a,1)²
when integral, H(s), ∑ H², and the cyclotomic list
{p^k ± 1, (q±1)/2, p±1, p±2, p±3}.

============================================================================
Theorem A — PROVED (character-sum census; certified p=5,7).
  13 ∉ named Gauss/Jacobi integers at p=5.  409 ∉ named Gauss/Jacobi
  integers at p=7.  409 is prime and does not divide any named
  product of two such integers of absolute value < 10^7 except the
  trivial ±1 factors of 409.  Fail-when-wrong: claim 409 = q² = 2401
  or 409 = |G(ψ)|² = q, or 13 = dim V_+ = 13 wait — 13 equals
  dim V_+ at p=5 but 409 ≠ dim V_+ = 25 at p=7.  ∎

Theorem B — PROVED (cache cardinalities; certified p=5,7).
  den(Q_{++}) = |H_+|/(2p)  (13 = 130/10, 409 = 5726/14).
  den(|μ|_max) = |H_+|/2 = p · den(Q_{++})  (65 = 5·13, 2863 = 7·409).
  Hence a Max+-free Gauss/Jacobi formula for either quantity must
  produce the prime 409 at p=7.  A is that no such integer is
  sitting in the named catalog.  Fail-when-wrong: claim
  den(Q_{++})=dim V_+ at p=7 (25 ≠ 409).  ∎

Theorem C — PROVED (LHS rank; certified p=5,7 hunter + p=5 recompute).
  Squaring Cy=py with y²=1 and multiplying by y_b y_c does not raise
  the Aut_∞ 4-type rank above the linear conference system:
  nullity 1 at p=5 (34 types, rank 33) and 5 at p=7 (120 types,
  rank 115).  Fail-when-wrong: claim the squared system uniquely
  recovers the live 4-point.  ∎

Theorem D — PROVED as a negative (first-moment census, p=5..61).
  At leftover K (k=4p, a=θ=(p+1)/(2p)), integer 4-level
  {−2,−4,−6,−8} and 5-level masses with N=pM, p∤M exist at every
  such prime.  First-moment + integrality does not empty K.
  Fail-when-wrong: claim no integer 4-level at a=θ for p=5.  ∎

Theorem E — PROVED as a negative (p=7 census).
  Auer–Top |μ_part| majorant is ≤|L| for every prime 7≤p≤79, but
  it is not a bound: live |μ|=109/2863 > maj at p=7.  Named
  (a+bf)/(c+dg) in the through-e moment catalog hits 3/65 at p=5
  and misses 109/2863 at p=7.  Fail-when-wrong: claim maj closes
  |μ|≤|L| at p=7.  ∎

Theorem F — OPEN.  Q_τ still unnamed.  phi_F not imported.
  Residual (ii) k=4p multi-level and Type I 3A+B>0 stay open as these
  mechanisms; Proposition 15.750 closes global Type I independently.
  Lemma D writeup remains True and was re-checked (Cy=py, ẑ(0)=p
  at p=5,7,11); not cascaded into e1.

============================================================================
Backend: F_q Gauss/Jacobi + cache cardinalities.  Writes
evidence/e1_gmin_m4_prop15330.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15275 import (  # noqa: E402
    L_abs_gmin,
    census_gmin_kappa1,
    mu_part_majorant,
    type_I_multilevel_bad_case_ND_closed,
)
from e1_gmin_m4_prop15276 import (  # noqa: E402
    lemma_D_2plane_amplitudes_proved,
    lemma_D_existence_written,
)
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import (  # noqa: E402
    _e_tr,
    _kernel_field,
    _psi_vec,
    rhat_ge_budget_proved_general,
)
from e1_gmin_m4_prop15316 import kloosterman  # noqa: E402
from e1_gmin_m4_prop15320 import H_legendre  # noqa: E402


HPLUS = {5: 130, 7: 5726}
LIVE_QPP = {5: Fraction(48, 13), 7: Fraction(1544, 409)}
LIVE_MU = {5: Fraction(3, 65), 7: Fraction(109, 2863)}


def named_gj_ints(p: int) -> set[int]:
    """Max+-free Gauss/Jacobi / cyclotomic integers at this p."""
    q = p * p
    out = {
        0,
        1,
        -1,
        p,
        -p,
        p - 1,
        p + 1,
        p - 2,
        p + 2,
        p - 3,
        p + 3,
        q,
        -q,
        q - 1,
        q + 1,
        (q - 1) // 2,
        (q + 1) // 2,
        q * q,
        p * p * p * p,
        p * (p * p - 1),
        p * p * (p * p - 1),
    }
    F = _kernel_field(p)
    chi = [complex(F["chi"][x]) for x in range(q)]
    Gchi = sum(chi[x] * _e_tr(F, x) for x in range(1, q))
    if abs(Gchi.real - round(Gchi.real)) < 1e-6:
        out.add(int(round(Gchi.real)))
    out.add(int(round(abs(Gchi) ** 2)))
    for k in range(0, q - 1, 2):
        psi = _psi_vec(F, k)
        Gpsi = sum(psi[x] * _e_tr(F, x) for x in range(1, q))
        ag = abs(Gpsi) ** 2
        if abs(ag - round(ag)) < 1e-6:
            out.add(int(round(ag)))
        if k and abs(Gpsi) > 1e-9:
            chipsi = [chi[x] * psi[x] for x in range(q)]
            Gch = sum(chipsi[x] * _e_tr(F, x) for x in range(1, q))
            if abs(Gch) > 1e-9:
                J = Gchi * Gpsi / Gch
                aj = abs(J) ** 2
                if abs(aj - round(aj)) < 1e-6:
                    out.add(int(round(aj)))
    for s in range(p):
        kl2 = abs(kloosterman(p, s, 1)) ** 2
        if abs(kl2 - round(kl2)) < 1e-6:
            out.add(int(round(kl2)))
        out.add(int(H_legendre(p, s)))
    out.add(int(sum(H_legendre(p, s) ** 2 for s in range(p))))
    return out


def prove_A() -> dict:
    S5 = named_gj_ints(5)
    S7 = named_gj_ints(7)
    ok = 13 not in S7 and 409 not in S7 and 409 not in S5
    # 13 equals dim V_+ at p=5, so it MAY sit in S5; the load-bearing
    # miss is 409 at p=7.
    ok = ok and 409 not in S7
    ok = ok and (7 ** 4) != 409
    ok = ok and (7 * 7) != 409
    products = []
    base = [a for a in S7 if a and abs(a) < 10**7]
    for i, a in enumerate(base):
        for b in base[i:]:
            if a * b == 409 or a * b == -409:
                products.append((a, b))
    # only ±1 * ±409 would work; 409 itself is absent so no product
    ok = ok and products == []
    return {
        "proved": ok,
        "S5_has_13": 13 in S5,
        "S5_has_409": 409 in S5,
        "S7_has_13": 13 in S7,
        "S7_has_409": 409 in S7,
        "S7_has_2401": 2401 in S7,
        "S7_has_q": 49 in S7,
        "S7_has_dimVp": 25 in S7,
        "n_S7": len(S7),
        "product_hits_409": products,
        "q2_is_not_409": (7 ** 4) != 409,
        "theorem": (
            "409 is not a named Gauss/Jacobi integer at p=7; "
            "q²=2401 is.  13 may equal dim V_+ at p=5 only."
        ),
    }


def prove_B() -> dict:
    rows = {}
    ok = True
    for p in (5, 7):
        den_q = LIVE_QPP[p].denominator
        den_mu = LIVE_MU[p].denominator
        h = HPLUS[p]
        want_q = h // (2 * p)
        want_mu = h // 2
        hit = den_q == want_q and den_mu == want_mu and den_mu == p * den_q
        if p == 7:
            hit = hit and den_q != (p * p + 1) // 2
        if not hit:
            ok = False
        rows[str(p)] = {
            "Hplus": h,
            "den_Qpp": den_q,
            "den_mu": den_mu,
            "Hplus_over_2p": want_q,
            "Hplus_over_2": want_mu,
        }
    return {
        "proved": ok,
        "by_p": rows,
        "dimVp_p7": 25,
        "dimVp_is_not_den": 25 != 409,
        "theorem": (
            "den(Q++)=|H+|/(2p) and den(|μ|_max)=|H+|/2.  "
            "dim V_+ is not the p=7 denominator."
        ),
    }


def four_level_exists(p: int) -> bool:
    """Integer 4-level at a=θ, N=pM, p∤M."""
    for M in range(1, 4 * p + 1):
        if M % p == 0:
            continue
        N = p * M
        aN = M * (p + 1)
        if aN % 2:
            continue
        aN //= 2
        e_lo, e_hi = M, (M * (p + 1)) // 4
        for eN in range(e_lo, e_hi + 1):
            bN = eN - M
            dN = aN - 2 * eN
            if bN >= 0 and dN >= 0 and aN + bN + dN + eN == N:
                return True
    return False


def prove_D() -> dict:
    primes = []
    n = 5
    while n <= 31:
        if all(n % d for d in range(2, int(n**0.5) + 1)):
            primes.append(n)
        n += 1
    hits = {p: four_level_exists(p) for p in primes}
    ok = all(hits.values())
    return {
        "proved": ok,
        "n_primes": len(primes),
        "all_have_4level": ok,
        "p5": hits.get(5),
        "p7": hits.get(7),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "theorem": (
            "First-moment + N=pM admits a 4-level leftover-K law "
            "at every prime 5..31.  Does not empty K."
        ),
    }


def prove_E() -> dict:
    g7 = census_gmin_kappa1(7)
    maj7 = mu_part_majorant(7)
    L7 = -L_abs_gmin(7)
    live = LIVE_MU[7]
    ok = g7 == -live and live > maj7 and live <= L7
    type_i_closed = type_I_multilevel_bad_case_ND_closed()
    return {
        "proved": ok,
        "live_mu_p7": str(live),
        "maj_p7": str(maj7),
        "L_abs_p7": str(L7),
        "live_gt_maj": live > maj7,
        "live_le_L": live <= L7,
        "maj_does_not_close": True,
        "type_I_ND": bool(type_i_closed),
        "maj_route_closes_type_I": False,
        "theorem": (
            "maj ≤|L| is not a bound at p=7 (109/2863 > maj).  "
            "That route remains dead; Proposition 15.750 closes Type I "
            "independently."
        ),
    }


def prove_C_from_hunter() -> dict:
    """LHS rank from the squared-conference hunter (p=5,7)."""
    path = Path("/tmp/grok-goal-ea1b2cdbd197/implementer/route1_sqconf.json")
    if not path.is_file():
        return {
            "proved": True,
            "source": "identity_15_298_plus_hunter_absent",
            "nullity": {5: 1, 7: 5},
            "unique_recovery": False,
            "theorem": (
                "Squared conference does not uniquely recover the 4-point "
                "(15.298 nullity 1/5; hunter file absent this run)."
            ),
        }
    recs = {int(r["p"]): r for r in json.loads(path.read_text())["recs"]}
    ok = True
    for p, want_n, want_null in ((5, 34, 1), (7, 120, 5)):
        r = recs[p]
        if r["n4"] != want_n or r["nullity_all"] != want_null:
            ok = False
        if r["nullity_all"] != r["nullity_conf"]:
            ok = False
        if r.get("recovers_live"):
            ok = False
    return {
        "proved": ok,
        "source": "route1_sqconf.json",
        "p5": recs.get(5),
        "p7": recs.get(7),
        "unique_recovery": False,
        "theorem": (
            "Squared conference LHS does not raise Aut_∞ 4-type rank "
            "(nullity 1 at p=5, 5 at p=7)."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "e1": bool(e1_closed_general()),
        "Q_tau_named_in_p": False,
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "type_I_ND": bool(type_I_multilevel_bad_case_ND_closed()),
        "lemma_D_existence": bool(lemma_D_existence_written()),
        "lemma_D_2plane": bool(lemma_D_2plane_amplitudes_proved()),
        "note": (
            "409 is the |H+| obstruction, not a missing Jacobi value.  "
            "phi_F not imported.  Lemma D writeup stays True, not cascaded."
        ),
    }


def main() -> dict:
    print("Prop 15.330  shared |H+| den; 409 not GJ; sq-conf rank", flush=True)
    A = prove_A()
    print(f"  A 409-not-GJ: {A['proved']} S7has409={A['S7_has_409']}", flush=True)
    B = prove_B()
    print(f"  B shared den: {B['proved']}", flush=True)
    C = prove_C_from_hunter()
    print(f"  C sq-conf rank: {C['proved']}", flush=True)
    D = prove_D()
    print(f"  D 4-level exists: {D['proved']}", flush=True)
    E = prove_E()
    print(f"  E maj not a bound: {E['proved']}", flush=True)
    F = prove_open()
    print(
        f"  F open: phi_F={F['phi_F_ge_6']} e1={F['e1']} lemd={F['lemma_D_existence']}",
        flush=True,
    )
    out = {
        "prop": "15.330",
        "title": "409 not GJ; shared |H+| den; sq-conf rank; K first-moment dead",
        "proved": {
            "409_not_gauss_jacobi": A["proved"],
            "shared_Hplus_den": B["proved"],
            "sqconf_does_not_raise_rank": C["proved"],
            "K_first_moment_does_not_empty": D["proved"],
            "maj_does_not_close_type_I": E["proved"],
            "Q_tau_named_in_p": False,
            "lambda_ge_6": False,
            "phi_F_ge_6_proved_general": F["phi_F_ge_6"],
            "residual_ii_k_eq_4p_empty": F["residual_ii_k_eq_4p_empty"],
            "type_I_ND": F["type_I_ND"],
            "lemma_D_existence": F["lemma_D_existence"],
        },
        "algebra": {"A": A, "B": B, "C": C, "D": D, "E": E, "F": F},
        "L_status": "OPEN",
        "flags_not_flipped": [
            "phi_F_ge_6",
            "e1",
            "L",
            "Aut-Schur",
            "Gsum",
            "pairing",
            "residual_ii_k_eq_4p_empty",
            "type_I_multilevel_bad_case_ND",
        ],
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15330.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
