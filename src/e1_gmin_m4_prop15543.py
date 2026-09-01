#!/usr/bin/env python3
"""
Prop 15.543 — Type+ 1D 3-point on |κ|=1 is κ/(p(p−2)).

Does **not** flip type_I / phi_F / e1 / L / Aut-Schur / Gsum /
pairing / residual_ii / 15.279–15.542 flags.  Does **not**
edit leftover-1 15.522–15.542 or residual-ii 15.528.

============================================================================
Setup.  Type+ 1D: ker a square F_p-line L, z=ε∘π_L, ε:F_p→{±1}
with #(+)=m=(p+1)/2 (15.538).  n_sq=(p+1)/2 square lines.
κ=χ(b−a)+χ(c−a)+χ(c−b) on an ∞-set {∞,a,b,c}.  Aut_e G>T
hinge (15.275 J) is G>T(p)=−(p−2)/(p(2p−1)) on |κ|=1.

============================================================================
Theorem A — PROVED Max+-free (Johnson 3-point + line count).
  t=1 (all three points in one coset of L) forces the three
  differences into L, hence κ=±3, so t=1 never occurs on |κ|=1.
  t=2 iff exactly one difference lies in L.  That L is square
  iff the difference is a square, so n₂ equals the number of
  +1 differences among the three: n₂=2 when κ=1 and n₂=1
  when κ=−1.  E_ε is 1/p on t=2 and μ₃=−3/(p(p−2)) on t=3
  (15.538).  Averaging over square lines:
      μ_{1d}=(n₂/p + n₃ μ₃)/n_sq = κ/(p(p−2)),
  with n₃=n_sq−n₂.  Certified all |κ|=1 ∞-sets at p=5,7.
  Fail: κ/p² (p=5: 1/15≠1/25).  Fail: μ_1d=T at p=7
  (1/35≠5/91).  ∎

Theorem B — PROVED as a corollary (does not close Aut_e).
  |μ_{1d}|=1/(p(p−2)) ≤ |T|=(p−2)/(p(2p−1)) for all primes
  p≥5, with equality iff p=5.  So the 1D slice is strictly
  inside G>T for p>5 and sits on the 3A+B=0 wall at p=5.
  Free/full-Ω μ is not this p-law (p=5: μ_full=±1/p²;
  p=7: μ_full∈{13,17,29}/(15 p²)).  Aut_e 3A+B on the
  full Max+ ensemble still G>T.  ∎

Theorem C — OPEN.  A_full on r-classes is still not a p-law
  (15.540: 0 vs −4p³/15).  Mixed A still den 13,409.
  This mechanism does not close Type I; Proposition 15.750 closes it independently.

============================================================================
Backend: Fraction identities + y_∞=+1 1D masks.  Writes
evidence/e1_gmin_m4_prop15543.json
"""
from __future__ import annotations

import json
import os
import sys
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    is_prime,
)
from e1_gmin_m4_prop15275 import (  # noqa: E402
    T_bitight,
    type_I_aut_e_3AB_positive_general,
    type_I_multilevel_bad_case_ND_closed,
)
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15534 import field_tables  # noqa: E402
from e1_gmin_m4_prop15536 import chi_Omega  # noqa: E402
from e1_gmin_m4_prop15538 import mu3_named  # noqa: E402

EV = ROOT / "evidence" / "e1_gmin_m4_prop15543.json"
YPLUS = {5: "/tmp/maxplus_p5.npy", 7: "/tmp/maxplus_p7.npy"}


def n2_of_kappa(kap: int) -> int:
    """Number of square lines with t=2.  Equals # of +1 differences."""
    return 2 if kap == 1 else 1


def mu_1d_named(p: int, kap: int) -> Fraction:
    return Fraction(kap, p * (p - 2))


def mu_1d_wrong_p2(p: int, kap: int) -> Fraction:
    return Fraction(kap, p * p)


def mu_1d_from_n2(p: int, kap: int) -> Fraction:
    n_sq = Fraction(p + 1, 2)
    n2 = Fraction(n2_of_kappa(kap), 1)
    n3 = n_sq - n2
    return (n2 / p + n3 * mu3_named(p)) / n_sq


def _omega_col(p, trt, mul_t, xi):
    om = np.exp(2j * np.pi / p)
    return om ** trt[mul_t[xi]].astype(np.int64)


@lru_cache(maxsize=1)
def live_mu_1d() -> dict:
    out = {}
    for p in (5, 7):
        F = field_tables(p)
        q, ch = F["q"], F["chi"]
        sub, trt, mul = F["sub"], F["tr"], F["mul"]
        Yp = np.sign(np.load(YPLUS[p]).astype(np.float64))
        Yp = Yp[Yp[:, 0] > 0]
        z = Yp[:, 1:]
        Omega = [xi for xi in range(1, q) if int(ch[xi]) == chi_Omega(p)]
        Psi = np.column_stack(
            [_omega_col(p, trt, mul, xi) for xi in range(q)]
        )
        nsupp = (np.abs((z @ Psi)[:, Omega]) > 1e-6).sum(axis=1)
        is_1d = nsupp == (p - 1)
        z1 = z[is_1d]
        n_ok = 0
        n_tot = 0
        worst = 0.0
        for a, b, c in combinations(range(q), 3):
            kap = int(ch[sub[b, a]]) + int(ch[sub[c, a]]) + int(ch[sub[c, b]])
            if abs(kap) != 1:
                continue
            n_tot += 1
            live = float(np.mean(z1[:, a] * z1[:, b] * z1[:, c]))
            want = float(mu_1d_named(p, kap))
            err = abs(live - want)
            if err > worst:
                worst = err
            if err < 1e-9:
                n_ok += 1
        out[p] = {
            "n_1d": int(is_1d.sum()),
            "n_kappa1": n_tot,
            "n_ok": n_ok,
            "max_err": worst,
            "ok": n_ok == n_tot and n_tot > 0,
        }
    return out


@lru_cache(maxsize=1)
def prove_A() -> dict:
    ok = True
    id_rows = {}
    for p in range(5, 32):
        if not is_prime(p):
            continue
        for kap in (1, -1):
            named = mu_1d_named(p, kap)
            from_n2 = mu_1d_from_n2(p, kap)
            if named != from_n2:
                ok = False
        if p in (5, 7, 11, 13):
            id_rows[str(p)] = {
                "plus": str(mu_1d_named(p, 1)),
                "minus": str(mu_1d_named(p, -1)),
                "from_n2": str(mu_1d_from_n2(p, 1)),
                "wrong_p2": str(mu_1d_wrong_p2(p, 1)),
            }
    live = live_mu_1d()
    ok = ok and live[5]["ok"] and live[7]["ok"]
    ok = ok and mu_1d_named(5, 1) == Fraction(1, 15)
    ok = ok and mu_1d_named(5, 1) != mu_1d_wrong_p2(5, 1)
    ok = ok and mu_1d_named(7, 1) != T_bitight(7)
    ok = ok and live[5]["n_1d"] == 30
    ok = ok and live[7]["n_1d"] == 140
    return {
        "proved": bool(ok),
        "identities": id_rows,
        "live": {str(k): v for k, v in live.items()},
        "theorem": (
            "μ_1d=κ/(p(p−2)) on |κ|=1.  Fail: κ/p²; fail T at p=7."
        ),
    }


@lru_cache(maxsize=1)
def prove_B() -> dict:
    ok = True
    rows = {}
    for p in range(5, 48):
        if not is_prime(p):
            continue
        abs_mu = Fraction(1, p * (p - 2))
        abs_T = -T_bitight(p)
        le = abs_mu <= abs_T
        eq = abs_mu == abs_T
        if not le:
            ok = False
        if eq != (p == 5):
            ok = False
        if p in (5, 7, 11, 13):
            rows[str(p)] = {
                "abs_mu_1d": str(abs_mu),
                "abs_T": str(abs_T),
                "le": le,
                "eq": eq,
            }
    ok = ok and rows["5"]["eq"] is True
    ok = ok and rows["7"]["eq"] is False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "|μ_1d|≤|T| for p≥5, equality iff p=5.  "
            "Does not close Aut_e (free μ unnamed)."
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
        "A_full_named_p_law": False,
        "Delta_conn_p_rational": False,
        "note": (
            "μ_1d named; free/full-Ω μ and A_full still unnamed as "
            "p-laws.  Aut_e G>T open.  Flags unflipped."
        ),
    }


def main() -> dict:
    print("Prop 15.543  μ_1d=κ/(p(p−2)) on |κ|=1", flush=True)
    A = prove_A()
    print(
        f"  A μ_1d: {A['proved']} "
        f"p5={A['identities']['5']['plus']} p7={A['identities']['7']['plus']}",
        flush=True,
    )
    B = prove_B()
    print(f"  B |μ_1d|≤|T|: {B['proved']}", flush=True)
    C = prove_open()
    print(
        f"  C open: ND={C['type_I_multilevel_bad_case_ND_closed']} "
        f"3AB={C['type_I_aut_e_3AB_positive_general']}",
        flush=True,
    )
    out = {
        "prop": "15.543",
        "title": "Type+ 1D 3-point on |κ|=1 is κ/(p(p−2))",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "C": C,
        "proved": {
            "mu_1d_named": A["proved"],
            "mu_1d_le_T": B["proved"],
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
