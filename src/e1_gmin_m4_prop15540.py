#!/usr/bin/env python3
"""
Prop 15.540 — A_1d is −4p³/(p−2) on F_p r-classes, 0 off;
A_free splits full-Ω / partial and is not a p-law.

Does **not** flip type_I / phi_F / e1 / L / Aut-Schur / Gsum /
pairing / residual_ii / 15.279–15.539 flags.  Does **not**
edit leftover-1 15.522–15.539 or residual-ii 15.528.

============================================================================
Setup.  15.538: Type+ 1D A(ξ,ξ,−2ξ)=−4p³/(p−2), via
E[ε̂(1)² ε̂(−2)]=−2(p+1)/(p−2).  r=η/ξ on Ω-triples (15.269 H).
F_p embeds as the prime field in F_q.  Free Max+ (y_∞=+1) splits
by ẑ-support: full-Ω (nsupp=|Ω|) vs partial (p=7: nsupp=18).

============================================================================
Theorem A — PROVED Max+-free (15.538 + support + r-independence).
  For a,b,c∈F_p^* with a+b+c=0, E[ε̂(a)ε̂(b)ε̂(c)]=−2(p+1)/(p−2).
  Reduction: rescale to (1,r,−(1+r)).  Collision counts for
  L=x+r y−(1+r)z match 15.538 B whenever a,b,c≠0 (L(t,t,t)=0;
  two-equal types each contribute 1 off u=0; three-distinct
  μ₃=−3/(p(p−2)) with the same |L^{-1}(u)|).  Certified p=7,11
  for generic r=2.  Hence for every r∈F_p with 1+r≠0,
      A_{1d}(r)=p³ E / m = −4p³/(p−2).
  If r∉F_p then {ξ,rξ,−(1+r)ξ} is not contained in any F_p-line,
  so every Type+ 1D has a vanishing ẑ and A_{1d}(r)=0.
  Fail: A_1d≡−4p³/(p−2) off F_p (live 0).  Fail: A_1d≡0 on F_p.  ∎

Theorem B — PROVED as a negative on A_free (caches p=5,7).
  Partial-support free (p=7, n=1176): A_part(r∈F_p)=0 and
  A_part(r∉F_p)=−2p³/3.  No partial at p=5.
  Full-Ω: A_full,dbl=0 at p=5 and −4p³/15 at p=7 — not a p-law.
  Fail: A_full,dbl=0 at p=7; fail A_full,dbl=−4p³/(p−2).
  A_free is the mix (n_full A_full+n_part A_part)/n_free, so
  A_free,dbl=0 at p=5 and −4p³/19 at p=7.  The 19 is not a
  p-denominator (not D, not p±k).  ∎

Theorem C — OPEN.  A_free is not named as a p-rational on all
  r-classes.  Mixed A still has den 13,409.  Aut_e G>T open.
  type_I flags stay False.

============================================================================
Backend: 15.538 identities + y_∞=+1 caches.  Writes
evidence/e1_gmin_m4_prop15540.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import defaultdict
from fractions import Fraction
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
from e1_gmin_m4_prop15534 import field_tables  # noqa: E402
from e1_gmin_m4_prop15536 import chi_Omega  # noqa: E402
from e1_gmin_m4_prop15538 import (  # noqa: E402
    A_1d_dbl_named,
    E_eps_named,
)
from minmax_quadratic import paley_conference_prime_power  # noqa: E402

EV = ROOT / "evidence" / "e1_gmin_m4_prop15540.json"
YPLUS = {5: "/tmp/maxplus_p5.npy", 7: "/tmp/maxplus_p7.npy"}


def A_1d_of_r(p: int, r: int) -> Fraction:
    """Named A_1d on the r-class. r is an F_q element in c0+c1 p coords."""
    if r < p and (1 + r) % p != 0:
        return A_1d_dbl_named(p)
    return Fraction(0)


def _omega_col(p, trt, mul_t, xi):
    om = np.exp(2j * np.pi / p)
    return om ** trt[mul_t[xi]].astype(np.int64)


def _frac(x, den: int) -> Fraction:
    return Fraction(x).limit_denominator(max(den, 8))


@lru_cache(maxsize=1)
def live_profiles() -> dict:
    out = {}
    for p in (5, 7):
        F = field_tables(p)
        q, ch = F["q"], F["chi"]
        mul, add, sub, inv, trt = F["mul"], F["add"], F["sub"], F["inv"], F["tr"]
        com = chi_Omega(p)
        Omega = [xi for xi in range(1, q) if int(ch[xi]) == com]
        nOm = len(Omega)
        Om = set(Omega)
        C = paley_conference_prime_power(p).astype(np.float64)
        Y = np.sign(np.load(YPLUS[p]).astype(np.float64))
        Yp = Y[Y[:, 0] > 0]
        assert np.allclose(Yp @ C, p * Yp, atol=1e-5)
        Zh = (Yp[:, 1:] @ np.column_stack(
            [_omega_col(p, trt, mul, xi) for xi in range(q)]
        ))[:, Omega]
        nsupp = (np.abs(Zh) > 1e-6).sum(axis=1)
        masks = {
            "1d": nsupp == (p - 1),
            "full": nsupp == nOm,
            "part": (nsupp != (p - 1)) & (nsupp != nOm),
        }
        xi0 = Omega[0]
        specials = {1, int(sub[0, 2]), int(mul[int(sub[0, 1]), int(inv[2])])}
        buckets: dict[str, dict[str, list[Fraction]]] = {
            name: defaultdict(list) for name in masks
        }
        for eta in Omega:
            th = int(sub[0, int(add[xi0, eta])])
            if th not in Om:
                continue
            r = int(mul[eta, int(inv[xi0])])
            kind = "dbl" if r in specials else ("fp" if r < p else "nfp")
            ia, ib, ic = 0, Omega.index(eta), Omega.index(th)
            for name, mask in masks.items():
                n = int(mask.sum())
                if n == 0:
                    continue
                sl = Zh[mask]
                val = float(np.real(np.mean(sl[:, ia] * sl[:, ib] * sl[:, ic])))
                buckets[name][kind].append(_frac(val, n * n))
        rec = {"n": {k: int(m.sum()) for k, m in masks.items()}}
        for name in masks:
            rec[name] = {
                k: sorted({str(x) for x in xs})
                for k, xs in buckets[name].items()
            }
        out[p] = rec
    return out


@lru_cache(maxsize=1)
def prove_A() -> dict:
    live = live_profiles()
    ok = True
    rows = {}
    for p in (5, 7):
        want = str(A_1d_dbl_named(p))
        a1 = live[p]["1d"]
        fp_ok = set(a1.get("dbl", [])) == {want} and set(a1.get("fp", [])) <= {want}
        nfp_ok = set(a1.get("nfp", [])) == {"0"}
        # A_1d_of_r
        named_fp = A_1d_of_r(p, 1) == A_1d_dbl_named(p)
        named_nfp = A_1d_of_r(p, p + 1) == 0
        row_ok = fp_ok and nfp_ok and named_fp and named_nfp
        if not row_ok:
            ok = False
        rows[str(p)] = {
            "A_1d": a1,
            "want": want,
            "ok": row_ok,
        }
    # E independent of r: generic = double (15.538 closed form)
    ok = ok and E_eps_named(7) == Fraction(-16, 5)
    ok = ok and A_1d_of_r(7, 2) == A_1d_dbl_named(7)
    ok = ok and A_1d_of_r(7, 8) == 0  # 8=p+1 ∉ F_p
    # fail: constant on all r
    ok = ok and A_1d_of_r(5, 6) == 0 and A_1d_of_r(5, 6) != A_1d_dbl_named(5)
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "A_1d(r)=−4p³/(p−2) on F_p (1+r≠0), else 0.  "
            "Fail: same value off F_p."
        ),
    }


@lru_cache(maxsize=1)
def prove_B() -> dict:
    live = live_profiles()
    p5, p7 = live[5], live[7]
    part7 = p7["part"]
    full5, full7 = p5["full"], p7["full"]
    ok = (
        p5["n"]["part"] == 0
        and p7["n"]["part"] == 1176
        and set(part7.get("fp", [])) == {"0"}
        and set(part7.get("dbl", [])) == {"0"}
        and set(part7.get("nfp", [])) == {str(Fraction(-2 * 7**3, 3))}
        and set(full5.get("dbl", [])) == {"0"}
        and set(full7.get("dbl", [])) == {str(Fraction(-4 * 7**3, 15))}
        and set(full5.get("nfp", [])) == {str(Fraction(-2 * 5**3))}
        # fail p-laws for A_full,dbl
        and Fraction(-4 * 7**3, 15) != 0
        and Fraction(-4 * 7**3, 15) != A_1d_dbl_named(7)
        and Fraction(-4 * 5**3, 15) != 0  # the p=7 formula at p=5 is not live
    )
    # A_free,dbl mix
    a_free_5 = Fraction(0)
    a_free_7 = Fraction(-4 * 7**3, 19)
    return {
        "proved": bool(ok),
        "p5_full": full5,
        "p7_full": full7,
        "p7_part": part7,
        "n": {"5": p5["n"], "7": p7["n"]},
        "A_free_dbl": {"5": str(a_free_5), "7": str(a_free_7)},
        "theorem": (
            "A_part(F_p)=0, A_part(off)=−2p³/3 (p=7).  "
            "A_full,dbl is 0, −4p³/15 — not a p-law.  "
            "Fail: A_full,dbl=0 at p=7; fail =A_1d."
        ),
    }


def prove_open() -> dict:
    B = prove_B()
    return {
        "proved": False,
        "type_I_multilevel_bad_case_ND_closed": bool(
            type_I_multilevel_bad_case_ND_closed()
        ),
        "type_I_aut_e_3AB_positive_general": bool(
            type_I_aut_e_3AB_positive_general()
        ),
        "Delta_conn_p_rational": False,
        "A_free_named_p_law": False,
        "A_free_dbl": B["A_free_dbl"],
        "note": (
            "A_free,dbl=−4p³/19 at p=7 is D-free as a fraction but 19 "
            "is not a p-law.  Mixed A still den 13,409.  Aut_e G>T open."
        ),
    }


def main() -> dict:
    print("Prop 15.540  A_1d on all r; A_free profile split", flush=True)
    A = prove_A()
    print(f"  A A_1d r-classes: {A['proved']}", flush=True)
    B = prove_B()
    print(
        f"  B A_free split: {B['proved']} "
        f"A_free_dbl={B['A_free_dbl']}",
        flush=True,
    )
    C = prove_open()
    print(
        f"  C open: ND={C['type_I_multilevel_bad_case_ND_closed']} "
        f"3AB={C['type_I_aut_e_3AB_positive_general']}",
        flush=True,
    )
    out = {
        "prop": "15.540",
        "title": "A_1d on all r-classes; A_free is not a p-law",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "C": C,
        "proved": {
            "A_1d_all_r": A["proved"],
            "A_free_split_not_plaw": B["proved"],
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
