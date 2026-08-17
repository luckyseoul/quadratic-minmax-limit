#!/usr/bin/env python3
"""
Prop 15.320 — χ-Kloosterman convolution is p times the Legendre
cubic sum: ∑_a χ(a) Kl(a)Kl(sa)=p ∑_x χ(x(x−1)(x−s)).
C(s) is then named in (p, H(s)).  A_□ / F_1(0) / Q_τ still
ensemble.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.319 flags.  Floor stays OPEN.

============================================================================
Setup.  Kl(n)=Kl(1,n).  χ=χ_p.  H(s)=∑_{x∈F_p} χ(x(x−1)(x−s)).
For s∉{0,1}, y²=x(x−1)(x−s) is the Legendre curve and
H(s)=−a_p(E_s).  ε=−χ_p(−1).  C(s)=(p+1)∑_{χ(a)=ε} Kl(a)Kl(sa).

============================================================================
Theorem A — PROVED (Gauss expansion; certified p=5,7,11,13).
  ∑_{a≠0} χ(a) Kl(a) Kl(sa) = p H(s) for every s∈F_p.
  Derivation: open both Kl, evaluate the a-sum as χ(λ)G, set
  t=y/x, the x-sum is a Gauss sum, and G²=χ(−1)p.  The remaining
  cubic character sum equals χ(−1) H(s) after u=−t, so the
  χ(−1) cancels.  At s=1, H(1)=−1 and the identity recovers
  ∑ χ Kl²=−p (15.319 C).  Fail-when-wrong: claim the convolution
  vanishes for s≠1; or claim H(s)=χ(s−1).  ∎

Theorem B — PROVED (15.319 B + Theorem A; certified p=5,7,11).
  C(s)=(p+1)/2 · ( κ(s) + ε p H(s) ) where
  κ(1)=p²−p−1 and κ(s≠1)=−(p+1).  In particular every C(s)
  is named once H(s) is.  Fail-when-wrong: drop the H term
  (C then constant off 1, false at p=11).  ∎

Theorem C — CERTIFIED; A_□ / F_1(0) still ensemble.
  |H(s)|≤2√p for s∉{0,1} (Hasse; checked p=5..31).  Substituting
  the named C into the 15.319 identity writes A_□ as a named-H
  linear form in the U-form factors F_{√c}(0).  Those F_c are
  still Max+ averages of Q, so A_□ and F_1(0) are not Max+-free
  numbers.  No poly-in-p formula for H on F_p^{×□} (values
  {0,±2,±4,…} with Hasse spread).  Q_τ unnamed.  phi_F not
  imported.  Fail-when-wrong: claim H≡0 on squares.  ∎

============================================================================
Backend: F_p Kloosterman / Legendre sums (Max+-free).
Writes evidence/e1_gmin_m4_prop15320.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15316 import chi_p_pm, kloosterman  # noqa: E402
from e1_gmin_m4_prop15319 import (  # noqa: E402
    C1_named,
    C_from_kl,
    chi_kl2_named,
    kl_conv_full_named,
    omega_chi_named,
    prove_open as prove_open_319,
)


def H_legendre(p: int, s: int) -> int:
    acc = 0
    for x in range(p):
        acc += chi_p_pm(p, (x * (x - 1) * (x - s)) % p)
    return int(acc)


def chiconv_live(p: int, s: int) -> float:
    acc = 0j
    for a in range(1, p):
        acc += (
            chi_p_pm(p, a)
            * kloosterman(p, 1, a)
            * kloosterman(p, 1, (s * a) % p)
        )
    return float(acc.real)


def chiconv_named(p: int, s: int) -> int:
    return p * H_legendre(p, s)


def C_named(p: int, s: int) -> float:
    eps = omega_chi_named(p)
    return (p + 1) / 2.0 * (kl_conv_full_named(p, s) + eps * chiconv_named(p, s))


def H_wrong_chi_sm1(p: int, s: int) -> int:
    return chi_p_pm(p, s - 1)


def prove_chiconv() -> dict:
    ok = True
    zero_hit = 0
    chi_sm1_hit = 0
    rows = {}
    for p in (5, 7, 11, 13):
        err = 0.0
        off_zero = True
        sm1_ok = True
        samples = {}
        for s in range(1, p):
            live = chiconv_live(p, s)
            named = chiconv_named(p, s)
            err = max(err, abs(live - named))
            if s != 1 and abs(live) > 1.0:
                off_zero = False
            if H_legendre(p, s) != H_wrong_chi_sm1(p, s):
                sm1_ok = False
            if s in (1, 2, p - 1):
                samples[str(s)] = {"live": live, "named": named, "H": H_legendre(p, s)}
        if err > 1e-6:
            ok = False
        if off_zero:
            zero_hit += 1
        if sm1_ok:
            chi_sm1_hit += 1
        rows[str(p)] = {
            "err": err,
            "H1": H_legendre(p, 1),
            "chiconv1": chiconv_named(p, 1),
            "off_vanishes": off_zero,
            "H_eq_chi_sm1": sm1_ok,
            "samples": samples,
        }
        if H_legendre(p, 1) != -1 or chiconv_named(p, 1) != chi_kl2_named(p):
            ok = False
    if zero_hit > 0:
        ok = False
    if chi_sm1_hit == 4:
        ok = False
    return {
        "proved": ok,
        "off1_vanishes_false": zero_hit == 0,
        "H_eq_chi_sm1_false": chi_sm1_hit < 4,
        "by_p": rows,
        "theorem": "∑ χ(a) Kl(a)Kl(sa)=p H(s), H=∑ χ(x(x−1)(x−s)).",
    }


def prove_C_named() -> dict:
    ok = True
    drop_H_hit = 0
    rows = {}
    for p in (5, 7, 11):
        err = 0.0
        drop_spread = 0.0
        cvals = {}
        dropped = {}
        for s in range(1, p):
            if chi_p_pm(p, s) != 1:
                continue
            live = C_from_kl(p, s)
            named = C_named(p, s)
            err = max(err, abs(live - named))
            dropped[s] = (p + 1) / 2.0 * kl_conv_full_named(p, s)
            cvals[s] = live
        others = [v for s, v in cvals.items() if s != 1]
        drop_others = [v for s, v in dropped.items() if s != 1]
        if others:
            drop_spread = max(drop_others) - min(drop_others) if drop_others else 0.0
            live_spread = max(others) - min(others)
            if drop_spread < 1.0 and live_spread > 1.0:
                drop_H_hit += 1
        if err > 1e-6:
            ok = False
        if abs(C_named(p, 1) - C1_named(p)) > 1e-6:
            ok = False
        rows[str(p)] = {
            "err": err,
            "C1_named": C_named(p, 1),
            "C1_319": C1_named(p),
            "live_spread_off1": (max(others) - min(others)) if others else 0.0,
            "drop_H_spread": drop_spread,
        }
    if drop_H_hit == 0:
        ok = False
    return {
        "proved": ok,
        "drop_H_fails": drop_H_hit > 0,
        "by_p": rows,
        "theorem": "C(s)=(p+1)/2 (κ(s)+ε p H(s)); drop H makes C constant off 1.",
    }


def prove_hasse_and_open() -> dict:
    ok = True
    zero_sq_hit = 0
    rows = {}
    for p in (5, 7, 11, 13, 17, 19, 23, 29, 31):
        bound = 2.0 * (p ** 0.5)
        mx = 0
        all_zero_sq = True
        for s in range(2, p):
            h = H_legendre(p, s)
            mx = max(mx, abs(h))
            if chi_p_pm(p, s) == 1 and h != 0:
                all_zero_sq = False
        if mx > bound + 1e-9:
            ok = False
        if all_zero_sq:
            zero_sq_hit += 1
        rows[str(p)] = {"max_abs_H": mx, "hasse": bound, "H_zero_on_squares": all_zero_sq}
    # fail-when-wrong: H≡0 on F_p^{×□} for every p — false at p=5,11,…
    if zero_sq_hit == 9:
        ok = False
    D319 = prove_open_319()
    return {
        "proved": ok,
        "H_zero_on_squares_false": zero_sq_hit < 9,
        "by_p": rows,
        "A_square_named": False,
        "F1_maxplus_free": False,
        "Q_tau_named_in_p": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "note": (
            "Hasse |H|≤2√p holds.  C is named in (p,H).  "
            "A_□ remains a linear form in the ensemble F_c.  "
            "phi_F not imported."
        ),
        "open_319": D319["note"],
    }


def prove_open() -> dict:
    C = prove_hasse_and_open()
    return {
        "proved": False,
        "phi_F_ge_6": C["phi_F_ge_6"],
        "rhat": C["rhat"],
        "e1": C["e1"],
        "A_square_named": False,
        "F1_maxplus_free": False,
        "Q_tau_named_in_p": False,
        "lambda_ge_6": False,
        "note": C["note"],
    }


def main() -> dict:
    print("Prop 15.320  χ-Kl conv = p H(s); C named; A_□ still ensemble", flush=True)
    A = prove_chiconv()
    print(
        f"  A chiconv: {A['proved']} off0={not A['off1_vanishes_false']} "
        f"H=chi_sm1={not A['H_eq_chi_sm1_false']}",
        flush=True,
    )
    B = prove_C_named()
    print(f"  B C-named: {B['proved']} dropH={B['drop_H_fails']}", flush=True)
    C = prove_hasse_and_open()
    print(
        f"  C Hasse/open: {C['proved']} H0sq={not C['H_zero_on_squares_false']} "
        f"phi_F={C['phi_F_ge_6']}",
        flush=True,
    )
    D = prove_open()
    out = {
        "prop": "15.320",
        "title": "χ-Kl conv = p H_Legendre; C named; A_□ / Q_τ open",
        "proved": {
            "chiconv_is_pH": A["proved"],
            "C_named_via_H": B["proved"],
            "hasse": C["proved"],
            "A_square_named": False,
            "F1_maxplus_free": False,
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
    path = ROOT / "evidence" / "e1_gmin_m4_prop15320.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
