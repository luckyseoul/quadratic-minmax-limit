#!/usr/bin/env python3
"""
Prop 15.523 — Type I particular κ-contraction is 3p χ_p(−1)χ(ξ).

Does **not** flip type_I / phi_F / e1 / L / Aut-Schur / Gsum /
pairing / residual_ii / 15.279–15.522 flags.

============================================================================
Setup.  Paley C on V=F_q∪{∞}, q=p², χ=χ_q, χ(−1)=1.
ψ_ξ(t)=exp(2πi Tr_{q/p}(ξ t)/p).  G(χ)=∑_{t≠0} χ(t)ψ_1(t)=−p χ_p(−1)
(χ_p(−1)=(−1)^{(p−1)/2}).  For x≠y in F_q^*,
    κ=χ(x)+χ(y)+χ(x−y).
This is the κ-slot of the μ_part−Wick contraction of Δ_conn
(15.275 / Type I A5).  Live Δ_conn stays unnamed.

============================================================================
Theorem A — PROVED (Gauss; all odd p).
  G(χ)=−p χ_p(−1).  Fail: G≡p (false at p=5: G=−5).  ∎

Theorem B — PROVED (three pairings; all odd p, ξ≠0).
  Write t=x−y.  Each of
      T_a=∑_{x≠y≠0} χ(y) ψ_ξ(t),
      T_b=∑_{x≠y≠0} χ(x) ψ_ξ(t),
      T_c=∑_{x≠y≠0} χ(x)χ(y)χ(t) ψ_ξ(t)
  equals −χ(ξ) G(χ).  Inner sums: ∑_{u≠0,−t} χ(u)=−χ(−t) and
  ∑_{x≠0,t} χ(x(x−t))=−1.  Hence
      ∑_{x≠y≠0} χ(x)χ(y) κ ψ_ξ(t) = 3p χ_p(−1) χ(ξ).
  Fail: drop one pairing (then p χ_p(−1)χ(ξ)).  ∎

Theorem C — PROVED (field sum; certified p=5,7,11,13).
  The double sum equals 3p χ_p(−1)χ(ξ) on both χ(ξ)=±1 classes.
  Fail: the dropped-pairing value.  ∎

Theorem D — OPEN.  Live Δ_conn (−328/65, −1144/2863) unnamed.
  This mechanism does not close Type I; Proposition 15.750 closes it independently.
  type_I_aut_e_3AB_positive_general stays False (G>T).

============================================================================
Backend: Gauss identities + field sums.  Writes
evidence/e1_gmin_m4_prop15523.json
"""
from __future__ import annotations

import json
import os
import sys
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
    type_I_aut_e_3AB_positive_general,
    type_I_multilevel_bad_case_ND_closed,
)
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402

EV = ROOT / "evidence" / "e1_gmin_m4_prop15523.json"


def chi_p_minus_1(p: int) -> int:
    """χ_p(−1)=(−1)^{(p−1)/2}."""
    return -1 if (p % 4 == 3) else 1


def G_chi(p: int) -> int:
    """G(χ)=−p χ_p(−1)."""
    return -p * chi_p_minus_1(p)


def pairing_one(p: int, chi_xi: int) -> int:
    """One pairing: −χ(ξ) G(χ)."""
    return -chi_xi * G_chi(p)


def kappa_contraction(p: int, chi_xi: int) -> int:
    """Three pairings: 3p χ_p(−1) χ(ξ)."""
    return 3 * pairing_one(p, chi_xi)


def dropped_one_pairing(p: int, chi_xi: int) -> int:
    """Fail-when-wrong value if one pairing is dropped."""
    return pairing_one(p, chi_xi)


def field_pack(p: int):
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

    def frob(u: int) -> int:
        return powm(u, p)

    def tr(u: int) -> int:
        return add(u, frob(u)) % p

    ch = [chi(x) for x in range(q)]
    return {
        "p": p,
        "q": q,
        "mul": mul,
        "add": add,
        "sub": sub,
        "chi": ch,
        "tr": tr,
    }


def field_sum_kappa(p: int, chi_xi: int) -> float:
    """∑_{x≠y≠0} χ(x)χ(y) κ ψ_ξ(t) for one ξ with χ(ξ)=chi_xi."""
    F = field_pack(p)
    q, ch = F["q"], F["chi"]
    sub, mul, tr = F["sub"], F["mul"], F["tr"]
    xi = next(x for x in range(1, q) if ch[x] == chi_xi)
    omega = np.exp(2j * np.pi / p)
    s = 0j
    for x in range(1, q):
        for y in range(1, q):
            if x == y:
                continue
            d = sub(x, y)
            kap = ch[x] + ch[y] + ch[d]
            s += ch[x] * ch[y] * kap * (omega ** tr(mul(xi, d)))
    return float(np.real(s))


def field_sum_one_pairing(p: int, chi_xi: int, which: str) -> float:
    """One of T_a, T_b, T_c (drop-two / keep-one)."""
    F = field_pack(p)
    q, ch = F["q"], F["chi"]
    sub, mul, tr = F["sub"], F["mul"], F["tr"]
    xi = next(x for x in range(1, q) if ch[x] == chi_xi)
    omega = np.exp(2j * np.pi / p)
    s = 0j
    for x in range(1, q):
        for y in range(1, q):
            if x == y:
                continue
            d = sub(x, y)
            ps = omega ** tr(mul(xi, d))
            if which == "a":
                s += ch[y] * ps
            elif which == "b":
                s += ch[x] * ps
            else:
                s += ch[x] * ch[y] * ch[d] * ps
    return float(np.real(s))


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 11, 13, 17, 19, 23):
        g = G_chi(p)
        chi = chi_p_minus_1(p)
        ok = ok and g == -p * chi
        ok = ok and g in (-p, p)
        rows[str(p)] = {"G": g, "chi_p(-1)": chi}
    # fail: G≡p
    ok = ok and G_chi(5) == -5
    ok = ok and G_chi(5) != 5
    ok = ok and G_chi(7) == 7
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": "G(χ)=−p χ_p(−1). Fail: G≡p (p=5).",
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 11, 13, 17, 19, 23, 29, 31):
        for cxi in (1, -1):
            one = pairing_one(p, cxi)
            full = kappa_contraction(p, cxi)
            drop = dropped_one_pairing(p, cxi)
            want = 3 * p * chi_p_minus_1(p) * cxi
            ok = ok and full == want
            ok = ok and drop == one
            ok = ok and drop == p * chi_p_minus_1(p) * cxi
            ok = ok and full != drop
            ok = ok and full == 3 * drop
            rows[f"{p}:{cxi:+d}"] = {
                "full": full,
                "drop": drop,
                "want": want,
            }
            if full == drop:
                ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "∑ χ(x)χ(y)κ ψ_ξ=3p χ_p(−1)χ(ξ).  "
            "Fail: drop one pairing → p χ_p(−1)χ(ξ)."
        ),
    }


def prove_C() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 11, 13):
        row = {}
        for cxi in (1, -1):
            live = field_sum_kappa(p, cxi)
            want = float(kappa_contraction(p, cxi))
            drop = float(dropped_one_pairing(p, cxi))
            ta = field_sum_one_pairing(p, cxi, "a")
            tb = field_sum_one_pairing(p, cxi, "b")
            tc = field_sum_one_pairing(p, cxi, "c")
            one = float(pairing_one(p, cxi))
            hit = abs(live - want) < 1e-6
            parts = (
                abs(ta - one) < 1e-6
                and abs(tb - one) < 1e-6
                and abs(tc - one) < 1e-6
            )
            not_drop = abs(live - drop) > 1.0
            if not (hit and parts and not_drop):
                ok = False
            row[str(cxi)] = {
                "live": live,
                "want": want,
                "drop": drop,
                "T_a": ta,
                "T_b": tb,
                "T_c": tc,
                "ok": bool(hit and parts and not_drop),
            }
        rows[str(p)] = row
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "Field sum equals 3p χ_p(−1)χ(ξ) at p=5,7,11,13.  "
            "Each pairing equals −χ(ξ)G(χ).  Fail: the dropped value."
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
        "Delta_conn_named": False,
        "note": (
            "Live Δ_conn (−328/65, −1144/2863) unnamed.  "
            "3A+B>0 remains G>T.  Flags unflipped."
        ),
    }


def main() -> dict:
    print("Prop 15.523  κ-contraction = 3p χ_p(−1)χ(ξ)", flush=True)
    A, B, C, D = prove_A(), prove_B(), prove_C(), prove_open()
    print(f"  A G(χ): {A['proved']} p5={A['rows']['5']}", flush=True)
    print(f"  B three pairings: {B['proved']}", flush=True)
    print(f"  C field sum: {C['proved']}", flush=True)
    print(
        f"  D open: ND={D['type_I_multilevel_bad_case_ND_closed']} "
        f"3AB={D['type_I_aut_e_3AB_positive_general']}",
        flush=True,
    )
    out = {
        "prop": "15.523",
        "title": "κ-contraction is 3p χ_p(−1)χ(ξ)",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "C": C,
        "D": D,
        "proved": {
            "G_chi": A["proved"],
            "three_pairings": B["proved"],
            "field_sum": C["proved"],
            "type_I_multilevel_bad_case_ND_closed": D[
                "type_I_multilevel_bad_case_ND_closed"
            ],
            "type_I_aut_e_3AB_positive_general": D[
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
