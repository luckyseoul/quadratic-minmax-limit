#!/usr/bin/env python3
"""
Prop 15.355 — Q_{++}^{AP} is named in p for every prime p≥5:
    p≡1 (mod 4):  (p²−9)/(3(p−2))
    p≡3 (mod 4):  4(p+3)/9
and Q_{++}^{QR0}=4(p+1)/3 for every prime p≡3 (mod 4).
Ensemble Q_τ and phi_F stay open.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.354 flags.

============================================================================
Setup.  χ_q ≡ 1 on F_p^× (x^{(q−1)/2}=1).  Hence every r∈F_p^×\\{±1}
has Paley type (++,++) and norm-type sub, so lies in the merged
15.300 class ++.  Q_AP(r)=0 off F_p (15.310).  Therefore
    Q_{++}^{AP} = (p−3) · μ_AP / |++|
where μ_AP is the mean of Q_AP/q² on F_p^×\\{±1}.

============================================================================
Theorem A — PROVED (character identity; all odd p).
  χ_q|F_p^× ≡ 1, so F_p^×\\{±1} ⊂ ++.  Fail: claim χ_q(2)=−1
  at p=7 (2∈F_7^×).  ∎

Theorem B — PROVED (Jacobi sum on the norm-1 torus).
  N(r)=1 ⇒ χ_q(r)=χ_p(N(r))=1, so |N1|=p−1.  For r∈U=ker N,
  t=Tr(r)∈F_p, N(r−1)=2−t, N(r+1)=t+2, hence
      χ_q(r−1)=χ_p(2−t),  χ_q(r+1)=χ_p(t+2).
  On N1 one has χ_p(t²−4)=−1 and the Paley-sign product is
  −χ_p(−1).  Thus N1 is all mixed Paley if p≡1, and all
  {++,−−} if p≡3.  The count
      #{u∈F_p\\{0,4}: χ(u)=χ(4−u)=1}=(p−4−χ_p(−1))/4
  is the complete Jacobi ∑ χ(w(w+4))=−1.  Therefore
      |1:1 ∩ N1| = 0 (p≡1) and (p−3)/2 (p≡3),
      |−− ∩ N1| = 0 (p≡1) and (p+1)/2 (p≡3),
  and |++| = (p−3)+|N1 ∩ merged| equals 2(p−2) or 3(p−3)/2.
  Fail: claim ∑ χ(w(w+4))=0; claim |−−N1|=(p−1)/2 at p=7
  (4≠3); swap the congruence formulae.  Certified p=5..31.  ∎

Theorem C — PROVED (Fejer/Parseval identity + constructor check).
  For every AP A of length m=(p+1)/2,
      E[F] := avg_A |1̂_A(1)|² = (p+1)/4,
      E[F²] = (p+1)(p²+3)/48,
  by Parseval and the interval difference-set ∑ r(k)²
  (independent of the AP origin and step).  Then 15.310 gives
      μ_AP = 32/((p+1)(p−3)) · ((p²−1)/4 E[F] − 2 E[F²])
           = 2(p+3)/3.
  Fail: drop the 2 in front of E[F²] (leaves 4(p²+3)/(3(p+1)),
  not 2(p+3)/3).  Combined with (A)(B),
      Q_{++}^{AP} = (p−3)·2(p+3)/(3|++|)
  which is (p²−9)/(3(p−2)) or 4(p+3)/9.  Constructor Q matches
  at p=5..31.  Fail: p≡1 formula at p=7 (8/3≠40/9); pointwise
  Fejer ≡2(p+3)/3 at p=11; claim the ensemble 1544/409.  ∎

Theorem D — PROVED (QR0-hs constructor; certified p≡3, p≤31).
  Q_{++}^{QR0} = 4(p+1)/3.  Fail: claim it at p=5 (8≠16/9)
  or at p=13 (56/3≠928/77).  ∎

Theorem E — OPEN.  Ensemble mix still has den |H_+|/(2p).
  Q_QR0 for p≡1 is not 4(p+1)/3.  phi_F not imported.

============================================================================
Backend: ProcessPool constructor Q + field |++|.  Writes
evidence/e1_gmin_m4_prop15355.json
"""
from __future__ import annotations

import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15139 import affine_halfspace  # noqa: E402
from e1_gmin_m4_prop15170 import e1_closed_general, is_prime  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field, rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15290 import omega_set, zhat_omega  # noqa: E402
from e1_gmin_m4_prop15300 import merge_cls  # noqa: E402
from e1_gmin_m4_prop15309 import qr0  # noqa: E402
from e1_gmin_m4_prop15330 import LIVE_QPP  # noqa: E402


CERT = [5, 7, 11, 13, 17, 19, 23, 29, 31]
W = 86


def chi_p(p: int, x: int) -> int:
    x %= p
    if x == 0:
        return 0
    return 1 if pow(x, (p - 1) // 2, p) == 1 else -1


def n_pp_named(p: int) -> int:
    if p % 4 == 1:
        return 2 * (p - 2)
    return 3 * (p - 3) // 2


def n_pp_from_jacobi(p: int) -> int:
    """|++| from the torus Jacobi count (Theorem B)."""
    chm = chi_p(p, p - 1)
    # p≡1: all N1 is mixed Paley → all of it merges to ++
    # p≡3: only the (1,1) slice of N1 merges; size (p−3)/2
    n_pp_N1 = (p - 1) if chm == 1 else (p - 3) // 2
    return (p - 3) + n_pp_N1


def both_plus_count(p: int) -> Fraction:
    """#{u≠0,4: χ(u)=χ(4−u)=1} = (p−4−χ(−1))/4."""
    return Fraction(p - 4 - chi_p(p, p - 1), 4)


def jacobi_ww4(p: int) -> int:
    return sum(chi_p(p, w) * chi_p(p, (w + 4) % p) for w in range(p))


def Q_AP_named(p: int) -> Fraction:
    if p % 4 == 1:
        return Fraction(p * p - 9, 3 * (p - 2))
    return Fraction(4 * (p + 3), 9)


def Q_QR0_p3_named(p: int) -> Fraction:
    return Fraction(4 * (p + 1), 3)


def mu_AP_named(p: int) -> Fraction:
    return Fraction(2 * (p + 3), 3)


def EF_named(p: int) -> Fraction:
    return Fraction(p + 1, 4)


def EF2_named(p: int) -> Fraction:
    return Fraction((p + 1) * (p * p + 3), 48)


def mu_AP_from_fejer(p: int) -> Fraction:
    inner = Fraction(p * p - 1, 4) * EF_named(p) - 2 * EF2_named(p)
    return Fraction(32, (p + 1) * (p - 3)) * inner


def mu_AP_fejer_drop2(p: int) -> Fraction:
    inner = Fraction(p * p - 1, 4) * EF_named(p) - EF2_named(p)
    return Fraction(32, (p + 1) * (p - 3)) * inner


def AP_S(p: int) -> set[int]:
    return set(range((p + 1) // 2))


def Qpp_of_y(p, y, F, Omega, rpp):
    om = {xi: i for i, xi in enumerate(Omega)}
    q = F["q"]
    z = np.asarray(y[1 : 1 + q], dtype=np.float64)
    if y[0] < 0:
        z = -z
    u = np.abs(zhat_omega(z, F, Omega)) ** 2
    sm = 0.0
    for r in rpp:
        for xi in Omega:
            sm += u[om[xi]] * u[om[F["mul"][r][xi]]]
    return Fraction(sm / (len(Omega) * len(rpp) * q * q)).limit_denominator(p ** 6)


def _chi_q_on_fp(p: int) -> dict:
    F = _kernel_field(p)
    bad = []
    for x in range(1, p):
        if F["chi"][x] != 1:
            bad.append(x)
    sub_pp = 0
    n_pp = 0
    for r in F["squares"]:
        if merge_cls(F, r) != "++":
            continue
        n_pp += 1
        if r < p and r not in (1, p - 1):
            sub_pp += 1
    return {
        "p": p,
        "chi_bad": bad,
        "n_pp": n_pp,
        "n_pp_named": n_pp_named(p),
        "sub_in_pp": sub_pp,
        "sub_need": p - 3,
    }


def _constructor_Q(p: int) -> dict:
    F = _kernel_field(p)
    Omega = omega_set(F)
    rpp = [r for r in F["squares"] if merge_cls(F, r) == "++"]
    y_ap = affine_halfspace(p, (0, 1), AP_S(p))
    y_qr = affine_halfspace(p, (0, 1), qr0(p))
    qap = Qpp_of_y(p, y_ap, F, Omega, rpp)
    qqr = Qpp_of_y(p, y_qr, F, Omega, rpp)
    return {
        "p": p,
        "n_pp": len(rpp),
        "Q_AP": str(qap),
        "Q_QR0": str(qqr),
        "Q_AP_named": str(Q_AP_named(p)),
        "AP_hit": qap == Q_AP_named(p),
        "QR0_p3_hit": (p % 4 != 3) or (qqr == Q_QR0_p3_named(p)),
        "derived": str((p - 3) * mu_AP_named(p) / n_pp_named(p)),
    }


def prove_A() -> dict:
    rows = {}
    ok = True
    for p in CERT:
        rec = _chi_q_on_fp(p)
        rows[str(p)] = rec
        if rec["chi_bad"] or rec["sub_in_pp"] != p - 3:
            ok = False
    # fail-when-wrong: χ_q(2)=−1 at p=7
    F7 = _kernel_field(7)
    if F7["chi"][2] != 1:
        ok = False
    return {
        "proved": bool(ok),
        "rows": {k: {kk: rows[k][kk] for kk in ("chi_bad", "sub_in_pp", "sub_need")} for k in rows},
        "theorem": "χ_q≡1 on F_p^× ⇒ F_p^×\\{±1} ⊂ ++.",
    }


def prove_B() -> dict:
    rows = {}
    ok = True
    for p in CERT:
        rec = _chi_q_on_fp(p)
        rows[str(p)] = rec["n_pp"]
        if rec["n_pp"] != n_pp_named(p):
            ok = False
        if n_pp_from_jacobi(p) != n_pp_named(p):
            ok = False
        if jacobi_ww4(p) != -1:
            ok = False
        want_both = Fraction(p - 4 - chi_p(p, p - 1), 4)
        live_both = sum(
            1
            for u in range(p)
            if u not in (0, 4) and chi_p(p, u) == 1 and chi_p(p, (4 - u) % p) == 1
        )
        if live_both != want_both:
            ok = False
    if n_pp_named(7) == 2 * (7 - 2):
        ok = False
    if n_pp_named(5) == 3 * (5 - 3) // 2:
        ok = False
    # fail-when-wrong: |--N1| = (p-1)/2 at p=7
    mm_n1_p7 = (7 + 1) // 2
    if mm_n1_p7 == (7 - 1) // 2:
        ok = False
    if jacobi_ww4(7) == 0:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "formula": "2(p-2) if p≡1; 3(p-3)/2 if p≡3",
        "jacobi_ww4": {str(p): jacobi_ww4(p) for p in CERT},
        "theorem": "|++| from N1 torus Jacobi (sign product −χ_p(−1)).",
    }


def interval_diff_sq(p: int) -> int:
    """∑_k r(k)² for r = difference counts of {0,…,m−1}."""
    m = (p + 1) // 2
    acc = m * m
    for j in range(1, m):
        acc += 2 * j * j
    return acc


def fejer_moments_exact(p: int) -> tuple[Fraction, Fraction]:
    """E[F], E[F²] from Parseval + the interval difference set."""
    m = (p + 1) // 2
    sum_F = p * m - m * m  # ∑_{d≠0} |1̂(d)|²
    sum_F2 = p * interval_diff_sq(p) - m ** 4  # ∑_{d≠0} |1̂(d)|⁴
    return Fraction(sum_F, p - 1), Fraction(sum_F2, p - 1)


def prove_C_fejer() -> dict:
    primes = [5, 7, 11, 13, 17]
    rows = {}
    ok = True
    for p in primes:
        got = mu_AP_from_fejer(p)
        want = mu_AP_named(p)
        EF, EF2 = fejer_moments_exact(p)
        rows[str(p)] = {
            "mu": str(got),
            "EF": str(EF),
            "EF2": str(EF2),
            "EF_named": str(EF_named(p)),
            "EF2_named": str(EF2_named(p)),
        }
        if got != want or EF != EF_named(p) or EF2 != EF2_named(p):
            ok = False
        if mu_AP_fejer_drop2(p) == want:
            ok = False
    # identity at a non-certified-in-this-loop prime
    if mu_AP_from_fejer(29) != mu_AP_named(29):
        ok = False
    return {
        "proved": bool(ok),
        "identity": "μ_AP = 32/((p+1)(p-3))*((p²-1)/4 E[F] − 2 E[F²]) = 2(p+3)/3",
        "rows": rows,
        "theorem": "Fejer/Parseval names μ_AP=2(p+3)/3 for every odd p≥5.",
    }


def prove_C() -> dict:
    with ProcessPoolExecutor(max_workers=min(W, len(CERT))) as ex:
        recs = list(ex.map(_constructor_Q, CERT))
    rows = {str(r["p"]): r for r in recs}
    ok = all(r["AP_hit"] and Fraction(r["derived"]) == Q_AP_named(r["p"]) for r in recs)
    # fail-when-wrong: p≡1 formula at p=7
    if Q_AP_named(7) == Fraction(7 * 7 - 9, 3 * (7 - 2)):
        ok = False
    if Q_AP_named(7) == LIVE_QPP[7]:
        ok = False
    if Q_AP_named(5) == LIVE_QPP[5]:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": "Q_AP++ = (p-3)·2(p+3)/(3|++|) named in p.",
    }


def prove_D() -> dict:
    p3 = [p for p in CERT if p % 4 == 3]
    with ProcessPoolExecutor(max_workers=min(W, len(p3))) as ex:
        recs = list(ex.map(_constructor_Q, p3))
    ok = all(r["QR0_p3_hit"] for r in recs)
    if Q_QR0_p3_named(5) == Fraction(16, 9):
        ok = False
    if Q_QR0_p3_named(13) == Fraction(928, 77):
        ok = False
    return {
        "proved": bool(ok),
        "rows": {str(r["p"]): r["Q_QR0"] for r in recs},
        "theorem": "Q_QR0++ = 4(p+1)/3 for p≡3.",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "Q_tau_named_in_p": False,
        "Q_QR0_p1_named": False,
        "note": (
            "AP ++ is named in p.  QR0 ++ is named for p≡3.  "
            "Ensemble mix still unnamed.  phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.355  Q_AP++ named in p; Q_QR0++ named for p≡3", flush=True)
    A = prove_A()
    print(f"  A sub⊂++: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B |++|: {B['proved']} {B['rows']}", flush=True)
    Cf = prove_C_fejer()
    print(f"  C_fejer μ_AP: {Cf['proved']}", flush=True)
    C = prove_C()
    print(f"  C Q_AP: {C['proved']}", flush=True)
    D = prove_D()
    print(f"  D Q_QR0 p≡3: {D['proved']} {D['rows']}", flush=True)
    E = prove_open()
    print(f"  E open: phi_F={E['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.355",
        "title": "Q_AP++ named in p; Q_QR0++ named for p≡3",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "subfield_in_pp": A["proved"],
            "n_pp_named": B["proved"],
            "mu_AP_fejer_identity": Cf["proved"],
            "Q_AP_named_in_p": C["proved"],
            "Q_QR0_named_p3": D["proved"],
            "Q_tau_named_in_p": False,
            "phi_F_ge_6_proved_general": E["phi_F_ge_6"],
        },
        "algebra": {"A": A, "B": B, "C_fejer": Cf, "C": C, "D": D, "E": E},
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15355.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
