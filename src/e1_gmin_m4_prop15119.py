#!/usr/bin/env python3
"""
Prop 15.119 — Residual budget dictionary (m4f / EW2 / ED4 closed forms);
weight-enumerator structure of Max+ dots; halfspace pin criterion.

Continues 15.117–15.118 toward N_ED4_SUF. Does **not** soft-close residual.
No class_key (F19).

============================================================================
Setup. Paley Max+, n=p²+1, d=n/2. ED4:=E[(y·z)⁴], W:=∑_{i<j} y_i y_j z_i z_j
so D=y·z satisfies D²=n+2W, E[W]=n/2, ED4=3n²+4 E[W²] (15.113).
ρ=m₄−κ/p²=ρ_min+δ, ED4=ED4_flat+24δ² (15.112). Path C primary:
  δ² ≤ room_hyp/24  ⇔  orth ≤ room_hyp  ⇔  ED4 ≤ ED4_bud.

============================================================================
Theorem A (simplified ED4 budgets) — PROVED Fraction, all primes p>√5.
  ED4_flat = 4(p²−3)(p²+1)(3p²+1)/(p²−5)
           = 4(p²−3)n(3p²+1)/(p²−5).
  ED4_bud  = 4(3p⁸+6p⁶−104p⁴+138p²−75)/((p²−5)(p²+1))
           = ED4_flat + room_hyp
  with room_hyp = 24·room_hyp/24 and room_hyp/24 as in Prop 15.117.
  Proof.  ED4_flat=wick_lo+proj=12n²−48n+64n(p²−3)/(p²−5); clear denoms.
  ED4_bud=ED4_flat+24·room_hyp/24; expand.  ∎

Theorem B (E[W²] residual channel) — PROVED.
  E[W²] = (ED4 − 3n²)/4.
  Write EW2_flat := (ED4_flat − 3n²)/4, EW2_bud := (ED4_bud − 3n²)/4.
  Closed forms (p>√5):
      EW2_flat = (p²+1)(9p⁴−20p²+3)/(4(p²−5)),
      EW2_bud  = (9p⁸+30p⁶−380p⁴+594p²−285)/(4(p²−5)(p²+1)).
  Path C residual ⇔ E[W²] ≤ EW2_bud.
  Moreover EW2_bud − EW2_flat = room_hyp/4 = 6·(room_hyp/24), matching
  E[W²]−EW2_flat = 6δ² (Prop 15.113.D).  ∎

Theorem C (⟨m₄,f_y⟩ ↔ ED4 dictionary) — PROVED.
  For every y∈Max+, with ED4(y):=E_z[(y·z)⁴],
      ⟨m₄, f_y⟩ = ED4(y)/24 + n(4−3n)/6 + n(n−2)/8.
  Proof.  ⟨m₄,f_y⟩=E_z e₄(y⊙z) and e₄(s)=s⁴/24+((4−3n)/12)s²+n(n−2)/8
  (Prop 15.116.A); use E[D²]=2n, E[D]=E[D³]=0.  ∎
  Average form: ⟨m₄,m₄⟩-style pairing gives the same with global ED4.
  Define
      m4f_flat := n(n−1)(n−2)/(8p²) + 4(p⁴−1)/(3(p²−5)),
      m4f_bud  := m4f_flat + room_hyp/24
  (the δ=0 and Path-C-budget values of ⟨m₄,f⟩ under bgf=0 / 15.118).
  Closed forms (p>√5):
      m4f_flat = (p−1)(p+1)(p²+1)(3p²+17)/(24(p²−5)),
      m4f_bud  = (p−1)(p+1)(3p²−5)(p⁴+20p²−61)/(24(p²−5)(p²+1)).
  Consistency: m4f_bud = ⟨m₄,f⟩ evaluated at ED4=ED4_bud via Thm C.  ∎

Theorem D (three-form residual equivalence) — PROVED.
  For all primes p≥5 the following are equivalent:
    (i)   δ² ≤ room_hyp/24;
    (ii)  ED4 ≤ ED4_bud;
    (iii) E[W²] ≤ EW2_bud.
  If moreover Q_δ(y):=⟨δ,f_y⟩ is constant on Max+ (certified p=3,5), then
  also equivalent:
    (iv)  ⟨m₄, f_y⟩ ≤ m4f_bud for every y (e.g. halfspace);
    (v)   ⟨m₄, f_hs⟩ ≤ m4f_bud for the halfspace boolean evec.
  Proof.  (i)⇔(ii) by ED4=ED4_flat+24δ² and ED4_bud=ED4_flat+room_hyp.
  (ii)⇔(iii) by Thm B.  Under Q_δ constant: δ²=Q_δ=⟨m₄,f⟩−m4f_flat
  (Prop 15.118.F with ⟨ρ_min,f⟩ closed), so (i)⇔(iv)⇔(v).  ∎

Theorem E (weight-enumerator structure) — PROVED structure + cert.
  For y,z∈Max+ on Paley n=p²+1:
    (a) y·z ≡ n ≡ 2 (mod 4) (n even, boolean products ±1);
    (b) y=±z ⇒ y·z=±n; else |y·z| ≤ D_max := p²−2p−1 (Prop 15.99);
    (c) the pair measure is antipodal: A(t)=A(−t).
  Certified full pair spectrum:
    p=3: dots ∈ {±10,±2} with mults (per fixed y) {1,5} on positive side;
    p=5: dots ∈ {±26,±14,±10,±6,±2} mults {1,13,20,36,60}.
  Global pair-average ED4 equals ED4_bud at p=3,5 (saturates Path C budget);
  at p=7 (prior census) ED4 < ED4_bud with δ²/room_hyp/24 = 124875/669124.
  The crude envelope ED4 ≤ 2n D_max² is strictly larger than ED4_bud for
  p≥5 (too weak — same obstruction as Prop 15.116.G).  ∎

Theorem F (orth·N at saturation) — CERTIFIED.
  orth=24δ².  At p=5 (δ²=room_hyp/24, N=|Max+|=260):
      orth · N = 147456 ∈ ℤ.
  At p=3: orth=0.  At p=7: orth·N = 12889497600/4499 (not integral).
  No general closed form for |Max+| or orth·N claimed.  ∎

============================================================================
OPEN residual (honest):
  Prove δ² ≤ room_hyp/24 for all primes p≥5 (equivalently any of Thm D
  (i)–(iii)).  Preferred: pin ⟨m₄,f_hs⟩ ≤ m4f_bud via Gauss sums / full
  Max+ weight enumerator closed form, or bound E[W²] independently.
  F19: no class_key.  Full E_4p energy of f_y too crude (15.115).

L OPEN. hyp_residual_for_all_p=false. ed4_le_suf_for_all_p=false.
F19/F20: Fraction algebra + Max+ cert p=3,5; GPU unused.
Writes evidence/e1_gmin_m4_prop15119.json
"""
from __future__ import annotations

import json
import sys
from collections import Counter
from fractions import Fraction
from itertools import combinations
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15100 import n_of, proj_kappa2, wick_lo  # noqa: E402
from e1_gmin_m4_prop15102 import delta2_budget_hyp, rho_min2  # noqa: E402
from e1_gmin_m4_prop15112 import ed4_bud, ed4_flat, ed4_suf  # noqa: E402
from e1_gmin_m4_prop15117 import room_hyp_over_24  # noqa: E402
from e1_gmin_m4_prop15118 import rho_min_dot_f_if_bgf_zero  # noqa: E402


def is_prime(p: int) -> bool:
    if p < 2:
        return False
    if p % 2 == 0:
        return p == 2
    k = 3
    while k * k <= p:
        if p % k == 0:
            return False
        k += 2
    return True


# ---------------------------------------------------------------------------
# Theorem A: ED4 budgets simplified
# ---------------------------------------------------------------------------

def ed4_flat_closed(p: int) -> Fraction:
    """4(p²−3)(p²+1)(3p²+1)/(p²−5)."""
    if p * p == 5:
        raise ValueError("p²≠5")
    return Fraction(
        4 * (p * p - 3) * (p * p + 1) * (3 * p * p + 1),
        p * p - 5,
    )


def ed4_bud_closed(p: int) -> Fraction:
    """4(3p⁸+6p⁶−104p⁴+138p²−75)/((p²−5)(p²+1))."""
    if p <= 3:
        return ed4_flat_closed(p) if p == 3 else ed4_flat(p)
    if p * p == 5:
        raise ValueError("p²≠5")
    num = 3 * p**8 + 6 * p**6 - 104 * p**4 + 138 * p**2 - 75
    return Fraction(4 * num, (p * p - 5) * (p * p + 1))


def prove_ed4_budgets(primes: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes:
        if p * p == 5:
            continue
        flat_c = ed4_flat_closed(p)
        flat_ref = ed4_flat(p)
        if p <= 3:
            bud_c = flat_c
            bud_ref = ed4_flat(p)
        else:
            bud_c = ed4_bud_closed(p)
            bud_ref = ed4_bud(p)
            # also flat + room_hyp
            room_hyp = 24 * room_hyp_over_24(p)
            bud_from_room = flat_c + room_hyp
            if bud_c != bud_from_room:
                ok = False
        match = flat_c == flat_ref and bud_c == bud_ref
        rows[str(p)] = {
            "flat_closed": str(flat_c),
            "flat_ref": str(flat_ref),
            "bud_closed": str(bud_c),
            "bud_ref": str(bud_ref),
            "match": match,
        }
        if not match:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "ED4_flat=4(p²−3)n(3p²+1)/(p²−5); "
            "ED4_bud=4(3p⁸+6p⁶−104p⁴+138p²−75)/((p²−5)(p²+1))."
        ),
        "by_p": rows,
    }


# ---------------------------------------------------------------------------
# Theorem B: EW2 channel
# ---------------------------------------------------------------------------

def ew2_of_ed4(p: int, ed4: Fraction) -> Fraction:
    n = n_of(p)
    return (Fraction(ed4) - 3 * n * n) / 4


def ew2_flat(p: int) -> Fraction:
    """(p²+1)(9p⁴−20p²+3)/(4(p²−5))."""
    if p * p == 5:
        raise ValueError("p²≠5")
    return Fraction(
        (p * p + 1) * (9 * p**4 - 20 * p * p + 3),
        4 * (p * p - 5),
    )


def ew2_bud(p: int) -> Fraction:
    """(9p⁸+30p⁶−380p⁴+594p²−285)/(4(p²−5)(p²+1))."""
    if p <= 3:
        return ew2_flat(p)
    if p * p == 5:
        raise ValueError("p²≠5")
    num = 9 * p**8 + 30 * p**6 - 380 * p**4 + 594 * p * p - 285
    return Fraction(num, 4 * (p * p - 5) * (p * p + 1))


def prove_ew2_channel(primes: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes:
        if p * p == 5:
            continue
        ef = ew2_flat(p)
        eb = ew2_bud(p)
        from_flat_ed4 = ew2_of_ed4(p, ed4_flat_closed(p))
        from_bud_ed4 = ew2_of_ed4(p, ed4_bud_closed(p) if p > 3 else ed4_flat_closed(p))
        if p > 3:
            gap = eb - ef
            gap_target = room_hyp_over_24(p) * 6  # room_hyp/4 = 6*(room/24)
            gap_alt = Fraction(24 * room_hyp_over_24(p), 4)  # room_hyp/4
        else:
            gap = Fraction(0)
            gap_target = Fraction(0)
            gap_alt = Fraction(0)
        match = (
            ef == from_flat_ed4
            and eb == from_bud_ed4
            and gap == gap_target == gap_alt
        )
        rows[str(p)] = {
            "EW2_flat": str(ef),
            "EW2_bud": str(eb),
            "gap": str(gap),
            "gap_eq_room_hyp_over_4": gap == gap_alt,
            "match": match,
        }
        if not match:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "E[W²]=(ED4−3n²)/4; residual ⇔ E[W²]≤EW2_bud; "
            "EW2_bud−EW2_flat=room_hyp/4=6·(room_hyp/24)."
        ),
        "by_p": rows,
    }


# ---------------------------------------------------------------------------
# Theorem C: m4f ↔ ED4
# ---------------------------------------------------------------------------

def m4f_of_ed4(p: int, ed4: Fraction) -> Fraction:
    """⟨m₄,f_y⟩ from ED4(y) via e₄ poly + design moments."""
    n = n_of(p)
    return (
        Fraction(ed4) / 24
        + Fraction(n * (4 - 3 * n), 6)
        + Fraction(n * (n - 2), 8)
    )


def m4f_flat(p: int) -> Fraction:
    """n(n−1)(n−2)/(8p²) + 4(p⁴−1)/(3(p²−5))."""
    n = n_of(p)
    return Fraction(n * (n - 1) * (n - 2), 8 * p * p) + rho_min_dot_f_if_bgf_zero(p)


def m4f_bud(p: int) -> Fraction:
    """m4f_flat + room_hyp/24."""
    if p <= 3:
        return m4f_flat(p)
    return m4f_flat(p) + room_hyp_over_24(p)


def m4f_flat_closed(p: int) -> Fraction:
    """(p−1)(p+1)(p²+1)(3p²+17)/(24(p²−5))."""
    if p * p == 5:
        raise ValueError("p²≠5")
    return Fraction(
        (p - 1) * (p + 1) * (p * p + 1) * (3 * p * p + 17),
        24 * (p * p - 5),
    )


def m4f_bud_closed(p: int) -> Fraction:
    """(p−1)(p+1)(3p²−5)(p⁴+20p²−61)/(24(p²−5)(p²+1))."""
    if p <= 3:
        return m4f_flat_closed(p)
    if p * p == 5:
        raise ValueError("p²≠5")
    return Fraction(
        (p - 1)
        * (p + 1)
        * (3 * p * p - 5)
        * (p**4 + 20 * p * p - 61),
        24 * (p * p - 5) * (p * p + 1),
    )


def prove_m4f_dictionary(primes: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes:
        if p * p == 5:
            continue
        mf = m4f_flat(p)
        mb = m4f_bud(p)
        mfc = m4f_flat_closed(p)
        mbc = m4f_bud_closed(p)
        # consistency with ED4 dictionary
        from_ed4_flat = m4f_of_ed4(p, ed4_flat_closed(p))
        from_ed4_bud = m4f_of_ed4(
            p, ed4_bud_closed(p) if p > 3 else ed4_flat_closed(p)
        )
        # e4 poly at design moments: E e4 = poly moments
        n = n_of(p)
        # E[e4(D)] with E[D^4]=ed4_flat, E[D^2]=2n
        e4_flat_check = (
            ed4_flat_closed(p) / 24
            + Fraction(4 - 3 * n, 12) * Fraction(2 * n)
            + Fraction(n * (n - 2), 8)
        )
        match = (
            mf == mfc == from_ed4_flat == e4_flat_check
            and mb == mbc == from_ed4_bud
        )
        # also m4f_flat parts
        flat_pair = Fraction(n * (n - 1) * (n - 2), 8 * p * p)
        rmf = rho_min_dot_f_if_bgf_zero(p)
        if mf != flat_pair + rmf:
            match = False
        rows[str(p)] = {
            "m4f_flat": str(mf),
            "m4f_bud": str(mb),
            "closed_flat": str(mfc),
            "closed_bud": str(mbc),
            "from_ed4_bud": str(from_ed4_bud),
            "match": match,
        }
        if not match:
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "⟨m₄,f_y⟩=ED4(y)/24+n(4−3n)/6+n(n−2)/8; "
            "m4f_bud=(p−1)(p+1)(3p²−5)(p⁴+20p²−61)/(24(p²−5)(p²+1))."
        ),
        "by_p": rows,
    }


# ---------------------------------------------------------------------------
# Theorem D: equivalence (structure)
# ---------------------------------------------------------------------------

def prove_equivalence() -> dict:
    return {
        "proved": True,
        "theorem": (
            "For p≥5: δ²≤room_hyp/24 ⇔ ED4≤ED4_bud ⇔ E[W²]≤EW2_bud. "
            "If Q_δ constant on Max+: also ⇔ ⟨m₄,f_y⟩≤m4f_bud for all y "
            "(halfspace pin)."
        ),
        "note": (
            "Q_δ constant certified at p=3,5 (15.118); at p=7 three ED4(y) "
            "types (15.114) so use pair-average ED4 / E[W²] forms."
        ),
    }


# ---------------------------------------------------------------------------
# Theorem E: weight enumerator structure + cert
# ---------------------------------------------------------------------------

def d_max(p: int) -> int:
    """p² − 2p − 1."""
    return p * p - 2 * p - 1


def crude_ed4_envelope(p: int) -> Fraction:
    """2n D_max² (extremal 2-point with |D|≤D_max, E[D²]=2n)."""
    n = n_of(p)
    return Fraction(2 * n * d_max(p) ** 2)


def prove_we_structure(primes: list[int]) -> dict:
    rows = {}
    ok = True
    for p in primes:
        n = n_of(p)
        env = crude_ed4_envelope(p)
        bud = ed4_bud(p) if p > 3 else ed4_flat(p)
        weak = env > bud  # too weak for closing
        rows[str(p)] = {
            "n": n,
            "n_mod_4": n % 4,
            "D_max": d_max(p),
            "crude_envelope_2n_Dmax2": str(env),
            "ED4_bud": str(bud),
            "envelope_too_weak": weak if p >= 5 else None,
            "parity_n_eq_2_mod_4": n % 4 == 2,
        }
        if n % 4 != 2:
            ok = False
        if p >= 5 and not weak:
            # envelope must be strictly above bud for p≥5
            ok = False
    return {
        "proved": ok,
        "theorem": (
            "Max+ dots: y·z ≡ 2 (mod 4); |y·z|≤p²−2p−1 off ±pairs; "
            "antipodal. Crude ED4≤2n D_max² too weak vs ED4_bud for p≥5."
        ),
        "by_p": rows,
    }


def _load_Y(p: int):
    try:
        from e1_gmin_cr_classify import load_maxplus

        Y = load_maxplus(p)
        if Y is not None:
            return np.asarray(Y, dtype=np.float64)
    except Exception:
        pass
    for path in (
        Path(f"/tmp/maxplus_p{p}.npy"),
        Path(f"/tmp/e1_p{p}/maxplus.npy"),
    ):
        if path.is_file():
            return np.load(path).astype(np.float64)
    return None


def certify_weight_enum(p: int) -> dict:
    Y = _load_Y(p)
    if Y is None:
        return {"p": p, "skipped": True}

    from minmax_quadratic import halfspace_boolean_vector, paley_conference_prime_power

    N, n = Y.shape
    C = paley_conference_prime_power(p).astype(np.float64)
    hs = halfspace_boolean_vector(p)
    # halfspace is Max+
    hs_err = float(np.max(np.abs(C @ hs - p * hs)))
    # dots from halfspace
    dots_hs = Y @ hs
    spec_hs = Counter(int(round(float(x))) for x in dots_hs)
    # global ED4
    # for N=260: full Gram OK; use pair mean of fourth powers via hs + full if small
    if N <= 300:
        G = Y @ Y.T
        ed4_global = Fraction(int(np.sum(np.rint(G).astype(np.int64) ** 4)), N * N)
        # full spectrum from first row (should match hs up to Aut)
        spec_row = Counter(int(round(float(x))) for x in G[0])
    else:
        # pair-average via all rows of fourth powers without full G if needed
        s4 = 0
        for i in range(N):
            d = Y @ Y[i]
            s4 += int(np.sum(np.rint(d).astype(np.int64) ** 4))
        ed4_global = Fraction(s4, N * N)
        spec_row = Counter(int(round(float(x))) for x in Y @ Y[0])

    ed4_hs = Fraction(
        int(np.sum(np.rint(dots_hs).astype(np.int64) ** 4)), N
    )
    d2 = (ed4_global - ed4_flat(p)) / 24
    hyp = room_hyp_over_24(p) if p > 3 else Fraction(0)
    m4f_hs = m4f_of_ed4(p, ed4_hs)
    m4f_avg = m4f_of_ed4(p, ed4_global)

    # parity / support checks
    all_dots = list(spec_hs.keys())
    parity_ok = all((t % 4) == (n % 4) for t in all_dots)
    max_off = max(abs(t) for t in all_dots if abs(t) != n) if N > 2 else 0
    dmax_ok = max_off <= d_max(p)

    # expected spectra
    expected = {
        3: {10: 1, 2: 5, -2: 5, -10: 1},
        5: {
            26: 1,
            14: 13,
            10: 20,
            6: 36,
            2: 60,
            -2: 60,
            -6: 36,
            -10: 20,
            -14: 13,
            -26: 1,
        },
    }
    spec_match = True
    if p in expected:
        spec_match = dict(spec_hs) == expected[p]

    saturates_bud = ed4_global == (ed4_bud(p) if p > 3 else ed4_flat(p))
    le_bud = ed4_global <= (ed4_bud(p) if p > 3 else ed4_flat(p)) + Fraction(0)

    return {
        "p": p,
        "N": N,
        "hs_Cy_py_err": hs_err,
        "hs_is_evec": hs_err < 1e-8,
        "dot_spectrum_hs": {str(k): v for k, v in sorted(spec_hs.items())},
        "dot_spectrum_row0": {str(k): v for k, v in sorted(spec_row.items())},
        "spectrum_matches_expected": spec_match,
        "parity_ok": parity_ok,
        "dmax_ok": dmax_ok,
        "max_off_antipode": max_off,
        "D_max": d_max(p),
        "ED4_global": str(ed4_global),
        "ED4_hs": str(ed4_hs),
        "ED4_flat": str(ed4_flat(p)),
        "ED4_bud": str(ed4_bud(p) if p > 3 else ed4_flat(p)),
        "delta2": str(d2),
        "room_hyp_24": str(hyp),
        "delta2_le_hyp": d2 <= hyp if p > 3 else d2 == 0,
        "saturates_bud": saturates_bud,
        "ED4_le_bud": le_bud,
        "m4f_hs": str(m4f_hs),
        "m4f_avg": str(m4f_avg),
        "m4f_bud": str(m4f_bud(p)),
        "m4f_hs_le_bud": m4f_hs <= m4f_bud(p),
        "m4f_hs_eq_bud": m4f_hs == m4f_bud(p),
        "EW2_global": str(ew2_of_ed4(p, ed4_global)),
        "EW2_bud": str(ew2_bud(p)),
        "EW2_le_bud": ew2_of_ed4(p, ed4_global) <= ew2_bud(p),
        "orth_times_N": str(24 * d2 * N),
    }


# ---------------------------------------------------------------------------
# Theorem F: orth·N
# ---------------------------------------------------------------------------

def prove_orth_N_note() -> dict:
    # p=5 equality
    d2_5 = room_hyp_over_24(5)
    orth_N_5 = 24 * d2_5 * 260
    # p=7 from prior census
    d2_7 = Fraction(19180800, 1840091)
    orth_N_7 = 24 * d2_7 * 11452
    r7 = d2_7 / room_hyp_over_24(7)
    return {
        "proved_structure": True,
        "p5_orth_N": str(orth_N_5),
        "p5_orth_N_int": orth_N_5.denominator == 1 and int(orth_N_5) == 147456,
        "p7_orth_N": str(orth_N_7),
        "p7_delta2_over_room": str(r7),
        "p7_ratio_float": float(r7),
        "note": (
            "orth·N=147456 at p=5 saturation; p=7 ratio "
            "δ²/(room_hyp/24)=124875/669124; no general N formula."
        ),
    }


def prove_general_status() -> dict:
    return {
        "ed4_budgets_proved": True,
        "ew2_channel_proved": True,
        "m4f_dictionary_proved": True,
        "equivalence_proved": True,
        "we_structure_proved": True,
        "hyp_residual_for_all_p": False,
        "ed4_le_suf_for_all_p": False,
        "delta_le_rho_min_for_all_p": False,
        "sixteen_N_proved_for_all_p": False,
        "open_residual": (
            "Prove δ²≤room_hyp/24 (⇔ ED4≤ED4_bud ⇔ E[W²]≤EW2_bud) for all "
            "primes p≥5. Have closed budgets + m4f dictionary + halfspace pin "
            "when Q_δ constant; need independent Max+ weight enumerator / "
            "Gauss-sum UB on E[W²] or ⟨m₄,f_hs⟩."
        ),
    }


def main() -> dict:
    primes = [p for p in range(3, 60) if is_prime(p)]
    print(
        "Prop 15.119: residual budgets m4f/EW2/ED4; weight enum; halfspace pin",
        flush=True,
    )

    a = prove_ed4_budgets(primes)
    print(f"  ED4 budgets proved={a['proved']}", flush=True)
    b = prove_ew2_channel(primes)
    print(f"  EW2 channel proved={b['proved']}", flush=True)
    c = prove_m4f_dictionary(primes)
    print(f"  m4f dictionary proved={c['proved']}", flush=True)
    d = prove_equivalence()
    e = prove_we_structure(primes)
    print(f"  WE structure proved={e['proved']}", flush=True)
    f = prove_orth_N_note()
    print(f"  orth·N p5={f['p5_orth_N']} int? {f['p5_orth_N_int']}", flush=True)

    cert = {}
    for p in (3, 5):
        cert[str(p)] = certify_weight_enum(p)
        if not cert[str(p)].get("skipped"):
            cp = cert[str(p)]
            print(
                f"  cert p={p}: ED4={cp['ED4_global']} sat_bud={cp['saturates_bud']} "
                f"m4f_hs={cp['m4f_hs']} spec_ok={cp['spectrum_matches_expected']}",
                flush=True,
            )

    status = prove_general_status()
    out = {
        "prop": "15.119",
        "title": (
            "Residual budget dictionary (m4f/EW2/ED4); "
            "weight-enumerator structure; halfspace pin criterion"
        ),
        "L_status": "OPEN",
        "settlement_note": (
            "No soft-close: δ²≤room_hyp/24 not proved for general p≥5. "
            "Shipped closed residual budgets + equivalence + WE structure."
        ),
        "proved": {
            "ed4_budgets": a["proved"],
            "ew2_channel": b["proved"],
            "m4f_dictionary": c["proved"],
            "equivalence": d["proved"],
            "we_structure": e["proved"],
            "hyp_residual_for_all_p": False,
            "ed4_le_suf_for_all_p": False,
            "sixteen_N_for_all_p": False,
        },
        "algebra": {
            "ed4_budgets": a,
            "ew2_channel": b,
            "m4f_dictionary": c,
            "equivalence": d,
            "we_structure": e,
            "orth_N": f,
        },
        "cert": cert,
        "status": status,
        "closed_forms": {
            "ED4_flat": "4(p^2-3)(p^2+1)(3p^2+1)/(p^2-5)",
            "ED4_bud": "4(3p^8+6p^6-104p^4+138p^2-75)/((p^2-5)(p^2+1))",
            "EW2_flat": "(p^2+1)(9p^4-20p^2+3)/(4(p^2-5))",
            "EW2_bud": "(9p^8+30p^6-380p^4+594p^2-285)/(4(p^2-5)(p^2+1))",
            "m4f_flat": "(p-1)(p+1)(p^2+1)(3p^2+17)/(24(p^2-5))",
            "m4f_bud": "(p-1)(p+1)(3p^2-5)(p^4+20p^2-61)/(24(p^2-5)(p^2+1))",
            "room_hyp_24": "4(p^2-9)(p^2-1)^2/(3(p^2-5)(p^2+1))",
        },
        "backend": "CPU Fraction algebra + Max+ cert p=3,5; GPU unused (F20)",
        "F19": "no class_key",
    }

    path = ROOT / "evidence" / "e1_gmin_m4_prop15119.json"
    path.write_text(json.dumps(out, indent=2) + "\n")
    print(f"  wrote {path}", flush=True)
    print(f"  L_status={out['L_status']} residual_closed={status['hyp_residual_for_all_p']}", flush=True)
    return out


if __name__ == "__main__":
    main()
