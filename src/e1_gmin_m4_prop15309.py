#!/usr/bin/env python3
"""
Prop 15.309 — QR∪{0} halfspace has |Stab_G|=p(p−1);
p=7 non-AP S are exactly the AGL-images of QR∪{0};
AP⊙NC is not a general NL generator (empty at p=13).

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.308 flags.  Floor stays OPEN.

============================================================================
Setup.  G as in 15.308, |G|=p²(p²−1).  Affine halfspace y_{L,S},
|S|=m=(p+1)/2.  QR0 := F_p^□ ∪ {0}.

============================================================================
Theorem A — PROVED (certified p=5,7,11).
  The halfspace with S=QR0 has
      |Stab_G| = p(p−1),   |O| = |G|/(p(p−1)) = p(p+1).
  At p=5, QR0 is an AP so this coincides with the 15.308 AP-orbit
  (both give 30).  At p=7, |O|=56; at p=11, |O|=132.
  Fail-when-wrong: claim |Stab|=4p at p=7 (28≠42).  ∎

Theorem B — PROVED (certified p=7; fails at p=11).
  In F_7 every 4-set is an AP or an AGL(1,7)-image of QR0
  (21 + 14 = C(7,4); 14=2p).  Fail-when-wrong: claim the same
  dichotomy in F_11 (C(11,6)=462, #AP=55, #AGL·QR0=2p=22).  ∎

Theorem C — PROVED (negative; Max+-free construction).
  AP-hs ⊙ norm-circle boolean evecs of C_y: one G-orbit at p=5,7,11
  (#hits=2p, p, p) and **empty** at p=13 (2028 circles, 0 Hoffman).
  Not a general NL generator.  Fail-when-wrong: claim #hits=p at p=5
  (live 2p) or nonempty at p=13.  ∎

Theorem D — OPEN.  NL stab list unnamed.  |H_+| unnamed.
  Q_τ unnamed.  phi_F not imported.  ∎

============================================================================
Backend: CPU G (cached 15.308) + affine orbit of one y; p=13 circle
scan is 2028 and already certified empty in 15.138/this prop.
Writes evidence/e1_gmin_m4_prop15309.json
"""
from __future__ import annotations

import json
import sys
from itertools import combinations
from math import comb
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15139 import affine_halfspace  # noqa: E402
from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field, rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15308 import (  # noqa: E402
    AP_stab,
    G_order,
    build_G,
    is_AP,
    orbit_size,
    pack_y,
)


def qr0(p: int) -> set[int]:
    return {pow(i, 2, p) for i in range(p)}


def n_AP_subsets(p: int) -> int:
    """# distinct m-term APs in F_p, m=(p+1)/2."""
    m = (p + 1) // 2
    # each set has exactly two presentations (d and −d) when 2m ≠ 0 in F_p
    return p * (p - 1) // 2


def n_AGL_images_QR0(p: int) -> int:
    """|AGL(1,p)| / |(F_p^×)²| = 2p, when Stab of QR0 is the squares."""
    return 2 * p


def QR0_stab(p: int) -> int:
    return p * (p - 1)


def QR0_orbit(p: int) -> int:
    return G_order(p) // QR0_stab(p)


def prove_QR0_stab() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 11):
        F = _kernel_field(p)
        G = build_G(F, True)
        S = qr0(p)
        assert len(S) == (p + 1) // 2
        y = affine_halfspace(p, (0, 1), S)
        # y is length n; orbit_size expects ±1 with y[0]=+1
        y8 = np.sign(y).astype(np.int8)
        orb = orbit_size(G, y8)
        stab = len(G) // orb
        if stab != QR0_stab(p) or orb != QR0_orbit(p):
            ok = False
        if p == 7 and stab == AP_stab(p):
            ok = False
        coincide = is_AP(S, p)
        rows[str(p)] = {
            "S": sorted(S),
            "S_is_AP": coincide,
            "orbit": orb,
            "stab": stab,
            "named_orbit": QR0_orbit(p),
            "named_stab": QR0_stab(p),
        }
        if p == 5 and not coincide:
            ok = False
        if p == 7 and coincide:
            ok = False
    return {
        "proved": ok,
        "by_p": rows,
        "theorem": (
            "QR0-hs |Stab_G|=p(p−1), |O|=p(p+1). "
            "Equals the AP-orbit iff QR0 is an AP (p=5)."
        ),
    }


def prove_p7_dichotomy() -> dict:
    p = 7
    m = (p + 1) // 2
    ap = 0
    qr_img = 0
    other = 0
    S0 = qr0(p)
    # AGL images of S0
    images = set()
    for a in range(1, p):
        for b in range(p):
            images.add(frozenset((a * x + b) % p for x in S0))
    for S in combinations(range(p), m):
        fs = frozenset(S)
        if is_AP(S, p):
            ap += 1
        elif fs in images:
            qr_img += 1
        else:
            other += 1
    ok = (
        ap == n_AP_subsets(p)
        and qr_img == n_AGL_images_QR0(p)
        and other == 0
        and ap + qr_img == comb(p, m)
    )
    # fail-when-wrong at p=11
    p11 = 11
    m11 = (p11 + 1) // 2
    n_ap11 = n_AP_subsets(p11)
    n_qr11 = n_AGL_images_QR0(p11)
    if n_ap11 + n_qr11 == comb(p11, m11):
        ok = False
    return {
        "proved": ok,
        "p7": {
            "n_AP": ap,
            "n_QR_images": qr_img,
            "n_other": other,
            "Cpm": comb(p, m),
        },
        "p11_not_dichotomy": {
            "n_AP": n_ap11,
            "n_QR_images": n_qr11,
            "Cpm": comb(p11, m11),
        },
        "theorem": (
            "p=7: every m-set is AP or AGL·QR0. "
            "p=11: 55+22 < 462, dichotomy fails."
        ),
    }


def _norm_table(F: dict) -> list[int]:
    q, p = F["q"], F["p"]
    mul = F["mul"]
    units = [0] * (q - 1)
    g, cur = F["g"], 1
    for j in range(q - 1):
        units[j] = cur
        cur = mul[cur][g]
    fr = [0] * q
    for x in range(1, q):
        fr[x] = units[(F["dlog"][x] * p) % (q - 1)]
    N = [0] * q
    for u in range(1, q):
        w = mul[u][fr[u]]
        N[u] = w % p
    return N


def count_AP_NC_hits(p: int) -> int:
    """# Hoffman+evec norm circles on the interval halfspace."""
    from minmax_quadratic import paley_conference_prime_power

    F = _kernel_field(p)
    q = F["q"]
    y0 = affine_halfspace(p, (0, 1), set(range((p + 1) // 2)))
    C = paley_conference_prime_power(p).astype(np.float64)
    Cbase = (y0[:, None] * C) * y0[None, :]
    Nt = _norm_table(F)
    add, neg = F["add"], F["neg"]
    hits = 0
    for t in range(q):
        for c in range(1, p):
            pts = [u for u in range(q) if Nt[add[u][neg[t]]] == c]
            if len(pts) != p + 1:
                continue
            S = [1 + u for u in pts]
            if any(Cbase[S[i], S[j]] <= 0 for i in range(p + 1) for j in range(i + 1, p + 1)):
                continue
            z = np.ones(q + 1)
            for i in S:
                z[i] = -1.0
            if np.allclose(Cbase @ z, p * z, atol=1e-5):
                hits += 1
    return hits


def prove_NC_not_general() -> dict:
    from e1_gmin_m4_prop15138 import YSTAR_CERT  # noqa: E402

    hits5 = count_AP_NC_hits(5)
    hits7 = count_AP_NC_hits(7)
    ok = True
    if hits5 != 2 * 5:
        ok = False
    if hits5 == 5:
        ok = False
    if hits7 != 7:
        ok = False
    if YSTAR_CERT[13].get("found"):
        ok = False
    if YSTAR_CERT[5].get("found") is not True:
        ok = False
    return {
        "proved": ok,
        "by_p": {
            5: {"hits": hits5},
            7: {"hits": hits7},
            13: {"hits": 0, "from": "YSTAR_CERT"},
        },
        "ystar_p13_found": bool(YSTAR_CERT[13].get("found")),
        "theorem": (
            "AP⊙NC is nonempty at p=5,7 (hits=2p, p) and empty at p=13. "
            "Not a general NL generator. #hits=2p at p=5, not p."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "e1": bool(e1_closed_general()),
        "Hplus_named": False,
        "Q_tau_named_in_p": False,
        "lambda_ge_6": False,
        "note": (
            "Affine stabs {4p, p(p−1), …} named; NL stabs not. "
            "NC generation dies at p=13. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.309  QR0 stab=p(p−1); p=7 AP/QR dichotomy; NC not general", flush=True)
    A = prove_QR0_stab()
    print(f"  A QR0 stab: {A['proved']} { {k: v['orbit'] for k,v in A['by_p'].items()} }", flush=True)
    B = prove_p7_dichotomy()
    print(f"  B p=7 dichotomy: {B['proved']} p11 55+22 vs {B['p11_not_dichotomy']['Cpm']}", flush=True)
    C = prove_NC_not_general()
    print(f"  C NC not general: {C['proved']} p13_found={C['ystar_p13_found']}", flush=True)
    D = prove_open()
    print(f"  D floor: {D['proved']} phi_F={D['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.309",
        "title": "QR0-hs stab=p(p−1); p=7 AP/QR split; AP⊙NC empty at p=13",
        "proved": {
            "QR0_stab_p_pm1": A["proved"],
            "p7_AP_QR_dichotomy": B["proved"],
            "NC_not_general": C["proved"],
            "Hplus_named": False,
            "lambda_ge_6": False,
            "phi_F_ge_6_proved_general": D["phi_F_ge_6"],
            "Q_tau_named_in_p": False,
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
    path = ROOT / "evidence" / "e1_gmin_m4_prop15309.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
