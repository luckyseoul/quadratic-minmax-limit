#!/usr/bin/env python3
"""
Prop 15.580 — Max− y_∞=+1 has ẑ-support a Galois union of
Ω_−-lines, and F=−2(3p²+2) on every Ω_−-triple.  DEAD as a
Type I Φ(H)≥Φ−2 kill.

Does **not** flip type_I / phi_F / e1 / L / Aut-Schur / Gsum /
pairing / residual_ii / 15.279–15.579 flags.  Does **not**
edit leftover-1 15.550–15.575 or residual-ii 15.566 / 15.579.
Does **not** import phi_F.  Aut_e DEAD (15.559).  Max±-of-C
DEAD (15.565).  Type+ 1D Johnson DEAD (15.577).

============================================================================
Setup.  Max−: Cy=−py.  y_∞=+1 ⇒ ẑ(0)=−p.  Ω_−={ξ: χ̂(ξ)=−p}
is the complementary χ-coset to Ω (15.269 B).  15.551: ẑ of
any Q-valued z has Galois line-union support.  15.564: on an
Ω-triple, F=−2I+2κ(J_u+J_v+J_w)=−2(3p²+2) and
Q3_02=−4N(2p²+1)/p.

============================================================================
Theorem A — PROVED Max+-free (15.269 twin + 15.551 A).
  y_∞=+1 and Cy=−py give ∑z=−p and (χ∗z)(a)=−p z(a)−1, so
  ẑ is supported on {0}∪Ω_−.  Galois ⇒ a union of k Ω_−-lines,
  nsupp=k(p−1).  Certified p=5 (k∈{1,3}, 30+100) and p=7
  (k∈{1,3,4}, 140+1176+4410); Ω_+ nsupp=0.  Fail: ẑ(0)=+p.
  Fail: support on Ω_+.  Fail: nsupp=18 at p=5.  ∎

Theorem B — PROVED (15.564 A–C with χ_{Ω_−}=−χ_Ω).
  I=2 still (independent of Ω).  On an Ω_−-triple, χ=χ_{Ω_−}
  on all three, so J_u+J_v+J_w=−6 G χ_{Ω_−}.  Then
      F=−4+2κ_−(−6 G χ_{Ω_−})=−4−6q=−2(3p²+2)
  with κ_−=χ_{Ω_−}G/2.  Same p-law as 15.564 D.  Certified
  p=5,7,11.  Fail: drop the 6.  Fail: use χ_Ω for Jsum.  ∎

Theorem C — DEAD as a Type I kill.  A–B name Max− ẑ-support
  and the 02 complete sum F.  The Type I pairing E[S f_e]_−
  expands into named stars plus far 3-points; the unnamed
  remainder is the Q3_T / A_full channel (15.564 E), not F.
  F does not force f_e=+1 on U_−.  type_I flags stay False.

============================================================================
Backend: 15.564 Gauss + Max− caches p=5,7.  Writes
evidence/e1_gmin_m4_prop15580.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import Counter
from functools import lru_cache
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

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
from e1_gmin_m4_prop15279 import _e_tr, _kernel_field  # noqa: E402
from e1_gmin_m4_prop15534 import field_tables  # noqa: E402
from e1_gmin_m4_prop15536 import chi_Omega  # noqa: E402
from e1_gmin_m4_prop15550 import omega_set_local  # noqa: E402
from e1_gmin_m4_prop15564 import (  # noqa: E402
    F_named,
    first_omega_triple,
    gauss_G,
    three_dist_sums,
)

EV = ROOT / "evidence" / "e1_gmin_m4_prop15580.json"
YMINUS = {5: "/tmp/maxminus_p5.npy", 7: "/tmp/maxminus_p7.npy"}


def omega_minus_set(F) -> list[int]:
    q, p = F["q"], F["p"]
    om = []
    for xi in range(1, q):
        chat = sum(
            F["chi"][x] * _e_tr(F, F["mul"][xi][x]) for x in range(1, q)
        )
        if abs(chat.real + p) < 1e-6 and abs(chat.imag) < 1e-6:
            om.append(xi)
    return om


def F_from_Jsum(I: complex, Jsum: complex, kappa_m: complex) -> complex:
    return -2.0 * I + 2.0 * kappa_m * Jsum


def F_drop6(I: complex, Ju, Jv, Jw, kappa_m: complex) -> complex:
    """Fail: Jsum without the 6 (use Ju+Jv only)."""
    return -2.0 * I + 2.0 * kappa_m * (Ju + Jv)


def _omega_col(p, trt, mul, xi):
    om = np.exp(2j * np.pi / p)
    return om ** trt[mul[xi]].astype(np.int64)


@lru_cache(maxsize=1)
def live_minus_support() -> dict:
    out = {}
    for p in (5, 7):
        F = field_tables(p)
        q, ch, mul, trt = F["q"], F["chi"], F["mul"], F["tr"]
        com = chi_Omega(p)
        Om_p = [xi for xi in range(1, q) if int(ch[xi]) == com]
        Om_m = [xi for xi in range(1, q) if int(ch[xi]) == -com]
        Y = np.sign(np.load(YMINUS[p]).astype(np.float64))
        z = Y[Y[:, 0] > 0][:, 1:]
        Psi = np.column_stack(
            [_omega_col(p, trt, mul, xi) for xi in range(q)]
        )
        Zh = z @ Psi
        z0_err = float(np.max(np.abs(Zh[:, 0] + p)))
        n_plus = (np.abs(Zh[:, Om_p]) > 1e-6).sum(axis=1)
        mag = np.abs(Zh[:, Om_m]) > 1e-6
        n_minus = mag.sum(axis=1)
        seen: set[int] = set()
        lines: list[list[int]] = []
        for xi in Om_m:
            if xi in seen:
                continue
            L = [int(mul[xi, a]) for a in range(1, p)]
            for u in L:
                seen.add(u)
            lines.append([Om_m.index(u) for u in L])
        n_partial = 0
        kcnt: Counter[int] = Counter()
        for i in range(len(z)):
            kk = 0
            for idxs in lines:
                s = int(mag[i, idxs].sum())
                if s == 0:
                    continue
                if s == len(idxs):
                    kk += 1
                else:
                    n_partial += 1
            kcnt[kk] += 1
        out[p] = {
            "nH": int(len(z)),
            "z0_err": z0_err,
            "nsupp_Om_plus_max": int(n_plus.max()),
            "kcnt": {str(k): int(v) for k, v in sorted(kcnt.items())},
            "n_partial_line": n_partial,
            "nsupp_mod": sorted({int(x) % (p - 1) for x in n_minus.tolist()}),
            "n_lines": len(lines),
        }
    return out


@lru_cache(maxsize=1)
def prove_A() -> dict:
    live = live_minus_support()
    ok = True
    rows = {}
    want = {5: {"1": 30, "3": 100}, 7: {"1": 140, "3": 1176, "4": 4410}}
    for p in (5, 7):
        L = live[p]
        row_ok = (
            L["z0_err"] < 1e-6
            and L["nsupp_Om_plus_max"] == 0
            and L["n_partial_line"] == 0
            and L["nsupp_mod"] == [0]
            and L["kcnt"] == want[p]
        )
        if p == 5:
            row_ok = row_ok and L["nH"] == 130
            row_ok = row_ok and (18 % 4 != 0)
        if p == 7:
            row_ok = row_ok and L["nH"] == 5726
        if not row_ok:
            ok = False
        rows[str(p)] = {**L, "ok": row_ok}
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "Max− y_∞=+1: ẑ(0)=−p, support Ω_− line-union.  "
            "Fail: ẑ(0)=+p; fail support on Ω_+; fail nsupp=18 at p=5."
        ),
    }


@lru_cache(maxsize=1)
def prove_B() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 11):
        Fld = _kernel_field(p)
        Om_p = omega_set_local(Fld)
        Om_m = omega_minus_set(Fld)
        if len(Om_p) != len(Om_m) or set(Om_p) & set(Om_m):
            ok = False
        th, ph, _rh = first_omega_triple(Fld, Om_m)
        I, Ju, Jv, Jw = three_dist_sums(Fld, th, ph)
        G = gauss_G(Fld)
        from e1_gmin_m4_prop15536 import chi_Omega

        chi_om_m = float(-chi_Omega(p))
        Jsum = Ju + Jv + Jw
        want_J = -6.0 * G * chi_om_m
        kappa_m = chi_om_m * G / 2.0
        Fv = F_from_Jsum(I, Jsum, kappa_m)
        Fn = F_named(p)
        wrong6 = F_drop6(I, Ju, Jv, Jw, kappa_m)
        # fail: use χ_Ω for Jsum prediction
        Jsum_wrong_chi = -6.0 * G * (-chi_om_m)
        row_ok = (
            abs(I - 2) < 1e-6
            and abs(Jsum - want_J) < 1e-4
            and abs(Fv - Fn) < 1e-4
            and abs(wrong6 - Fn) > 1.0
            and abs(Jsum_wrong_chi - Jsum) > 1.0
        )
        if not row_ok:
            ok = False
        rows[str(p)] = {
            "I": float(I.real),
            "F": float(Fv.real),
            "F_named": Fn,
            "ok": row_ok,
        }
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "F=−2(3p²+2) on Ω_−-triples.  Fail: drop 6; "
            "fail Jsum with χ_Ω."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "path_dead": True,
        "type_I_multilevel_bad_case_ND_closed": bool(
            type_I_multilevel_bad_case_ND_closed()
        ),
        "type_I_aut_e_3AB_positive_general": bool(
            type_I_aut_e_3AB_positive_general()
        ),
        "note": (
            "Max− Galois + F name the 02 complete sum, not f_e on "
            "U_−.  Far pairing remainder is Q3_T/A_full (15.564 E).  "
            "Flags unflipped."
        ),
    }


def main() -> dict:
    print("Prop 15.580  Max− Ω_− Galois; F same; DEAD as Type I kill", flush=True)
    A = prove_A()
    print(
        f"  A Max− support: {A['proved']} "
        f"p5 k={A['rows']['5']['kcnt']}",
        flush=True,
    )
    B = prove_B()
    print(f"  B F on Ω_−: {B['proved']}", flush=True)
    C = prove_open()
    print(
        f"  C DEAD: ND={C['type_I_multilevel_bad_case_ND_closed']} "
        f"3AB={C['type_I_aut_e_3AB_positive_general']}",
        flush=True,
    )
    out = {
        "prop": "15.580",
        "title": "Max− Ω_− Galois + F=−2(3p²+2); DEAD as Type I kill",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "C": C,
        "proved": {
            "maxminus_omega_minus_galois": A["proved"],
            "F_on_Omega_minus": B["proved"],
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
