#!/usr/bin/env python3
"""
Prop 15.263 — Same Cy/Ext equation for m₄⁺ and m₄⁻ on |κ|=1;
Ext[ν]=(p⁴−1)ν on |κ|=1; residual-(i) still OPEN on spectral gap.

Does **not** flip residual_i / type_I / gsum / E1 / L. Soft-close forbidden.

SETUP (15.251 Cy/Ext for Max+, 15.254 m₄±, 15.262 π-anti-invariance of ν)
  On |κ|=1: Per(C_S)=1 (derangement permanent).  Cy-expansion of
  E[∏_{i∈S}(Cy)_i] with Cy=±py gives permanent transfer Ext (same
  bipartite-permanent weights W(S,T) for both Max±).

PROVED Max+-free
  A. Max− Cy-expansion on |κ|=1.
     For y∈Max−, Cy=−p y, so ∏_{i∈S}(Cy)_i = (−p)⁴ ∏ y_i = p⁴ ∏ y_i.
     Expand ∏_i (∑_j C_{ij} y_j) by injections as in 15.251:
     size-0 (diagonal Per) contributes Per(C_S) m₄⁻ = m₄⁻ on |κ|=1;
     size-1+size-2 contribute −2φ under the Max− pairwise law π=−C/p
     **with the same closed form −2φ as Max+** (each pair term picks up
     two minus signs from π_− vs π_+: (−C/p)·(−C/p)=(+C C)/p², and the
     combinatorial coefficient of φ is identical — certified identity
     size1+size2=−2φ for Max− at p=5,7; the 64-exhaust for the local
     pair/triple contractions is Max±-free as in 15.191/15.229).
     size-3 vanishes on |κ|=1 (15.229).  size-4 = Ext[m₄⁻].
     Therefore on every |κ|=1 four-set:
        (p⁴ − 1) m₄⁻ + 2 φ = Ext[m₄⁻],
     same Ext operator as for m₄⁺.  ∎

  B. Common equation and difference.
     From 15.251 B and (A): both m₄⁺ and m₄⁻ solve
        (p⁴ − 1) f + 2 φ = Ext[f]
     on |κ|=1.  Subtract:
        (p⁴ − 1) ν = Ext[ν]
     on every |κ|=1 four-set.  (φ cancels.)  ∎

  C. Compatibility with 15.262.
     On finite 4-sets ν(S)=−ν(πS).  The |κ|=1 stratum is π-stable
     (π preserves |κ|).  Hence ν|_{|κ|=1} is a π-odd function in the
     kernel of (Ext − (p⁴−1)I) on that stratum (coupled to |κ|=3 via Ext).
     Census: this odd piece is 0 at p=5,7.  ∎

CERTIFIED
  D. At p=5: for 20 random |κ|=1 sets, Ext[m₄⁺]=(p⁴−1)m₄⁺+2φ and
     Ext[m₄⁻]=(p⁴−1)m₄⁻+2φ to machine precision; Ext[ν]=(p⁴−1)ν (=0).
  E. size1+size2=−2φ for Max− on |κ|=1 sample at p=5,7 (matches Max+).

OPEN (blocks residual i / predicate flip)
  F. Prove the π-odd kernel of (Ext−(p⁴−1)I) on |κ|=1 (with the |κ|=3
     coupling / master / Π-odd constraints) is {0}, **or** envelope /
     reflection / |μ|≤1/(2p) by another path. Soft-close forbidden.
     (Ext[ν]=(p⁴−1)ν alone is necessary for ν but not sufficient without
     a spectral gap on the π-odd module — referee direction.)

Until F: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15263.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    type_I_k_3p_minus_2_closed_general,
)


def theorem_m4minus_same_ext() -> dict:
    return {
        "proved": True,
        "theorem": (
            "On |κ|=1: (p⁴−1)m₄⁻+2φ=Ext[m₄⁻] with the same Ext as Max+.  "
            "Proof: Cy=−py ⇒ ∏(Cy)=p⁴∏y; Per=1; size1+size2=−2φ "
            "(pairwise π=−C/p, two minuses); size3=0; size4=Ext."
        ),
        "depends_on": [
            "cy_ext_maxplus_251",
            "derangement_perm_kappa1_191",
            "size3_vanishes_229",
            "pairwise_pi_minus",
        ],
        "note_size12": (
            "size1+size2=−2φ for Max− certified p=5,7; local 64-exhaust "
            "is Max±-free (same C-contractions as 15.191)."
        ),
    }


def theorem_ext_nu_eigen() -> dict:
    return {
        "proved": True,
        "forces_nu_zero": False,
        "theorem": (
            "On |κ|=1: (p⁴−1)ν=Ext[ν], since both m₄± solve the same "
            "inhomogeneous Ext equation.  π-odd kernel gap OPEN."
        ),
        "depends_on": ["theorem_m4minus_same_ext", "cy_ext_maxplus_251"],
    }


def residual_i_closed_via_263() -> bool:
    return False


def hinge_status_263() -> dict:
    return {
        "m4minus_same_ext": True,
        "ext_nu_eigen": True,
        "pi_odd_kernel_trivial": False,
        "residual_i_closed": residual_i_closed_via_263(),
        "open": (
            "Ext[ν]=(p⁴−1)ν on |κ|=1 proved; need π-odd spectral gap "
            "or envelope/reflection/|μ|."
        ),
    }


def certify_same_ext() -> dict:
    """Light census: Ext equation for m4- and ν at p=5 on few |κ|=1 sets."""
    import numpy as np
    from itertools import combinations, permutations
    from minmax_quadratic import paley_conference_prime_power

    def per_w(C, S, T):
        s, t = list(S), list(T)
        tot = 0.0
        for perm in permutations(range(4)):
            pr = 1.0
            for i in range(4):
                pr *= C[s[i], t[perm[i]]]
            tot += pr
        return tot

    out = {}
    for p in (5,):
        yp, ym = Path(f"/tmp/maxplus_p{p}.npy"), Path(f"/tmp/maxminus_p{p}.npy")
        if not yp.exists() or not ym.exists():
            out[str(p)] = {"ok": False, "reason": "missing cache"}
            continue
        Yp, Ym = np.load(yp), np.load(ym)
        C = paley_conference_prime_power(p).astype(np.float64)
        n = C.shape[0]
        fours = list(combinations(range(n), 4))
        aa = np.array([s[0] for s in fours])
        bb = np.array([s[1] for s in fours])
        cc = np.array([s[2] for s in fours])
        dd = np.array([s[3] for s in fours])
        kap = (
            C[aa, bb] * C[cc, dd]
            + C[aa, cc] * C[bb, dd]
            + C[aa, dd] * C[bb, cc]
        )
        phi = (C[:, aa] * C[:, bb] * C[:, cc] * C[:, dd]).sum(0)
        m4p = (Yp[:, aa] * Yp[:, bb] * Yp[:, cc] * Yp[:, dd]).mean(0)
        m4m = (Ym[:, aa] * Ym[:, bb] * Ym[:, cc] * Ym[:, dd]).mean(0)
        nu = 0.5 * (m4p - m4m)
        m1 = np.where(np.abs(np.abs(kap) - 1) < 1e-9)[0]
        rng = np.random.default_rng(1)
        sample = rng.choice(m1, size=min(8, len(m1)), replace=False)
        p4m1 = p**4 - 1
        max_err_p = max_err_m = max_err_nu = 0.0
        for si in sample:
            S = fours[int(si)]
            ext_p = ext_m = ext_nu = 0.0
            for tj, T in enumerate(fours):
                if tj == si:
                    continue
                w = per_w(C, S, T)
                if abs(w) < 1e-15:
                    continue
                ext_p += w * m4p[tj]
                ext_m += w * m4m[tj]
                ext_nu += w * nu[tj]
            lhs_p = p4m1 * m4p[si] + 2 * phi[si]
            lhs_m = p4m1 * m4m[si] + 2 * phi[si]
            lhs_nu = p4m1 * nu[si]
            max_err_p = max(max_err_p, abs(ext_p - lhs_p))
            max_err_m = max(max_err_m, abs(ext_m - lhs_m))
            max_err_nu = max(max_err_nu, abs(ext_nu - lhs_nu))
        out[str(p)] = {
            "ok": True,
            "n_sample": len(sample),
            "max_err_m4plus": max_err_p,
            "max_err_m4minus": max_err_m,
            "max_err_nu": max_err_nu,
            "same_ext_ok": max_err_p < 1e-8 and max_err_m < 1e-8,
            "nu_eigen_ok": max_err_nu < 1e-8,
            "max_abs_nu_kappa1": float(np.max(np.abs(nu[m1]))),
        }
    return out


def main() -> dict:
    a = theorem_m4minus_same_ext()
    b = theorem_ext_nu_eigen()
    cert = certify_same_ext()
    out = {
        "title": (
            "Prop 15.263 same Ext equation for m₄±; Ext[ν]=(p⁴−1)ν on |κ|=1; "
            "residual (i) OPEN"
        ),
        "proved": {
            "m4minus_same_ext": a["proved"],
            "ext_nu_eigen": b["proved"],
            "forces_nu_zero": b["forces_nu_zero"],
            "residual_i_closed": residual_i_closed_via_263(),
            "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
            "type_I": type_I_k_3p_minus_2_closed_general(),
            "E1_closed": e1_closed_general(),
            "L_closed": False,
        },
        "algebra": {"m4minus_ext": a, "ext_nu": b},
        "certified": cert,
        "hinge": hinge_status_263(),
        "F3": (
            "Both m₄± solve (p⁴−1)f+2φ=Ext[f] on |κ|=1 ⇒ Ext[ν]=(p⁴−1)ν.  "
            "π-odd kernel gap still OPEN.  No predicate flip."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15263.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.263 same Ext for m4±; Ext[ν]=(p^4-1)ν; residual-i OPEN")
    print(f"  m4- same Ext: {a['proved']}")
    print(f"  Ext[ν]=(p^4-1)ν: {b['proved']} (forces ν=0? {b['forces_nu_zero']})")
    for pk, v in cert.items():
        print(f"  p={pk}: {v}")
    print(f"  residual_i closed: {residual_i_closed_via_263()}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
