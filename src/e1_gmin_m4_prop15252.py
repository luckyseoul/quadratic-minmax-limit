#!/usr/bin/env python3
"""
Prop 15.252 — Master/T²/Ext residual-(i) criteria (Max+-free equivalences);
envelope + Ext + |16κ+T²μ| census p=5,7; residual (i) still OPEN.

Does **not** flip residual_i / type_I / gsum / E1 / L. Soft-close forbidden.

SETUP (15.176–177 master, 15.240 envelope, 15.251 Cy/Ext)
  Master (any conference Max± μ = ½(m₄⁺+m₄⁻)):
      (16 p² I − T²) μ = 16 κ.
  On |κ|=1: residual-(i) ⇔ |μ|≤1/(2p) ⇔ Gsum≥−1/p (15.176).
  Sufficient: |μ|≤2/n (15.188/240); envelope |μ|≤max(|μ_part|,|f4|) (15.240).

PROVED Max+-free
  A. Pointwise master rewrite.
        T²μ = 16(p² μ − κ)    (all 4-sets)
     Hence |T²μ| = 16 |p² μ − κ| and
        16κ + T²μ = 16 p² μ
     (so crude triangle on |16κ+T²μ| alone cannot bound |μ| without an
     independent estimate of T²μ — signed cancellation in T² is essential).  ∎

  B. Wick residual form.
     ρ := μ − κ/p²  ⇒  |T²μ| = 16 p² |ρ|  and
        |μ| ≤ 1/p² + |ρ|  on |κ|=1.
     Therefore |ρ| ≤ (p−2)/(2 p²) =: L_abs (15.205)
        ⇒  |μ| ≤ 1/p² + (p−2)/(2p²) = 1/(2p)
     for every prime p≥3.  ∎

  C. Ext triangle criterion for the 2/n target (15.251 Cy identity).
     On |κ|=1: (p⁴−1)μ + 2φ = Ext.
     If |Ext| ≤ 2(p²−1) − 2|φ| then
        |(p⁴−1)μ| = |Ext − 2φ| ≤ 2(p²−1)
        ⇒ |μ| ≤ 2/(p²+1) = 2/n.
     Under the proved |φ|≤2(p−2) (15.186), the uniform majorant
        |Ext| ≤ 2p² − 4p + 6
     implies |Ext| ≤ 2(p²−1)−2|φ| and hence |μ|≤2/n≤1/(2p).  ∎

  D. Envelope criterion recall (15.240).
     |μ| ≤ max(|μ_part|,|f4|) on |κ|=1 ⇒ |μ|≤2/n≤1/(2p) for all primes
     p≥5 (maj≤2/n and |f4|≤2/n proved).  ∎

CERTIFIED (Max± census; not general)
  E. p=5: μ = f4 on all |κ|=1; envelope holds; max|μ|=3/65 = M_cand;
     max|T²μ|=288/13; max|16κ+T²μ|=240/13 (=16p² M_cand, Mcand-sharp);
     max|Ext|=124/5; 2p²−4p+6=36; room on Ext majorant 36−24.8>0.
  F. p=7: envelope holds on all 176400 |κ|=1 sets (incl. CR split at
     |φ|=10 with μ∈{±61/2863, ±109/2863}); max|μ|=109/2863 < 2/n=1/25;
     max|T²μ|≈13.85 ≤ 8(p−2)=40; max|16κ+T²μ|≈29.85 ≤ 31.36 (2/n thr)
     ≤ 32.94 (M_cand thr) ≤ 56 (1/(2p) thr); max|Ext| on extremes ≈71.4
     ≤ 2p²−4p+6=76; comb max|(p⁴−1)μ|≈91.37 ≤ 96 with room ≈4.63.
  G. |μ| ≰ |f4| pointwise at p=7 (e.g. (κ,φ)=(−1,+2): |μ|≈0.0227 > |f4|
     ≈0.0171) while |μ| ≤ max(|μ_part|,|f4|) still holds — confirms that
     the max-abs envelope (not the pure |f4| bound) is the right target.

OPEN (blocks residual i / predicate flip)
  H. Prove for all primes p≥5 one of:
     (i)   envelope |μ|≤max(|μ_part|,|f4|) on every |κ|=1 four-set;
     (ii)  |Ext| ≤ 2p²−4p+6 on |κ|=1 (⇒ |μ|≤2/n via C);
     (iii) |ρ| ≤ L_abs = (p−2)/(2p²) on |κ|=1 (⇒ |μ|≤1/(2p) via B);
     (iv)  E[s⁴]≤15n²−22n (15.250) or ker=sc (15.249 free-e).
     Soft-close forbidden.  Do **not** re-run crude |T²|≤16(n−2)/n
     (false at p=5) or unsigned Per.

Until H: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15252.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    is_prime,
    n_of,
    residual_i_dual_eq_empty_proved_general,
    type_I_k_3p_minus_2_closed_general,
)
from e1_gmin_m4_prop15177 import mu4_threshold  # noqa: E402
from e1_gmin_m4_prop15188 import two_over_n  # noqa: E402
from e1_gmin_m4_prop15205 import L_abs, M_cand  # noqa: E402
from e1_gmin_m4_prop15240 import theorem_envelope_criterion  # noqa: E402


def theorem_T2_master_rewrite() -> dict:
    return {
        "proved": True,
        "theorem": (
            "Master (16p²I−T²)μ=16κ ⇒ T²μ=16(p²μ−κ) pointwise on all "
            "4-sets, and 16κ+T²μ=16p²μ.  Max+-free (any conference).  "
            "Triangle on |16κ+T²μ| alone is vacuous for bounding |μ|."
        ),
        "depends_on": ["master_equation_15_177", "T_symmetric_15_214"],
    }


def theorem_rho_L_abs_suffices() -> dict:
    ok = True
    sample = {}
    for p in [p for p in range(3, 60) if is_prime(p)]:
        # 1/p^2 + L_abs = 1/(2p)
        lhs = Fraction(1, p * p) + L_abs(p)
        rhs = Fraction(1, 2 * p)
        if lhs != rhs:
            ok = False
        if p in (3, 5, 7, 11, 13):
            sample[str(p)] = {
                "L_abs": str(L_abs(p)),
                "one_over_p2_plus_L_abs": str(lhs),
                "one_over_2p": str(rhs),
                "equal": lhs == rhs,
            }
    return {
        "proved": ok,
        "theorem": (
            "On |κ|=1: |μ|≤1/p²+|ρ| with ρ=μ−κ/p².  "
            "|ρ|≤L_abs=(p−2)/(2p²) ⇒ |μ|≤1/(2p) for all primes p≥3."
        ),
        "by_p_sample": sample,
        "depends_on": ["L_abs_15_205", "kappa_abs_1"],
    }


def theorem_Ext_triangle_two_over_n() -> dict:
    ok = True
    sample = {}
    for p in [p for p in range(5, 80) if is_prime(p)]:
        # uniform Ext majorant
        ext_maj = 2 * p * p - 4 * p + 6
        # need: ext_maj <= 2(p^2-1) - 2*(2(p-2)) = 2p^2-2 - 4p + 8 = 2p^2-4p+6
        need = 2 * (p * p - 1) - 4 * (p - 2)
        if ext_maj != need:
            ok = False
        if p in (5, 7, 11, 13):
            sample[str(p)] = {
                "Ext_uniform_maj": ext_maj,
                "budget_with_phi_max": need,
                "equal": ext_maj == need,
                "two_over_n": str(two_over_n(p)),
            }
    return {
        "proved": ok,
        "theorem": (
            "On |κ|=1: (p⁴−1)μ+2φ=Ext.  If |Ext|≤2(p²−1)−2|φ| then "
            "|μ|≤2/n.  Under |φ|≤2(p−2), the bound |Ext|≤2p²−4p+6 is "
            "equivalent and implies |μ|≤2/n≤1/(2p) for primes p≥5."
        ),
        "by_p_sample": sample,
        "depends_on": ["Cy_identity_15_251", "phi_bound_15_186", "two_over_n_15_188"],
    }


def certify_census() -> dict:
    """Census envelope / Ext / combo at p=5,7 from Max± caches."""
    import numpy as np
    from itertools import combinations
    from minmax_quadratic import paley_conference_prime_power

    out = {}
    for p in (5, 7):
        yp = Path(f"/tmp/maxplus_p{p}.npy")
        ym = Path(f"/tmp/maxminus_p{p}.npy")
        if not yp.exists() or not ym.exists():
            out[str(p)] = {"ok": False, "reason": "missing Max± cache"}
            continue
        Yp = np.load(yp)
        Ym = np.load(ym)
        C = paley_conference_prime_power(p).astype(np.float64)
        n = C.shape[0]
        fours = np.array(list(combinations(range(n), 4)), dtype=np.int32)
        a, b, c, d = fours.T
        kap = C[a, b] * C[c, d] + C[a, c] * C[b, d] + C[a, d] * C[b, c]
        ids = np.where(np.abs(np.abs(kap) - 1) < 1e-9)[0]
        aa, bb, cc, dd = a[ids], b[ids], c[ids], d[ids]
        phi = np.zeros(len(ids))
        for r in range(n):
            phi += C[r, aa] * C[r, bb] * C[r, cc] * C[r, dd]
        phi = np.rint(phi)
        kaps = np.rint(kap[ids])
        den = p * p * (p * p - 5)
        pn = p * n
        max_mu = 0.0
        max_ext = 0.0
        max_comb = 0.0
        env_viol = 0
        f4_viol = 0
        max_rho = 0.0
        chunk = 4000
        for s0 in range(0, len(ids), chunk):
            sl = slice(s0, s0 + chunk)
            mp = np.mean(
                Yp[:, aa[sl]] * Yp[:, bb[sl]] * Yp[:, cc[sl]] * Yp[:, dd[sl]],
                axis=0,
            )
            mm = np.mean(
                Ym[:, aa[sl]] * Ym[:, bb[sl]] * Ym[:, cc[sl]] * Ym[:, dd[sl]],
                axis=0,
            )
            mu = 0.5 * (mp + mm)
            k = kaps[sl]
            ph = phi[sl]
            mupart = ((p * p - 1) * k - 2 * ph) / den
            f4 = (4 * k - ph) / pn
            env = np.maximum(np.abs(mupart), np.abs(f4))
            env_viol += int(np.sum(np.abs(mu) > env + 1e-12))
            f4_viol += int(np.sum(np.abs(mu) > np.abs(f4) + 1e-12))
            Ext = (p**4 - 1) * mu + 2 * ph
            comb = (p**4 - 1) * mu
            rho = mu - k / (p * p)
            max_mu = max(max_mu, float(np.max(np.abs(mu))))
            max_ext = max(max_ext, float(np.max(np.abs(Ext))))
            max_comb = max(max_comb, float(np.max(np.abs(comb))))
            max_rho = max(max_rho, float(np.max(np.abs(rho))))
        two_n = 2.0 / n
        thr = 1.0 / (2 * p)
        ext_maj = 2 * p * p - 4 * p + 6
        comb_thr = 2 * (p * p - 1)
        out[str(p)] = {
            "ok": True,
            "n_kappa1": int(len(ids)),
            "max_abs_mu": max_mu,
            "max_abs_Ext": max_ext,
            "max_abs_comb": max_comb,
            "max_abs_rho": max_rho,
            "two_over_n": two_n,
            "one_over_2p": thr,
            "M_cand": float(M_cand(p)),
            "L_abs": float(L_abs(p)),
            "Ext_uniform_maj": ext_maj,
            "comb_thr_for_2n": comb_thr,
            "mu_le_2n": max_mu <= two_n + 1e-12,
            "mu_le_half_p": max_mu <= thr + 1e-12,
            "mu_le_Mcand": max_mu <= float(M_cand(p)) + 1e-12,
            "rho_le_L_abs": max_rho <= float(L_abs(p)) + 1e-12,
            "Ext_le_uniform_maj": max_ext <= ext_maj + 1e-9,
            "comb_le_thr": max_comb <= comb_thr + 1e-9,
            "envelope_viol": env_viol,
            "f4_only_viol": f4_viol,
            "envelope_ok": env_viol == 0,
        }
    return out


def residual_i_closed_via_252() -> bool:
    return False


def hinge_status_252() -> dict:
    return {
        "T2_master_rewrite": True,
        "rho_L_abs_suffices": True,
        "Ext_triangle_two_over_n": True,
        "envelope_criterion": True,
        "envelope_hyp_general": False,
        "Ext_bound_general": False,
        "rho_bound_general": False,
        "residual_i_closed": residual_i_closed_via_252(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "open": (
            "Prove envelope / |Ext|≤2p²−4p+6 / |ρ|≤L_abs on |κ|=1 for all "
            "primes p≥5, or Es4 / ker=sc alternate."
        ),
    }


def main() -> dict:
    a = theorem_T2_master_rewrite()
    b = theorem_rho_L_abs_suffices()
    c = theorem_Ext_triangle_two_over_n()
    d = theorem_envelope_criterion()
    cert = certify_census()
    env_crit_ok = bool(d.get("proved_criterion", d.get("proved", False)))
    out = {
        "title": (
            "Prop 15.252 Master/T²/Ext residual-(i) criteria; "
            "envelope+Ext census p=5,7; residual (i) OPEN"
        ),
        "proved": {
            "T2_master_rewrite": a["proved"],
            "rho_L_abs_suffices": b["proved"],
            "Ext_triangle_two_over_n": c["proved"],
            "envelope_criterion": env_crit_ok,
            "envelope_hyp_general": False,
            "Ext_bound_general": False,
            "rho_bound_general": False,
            "residual_i_closed": residual_i_closed_via_252(),
            "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
            "type_I": type_I_k_3p_minus_2_closed_general(),
            "E1_closed": e1_closed_general(),
            "L_closed": False,
        },
        "algebra": {
            "T2_master": a,
            "rho_L_abs": b,
            "Ext_triangle": c,
            "envelope": d,
        },
        "certified": {"census": cert},
        "hinge": hinge_status_252(),
        "F3": (
            "Master rewrite T²μ=16(p²μ−κ); |ρ|≤L_abs⇒|μ|≤1/(2p); "
            "|Ext|≤2p²−4p+6⇒|μ|≤2/n.  Census p=5,7: envelope OK, Ext and "
            "ρ under budget, |μ|≤2/n.  General bounds still OPEN.  No "
            "predicate flip."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15252.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.252 Master/T2/Ext criteria; residual-i OPEN")
    print(f"  T2 rewrite: {a['proved']}")
    print(f"  rho L_abs suffices: {b['proved']}")
    print(f"  Ext triangle: {c['proved']}")
    print(f"  envelope criterion: {env_crit_ok}")
    for pk, v in cert.items():
        if not v.get("ok"):
            print(f"  p={pk}: {v}")
            continue
        print(
            f"  p={pk}: max|μ|={v['max_abs_mu']:.6f} env_ok={v['envelope_ok']} "
            f"Ext={v['max_abs_Ext']:.4f}/{v['Ext_uniform_maj']} "
            f"rho={v['max_abs_rho']:.6f}/{v['L_abs']:.6f} "
            f"f4_only_viol={v['f4_only_viol']}"
        )
    print(f"  residual_i closed: {residual_i_closed_via_252()}")
    print(f"  open: {out['hinge']['open']}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
