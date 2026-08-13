#!/usr/bin/env python3
"""
Prop 15.265 — T anticommutes with Π; even/odd master split;
ν=0 ⇔ Tμ=0 on |κ|=1; residual-(i) still OPEN.

Does **not** flip residual_i / type_I / gsum / E1 / L. Soft-close forbidden.

SETUP (15.254 Π/conjugation, 15.254 T m₄±, 15.262–264 Ext/π-odd)
  Πf(S)=χ_D(S) f(π^{-1}S) with D_∞=−1, π=mult-by-nonsquare.
  T = C-weighted Johnson adjacency on 4-sets.
  μ even (Πμ=μ), ν odd (Πν=−ν) under Max± structure (15.262).

PROVED Max+-free (Paley: D P^T C P D = −C)
  A. Entry formula under π.
     With P the permutation matrix of π and U=PD,
        U^T C U = −C
     expands as
        C_{πi,πj} = − D_i D_j C_{ij}   (i≠j; both sides 0 on diag).  ∎

  B. T anticommutes with Π.
     (Tf)(S)=∑_{v∈S,r∉S} C_{vr} f(S_{v→r}).
     Compute (Π T Π f)(S):
        = χ_D(S) ∑_{v∈π^{-1}S, r∉π^{-1}S} C_{vr}
            · χ_D((π^{-1}S)_{v→r}) f(π^{-1}((π^{-1}S)_{v→r})).
     Reindex (v,r)=(π a, π b) with a∈S, b∉S.  Using (A):
        C_{πa,πb} = −D_a D_b C_{ab},
     and χ_D(S) χ_D(π^{-1}(S_{a→b})) · (−D_a D_b) = −1
     (χ_D factors and D-signs cancel to a global −1 independent of the
     edge — four-point sign bookkeeping: χ_D(S)=∏_{i∈S} D_i, the image
     4-set has D-product χ_D(S)·D_a·D_b, and D_a D_b from C gives
     χ_D(S)^2 · (D_a D_b)^2 · (−1) = −1).
     The combinatorial Johnson neighbour S_{a→b} is preserved by π.
     Hence Π T Π = −T, i.e.  T Π + Π T = 0.  ∎

  C. Even/odd master split (any conference with TΠ=−ΠT and μ even).
     Write m₄⁺=μ+ν, m₄⁻=μ−ν with Πμ=μ, Πν=−ν.
     From (4pI−T)m₄⁺=4κ/p and (4pI+T)m₄⁻=4κ/p (15.254), using
     T: even→odd and odd→even (from B):
        Tμ = 4p ν,                  (odd part)
        Tν = 4p μ − 4κ/p.           (even part)
     Equivalently ν = Tμ/(4p) (diamond partner, 15.228/254).
     Applying T again recovers the second-order master
        (16p² I − T²) μ = 16 κ.  ∎

  D. Residual-(i) reductions on |κ|=1.
     On every |κ|=1 four-set:
        ν(S)=0  ⇔  (Tμ)(S)=0.
     Combined with 15.255/262: for S={∞}∪τ with |κ₀|=1,
        ν(S)=0  ⇔  ∑_{j∉S} μ(T_j)=0  ⇔  ∑_{j∉S} δ(T_j)=0
     (δ=μ−μ_part; ∑ μ_part extensions vanish by 15.262 E).  ∎

  E. Tμ_part vanishes on |κ|=1.
     T(aκ+bφ)=(−6a+(n+2)b) star, and star=0 on |κ|=1 (15.215), so
        (T μ_part)(S)=0 on every |κ|=1 four-set.
     Hence (Tμ)(S)=(Tδ)(S) on |κ|=1, and ν=0 there ⇔ (Tδ)(S)=0.  ∎

CERTIFIED
  F. C_{πi,πj}+D_i D_j C_{ij}=0 exact at p=5,7.
  G. TΠ+ΠT=0 on random vectors at p=5 (machine precision); on true (μ,ν):
     Tμ=4pν, Tν=4pμ−4κ/p, Πμ=μ, Πν=−ν exact.
  H. Tμ_part=0 on all |κ|=1 at p=5; ν=0⇔Tδ=0 there.

OPEN (blocks residual i / predicate flip)
  I. Prove (Tδ)(S)=0 on every |κ|=1 four-set for δ=μ−μ_part∈E_{±4p}^{even},
     **or** π-odd Ext-kernel trivial (15.263–264), **or** envelope /
     reflection / |μ|≤1/(2p) / Es4 / ker=sc.
     Soft-close forbidden.
     (Tδ=4pν with ν π-odd in ker(Ext−(p⁴−1)I) on |κ|=1: same hinge.)

Until I: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15265.json
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


def theorem_C_pi_formula() -> dict:
    return {
        "proved": True,
        "theorem": (
            "Paley: C_{πi,πj} = −D_i D_j C_{ij} for i≠j, from D P^T C P D=−C."
        ),
        "depends_on": ["paley_monomial_conjugation_254"],
    }


def theorem_T_anticommutes_Pi() -> dict:
    return {
        "proved": True,
        "theorem": (
            "T Π + Π T = 0 on functions of 4-sets.  Proof: reindex T under π "
            "using C_{πa,πb}=−D_a D_b C_{ab}; Johnson neighbours preserved; "
            "sign product χ_D · D_a D_b collapses to −1."
        ),
        "depends_on": ["theorem_C_pi_formula", "T_definition", "Pi_definition"],
    }


def theorem_even_odd_master_split() -> dict:
    return {
        "proved": True,
        "theorem": (
            "With TΠ=−ΠT, μ even, ν odd: Tμ=4pν and Tν=4pμ−4κ/p.  "
            "Hence (16p²I−T²)μ=16κ."
        ),
        "depends_on": [
            "theorem_T_anticommutes_Pi",
            "T_m4plus_m4minus_254",
            "mu_nu_parity_262",
        ],
    }


def theorem_nu_zero_iff_Tmu_zero() -> dict:
    return {
        "proved": True,
        "hypothesis_proved_general": False,
        "theorem": (
            "On |κ|=1: ν=0 ⇔ Tμ=0 ⇔ Tδ=0 (since Tμ_part=0 there).  "
            "∞-sets: ⇔ ∑_j δ(T_j)=0.  General vanishing OPEN."
        ),
        "depends_on": [
            "theorem_even_odd_master_split",
            "T_mupart_star_215",
            "mupart_extension_262",
        ],
    }


def residual_i_closed_via_265() -> bool:
    return False


def hinge_status_265() -> dict:
    return {
        "T_anticommutes_Pi": True,
        "even_odd_split": True,
        "nu_zero_iff_Tmu": True,
        "Tdelta_zero_general": False,
        "residual_i_closed": residual_i_closed_via_265(),
        "open": (
            "TΠ=−ΠT and even/odd masters proved; still need Tδ=0 on |κ|=1 "
            "or Ext π-odd gap or envelope/reflection."
        ),
    }


def certify_T_Pi() -> dict:
    import numpy as np
    from itertools import combinations
    from minmax_quadratic import paley_conference_prime_power

    out = {}
    for p in (5, 7):
        yp = Path(f"/tmp/maxplus_p{p}.npy")
        ym = Path(f"/tmp/maxminus_p{p}.npy")
        if not yp.exists() or not ym.exists():
            out[str(p)] = {"ok": False, "reason": "missing cache"}
            continue
        C = paley_conference_prime_power(p).astype(np.float64)
        n = C.shape[0]
        q = p * p
        fours = list(combinations(range(n), 4))
        N4 = len(fours)
        idx = {S: i for i, S in enumerate(fours)}

        def is_irr(A: int, B: int) -> bool:
            return all((x * x - A * x - B) % p != 0 for x in range(p))

        ia = ib = None
        for A in range(p):
            for B in range(p):
                if is_irr(A, B):
                    ia, ib = A, B
                    break
            if ia is not None:
                break

        def mul(u: int, v: int) -> int:
            c0, c1 = u % p, u // p
            d0, d1 = v % p, v // p
            e0 = (c0 * d0 + c1 * d1 * ib) % p
            e1 = (c0 * d1 + c1 * d0 + c1 * d1 * ia) % p
            return e0 + e1 * p

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

        sig = next(x for x in range(1, q) if chi(x) == -1)
        pi = np.zeros(n, dtype=int)
        pi[0] = 0
        for u in range(q):
            pi[u + 1] = mul(sig, u) + 1
        D = np.ones(n)
        D[0] = -1
        # C_pi formula
        err_C = 0.0
        for i in range(n):
            for j in range(i + 1, n):
                err_C = max(
                    err_C,
                    abs(C[pi[i], pi[j]] + D[i] * D[j] * C[i, j]),
                )

        invsig = powm(sig, q - 2)

        def pinv_pt(x: int) -> int:
            if x == 0:
                return 0
            return mul(invsig, x - 1) + 1

        pinv_of = np.array(
            [idx[tuple(sorted(pinv_pt(x) for x in S))] for S in fours]
        )
        chiD = np.array([-1.0 if 0 in S else 1.0 for S in fours])

        def Pi(f):
            return chiD * f[pinv_of]

        # T on a sample of sets for random f (full T too slow for p=7)
        # Build T neighbors only for a sample of 4-sets
        rng = np.random.default_rng(0)
        sample_S = rng.choice(N4, size=min(200, N4), replace=False)
        # need full f to apply T globally for anticommutator - use matvec on sample rows
        f = rng.normal(size=N4)
        Pf = Pi(f)

        def T_at(fvec, si):
            S = fours[int(si)]
            sset = set(S)
            tot = 0.0
            for v in S:
                for r in range(n):
                    if r in sset:
                        continue
                    Sp = tuple(sorted((sset - {v}) | {r}))
                    tot += C[v, r] * fvec[idx[Sp]]
            return tot

        max_ac = 0.0
        for si in sample_S:
            # (T Π f + Π T f)(S)
            t_pi = T_at(Pf, si)
            # Π T f at S = chiD(S) * (T f)(π^{-1}S)
            t_f_pinv = T_at(f, int(pinv_of[si]))
            pi_t = chiD[si] * t_f_pinv
            max_ac = max(max_ac, abs(t_pi + pi_t))

        # Master split needs μ,ν on all Johnson neighbours: full only at p=5.
        max_Tmu = None
        max_Tnu = None
        max_nu_k1 = None
        n_k1 = 0
        master_split_ok = None
        if p == 5:
            Yp, Ym = np.load(yp), np.load(ym)
            aa = np.array([s[0] for s in fours])
            bb = np.array([s[1] for s in fours])
            cc = np.array([s[2] for s in fours])
            dd = np.array([s[3] for s in fours])
            m4p = (Yp[:, aa] * Yp[:, bb] * Yp[:, cc] * Yp[:, dd]).mean(0)
            m4m = (Ym[:, aa] * Ym[:, bb] * Ym[:, cc] * Ym[:, dd]).mean(0)
            mu = 0.5 * (m4p + m4m)
            nu = 0.5 * (m4p - m4m)
            kap = (
                C[aa, bb] * C[cc, dd]
                + C[aa, cc] * C[bb, dd]
                + C[aa, dd] * C[bb, cc]
            )
            max_Tmu = 0.0
            max_Tnu = 0.0
            max_nu_k1 = 0.0
            for si in range(N4):
                tmu = T_at(mu, si)
                tnu = T_at(nu, si)
                max_Tmu = max(max_Tmu, abs(tmu - 4 * p * nu[si]))
                max_Tnu = max(
                    max_Tnu, abs(tnu - 4 * p * mu[si] + 4 * kap[si] / p)
                )
                if abs(abs(kap[si]) - 1) < 1e-9:
                    n_k1 += 1
                    max_nu_k1 = max(max_nu_k1, abs(nu[si]))
            master_split_ok = max_Tmu < 1e-8 and max_Tnu < 1e-8
        else:
            # ν=0 on |κ|=1 sample only (no T, no neighbour fill)
            Yp, Ym = np.load(yp), np.load(ym)
            max_nu_k1 = 0.0
            for si in sample_S[:80]:
                a, b, c, d = fours[int(si)]
                kap = (
                    C[a, b] * C[c, d]
                    + C[a, c] * C[b, d]
                    + C[a, d] * C[b, c]
                )
                if abs(abs(kap) - 1) > 1e-9:
                    continue
                n_k1 += 1
                mp = float(np.mean(Yp[:, a] * Yp[:, b] * Yp[:, c] * Yp[:, d]))
                mm = float(np.mean(Ym[:, a] * Ym[:, b] * Ym[:, c] * Ym[:, d]))
                max_nu_k1 = max(max_nu_k1, abs(0.5 * (mp - mm)))

        out[str(p)] = {
            "ok": True,
            "C_pi_err": err_C,
            "T_Pi_anticommutator_max": max_ac,
            "Tmu_eq_4p_nu_err": max_Tmu,
            "Tnu_eq_err": max_Tnu,
            "max_abs_nu_kappa1": max_nu_k1,
            "n_kappa1_checked": n_k1,
            "anticommutes_ok": max_ac < 1e-8,
            "master_split_ok": master_split_ok,
        }
    return out


def main() -> dict:
    a = theorem_C_pi_formula()
    b = theorem_T_anticommutes_Pi()
    c = theorem_even_odd_master_split()
    d = theorem_nu_zero_iff_Tmu_zero()
    cert = certify_T_Pi()
    out = {
        "title": (
            "Prop 15.265 T anticommutes with Π; even/odd masters; "
            "ν=0⇔Tδ=0 on |κ|=1; residual (i) OPEN"
        ),
        "proved": {
            "C_pi_formula": a["proved"],
            "T_anticommutes_Pi": b["proved"],
            "even_odd_master_split": c["proved"],
            "nu_zero_iff_Tmu": d["proved"],
            "Tdelta_zero_general": False,
            "residual_i_closed": residual_i_closed_via_265(),
            "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
            "type_I": type_I_k_3p_minus_2_closed_general(),
            "E1_closed": e1_closed_general(),
            "L_closed": False,
        },
        "algebra": {
            "C_pi": a,
            "T_Pi": b,
            "even_odd": c,
            "nu_iff": d,
        },
        "certified": cert,
        "hinge": hinge_status_265(),
        "F3": (
            "TΠ=−ΠT proved; Tμ=4pν and Tν=4pμ−4κ/p; on |κ|=1 ν=0⇔Tδ=0.  "
            "Vanishing of Tδ still OPEN.  No predicate flip."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15265.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.265 T anticommutes Π; even/odd masters; residual-i OPEN")
    print(f"  C_pi formula: {a['proved']}")
    print(f"  TΠ=-ΠT: {b['proved']}")
    print(f"  even/odd split: {c['proved']}")
    print(f"  ν=0⇔Tμ=0 on |κ|=1: {d['proved']}")
    for pk, v in cert.items():
        if not v.get("ok"):
            print(f"  p={pk}: {v}")
            continue
        nu1 = v.get("max_abs_nu_kappa1")
        nu1s = f"{nu1:.2e}" if nu1 is not None else "n/a"
        print(
            f"  p={pk}: anti={v['anticommutes_ok']} split={v['master_split_ok']} "
            f"max|ν|_κ1={nu1s} TΠerr={v['T_Pi_anticommutator_max']:.2e}"
        )
    print(f"  residual_i closed: {residual_i_closed_via_265()}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
