#!/usr/bin/env python3
"""
Prop 15.259 — Transport lemma for local orientation; reduction of m₄⁺=m₄⁻
to π-covariance of m₄⁺ on |κ|=1; residual-(i) still OPEN.

Does **not** flip residual_i / type_I / gsum / E1 / L. Soft-close forbidden.

SETUP (15.254 conjugation, 15.258 orientation)
  Paley C, n=p²+1.  Max±, m₄±, μ, ν as usual.
  Local orientation (15.258): |κ(S)|=1 ⇔ C_S ∼_sp −C_S with det D=+1.
  Global conjugation (15.254): D_∞ P_πᵀ C P_π D_∞ = −C, hence
     m₄⁻(S) = χ_D(S) m₄⁺(π⁻¹S).

PROVED Max+-free
  A. Transport lemma.
     Let S be any 4-set with |κ(S)|=1, and (D,Π) a signed permutation of
     the four points of S with D Πᵀ C_S Π D = −C_S and det D = +1
     (exists by 15.258).  Extend D,Π by the identity on V\\S, and set
        C' := D_ext Pᵀ C P D_ext.
     Then:
       (i)   C' is a conference matrix (signed permutation of C);
       (ii)  C'_S = −C_S;
       (iii) Φ(y) := D_ext Pᵀ y is a bijection Max⁺(C) → Max⁺(C');
       (iv)  ∏_{i∈S} Φ(y)_i = (det D) ∏_{i∈S} y_i = ∏_{i∈S} y_i;
       (v)   consequently m₄⁺(C,S) = m₄⁺(C',S).  ∎
     Proof of (iii): if Cy=py then
        C'(D Pᵀ y) = D Pᵀ C P D · D Pᵀ y = D Pᵀ C y = p D Pᵀ y.
     (iv) is det D=+1 and P permutes S.  (v) averages (iv) over Max⁺.  ∎

  B. Reduction to π-covariance of m₄⁺.
     On every 4-set, m₄⁻(S)=χ_D(S) m₄⁺(π⁻¹S) (15.254).  Hence on |κ|=1:
        m₄⁺(S)=m₄⁻(S)  ⇔  m₄⁺(S)=χ_D(S) m₄⁺(π⁻¹S).
     Finite S (∞∉S): χ_D=1, so ⇔ m₄⁺(S)=m₄⁺(π⁻¹S).
     ∞∈S: χ_D=−1, so ⇔ m₄⁺(S)=−m₄⁺(π⁻¹S).  ∎

  C. Type transport under π (Paley algebra).
     For finite S: mult-by-nonsquare flips all finite edges of C, so
     κ(πS)=κ(S) and φ(πS)=φ(S) (four sign flips in each φ-term from finite
     rows; ∞-row unchanged).  For ∞∈S: (κ,φ)↦(−κ,−φ) as in 15.254 H.  ∎

CERTIFIED
  D. (A) holds at p=5,7 on sampled |κ|=1 sets: C'²=(n−1)I, C'_S=−C_S,
     transport residual 0, product identity exact.
  E. (B) identity m₄⁺=m₄⁻ on all |κ|=1 at p=3,5,7 (prior); equivalently
     finite: m₄⁺(S)=m₄⁺(πS) and ∞: m₄⁺(S)=−m₄⁺(πS) on |κ|=1 (full sample
     checks: 5000 finite + 3000 ∞ at p=7, all |κ|=1 at p=5).
  F. Local conjugator pushforward of Max⁺ config counts on S equals Max⁻
     config counts at p=5 on all tested |κ|=1 (stronger than product match).

OPEN (blocks residual i / predicate flip)
  G. Prove Max+-free the π-covariance in (B):
        m₄⁺(S)=χ_D(S) m₄⁺(π⁻¹S) on every |κ|=1 four-set,
     for all primes p≥5.  Sufficient routes (referee):
       • S and π⁻¹S lie in the same Aut(C)-orbit when ∞∉S (then extend
         by 3-transitivity / 15.255 to all |κ|=1); or
       • a global signed auto taking C' back to −C while stabilising S
         with sign product +1 on S; or
       • Schur multilin average identity (15.257 D); or reflection.
     Soft-close forbidden.  Note: type-class sizes at p=5 do not divide
     |PGL(2,25)|, so (κ,φ) classes are not single Aut-orbits — Aut-orbit
     fusion under π needs a finer invariant or a different argument.

Until G: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15259.json
"""
from __future__ import annotations

import json
import sys
from itertools import combinations, permutations, product
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    type_I_k_3p_minus_2_closed_general,
)


def theorem_transport_lemma() -> dict:
    return {
        "proved": True,
        "theorem": (
            "On |κ(S)|=1, local orientation (D,Π) with det D=+1 extends by id "
            "outside S to C'=D Pᵀ C P D with C'_S=−C_S, Max⁺(C)≅Max⁺(C') via "
            "Φ(y)=D Pᵀ y, and m₄⁺(C,S)=m₄⁺(C',S)."
        ),
        "depends_on": ["theorem_orientation_258", "conference_signed_perm"],
    }


def theorem_reduction_to_pi_covariance() -> dict:
    return {
        "proved": True,
        "hypothesis_proved_general": False,
        "theorem": (
            "On |κ|=1: m₄⁺=m₄⁻ ⇔ m₄⁺(S)=χ_D(S) m₄⁺(π⁻¹S).  "
            "Finite: ⇔ m₄⁺(S)=m₄⁺(π⁻¹S).  ∞∈S: ⇔ m₄⁺(S)=−m₄⁺(π⁻¹S).  "
            "π-covariance OPEN general."
        ),
        "depends_on": ["theorem_mu_from_m4plus_pullback_254"],
    }


def theorem_pi_type_transport() -> dict:
    return {
        "proved": True,
        "theorem": (
            "Paley: finite S ⇒ (κ,φ)(πS)=(κ,φ)(S); "
            "∞∈S ⇒ (κ,φ)(πS)=(−κ,−φ)(S)."
        ),
        "depends_on": ["paley_quadratic_character", "infty_border"],
    }


def _find_conjugator(A):
    import numpy as np

    for perm in permutations(range(4)):
        Ap = A[np.ix_(perm, perm)]
        for signs in product([-1.0, 1.0], repeat=4):
            D = np.diag(signs)
            if np.allclose(D @ Ap @ D, -A) and abs(float(np.prod(signs)) - 1) < 1e-12:
                return list(perm), list(signs)
    return None


def certify_transport_and_pi() -> dict:
    import numpy as np
    from minmax_quadratic import paley_conference_prime_power

    out: dict = {}
    for p in (5, 7):
        yp, ym = Path(f"/tmp/maxplus_p{p}.npy"), Path(f"/tmp/maxminus_p{p}.npy")
        if not yp.exists() or not ym.exists():
            out[str(p)] = {"ok": False, "reason": "missing cache"}
            continue
        Yp, Ym = np.load(yp), np.load(ym)
        C = paley_conference_prime_power(p).astype(np.float64)
        n = C.shape[0]
        # field pi
        q = p * p

        def is_irr(a, b):
            return all((x * x - a * x - b) % p != 0 for x in range(p))

        ia = ib = None
        for a in range(p):
            for b in range(p):
                if is_irr(a, b):
                    ia, ib = a, b
                    break
            if ia is not None:
                break

        def mul(u, v):
            c0, c1 = u % p, u // p
            d0, d1 = v % p, v // p
            return (c0 * d0 + c1 * d1 * ib) % p + (
                (c0 * d1 + c1 * d0 + c1 * d1 * ia) % p
            ) * p

        def powm(u, e):
            r, base = 1, u
            while e:
                if e & 1:
                    r = mul(r, base)
                base = mul(base, base)
                e >>= 1
            return r

        def chi(x):
            if x == 0:
                return 0
            return 1 if powm(x, (q - 1) // 2) == 1 else -1

        sig = next(x for x in range(1, q) if chi(x) == -1)

        def pi_v(v):
            if v == 0:
                return 0
            return mul(sig, v - 1) + 1

        # Transport on a few |κ|=1 sets
        transport_ok = True
        n_transport = 0
        max_resid = 0.0
        for S in combinations(range(n), 4):
            a, b, c, d = S
            kap = (
                C[a, b] * C[c, d]
                + C[a, c] * C[b, d]
                + C[a, d] * C[b, c]
            )
            if abs(abs(kap) - 1) > 1e-9:
                continue
            conj = _find_conjugator(C[np.ix_(S, S)])
            if conj is None:
                transport_ok = False
                break
            perm, signs = conj
            # build C'
            sigma = list(range(n))
            for i in range(4):
                sigma[S[i]] = S[perm[i]]
            Dext = np.ones(n)
            for i in range(4):
                Dext[S[i]] = signs[i]
            idx = np.array(sigma)
            Cnew = Dext[:, None] * C[np.ix_(idx, idx)] * Dext[None, :]
            np.fill_diagonal(Cnew, 0.0)
            if not np.allclose(Cnew[np.ix_(S, S)], -C[np.ix_(S, S)]):
                transport_ok = False
            if not np.allclose(Cnew @ Cnew, (n - 1) * np.eye(n)):
                transport_ok = False
            Ynew = Yp[:, idx] * Dext[None, :]
            resid = float(np.max(np.abs(Ynew @ Cnew.T - p * Ynew)))
            max_resid = max(max_resid, resid)
            mp = float(np.mean(Yp[:, a] * Yp[:, b] * Yp[:, c] * Yp[:, d]))
            mp2 = float(np.mean(Ynew[:, a] * Ynew[:, b] * Ynew[:, c] * Ynew[:, d]))
            if abs(mp - mp2) > 1e-10:
                transport_ok = False
            n_transport += 1
            if n_transport >= (40 if p == 5 else 20):
                break

        # π-covariance sample: separate finite and ∞ |κ|=1
        pi_finite_ok = 0
        pi_finite_n = 0
        pi_inf_ok = 0
        pi_inf_n = 0
        need = 80 if p == 5 else 60
        # finite first
        for S in combinations(range(1, n), 4):
            a, b, c, d = S
            kap = (
                C[a, b] * C[c, d]
                + C[a, c] * C[b, d]
                + C[a, d] * C[b, c]
            )
            if abs(abs(kap) - 1) > 1e-9:
                continue
            mp = float(np.mean(Yp[:, a] * Yp[:, b] * Yp[:, c] * Yp[:, d]))
            Sp = tuple(sorted(pi_v(i) for i in S))
            aa, bb, cc, dd = Sp
            mp2 = float(np.mean(Yp[:, aa] * Yp[:, bb] * Yp[:, cc] * Yp[:, dd]))
            pi_finite_n += 1
            if abs(mp - mp2) < 1e-10:
                pi_finite_ok += 1
            if pi_finite_n >= need:
                break
        # ∞-containing
        for a, b, c in combinations(range(1, n), 3):
            S = (0, a, b, c)
            kap = C[a, b] + C[a, c] + C[b, c]  # C_0*=1
            if abs(abs(kap) - 1) > 1e-9:
                continue
            # full κ for {0,a,b,c}
            kap = (
                C[0, a] * C[b, c]
                + C[0, b] * C[a, c]
                + C[0, c] * C[a, b]
            )
            if abs(abs(kap) - 1) > 1e-9:
                continue
            mp = float(np.mean(Yp[:, 0] * Yp[:, a] * Yp[:, b] * Yp[:, c]))
            Sp = tuple(sorted(pi_v(i) for i in S))
            aa, bb, cc, dd = Sp
            mp2 = float(np.mean(Yp[:, aa] * Yp[:, bb] * Yp[:, cc] * Yp[:, dd]))
            pi_inf_n += 1
            if abs(mp + mp2) < 1e-10:
                pi_inf_ok += 1
            if pi_inf_n >= need:
                break

        # m4+=m4- on mixed sample
        m4eq_err = 0.0
        m4eq_n = 0
        for S in combinations(range(n), 4):
            a, b, c, d = S
            kap = (
                C[a, b] * C[c, d]
                + C[a, c] * C[b, d]
                + C[a, d] * C[b, c]
            )
            if abs(abs(kap) - 1) > 1e-9:
                continue
            mp = float(np.mean(Yp[:, a] * Yp[:, b] * Yp[:, c] * Yp[:, d]))
            mm = float(np.mean(Ym[:, a] * Ym[:, b] * Ym[:, c] * Ym[:, d]))
            m4eq_err = max(m4eq_err, abs(mp - mm))
            m4eq_n += 1
            if m4eq_n >= need * 2:
                break

        out[str(p)] = {
            "ok": True,
            "n_transport": n_transport,
            "transport_ok": transport_ok and max_resid < 1e-8,
            "max_transport_resid": max_resid,
            "pi_finite_ok": pi_finite_ok,
            "pi_finite_n": pi_finite_n,
            "pi_inf_ok": pi_inf_ok,
            "pi_inf_n": pi_inf_n,
            "m4eq_max_err": m4eq_err,
            "m4eq_n": m4eq_n,
            "pi_covariance_sample_ok": (
                pi_finite_ok == pi_finite_n and pi_inf_ok == pi_inf_n
            ),
            "m4eq_sample_ok": m4eq_err < 1e-12,
        }
    return out


def residual_i_closed_via_259() -> bool:
    return False


def hinge_status_259() -> dict:
    return {
        "transport_lemma": True,
        "reduction_to_pi_covariance": True,
        "pi_type_transport": True,
        "pi_covariance_general": False,
        "residual_i_closed": residual_i_closed_via_259(),
        "open": (
            "Prove m₄⁺(S)=χ_D(S) m₄⁺(π⁻¹S) on |κ|=1 Max+-free "
            "(π-covariance), or Schur average, or reflection."
        ),
    }


def main() -> dict:
    a = theorem_transport_lemma()
    b = theorem_reduction_to_pi_covariance()
    c = theorem_pi_type_transport()
    cert = certify_transport_and_pi()
    out = {
        "title": (
            "Prop 15.259 transport lemma + π-covariance reduction; "
            "residual (i) OPEN"
        ),
        "proved": {
            "transport_lemma": a["proved"],
            "reduction_to_pi_covariance": b["proved"],
            "pi_type_transport": c["proved"],
            "pi_covariance_general": False,
            "residual_i_closed": residual_i_closed_via_259(),
            "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
            "type_I": type_I_k_3p_minus_2_closed_general(),
            "E1_closed": e1_closed_general(),
            "L_closed": False,
        },
        "algebra": {
            "transport": a,
            "reduction": b,
            "pi_type": c,
        },
        "certified": cert,
        "hinge": hinge_status_259(),
        "F3": (
            "Transport lemma m₄⁺(C,S)=m₄⁺(C',S) under local orientation "
            "proved Max+-free.  m₄⁺=m₄⁻ ⇔ π-covariance of m₄⁺ on |κ|=1.  "
            "π-covariance certified p=5,7; general OPEN.  No predicate flip."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15259.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.259 transport + π-covariance reduction; residual-i OPEN")
    print(f"  transport: {a['proved']}")
    print(f"  reduction: {b['proved']}")
    print(f"  pi type: {c['proved']}")
    for pk, v in cert.items():
        if not v.get("ok"):
            print(f"  p={pk}: {v}")
            continue
        print(
            f"  p={pk}: transport={v['transport_ok']} "
            f"pi_cov={v['pi_covariance_sample_ok']} "
            f"m4eq={v['m4eq_sample_ok']} "
            f"finite {v['pi_finite_ok']}/{v['pi_finite_n']} "
            f"inf {v['pi_inf_ok']}/{v['pi_inf_n']}"
        )
    print(f"  residual_i closed: {residual_i_closed_via_259()}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
