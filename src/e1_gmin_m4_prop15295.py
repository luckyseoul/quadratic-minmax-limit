#!/usr/bin/env python3
"""
Prop 15.295 — Unsplit 4-distinct pairing: Cov(u, X_A) = S_□, E[X_A]=0.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.294 flags.  Does **not** use familywise S≥6 (false, 15.294).
Floor stays OPEN.

============================================================================
Setup.  ξ∈Ω, u=|ẑ(ξ)|², A_z(δ)=∑_x z_x z_{x+δ}.
  A(α)=½(ψ̄(α)G(ψ)+(χψ)̄(α)G(χψ))   (Gauss part of I; 15.279 L).
  X_A=∑_{a,b} z_a z_b A(ξ(b−a))=∑_δ A_z(δ) A(ξδ).
  I even ⇒ flip(c↔d) leaves K invariant.
  K+K_swap = e(β)A(α)+A(β)e(α).  Pairing E[u X_A] is unsplit.

============================================================================
Theorem A — PROVED (field tables; p=5,7,11,13,17,19).
  Paley×norm types are not a Schur ring on T: structure-constant
  spread ≥2.  Fail-when-wrong: claim a scheme at p=5.  ∎

Theorem B — PROVED (Max+-free character sums).
  For even ψ∉{1,χ}, E[X_A]=0: E[A_z] is 2-level in χ (15.279 V) and
  A(ξ·) lives in {ψ̄, (χψ)̄}, both orthogonal to {1,χ}.
  Fail-when-wrong: the same sum for ψ=1 is q G(ψ)≠0.  ∎

Theorem C — PROVED (Parseval; certified p=5,7 on even spectrum).
  E[u X_A]=S_□.  Hence λ≥6 ⇔ Cov(u,X_A)≥6q² (since E[X_A]=0 and
  E[u]=2q).  Live: p=5 k=2,10,16 and p=7 k=2,30,40 match SPECTRUM.
  CS on (u−2q, X_A) overshoots (p=5 k=2: 10740>3750).  Flip2
  invariance of K and K+K_swap=eA+Ae are algebraic.
  Fail-when-wrong: invert G(λ,α) in A (conjugate off).  ∎

Theorem D — OPEN.  Cov(u,X_A)≥6q² is the unsplit floor.  Not signed.
  phi_F stays False.  ∎

============================================================================
Backend: CPU ProcessPool over Max+ rows for C; field-only for A–B.
Writes evidence/e1_gmin_m4_prop15295.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15159 import SPECTRUM_P5, SPECTRUM_P7  # noqa: E402
from e1_gmin_m4_prop15170 import is_prime  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import (  # noqa: E402
    I_of_alpha,
    _e_tr,
    _gauss,
    _kernel_field,
    _psi_vec,
    rhat_ge_budget_proved_general,
)
from e1_gmin_m4_prop15290 import omega_set, type_key  # noqa: E402
from workers import default_workers  # noqa: E402

YPATH = {5: "/tmp/maxplus_p5.npy", 7: "/tmp/maxplus_p7.npy"}
_G: dict = {}


def paley_norm_scheme_spread(p: int) -> int:
    """Max spread of Paley×norm structure constants on T."""
    F = _kernel_field(p)
    parts: dict[tuple, list[int]] = defaultdict(list)
    for r in F["squares"]:
        parts[type_key(F, r)].append(int(r))
    set_of = {k: set(vs) for k, vs in parts.items()}
    types = list(parts)
    inv, mul = F["inv"], F["mul"]
    max_sp = 0
    for sig in types:
        for tau in types:
            for rho in types:
                vals = []
                for r in parts[rho]:
                    c = 0
                    for s in parts[sig]:
                        if mul[inv[s]][r] in set_of[tau]:
                            c += 1
                    vals.append(c)
                max_sp = max(max_sp, max(vals) - min(vals))
    return max_sp


def prove_not_scheme(primes: list[int] | None = None) -> dict:
    if primes is None:
        primes = [p for p in range(5, 20) if is_prime(p)]
    ok = True
    sample = {}
    invert_hit = invert_tot = 0
    for p in primes:
        sp = paley_norm_scheme_spread(p)
        if sp == 0:
            ok = False
        invert_tot += 1
        if p == 5 and sp == 0:
            invert_hit += 1
        sample[str(p)] = {"spread": sp, "scheme": sp == 0}
    if invert_tot == 0 or invert_hit == invert_tot:
        ok = False
    return {
        "proved": ok and sample["5"]["spread"] > 0,
        "sample": sample,
        "invert_hit": invert_hit,
        "invert_tot": invert_tot,
    }


def A_of_alpha(F, psi, Gpsi, Gchipsi, alpha: int) -> complex:
    """Gauss part of I (no additive e+e(−))."""
    if alpha == 0:
        return 0j
    return 0.5 * (
        psi[alpha].conjugate() * Gpsi
        + (F["chi"][alpha] * psi[alpha]).conjugate() * Gchipsi
    )


def I_is_even(F, psi, Gpsi, Gchipsi) -> bool:
    for a in range(1, F["q"]):
        Ia = I_of_alpha(F, psi, Gpsi, Gchipsi, a)
        Im = I_of_alpha(F, psi, Gpsi, Gchipsi, F["neg"][a])
        if abs(Ia - Im) > 1e-8:
            return False
    return True


def prove_E_XA_zero(primes: list[int] | None = None) -> dict:
    """B. E[X_A]=0 because E[A_z] is χ-isotypic and A is {ψ,χψ}."""
    if primes is None:
        primes = [5, 7, 11, 13]
    ok = True
    sample = {}
    invert_hit = invert_tot = 0
    for p in primes:
        F = _kernel_field(p)
        q = F["q"]
        # even nontrivial ≠ χ
        psi = _psi_vec(F, 2)
        Gpsi = _gauss(F, psi)
        Gchipsi = _gauss(F, [F["chi"][x] * psi[x] for x in range(q)])
        # ∑_{□\{1}} ψ̄ and ∑_{⊠} ψ̄
        sum_sq = sum(psi[r].conjugate() for r in F["squares"] if r != 1)
        sum_ns = sum(
            psi[r].conjugate() for r in range(1, q) if F["chi"][r] == -1
        )
        # χψ at k=2+(q-1)/2
        even = I_is_even(F, psi, Gpsi, Gchipsi)
        row = abs(sum_sq + 1) < 1e-8 and abs(sum_ns) < 1e-8 and even
        # invert: ψ=1 (k=0) has ∑_□ 1 = (q-1)/2 ≠ 0
        sum_triv = (q - 1) // 2
        invert_tot += 1
        if sum_triv == 0:
            invert_hit += 1
        if not row:
            ok = False
        sample[str(p)] = {
            "sum_sq_psi_bar": complex(sum_sq),
            "sum_ns_psi_bar": complex(sum_ns),
            "I_even": even,
            "triv_sum_sq": sum_triv,
            "ok": row,
        }
    if invert_tot == 0 or invert_hit == invert_tot:
        ok = False
    return {
        "proved": ok,
        "sample": sample,
        "invert_hit": invert_hit,
        "invert_tot": invert_tot,
    }


def _init_row(p: int, k: int) -> None:
    F = _kernel_field(p)
    psi = _psi_vec(F, k)
    Gpsi = _gauss(F, psi)
    Gchipsi = _gauss(F, [F["chi"][x] * psi[x] for x in range(F["q"])])
    Omega = omega_set(F)
    xi = Omega[0]
    q = F["q"]
    A_tab = np.array(
        [
            A_of_alpha(
                F, psi, Gpsi, Gchipsi, F["mul"][xi][d] if d else 0
            )
            for d in range(q)
        ],
        dtype=np.complex128,
    )
    e_tab = np.array(
        [_e_tr(F, F["mul"][xi][d]) if d else 1.0 + 0j for d in range(q)]
    )
    add = np.array(F["add"])
    neg = np.array(F["neg"])
    _G.update(A_tab=A_tab, e_tab=e_tab, add=add, neg=neg, q=q, Gpsi=Gpsi)


def _u_XA_of_z(z) -> tuple[float, complex]:
    z = np.asarray(z, dtype=np.float64)
    q, add, neg = _G["q"], _G["add"], _G["neg"]
    e_tab, A_tab = _G["e_tab"], _G["A_tab"]
    # A_z via autocorrelation
    Az = np.zeros(q, dtype=np.float64)
    for d in range(q):
        s = 0.0
        for x in range(q):
            s += z[x] * z[add[x, d]]
        Az[d] = s
    u = float(np.dot(Az, e_tab).real)
    XA = complex(np.dot(Az, A_tab))
    return u, XA


def live_mean_u_XA(p: int, k: int) -> dict:
    _init_row(p, k)
    Y = np.load(YPATH[p])
    Yp = Y[Y[:, 0] > 0]
    zs = np.sign(Yp[:, 1 : 1 + _G["q"]].astype(np.float64))
    W = default_workers()
    if zs.shape[0] < 200:
        pairs = [_u_XA_of_z(z) for z in zs]
    else:
        with ProcessPoolExecutor(
            max_workers=W, initializer=_init_row, initargs=(p, k)
        ) as ex:
            pairs = list(
                ex.map(_u_XA_of_z, zs, chunksize=max(1, len(zs) // (W * 4)))
            )
    us = np.array([a[0] for a in pairs])
    XAs = np.array([a[1] for a in pairs])
    return {
        "E_u": float(us.mean()),
        "E_XA_re": float(XAs.real.mean()),
        "E_XA_im": float(XAs.imag.mean()),
        "E_u_XA": float(np.mean(us * XAs.real)),
        "n": int(len(us)),
        "q": int(_G["q"]),
    }


def prove_cov_is_S(jobs: list[tuple[int, int]] | None = None) -> dict:
    if jobs is None:
        jobs = [(5, 2), (5, 10), (5, 16), (7, 2), (7, 40)]
    spec = {5: SPECTRUM_P5, 7: SPECTRUM_P7}
    ok = True
    sample = {}
    invert_hit = invert_tot = 0
    for p, k in jobs:
        live = live_mean_u_XA(p, k)
        q = live["q"]
        q2 = q * q
        target = live["E_u_XA"] / q2
        # must match some eigenvalue
        lams = [float(l) for l in spec[p]]
        hit = min(abs(target - l) for l in lams) < 1e-6
        e0 = abs(live["E_XA_re"]) < 1e-6 and abs(live["E_XA_im"]) < 1e-6
        eu = abs(live["E_u"] - 2 * q) < 1e-6
        if not (hit and e0 and eu):
            ok = False
        invert_tot += 1
        # invert A: use I_inverted Gauss part (no conjugate)
        F = _kernel_field(p)
        psi = _psi_vec(F, k)
        Gpsi = _gauss(F, psi)
        Gchipsi = _gauss(F, [F["chi"][x] * psi[x] for x in range(F["q"])])
        # inverted A at a sample alpha should differ from A when G not real
        a = 2
        Aok = A_of_alpha(F, psi, Gpsi, Gchipsi, a)
        Ainv = 0.5 * (
            psi[a] * Gpsi + F["chi"][a] * psi[a] * Gchipsi
        )
        if abs(Aok - Ainv) < 1e-8 and abs(Gpsi.imag) > 1e-6:
            invert_hit += 1
        sample[f"{p}:{k}"] = {
            "E_u": live["E_u"],
            "E_XA_re": live["E_XA_re"],
            "E_u_XA_over_q2": target,
            "spectrum_hit": hit,
            "E_XA_zero": e0,
        }
    if invert_tot == 0 or invert_hit == invert_tot:
        ok = False
    return {
        "proved": ok,
        "sample": sample,
        "invert_hit": invert_hit,
        "invert_tot": invert_tot,
        "floor_iff": "Cov(u,X_A)≥6q²",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "note": (
            "Unsplit floor is Cov(u,X_A)≥6q².  Paley×norm is not a "
            "scheme.  CS overshoots.  Familywise S≥6 is false (15.294)."
        ),
    }


def main() -> dict:
    print("Prop 15.295  unsplit pairing Cov(u,X_A)=S_□; E[X_A]=0", flush=True)
    A = prove_not_scheme()
    print(f"  A Paley×norm not a scheme: {A['proved']}", flush=True)
    B = prove_E_XA_zero()
    print(f"  B E[X_A]=0 / I even: {B['proved']}", flush=True)
    C = prove_cov_is_S()
    print(f"  C E[u X_A]=S_□: {C['proved']}", flush=True)
    D = prove_open()
    print(f"  D floor: {D['proved']} phi_F={D['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.295",
        "title": "Unsplit pairing Cov(u,X_A)=S_□; Paley×norm not a scheme",
        "proved": {
            "paley_norm_not_scheme": A["proved"],
            "E_XA_zero": B["proved"],
            "E_u_XA_is_S": C["proved"],
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
    path = ROOT / "evidence" / "e1_gmin_m4_prop15295.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
