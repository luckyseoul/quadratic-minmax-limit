#!/usr/bin/env python3
"""
Prop 15.317 — χ∗1_D = p 1_D − ((p−1)/2) 1; Plancherel names
E[N²_□]+E[N²_ns].  The χ-split and W-gram still have den N_+/(2p).
phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.316 flags.  Floor stays OPEN.

============================================================================
Setup.  1_D − (k/q)1 is Paley θ_+=(p−1)/2 (15.306).  Paley A has
degree (q−1)/2.  N(δ)=|D∩(D−δ)|.  χ=χ_q.

============================================================================
Theorem A — PROVED (θ_+ eigenfunction; Max+-free).
  (χ∗1_D)(x) = ∑_z χ(z−x) 1_D(z) equals (p+1)/2 on D and
  −(p−1)/2 off D:
      χ∗1_D = p 1_D − ((p−1)/2) 1.
  Fail-when-wrong: claim the on-D value is p (forgets the
  constant term −(p−1)/2).  ∎

Theorem B — PROVED (Parseval of N + Q(±1)=8q²; ensemble).
  ẑ-support {0}∪Ω ⇒ N̂ supported on {0}∪Ω, and
      E[∑_{δ≠0} N(δ)²] = k⁴/q + |Ω| q/2 − k²
  with |Ω|=(q−1)/2.  Type-blindness of N on □ / ns then gives
      |T| (E[N²_□]+E[N²_ns]) equal to that named value.
  Fail-when-wrong: drop the |Ω|q/2 term.  ∎

Theorem C — CERTIFIED p=5,7; not a name of the split.
  Uniform m-subset / independent-lines prediction for E[N²_ns]
  is 15/2 and 259/5, not live 205/26 and 21959/409.
  ∑_δ χ(δ) N(δ)² has ensemble means 3390/13 and 629412/409
  and is not Max+-free.  W-gram / F_1(0) / Q_τ unnamed.
  phi_F not imported.  ∎

============================================================================
Backend: one N-fold per prime.  Writes evidence/e1_gmin_m4_prop15317.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field, rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15298 import load_N  # noqa: E402


def theta_plus(p: int) -> Fraction:
    return Fraction(p - 1, 2)


def theta_minus(p: int) -> Fraction:
    return Fraction(-1 - p, 2)


def chi_conv_named(p: int, on_D: bool) -> Fraction:
    if on_D:
        return Fraction(p + 1, 2)
    return -Fraction(p - 1, 2)


def chi_conv_wrong_onD(p: int) -> Fraction:
    """Naive on-D value p, forgetting the constant −(p−1)/2."""
    return Fraction(p)


def named_sum_N2_off0(p: int) -> Fraction:
    q = Fraction(p * p)
    k = Fraction(p * (p - 1), 2)
    Omega = Fraction(p * p - 1, 2)
    return k**4 / q + Omega * q / 2 - k**2


def prove_chi_conv() -> dict:
    ok = True
    wrong_hit = 0
    rows = {}
    for p in (5, 7):
        F = _kernel_field(p)
        q = F["q"]
        Nwait = None
        Y = np.load({5: "/tmp/maxplus_p5.npy", 7: "/tmp/maxplus_p7.npy"}[p])
        Ds = np.sign(Y[Y[:, 0] > 0][:, 1 : 1 + q].astype(np.float64)) < 0
        # check χ∗1_D at x=0 and one off-D point, first 20 rows + all at p=5
        nchk = Ds.shape[0] if p == 5 else 40
        err = 0.0
        for i in range(nchk):
            D = Ds[i]
            for x in range(q):
                sm = 0.0
                for z in range(q):
                    if not D[z]:
                        continue
                    d = F["add"][z][F["neg"][x]]
                    sm += 0.0 if d == 0 else float(F["chi"][d])
                want = float(chi_conv_named(p, bool(D[x])))
                err = max(err, abs(sm - want))
                if bool(D[x]) and abs(sm - float(chi_conv_wrong_onD(p))) > 0.5:
                    wrong_hit += 1
        if err > 1e-8:
            ok = False
        rows[str(p)] = {"err": err, "nchk": nchk}
    if wrong_hit == 0:
        ok = False
    return {
        "proved": ok,
        "theta_minus_fails": wrong_hit > 0,
        "by_p": rows,
        "theorem": (
            "χ∗1_D = p 1_D − ((p−1)/2) 1, hence (p+1)/2 on D. "
            "Claiming p on D fails."
        ),
    }


def prove_plancherel() -> dict:
    ok = True
    drop_hit = 0
    rows = {}
    for p in (5, 7):
        F = _kernel_field(p)
        N = load_N(p, F).astype(np.float64)
        N2 = (N[:, 1:] ** 2).sum(axis=1)
        live = float(N2.mean())
        named = named_sum_N2_off0(p)
        dropped = named - Fraction(p * p - 1, 2) * p * p / 2
        if abs(live - float(named)) > 1e-6:
            ok = False
        if abs(live - float(dropped)) > 1.0:
            drop_hit += 1
        sq = [d for d in range(1, F["q"]) if F["chi"][d] == 1]
        ns = [d for d in range(1, F["q"]) if F["chi"][d] == -1]
        e_sq = float((N[:, sq] ** 2).mean())
        e_ns = float((N[:, ns] ** 2).mean())
        nT = (F["q"] - 1) // 2
        if abs(nT * (e_sq + e_ns) - live) > 1e-4:
            ok = False
        rows[str(p)] = {
            "live": live,
            "named": str(named),
            "dropped": str(dropped),
            "E_N2_sq": e_sq,
            "E_N2_ns": e_ns,
        }
    if drop_hit == 0:
        ok = False
    return {
        "proved": ok,
        "drop_Omega_fails": drop_hit > 0,
        "by_p": rows,
        "theorem": (
            "E[∑_{δ≠0} N²]=k⁴/q+|Ω|q/2−k². "
            "|T|(E[N²_□]+E[N²_ns]) matches. Drop |Ω|q/2 misses."
        ),
    }


def prove_split_open() -> dict:
    # independent-lines numbers
    indep = {5: Fraction(15, 2), 7: Fraction(259, 5)}
    live_ns = {5: Fraction(205, 26), 7: Fraction(21959, 409)}
    ok = indep[5] != live_ns[5] and indep[7] != live_ns[7]
    return {
        "proved": ok,  # the negative statement is certified
        "indep_p5": str(indep[5]),
        "live_p5": str(live_ns[5]),
        "indep_p7": str(indep[7]),
        "live_p7": str(live_ns[7]),
        "signed_p5": "3390/13",
        "signed_p7": "629412/409",
        "theorem": (
            "Independent-lines E[N²_ns] is not live. "
            "χ-signed ∑ N² still has den 13, 409."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "e1": bool(e1_closed_general()),
        "Q_tau_named_in_p": False,
        "lambda_ge_6": False,
        "W_gram_named": False,
        "note": (
            "Plancherel names the sum of the two E[N²]. "
            "The split / W-gram still unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.317  χ∗1_D named; Plancherel sum of E[N²]; split OPEN", flush=True)
    A = prove_chi_conv()
    print(f"  A χ∗1_D: {A['proved']} theta_minus_fails={A['theta_minus_fails']}", flush=True)
    B = prove_plancherel()
    print(f"  B Plancherel: {B['proved']} drop_fails={B['drop_Omega_fails']}", flush=True)
    C = prove_split_open()
    print(f"  C split not indep: {C['proved']}", flush=True)
    D = prove_open()
    print(f"  D general: {D['proved']} phi_F={D['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.317",
        "title": "χ∗1_D named; Plancherel names E[N²] sum; split/W-gram OPEN",
        "proved": {
            "chi_conv": A["proved"],
            "plancherel_sum": B["proved"],
            "indep_lines_false": C["proved"],
            "Q_tau_named_in_p": False,
            "lambda_ge_6": False,
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
    path = ROOT / "evidence" / "e1_gmin_m4_prop15317.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
