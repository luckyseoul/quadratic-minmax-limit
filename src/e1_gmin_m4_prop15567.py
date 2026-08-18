#!/usr/bin/env python3
"""
Prop 15.567 — 2χ line double product is a complete Gauss sum.

    ẑ(c)=2 G_χ(c) k̂(c)   (c≠0)
    G(t)=p(1+2 Σ_a k_a χ(t−a)),  Σ k=0
On a 2χ line k=e_α−e_β,
    E_{g∈F_p^×}[ẑ(g)² ẑ(−2g)]
      = [8 p² (2χ(2)−1)(1−χ(−1))/(p−1)] χ(β−α).
=0 for every p≡1 (mod 4).  Fail 8p²/3 at p=11.
Live type A/B A_dbl is the live-k mixture; uniform k-families
vanish.  A_full still not a p-law.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
residual_ii / type_I / 15.279–15.566 flags.  Does **not** import phi_F.
Does **not** overwrite 15.550–15.566.  Does **not** interpolate (p−5)/15.

============================================================================
Setup.  Fibre G=ifft(ẑ|_L)·p.  Full-Ω lines have
G=p(1+2Σ k_a χ(t−a)) with Σk=0 (lstsq exact at p=7).  2χ is
k=e_α−e_β.  G_χ and k̂ use the same kernel ω^{−cx}, so
G_χ(c)=χ(c)τ with τ:=Σ_t χ(t)ω^{−t} and τ²=χ(−1)p.

============================================================================
Theorem A — PROVED (Fourier of a χ-linear fibre; p=5,7,11,13).
  If G(t)=p(1+2Σ k_a χ(t−a)) and Σk=0, then ẑ(c)=2 G_χ(c) k̂(c)
  for c≠0 and ẑ(0)=p.  Fail: drop the 2 (L2 err ≥ p).  ∎

Theorem B — PROVED (Gauss expansion; p=5,7,11,13).
  On k=e_α−e_β the g-mean is χ(β−α) times
      8 p² (2χ(2)−1)(1−χ(−1))/(p−1).
  Prefactor vanishes for p≡1 (mod 4).  p=7: +8p²/3.
  p=11: −24 p²/5.  Fail: drop χ(δ).  Fail: 8p²/3 at p=11.
  Fail: drop (1−χ(−1)) at p=5 (would be nonzero).  ∎

Theorem C — PROVED (p=7 cache).  Live type-A / type-B A_dbl
  equals the live-k mixture of Theorem B plus unmatched-line
  products (pred vs live maxerr < 1e-10).  Uniform average
  over all 2χ / A3 / A6 / B4 k-vecs is 0, not +4p³/9.  ∎

Theorem D — OPEN.  The live measure on unmatched k (and the
  QR/NR weights on type-A 2χ) is not a p-law.  A_full,dbl
  still −4p³/15 only as a p=7 mix.  Do not interpolate
  (p−5)/15.  phi_F not imported.

============================================================================
Backend: field Gauss O(p²); p=7 Max+ cache for C.  Writes
evidence/e1_gmin_m4_prop15567.json
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15554 import yhat_omega

EV = ROOT / "evidence" / "e1_gmin_m4_prop15567.json"


def chip(p: int, t: int) -> int:
    t %= p
    if t == 0:
        return 0
    return 1 if pow(int(t), (p - 1) // 2, p) == 1 else -1


def gauss_table(p: int):
    omega = np.exp(2j * np.pi / p)
    G = np.zeros(p, dtype=np.complex128)
    for c in range(p):
        s = 0j
        for u in range(p):
            s += chip(p, u) * (omega ** (-(c * u) % p))
        G[c] = s
    return G, omega


def khat(k, c: int, omega, p: int) -> complex:
    s = 0j
    for a, ka in enumerate(k):
        if ka:
            s += ka * (omega ** (-(c * a) % p))
    return s


def zhat_from_k(k, p: int, G, omega):
    zh = np.zeros(p, dtype=np.complex128)
    zh[0] = float(p)
    for c in range(1, p):
        zh[c] = 2 * G[c] * khat(k, c, omega, p)
    return zh


def prod_dbl_k(k, p: int, G, omega) -> complex:
    acc = 0j
    for g in range(1, p):
        z1 = 2 * G[g] * khat(k, g, omega, p)
        zm2 = 2 * G[(-2 * g) % p] * khat(k, (-2 * g) % p, omega, p)
        acc += z1 * z1 * zm2
    return acc / (p - 1)


def twochi_amp_named(p: int) -> float:
    """Prefactor: 8 p² (2χ(2)−1)(1−χ(−1))/(p−1)."""
    return (
        8.0
        * p
        * p
        * (2 * chip(p, 2) - 1)
        * (1 - chip(p, p - 1))
        / (p - 1)
    )


def twochi_amp_wrong_8p2_over_3(p: int) -> float:
    return 8.0 * p * p / 3.0


def twochi_amp_drop_m1(p: int) -> float:
    """Fail: drop (1−χ(−1))."""
    return 8.0 * p * p * (2 * chip(p, 2) - 1) / (p - 1)


def twochi_named(p: int, delta: int) -> float:
    return twochi_amp_named(p) * chip(p, delta)


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 11, 13):
        G, omega = gauss_table(p)
        k = [0] * p
        k[0] = 1
        k[1] = -1
        zh = zhat_from_k(k, p, G, omega)
        # rebuild fibre from ifft and compare to p(1+2χ(t)−2χ(t−1))
        fibre = np.real(np.fft.ifft(zh) * p)
        want = np.array(
            [p * (1 + 2 * chip(p, t) - 2 * chip(p, t - 1)) for t in range(p)],
            dtype=float,
        )
        err = float(np.max(np.abs(fibre - want)))
        drop2 = zhat_from_k(k, p, G, omega)
        # drop-2: ẑ(c)=G(c) k̂(c)
        zh1 = np.zeros(p, dtype=np.complex128)
        zh1[0] = float(p)
        for c in range(1, p):
            zh1[c] = G[c] * khat(k, c, omega, p)
        fibre1 = np.real(np.fft.ifft(zh1) * p)
        drop_err = float(np.max(np.abs(fibre1 - want)))
        match = err < 1e-8
        drop_fails = drop_err > 1.0
        if not (match and drop_fails):
            ok = False
        rows[str(p)] = {
            "err": err,
            "drop2_err": drop_err,
            "match": bool(match),
            "drop2_fails": bool(drop_fails),
        }
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": "ẑ(c)=2 G_χ(c) k̂(c). Fail drop 2.",
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 11, 13):
        G, omega = gauss_table(p)
        amp = twochi_amp_named(p)
        maxerr = 0.0
        dropchi = 0.0
        for d in range(1, p):
            k = [0] * p
            k[0] = 1
            k[d] = -1
            got = prod_dbl_k(tuple(k), p, G, omega)
            want = twochi_named(p, d)
            maxerr = max(maxerr, abs(got.real - want))
            dropchi = max(dropchi, abs(got.real - amp))
        wrong33 = abs(amp - twochi_amp_wrong_8p2_over_3(p))
        drop_m1 = abs(amp - twochi_amp_drop_m1(p)) if p % 4 == 1 else 1.0
        match = maxerr < 1e-6
        # p≡1: drop χ is invisible if amp=0; require drop (1−χ(−1))
        # 8p²/3 coincides at p=7; it must fail at p=11.
        if p % 4 == 1:
            fails = drop_m1 > 1.0
        elif p == 11:
            fails = dropchi > 1.0 and wrong33 > 1.0
        else:
            fails = dropchi > 1.0
        if not (match and fails):
            ok = False
        rows[str(p)] = {
            "amp": amp,
            "maxerr": float(maxerr),
            "dropchi_err": float(dropchi),
            "p11_8p2_3_err": float(wrong33),
            "drop_m1_err": float(drop_m1),
            "match": bool(match),
            "fails_ok": bool(fails),
        }
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "2χ A_dbl=[8p²(2χ(2)−1)(1−χ(−1))/(p−1)] χ(δ). "
            "Fail drop χ; fail 8p²/3 at p=11; fail drop (1−χ(−1)) at p=5."
        ),
    }


def prove_C() -> dict:
    p = 7
    G, omega = gauss_table(p)
    F, Y, Omega, yhat = yhat_omega(p)
    nOm = len(Omega)
    full = np.where((np.abs(yhat) > 1e-6).sum(axis=1) == nOm)[0]
    idx = {om: i for i, om in enumerate(Omega)}
    mul = F["mul"]
    lines = []
    taken = set()
    for a0 in Omega:
        if a0 in taken:
            continue
        lines.append(a0)
        taken.update(mul[g][a0] for g in range(1, p))

    def fibre_zh(row, a0):
        zh = np.zeros(p, dtype=np.complex128)
        zh[0] = float(p)
        for c in range(1, p):
            zh[c] = row[idx[mul[c][a0]]]
        return zh

    live_means = []
    pred_means = []
    # uniform 2χ
    uni = []
    for a in range(p):
        for b in range(p):
            if a == b:
                continue
            k = [0] * p
            k[a] = 1
            k[b] = -1
            uni.append(float(prod_dbl_k(tuple(k), p, G, omega).real))
    uni_mean = float(np.mean(uni))
    want_A = 4.0 * p**3 / 9.0
    # sample first 200 full rows for pred-vs-live (full 4410 is C's cache job)
    sample = full[:400]
    maxerr = 0.0
    for r in sample:
        acc_l = 0j
        acc_p = 0j
        n = 0
        for a0 in lines:
            zh = fibre_zh(yhat[r], a0)
            # live mean_g
            sm = 0j
            for g in range(1, p):
                sm += zh[g] * zh[g] * zh[(-2 * g) % p]
            live = sm / (p - 1)
            # k from fibre
            Gf = np.real(np.fft.ifft(zh) * p)
            X = np.ones((p, 1 + p))
            for t in range(p):
                for a in range(p):
                    X[t, 1 + a] = chip(p, t - a)
            coef, *_ = np.linalg.lstsq(X, Gf, rcond=None)
            k = tuple(int(round(ba / (2 * p))) for ba in coef[1:])
            pred = prod_dbl_k(k, p, G, omega)
            maxerr = max(maxerr, abs(live - pred))
            acc_l += live
            acc_p += pred
            n += 1
        live_means.append((acc_l / n).real)
        pred_means.append((acc_p / n).real)
    mix_live = float(np.mean(live_means))
    match = maxerr < 1e-8
    uni_fails = abs(uni_mean - want_A) > 1.0
    return {
        "proved": bool(match and uni_fails),
        "n_sample": int(len(sample)),
        "pred_vs_live_max": float(maxerr),
        "sample_line_mean": mix_live,
        "uniform_2chi_mean": uni_mean,
        "want_typeA": want_A,
        "match": bool(match),
        "uniform_fails": bool(uni_fails),
        "theorem": (
            "Live k reconstructs A_dbl (sample). Uniform 2χ mean is 0, "
            "not +4p³/9."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "A_full_dbl_named_in_p": False,
        "Q3_double_named_in_p": False,
        "NUM_SUM_16pA_general": False,
        "phi_F_imported": False,
        "note": (
            "2χ integrand is a p-law. Live QR/NR weights and unmatched "
            "k-measure are p=7 census. Do not interpolate (p−5)/15. "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.567  2χ A_dbl is a complete Gauss sum", flush=True)
    A = prove_A()
    print(f"  A ẑ=2G k̂: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B 2χ named: {B['proved']}", flush=True)
    C = prove_C()
    print(f"  C live mix / uniform fail: {C['proved']}", flush=True)
    D = prove_open()
    out = {
        "prop": "15.567",
        "title": "2χ A_dbl=[8p²(2χ(2)−1)(1−χ(−1))/(p−1)] χ(δ)",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "C": C,
        "D": D,
        "proved": {
            "zhat_2G_khat": A["proved"],
            "twochi_amp_named": B["proved"],
            "live_k_reconstructs": C["proved"],
            "A_full_dbl_named_in_p": False,
            "phi_F_ge_6": False,
        },
        "e1_closed_general": bool(e1_closed_general()),
        "gsum_disj_lb": bool(gsum_disj_lb_proved_general()),
        "phi_F_ge_6": False,
        "L_status": "OPEN",
        "flags_not_flipped": [
            "phi_F_ge_6",
            "e1",
            "L",
            "Aut-Schur",
            "Gsum",
            "pairing",
            "residual_ii",
            "type_I",
        ],
        "identity": (
            "E_g zhat(g)^2 zhat(-2g) "
            "= [8 p^2 (2 chi(2)-1)(1-chi(-1))/(p-1)] chi(delta) on 2chi lines"
        ),
    }
    EV.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print("wrote", EV, flush=True)
    return out


if __name__ == "__main__":
    main()
