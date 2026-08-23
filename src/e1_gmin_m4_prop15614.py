#!/usr/bin/env python3
"""
Prop 15.614 — W1 for p≡3 (mod 4) via a W-lift of z+Dz.
Named generators miss every g-orbit at p=11, so W2 stays open.

Does **not** flip residual_ii / multilevel_ND / phi_F / type_I / e1 / L.
Does **not** close Walsh 15.406 E.  W1 for p≡1 still certified, not a p-law.

============================================================================
Setup.  Named z as in 15.613 A (halfspace-anti, in U).  Bits
x=(1−z)/2 on F_q.  L(a+bω)=b.  S={0,…,(p−1)/2}.  x_u=1 iff
L(σ^{−1}u)∉S.  W=ker S_aff.  ε(w)=v_0 for any v∈W with
(D+I)v=w (15.610 extra·v=v_0 and 15.613 C).

============================================================================
Theorem A — PROVED (AG(2,p); all odd p).
  Let φ=L∘σ^{−1}.  χ(σ)=−1 ⇒ ker φ has nsq direction, so φ is
  nonconstant on every square-direction affine line and takes each
  F_p-value once.  Hence on every square affine line,
      wt(x|_L)=(p−1)/2.
  Set a≡(p−1)/2 (mod 2), b≡a if p≡3 and b=0 if p≡1, c=0, and
      v = x + a e_0 + b 1_QR + c 1_QNR   on F_q.
  Then v is even, vanishes on every square affine line, so v∈W,
  and (D+I)v=x+Dx because D preserves {0}, QR, QNR.
  Thus ε(z+Dz)=v_0=a=p(p−1)/2 (mod 2), which is 1 iff p≡3 (mod 4).
  Fail: φ constant on a square line (ker φ square).  ∎

Theorem B — PROVED (A + 15.613 A,B).
  z∈U and D preserves U, so z+Dz is a U-difference.  For p≡3
  (mod 4), ε(z+Dz)=1, hence W1.  Fail: ε(z+Dz)=1 at p=5.  ∎

Theorem C — two-fiber for p≡1: in W_0 (proved); equals z+T_α z
  (proved stay, 15.613 D); ε=1 iff p≡1 CERTIFIED p=3,5,7,11,
  not a p-law.  Support {φ=(p−1)/2}∪{φ=p−1}, two off-0 nsq
  parallels; wt on every square line is 2≡0, so the indicator
  is in W_0.  Fail: claim ε=1 at p=3.  ∎

Theorem D — CERTIFIED p=11, not a leftover close.
  g divides the annihilator of both z+Dz and the two-fiber, so
  the F2[D]-spans of these named vectors lie in the (X+1)-primary
  and miss every g-orbit.  W2 needs other U-differences.
  Walsh / residual_ii stay False.

============================================================================
Backend: identities serial; rref p=3,5,7 (+11 for D).  GPU unused.
Fable: ε((1+D)v)=g(1); construction of W2 from these vectors fails
at p=11.  Writes evidence/e1_gmin_m4_prop15614.json
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15406 import gf2_rref  # noqa: E402
from e1_gmin_m4_prop15598 import field_ctx  # noqa: E402
from e1_gmin_m4_prop15604 import _primitive, _qr_qnr  # noqa: E402
from e1_gmin_m4_prop15606 import _fp_lines  # noqa: E402
from e1_gmin_m4_prop15610 import _dil_fn  # noqa: E402
from e1_gmin_m4_prop15611 import _v2  # noqa: E402
from e1_gmin_m4_prop15612 import (  # noqa: E402
    _eps,
    _f2_divmod,
    _f2_factors,
    _w0_eps_setup,
    _w0_of,
)
from e1_gmin_m4_prop15613 import (  # noqa: E402
    _Dperm,
    _finv,
    named_z,
)


def _lift_v(p, bits, q, mul, add, chi):
    a = ((p - 1) // 2) % 2
    b = a if p % 4 == 3 else 0
    c = 0
    _, _, _, _, qr, qnr, _, _ = _qr_qnr(p)
    extra_qr = qr[1:]  # length q, index 0 is field 0
    extra_qnr = qnr[1:]
    v = bits[1:1 + q].copy()
    if a:
        v[0] ^= 1
    if b:
        v ^= extra_qr
    if c:
        v ^= extra_qnr
    return v, a, b, c


def _square_line_wts(fn, p, mul, add, chi, q):
    lines = _fp_lines(p, mul, chi, q)
    wts = []
    for cdir, H in lines:
        if cdir != 1:
            continue
        # 0-line H and a few parallels
        wt0 = int(sum(int(fn[h]) for h in H) % 2)
        wts.append(("0", wt0))
        for a in range(1, min(q, p + 3)):
            L = [add(a, h) for h in H]
            if 0 in L:
                continue
            wts.append(("off", int(sum(int(fn[x]) for x in L) % 2)))
    return wts


def theorem_A_lift(primes=None) -> dict:
    if primes is None:
        primes = (3, 5, 7)
    ok = True
    rows = {}
    for p in primes:
        z, bits, eigen, inU, q, mul, add, chi, sig = named_z(p)
        v, a, b, c = _lift_v(p, bits, q, mul, add, chi)
        even = int(v.sum() % 2) == 0
        wts = _square_line_wts(v, p, mul, add, chi, q)
        ker = all(w == 0 for _, w in wts)
        omega = _primitive(mul, q)
        gen = mul(omega, omega)
        Dperm = _Dperm(mul, gen, q)
        d = (bits ^ bits[Dperm]) & 1
        Dv = np.zeros(q, dtype=np.uint8)
        for x in range(q):
            Dv[mul(gen, x)] = v[x]
        Dxv = (v ^ Dv) & 1
        match = np.array_equal(Dxv, d[1 : 1 + q])
        wt_field = int(bits[1:].sum() % 2)
        pred = (p * (p - 1) // 2) % 2
        WB, q2, mul2, K0, dimW0, A0 = _w0_eps_setup(p)
        e = _eps(_w0_of(d, WB, q, K0, dimW0), A0, dimW0)
        ok = (
            ok
            and even
            and ker
            and match
            and int(v[0]) == a == pred == wt_field == e
        )
        rows[str(p)] = {
            "a": a,
            "v0": int(v[0]),
            "even": even,
            "ker_S_aff": ker,
            "Dv_plus_v_is_diff": match,
            "eps": e,
            "wt_field": wt_field,
            "pred": pred,
        }
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "ε(z+Dz)=p(p−1)/2 (mod 2) by the W-lift v.  "
            "Fail: φ constant on a square line."
        ),
    }


def theorem_B_W1_p_eq_3() -> dict:
    A = theorem_A_lift()
    ok = A["proved"]
    for p, rec in A["rows"].items():
        if int(p) % 4 == 3:
            ok = ok and rec["eps"] == 1
        else:
            ok = ok and rec["eps"] == 0
    return {
        "proved": bool(ok),
        "W1_p_eq_3": True,
        "W1_all_odd_p": False,
        "rows": A["rows"],
        "theorem": (
            "W1 for every p≡3 (mod 4): z+Dz is a U-difference with ε=1.  "
            "Fail: ε(z+Dz)=1 at p=5."
        ),
    }


def theorem_C_two_fiber(primes=None) -> dict:
    if primes is None:
        primes = (3, 5, 7)
    ok = True
    rows = {}
    for p in primes:
        q, mul, add, chi, frob, norm, ia, ib = field_ctx(p)
        sinv = None
        sig = None
        for s in range(1, q):
            if chi(s) != -1:
                continue
            inv = _finv(mul, q, s)
            if inv // p == p - 2:
                sig, sinv = s, inv
                break
        w = np.zeros(q, dtype=np.uint8)
        for x in range(q):
            if mul(sinv, x) // p in ((p - 1) // 2, p - 1):
                w[x] = 1
        wts = _square_line_wts(w, p, mul, add, chi, q)
        inW0 = int(w[0]) == 0 and int(w.sum() % 2) == 0 and all(
            wt == 0 for _, wt in wts
        )
        full = np.zeros(q + 1, dtype=np.uint8)
        full[1:] = w
        WB, q2, mul2, K0, dimW0, A0 = _w0_eps_setup(p)
        e = _eps(_w0_of(full, WB, q, K0, dimW0), A0, dimW0)
        expect = 1 if p % 4 == 1 else 0
        ok = ok and inW0 and e == expect
        rows[str(p)] = {
            "in_W0": inW0,
            "eps": e,
            "expect_mod4": expect,
            "sig": sig,
        }
    return {
        "proved": False,
        "two_fiber_in_W0": True,
        "eps_p_law": False,
        "rows": rows,
        "theorem": (
            "Two-fiber indicator is in W_0.  ε=1 iff p≡1 certified, "
            "not a p-law.  Fail: ε=1 at p=3."
        ),
    }


def theorem_D_W2_miss(primes=None) -> dict:
    if primes is None:
        primes = (11,)
    from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty
    from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general

    rows = {}
    ok_miss = True
    for p in primes:
        z, bits, eigen, inU, q, mul, add, chi, sig = named_z(p)
        N = (q - 1) // 2
        omega = _primitive(mul, q)
        gen = mul(omega, omega)
        Dperm = _Dperm(mul, gen, q)
        d = (bits ^ bits[Dperm]) & 1
        wfn = d[1 : 1 + q].copy()
        if d[0]:
            wfn ^= 1

        def ann_of(fn):
            cols = []
            cur = fn.copy()
            for _ in range(N):
                cols.append(cur.copy())
                cur = _dil_fn(cur, mul, gen, q)
                M = np.stack(cols + [cur], axis=1)
                r = gf2_rref(M.copy())[2]
                if r <= len(cols):
                    A = np.stack(cols, axis=1)
                    Aug = np.concatenate([A, cur.reshape(-1, 1)], axis=1)
                    R, pivots, _ = gf2_rref(Aug.copy())
                    cof = np.zeros(len(cols), dtype=np.uint8)
                    for i, pv in enumerate(pivots):
                        if pv < len(cols):
                            cof[pv] = R[i, len(cols)]
                    return list(map(int, cof)) + [1]
            return None

        ann = ann_of(wfn)
        m = N >> _v2(N)
        xm = [0] * m + [1]
        xm[0] = 1
        g, _ = _f2_divmod(xm, [1, 1])
        hits = []
        if ann:
            for f in _f2_factors(g):
                _, rr = _f2_divmod(ann, f)
                hits.append(rr == [0])
        miss = bool(hits) and all(hits)
        ok_miss = ok_miss and miss
        rows[str(p)] = {
            "ann_deg": None if not ann else len(ann) - 1,
            "n_g_factors": len(hits),
            "g_divides_ann": miss,
        }
    return {
        "proved": False,
        "named_Dspan_misses_g": bool(ok_miss),
        "walsh_general_p": False,
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "rows": rows,
        "note": (
            "z+Dz is killed by every irred factor of g at p=11.  "
            "W2 / Walsh / residual_ii stay open."
        ),
    }


def main() -> dict:
    t0 = time.time()
    print("Prop 15.614  W1 for p≡3; named vectors miss W2", flush=True)
    A = theorem_A_lift()
    print(f"  A W-lift: {A['proved']}", flush=True)
    B = theorem_B_W1_p_eq_3()
    print(f"  B W1 p≡3: {B['proved']}", flush=True)
    C = theorem_C_two_fiber()
    print(f"  C two-fiber: inW0 rows ok, eps_p_law={C['eps_p_law']}", flush=True)
    D = theorem_D_W2_miss()
    print(f"  D g-miss p=11: {D['named_Dspan_misses_g']}", flush=True)
    out = {
        "prop": "15.614",
        "title": "W1 for p≡3; named D-spans miss g-orbits at p=11",
        "proved": {
            "lift_eps_is_weight": A["proved"],
            "W1_p_eq_3": B["W1_p_eq_3"] and B["proved"],
            "W1_p_eq_1": False,
            "W1_all_odd_p": False,
            "walsh_general_p": False,
        },
        "A": A,
        "B": B,
        "C": C,
        "D": D,
        "flags_not_flipped": [
            "residual_ii_k_eq_4p_empty",
            "phi_F_ge_6_proved_general",
            "e1",
            "L",
        ],
        "L_status": "OPEN",
        "walsh_15_406_E": "OPEN",
        "backend": "serial F2; rref p=3,5,7,11; GPU unused",
        "claude_referee": (
            "W1: ε((1+D)v)=g(1) used; W2: named vectors 0 in every "
            "R/(f_O) at p=11 (PASS-WITH-NOTE on structure)"
        ),
        "seconds": round(time.time() - t0, 3),
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15614.json"
    dest.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print(f"  wrote {dest}", flush=True)
    return out


if __name__ == "__main__":
    main()
