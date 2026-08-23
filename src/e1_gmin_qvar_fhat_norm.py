#!/usr/bin/env python3
"""
F̂(ψ) is not a Paley-field square or field-norm.

Unnumbered kill.  Does **not** prove F̂(ψ)≥0.  Does **not** flip leftover
1/2/3, phi_F, L, Aut-Schur, Gsum, pairing.

SETUP (15.279 D/L/M)
  Even ψ∉{1,χ} of F_q^×, q=p².  λ(ψ)=S_□/q²=8+R̂_rest/q².
  Floor λ≥6 ⇔ F̂(ψ)≥0 with F̂=S_□−6q²=(λ−6)q² (15.279 M).
  Paley character field K = Q(√p*) Q(ζ_{q−1}), p*=(−1)^{(p−1)/2} p.

KILL (exact Fractions from the named Φ|_F spectrum)
  At p=5, λ ∈ {80,144,176}/13 and F̂=(λ−6)·625 has reduced
  denominator 13.  13²≡1 (mod 24), so residue degree of 13 in
  Q(ζ_{24}) is 2; (5/13)=−1 so 13 is inert in Q(√5).  v_13(F̂)=−1
  is odd, hence F̂ is not a norm from K (nor a square in Q).
  Fail: claim 1250/13 is a square; claim v_13 is even.

  At p=7, λ ∈ {3072,3360,3648,4032,4320}/409 and every F̂ has
  reduced denominator 409.  25²≡1 (mod 48) so f_{409}(Q(ζ_{48}))=2;
  (−7/409)=−1 so 409 is inert in Q(√−7).  v_409(F̂)=−1 odd.
  409 does not divide p(p²−1)(p²+1).  Fail: claim F̂ is a Gauss /
  Jacobi monomial-norm in p only; claim v_409 even.

  Therefore a Bochner / convolution-square certificate F̂=|A|² with
  A in the Paley character field is false.  Positivity of F̂ remains
  OPEN; it cannot be a completed square in that field.

Backend: serial exact arithmetic (a handful of Fractions).  GPU unused.
"""
from __future__ import annotations

import json
import math
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import is_prime  # noqa: E402
from e1_gmin_m4_prop15589 import q_of  # noqa: E402


# Named Φ|_F eigenvalues (15.242 / scripts/lambda_as_char_moment.py).
# Keys are even character indices k with α(x)=exp(2πi k dlog x /(q−1)).
LAMBDA = {
    5: {
        4: Fraction(80, 13),
        8: Fraction(144, 13),
        12: Fraction(176, 13),  # quadratic, k=(q-1)/2
    },
    7: {
        4: Fraction(3360, 409),
        8: Fraction(4032, 409),
        12: Fraction(3648, 409),
        16: Fraction(3072, 409),  # binding min
        20: Fraction(3360, 409),
        24: Fraction(4320, 409),  # quadratic
    },
}


def fhat_from_lambda(p: int, lam: Fraction) -> Fraction:
    q = q_of(p)
    return (lam - 6) * q * q


def trial_factor(n: int) -> dict[int, int]:
    """Prime factorization of |n|.  Empty for 0."""
    n = abs(int(n))
    if n <= 1:
        return {}
    out: dict[int, int] = {}
    while n % 2 == 0:
        out[2] = out.get(2, 0) + 1
        n //= 2
    f = 3
    while f * f <= n:
        while n % f == 0:
            out[f] = out.get(f, 0) + 1
            n //= f
        f += 2
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def is_q_square(x: Fraction) -> bool:
    if x < 0:
        return False
    n, d = x.numerator, x.denominator
    sn = int(math.isqrt(abs(n)))
    sd = int(math.isqrt(d))
    return sn * sn == abs(n) and sd * sd == d


def kronecker(a: int, n: int) -> int:
    """Kronecker symbol (a/n)."""
    if n == 0:
        return 1 if abs(a) == 1 else 0
    if n < 0:
        s = -1 if a < 0 else 1
        return s * kronecker(a, -n)
    a %= n if n else a
    # handle 2
    t = 1
    while n % 2 == 0:
        if a % 2 == 0:
            return 0
        n //= 2
        a8 = a % 8
        if a8 in (3, 5):
            t = -t
    if n == 1:
        return t
    # odd part: quadratic reciprocity
    a %= n
    while a != 0:
        while a % 2 == 0:
            a //= 2
            n8 = n % 8
            if n8 in (3, 5):
                t = -t
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            t = -t
        a %= n
    return 0 if n > 1 else t


def paley_star(p: int) -> int:
    return p if (p % 4 == 1) else -p


def quadratic_f(prime: int, d: int) -> int | None:
    """Residue degree of odd prime in Q(√d), or None if ramified/divides 2d."""
    if prime == 2 or d % prime == 0:
        return None
    chi = kronecker(d, prime)
    if chi == 0:
        return None
    return 1 if chi == 1 else 2


def mul_order_mod(a: int, m: int) -> int:
    """Order of a in (Z/mZ)*; m coprime to a."""
    a %= m
    if math.gcd(a, m) != 1:
        raise ValueError(f"{a} not a unit mod {m}")
    x, f = a, 1
    while x != 1:
        x = (x * a) % m
        f += 1
        if f > m:
            raise RuntimeError(f"no order for {a} mod {m}")
    return f


def cyclotomic_f(prime: int, m: int) -> int | None:
    """Residue degree of prime in Q(ζ_m).  None if prime | m (ramified)."""
    if m % prime == 0:
        return None
    return mul_order_mod(prime, m)


def valuation(x: Fraction, prime: int) -> int:
    n, d = x.numerator, x.denominator
    v = 0
    while n % prime == 0:
        n //= prime
        v += 1
    while d % prime == 0:
        d //= prime
        v -= 1
    return v


def paley_field_primes(p: int) -> set[int]:
    """Primes that may ramify in Q(√p*) Q(ζ_{q-1}), plus n=p²+1."""
    q = q_of(p)
    n = q + 1
    fac: dict[int, int] = {}
    for v in (p, q - 1, n, 2):
        for r, e in trial_factor(v).items():
            fac[r] = fac.get(r, 0) + e
    return set(fac)


def row_for_lambda(p: int, k: int, lam: Fraction) -> dict:
    q = q_of(p)
    fh = fhat_from_lambda(p, lam)
    fn = trial_factor(fh.numerator)
    fd = trial_factor(fh.denominator)
    pstar = paley_star(p)
    den_primes = sorted(fd)
    cyc_m = q - 1
    cyc = {}
    quad = {}
    odd_block = False
    for r in den_primes:
        fc = cyclotomic_f(r, cyc_m)
        fq = quadratic_f(r, pstar)
        v = valuation(fh, r)
        cyc[str(r)] = fc
        quad[str(r)] = fq
        if fc is not None and v % fc != 0:
            odd_block = True
        if fq is not None and v % fq != 0:
            odd_block = True
    return {
        "p": p,
        "k": k,
        "lambda": str(lam),
        "Fhat": str(fh),
        "Fhat_num_fac": {str(a): b for a, b in sorted(fn.items())},
        "Fhat_den_fac": {str(a): b for a, b in sorted(fd.items())},
        "is_Q_square": is_q_square(fh),
        "pstar": pstar,
        "cyclotomic_m": cyc_m,
        "cyclotomic_f_at_den": cyc,
        "quadratic_f_at_den": quad,
        "odd_val_blocks_norm": odd_block,
        "den_outside_paley_ramification": sorted(
            r for r in den_primes if r not in paley_field_primes(p)
        ),
    }


def all_rows() -> dict[int, list[dict]]:
    out: dict[int, list[dict]] = {}
    for p, tab in LAMBDA.items():
        out[p] = [row_for_lambda(p, k, lam) for k, lam in sorted(tab.items())]
    return out


def theorem_fhat_not_paley_field_norm() -> dict:
    """Proved kill of Paley-field square/norm for F̂.  Sign of F̂ still OPEN."""
    rows = all_rows()
    r5 = rows[5]
    r7 = rows[7]
    min5 = next(r for r in r5 if r["k"] == 4)
    min7 = next(r for r in r7 if r["k"] == 16)
    ok = True
    # p=5: 1250/13, not a square, v_13=-1, f_cyc=2, inert in Q(√5)
    fh5 = Fraction(min5["Fhat"])
    ok = ok and fh5 == Fraction(1250, 13)
    ok = ok and min5["is_Q_square"] is False
    ok = ok and min5["cyclotomic_f_at_den"]["13"] == 2
    ok = ok and min5["quadratic_f_at_den"]["13"] == 2
    ok = ok and min5["odd_val_blocks_norm"] is True
    ok = ok and all(r["odd_val_blocks_norm"] for r in r5)
    ok = ok and all(not r["is_Q_square"] for r in r5)
    # p=7: den 409, f=2, inert in Q(√-7), 409 ∤ p(q-1)n
    fh7 = Fraction(min7["Fhat"])
    q7 = q_of(7)
    ok = ok and fh7 == (Fraction(3072, 409) - 6) * q7 * q7
    ok = ok and min7["is_Q_square"] is False
    ok = ok and min7["cyclotomic_f_at_den"]["409"] == 2
    ok = ok and min7["quadratic_f_at_den"]["409"] == 2
    ok = ok and min7["odd_val_blocks_norm"] is True
    ok = ok and 409 in min7["den_outside_paley_ramification"]
    ok = ok and all(r["odd_val_blocks_norm"] for r in r7)
    ok = ok and all(not r["is_Q_square"] for r in r7)
    # sanity of the Kronecker / order helpers
    ok = ok and kronecker(5, 13) == -1
    ok = ok and kronecker(-7, 409) == -1
    ok = ok and mul_order_mod(13, 24) == 2
    ok = ok and mul_order_mod(25, 48) == 2
    ok = ok and is_prime(13) and is_prime(409)
    return {
        "proved": bool(ok),
        "inequality_proved": False,
        "method_killed": "Fhat = Paley-field square or field-norm",
        "claim_Fhat_is_Q_square": False,
        "claim_Fhat_is_paley_field_norm": False,
        "rows": {str(p): v for p, v in rows.items()},
        "p5_min_Fhat": str(fh5),
        "p7_min_Fhat": str(fh7),
        "theorem": (
            "F̂=(λ−6)q² at the named even characters.  At p=5 every "
            "value has v_13=−1 with f=2 in Q(ζ_{24}) and Q(√5); 1250/13 "
            "is not a square in Q.  At p=7 every value has v_409=−1 "
            "with f=2 in Q(ζ_{48}) and Q(√−7), and 409 does not divide "
            "p(p²−1)(p²+1).  Fail: Q-square; fail: Paley-field norm "
            "(valuations would be multiples of the residue degree).  "
            "F̂≥0 stays OPEN."
        ),
    }


def main() -> dict:
    from io_atomic import write_json_atomic

    T = theorem_fhat_not_paley_field_norm()
    out = {
        "title": "Fhat is not a Paley-field square or field-norm",
        "numbered": False,
        "theorem": T,
        "global_qvar_not_claimed": True,
        "backend": "serial exact Fractions (inherently sequential; GPU unused)",
    }
    path = ROOT / "evidence" / "e1_gmin_qvar_fhat_norm.json"
    write_json_atomic(path, out)
    print("Fhat Paley-field square/norm kill", flush=True)
    print(f"  proved_kill={T['proved']} inequality={T['inequality_proved']}", flush=True)
    print(f"  p5 min Fhat={T['p5_min_Fhat']}", flush=True)
    print(f"  p7 min Fhat={T['p7_min_Fhat']}", flush=True)
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
