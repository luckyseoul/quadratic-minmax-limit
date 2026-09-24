"""Focused remote check of the new two-phase scalar argument.

Requires SymPy. Run serially with: python3 check.py
Exact identities and rational certificates are symbolic; decimal evaluations
use 60-digit working precision, not certified interval arithmetic. This does
not verify Gaussianization, cap removal, or convergence of alpha_n.
"""

try:
    import sympy as sp
except ImportError as exc:
    raise SystemExit("SymPy is required; run this on a host with SymPy installed.") from exc


checks = 0


def check(label, condition):
    global checks
    if not bool(condition):
        raise AssertionError(label)
    checks += 1
    print(f"PASS: {label}", flush=True)


q, t, kappa, aq = sp.symbols("q t kappa aq", positive=True)
r = sp.symbols("r", real=True)
d = 1 + t**2 * q
b = kappa * t**2 / d
for s in (-1, 1):
    cov = sp.sqrt(kappa) * (2*t*q + s*t**2*r) / d
    variance_lower = aq + b*r**2/q + s*2*kappa*t*r/d
    residual_lower = aq - kappa*q + kappa*(t*r + s*q*(1-t**2*q))**2/(q*d**2)
    check(f"exact residual identity, phase {s:+d}",
          sp.factor(variance_lower - cov**2 - residual_lower) == 0)

cov_plus = sp.sqrt(kappa) * (2*t*q + t**2*r) / d
cov_minus = sp.sqrt(kappa) * (2*t*q - t**2*r) / d
check("exact mean input/field covariance",
      sp.simplify((cov_plus + cov_minus)/2 - 2*sp.sqrt(kappa)*t*q/d) == 0)

x, beta = sp.symbols("x beta", positive=True)
phi = sp.atan(sp.sqrt(beta)/x)/sp.pi
check("exact convexity formula",
      sp.simplify(sp.diff(phi, x, 2) - 2*sp.sqrt(beta)*x/(sp.pi*(beta+x**2)**2)) == 0)
check("angle at zero covariance", sp.limit(phi, x, 0, dir="+") == sp.Rational(1, 2))
check("C squared at most kappa",
      sp.factor(kappa - 4*kappa*t**2/(1+t**2)**2
                - kappa*(1-t**2)**2/(1+t**2)**2) == 0)

T = sp.Rational(993, 1000)
Z = T**2/(1+T**2)
check("z exceeds 49/100", Z > sp.Rational(49, 100))
check("rational residual-floor boundary is 1/75",
      (sp.Rational(2, 3)-sp.Rational(16, 25))/2 == sp.Rational(1, 75))
rb = sp.Rational(16, 25)/(sp.Rational(2, 3)+sp.Rational(7, 11)*sp.Rational(49, 100))
check("exact rational angle bound", rb == sp.Rational(2112, 3229))
check("rational angle comparison", rb < sp.Rational(1309, 2000))
check("strict rational lower bound for sqrt(5)", sp.Rational(559, 250)**2 < 5)
check("cosine-square certificate link",
      (3+sp.Rational(559, 250))/8 == sp.Rational(1309, 2000))
check("cubing the changed fraction", sp.Rational(1, 5)**3 == sp.Rational(1, 125))
check("eightfold eta improvement", sp.Rational(1, 125)/sp.Rational(1, 1000) == 8)

B, p = sp.symbols("B p", positive=True)
new_bound = (B + sp.sqrt(B**2 + p**2/(750*(1+p**2))))/2
check("exact new lower-bound quadratic",
      sp.simplify(new_bound**2-B*new_bound-p**2/(3000*(1+p**2))) == 0)

K = 2/sp.pi
A = 1-K*sp.asin(Z)+K*Z
E = K*T/(1+T**2)
F = sp.sqrt(K*A)
C = 2*sp.sqrt(K)*T/(1+T**2)
H = sp.Min(sp.Rational(1, 4), sp.atan(sp.sqrt(A-K)/C)/sp.pi)
P = (F-2*E)/(F+sp.sqrt(F**2+(F-2*E)**2))
B0 = E-F/2+sp.sqrt(F**2+(F-2*E)**2)/2
B1 = (B0+sp.sqrt(B0**2+P**2/(6000*(1+P**2))))/2
B2 = new_bound.subs({B: B0, p: P})
check("numerical paired fraction exceeds 1/5", sp.N(H, 60) > sp.Rational(1, 5))
check("numerical update weight admissible", 0 < sp.N(P, 60) < sp.Rational(1, 2))
check("numerical strict lower-bound gain", sp.N(B2-B1, 60) > 0)

for label, expression in (("paired_fraction_floor", H), ("p", P),
                          ("B_tilt", B0), ("B_int", B1),
                          ("B_both", B2), ("gain_over_B_int", B2-B1)):
    print(f"{label} = {sp.N(expression, 60)}", flush=True)

print(f"ALL {checks} CHECKS PASSED; original convergence OPEN", flush=True)
