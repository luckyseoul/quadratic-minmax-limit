"""Focused remote scalar checks for the actual successor flip-noise proof.

Requires SymPy; run serially with: python3 check.py
Exact algebra and rational comparisons are separate from 60-digit numerical
evaluations. The latter are not certified interval arithmetic. This does
not formally verify probabilistic implications, cap removal, or convergence.
No matrix census or simulation is run.
"""

try:
    import sympy as sp
except ImportError as exc:
    raise SystemExit("SymPy is required; use a host with SymPy installed.") from exc


checks = 0


def check(label, condition):
    global checks
    if not bool(condition):
        raise AssertionError(label)
    checks += 1
    print(f"PASS: {label}", flush=True)


p, b, u, a, variance = sp.symbols("p b u a variance", positive=True)
for power, expected in (
    (2, p*(1-p)),
    (3, p*(1-p)*(1-2*p)),
    (4, p*(1-p)*(1-3*p*(1-p))),
):
    moment = p*(1-p)**power + (1-p)*(-p)**power
    check(f"centered Bernoulli moment {power}",
          sp.factor(moment-expected) == 0)

v = a**2+variance
check("shifted fourth-moment quadratic remainder",
      sp.expand(3*v**2-(a**4+6*a**2*variance+3*variance**2)) == 2*a**4)
check("positive square majorant remainder",
      sp.expand((sp.sqrt(3)*u**2+2*b*u)**2
                -(3*u**4+4*b*u**3+b**2*u**2)
                -(4*(sp.sqrt(3)-1)*b*u**3+3*b**2*u**2)) == 0)
check("Holder lower-floor remainder",
      sp.simplify(u**2/(sp.sqrt(3)*u+2*b)
                  -(u/sp.sqrt(3)-2*b/3)
                  -4*b**2/(3*(sp.sqrt(3)*u+2*b))) == 0)

dp = 2*sp.sqrt(p*(1-p)/3)
check("successor field coefficient squared", sp.expand(dp**2-4*p*(1-p)/3) == 0)
G_expr = 2*p*dp/sp.sqrt(125)
G_squared = 16*p**3*(1-p)/375
check("source-seeded gap squared", sp.simplify(G_expr**2-G_squared) == 0)
check("single-bad-coordinate endpoint", 1-sp.Rational(1, 2)-sp.Rational(1, 2)**2 > 0)

e, g, alpha = sp.symbols("e g alpha", real=True)
check("energy recursion normalization",
      sp.expand((1-p)**2*e+p*(1-p)*(2*e+g)-p**2*alpha
                -((1-p**2)*e+p*(1-p)*g-p**2*alpha)) == 0)

alpha, g = sp.symbols("alpha g", positive=True)
D = alpha+g/2-sp.sqrt(alpha**2+alpha*g)
check("sharp deficit inversion polynomial",
      sp.simplify(D**2-(g+2*alpha)*D+g**2/4) == 0)
check("sharp deficit rationalization",
      sp.simplify(4*D*(alpha+g/2+sp.sqrt(alpha**2+alpha*g))-g**2) == 0)
check("relaxed deficit floor increases with gap",
      sp.factor(sp.diff(g**2/(8*alpha+4*g), g)
                -g*(4*alpha+g)/(4*(2*alpha+g)**2)) == 0)

B, G, x = sp.symbols("B G x", positive=True)
noise_root = (B-G/2+sp.sqrt((B+G/2)**2+G**2/(2*(1+p**2))))/2
check("exact noise-root quadratic",
      sp.simplify((noise_root-B)*(noise_root+G/2)-G**2/(8*(1+p**2))) == 0)
comparison = p**2*(8*x+4*G)/(3000*x)-G_squared
check("comparison with interior-only bound",
      sp.factor(comparison-p**2/(750*x)*(G-(32*p*(1-p)-2)*x)) == 0)

rb = sp.Rational(2112, 3229)
check("exact source comparison for p above 1/12", rb < sp.Rational(119, 143)**2)
check("update-weight inverse at 1/12",
      2*sp.Rational(1, 12)/(1-sp.Rational(1, 12)**2) == sp.Rational(24, 143))
check("source ratio threshold complement", 1-sp.Rational(119, 143) == sp.Rational(24, 143))
check("strict-comparison coefficient at 1/12",
      32*sp.Rational(1, 12)*sp.Rational(11, 12)-2 == sp.Rational(4, 9))
check("source energy exceeds 1/4 by rational bounds",
      sp.Rational(7, 11)*sp.Rational(993, 1000)/2 > sp.Rational(1, 4))
check("gap-square upper-bound decomposition",
      sp.factor(sp.Rational(1, 16)-p**3*(1-p)
                -(sp.Rational(1, 4)-p**2)/4
                -p**2*(sp.Rational(1, 4)-p*(1-p))) == 0)
check("G squared conservative ceiling", sp.Rational(16, 375)/16 == sp.Rational(1, 375))
check("strict noise-versus-interior certificate", sp.Rational(1, 375) < sp.Rational(1, 81))
check("iteration: count factor cubed", (p**2)**3 == p**6)

T = sp.Rational(993, 1000)
K = 2/sp.pi
Z = T**2/(1+T**2)
A = 1-K*sp.asin(Z)+K*Z
E = K*T/(1+T**2)
F = sp.sqrt(K*A)
P = (F-2*E)/(F+sp.sqrt(F**2+(F-2*E)**2))
B0 = E-F/2+sp.sqrt(F**2+(F-2*E)**2)/2
B1 = (B0+sp.sqrt(B0**2+P**2/(6000*(1+P**2))))/2
B2 = (B0+sp.sqrt(B0**2+P**2/(750*(1+P**2))))/2
GG = sp.sqrt(G_squared.subs(p, P))
B3 = noise_root.subs({B: B0, G: GG, p: P})
check("numerical admissible weight exceeds 1/12",
      sp.Rational(1, 12) < sp.N(P, 60) < sp.Rational(1, 2))
check("numerical noise bound exceeds both-phase bound", sp.N(B3-B2, 60) > 0)

for label, expression in (
    ("p", P), ("actual_successor_gap_floor_G", GG),
    ("B_tilt", B0), ("B_int", B1), ("B_both", B2),
    ("B_noise", B3), ("gain_over_B_both", B3-B2),
    ("gain_over_B_int", B3-B1),
):
    print(f"{label} = {sp.N(expression, 60)}", flush=True)

print(f"ALL {checks} CHECKS PASSED; original convergence OPEN", flush=True)
