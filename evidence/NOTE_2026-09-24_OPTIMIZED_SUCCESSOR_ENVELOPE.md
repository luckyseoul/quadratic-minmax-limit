# Optimizing the actual-successor lower-bound envelope

2026-09-24. Positive strengthening of the September 24 lower bound.
The original convergence problem remains OPEN. No independent human or
proof-assistant review is claimed.

The previous successor-noise theorem fixed two parameters for historical
reasons: the Gaussian tilt was held at t=993/1000, and the partial-flip
probability was kept equal to the optimizer of the PRE-successor mean-update
bound. Neither restriction is needed by the actual-successor proof.

This note retains the explicit both-phase disagreement angle instead of
rounding it down to 1/5, then frees the update probability. The result is
a strictly stronger unconditional lower bound.

## 1. A two-parameter source family

Let

    kappa = 2/pi,
    99/100 <= t <= 1,
    z(t) = t^2/(1+t^2),
    a(t) = 1-kappa asin(z(t))+kappa z(t),
    e(t) = kappa t/(1+t^2),
    f(t) = sqrt(kappa a(t)),
    beta(t) = a(t)-kappa,
    C(t) = 2 sqrt(kappa)t/(1+t^2).

The tilted paired Gaussian construction and the both-phase residual identity
are valid for every fixed t in this interval. Define

    H(t)=min{1/4, atan(sqrt(beta(t))/C(t))/pi}.             (1)

Then the phase-averaged initial changed-coordinate fraction satisfies

    mu_0 >= H(t)-o_L(1).                                   (2)

Moreover H(t)>1/5 throughout 99/100<=t<=1. Indeed,
z(t)>=9801/19801>49/100, while z(t)<=1/2. The rational certificate already
used in NOTE_2026-09-24_BOTH_PHASE_DISAGREEMENT therefore applies unchanged:

    C(t)^2/(C(t)^2+beta(t))
      <= kappa/a(t)
      < (16/25)/[2/3+(7/11)(49/100)]
      = 2112/3229
      < 1309/2000
      < cos(pi/5)^2.

Hence the angle branch in (1) is >1/5, and the other branch is 1/4.

No uniform-in-t Gaussianization statement is needed. The conclusion is
applied at one fixed t at a time; the final supremum is taken only after
the fixed-parameter theorem is established.

## 2. Free the actual partial-flip probability

Let 0<u<=1/2 be the probability used for the ACTUAL independent partial
best-response update. Put

    d(u)=2 sqrt(u(1-u)/3),
    G(t,u)=2u d(u) H(t)^(3/2),                              (3)
    A(t,u)=(1-u)^2 e(t)+u(1-u) f(t),                       (4)
    c(u)=1+u^2.

The successor-noise theorem and (2) give

    g_1 >= G(t,u)-o_L(1),                                  (5)
    e_1 >= A(t,u)-u^2 alpha-o_L(1).                        (6)

The sharp pointwise deficit-disagreement closure, averaged over the two
actual successor laws, gives for D_1=alpha-e_1

    D_1 >= g_1^2/(8 alpha+4 g_1).                          (7)

Combining (5)--(7) and taking n->infinity at fixed source cap yields

    [c(u) alpha-A(t,u)] [alpha+G(t,u)/2] >= G(t,u)^2/8.    (8)

The positive root is

    L(t,u)=
      { A(t,u)-c(u)G(t,u)/2
        +sqrt([A(t,u)+c(u)G(t,u)/2]^2
              +c(u)G(t,u)^2/2) }
      /[2c(u)].                                            (9)

Same-order regularization removes the initial fixed cap exactly as in the
preceding tilted/source/successor proofs. Therefore

    liminf alpha_n >= L(t,u)                               (10)

for every fixed pair in the displayed parameter rectangle.

Consequently the exact optimized constant

    B_opt =
      sup_{99/100<=t<=1, 0<u<=1/2} L(t,u)                  (11)

is an unconditional lower bound for the original minimax problem.

## 3. Strict improvement over B_noise without numerical optimization

Let

    t0=993/1000

and let u0 be the old optimizer of the PRE-successor mean-update ratio

    R_t0(u)=A(t0,u)/c(u).

Thus R_t0(u0)=B_tilt and, because u0 is an interior optimizer,

    A_u(t0,u0)=c'(u0) B_tilt.                              (12)

The previous B_noise used only the coarse floor H>=1/5, so its
successor gap was

    G0(u)=2u d(u)/5^(3/2).

But H(t0)>1/5 strictly. For fixed A,c, the right side of

    c alpha-A >= G^2/(8alpha+4G)

is strictly increasing in G>0. Hence the exact-angle value already gives

    L(t0,u0) > B_noise.                                    (13)

There is a second strict gain: u0 is not stationary for the new
successor-corrected envelope.

Fix t=t0, abbreviate G=G(t0,u), c=c(u), A=A(t0,u), and let x=L(t0,u).
The defining equation is

    F(x,u)=(cx-A)(x+G/2)-G^2/8=0.                          (14)

At u=u0 write B=B_tilt, X=c(x-B)>0, Y=x+G/2. Since XY=G^2/8,
differentiating (14) at fixed x and using (12) gives

    F_u =
      (c'/c) G^2/8
      - G G'(4x+G)/(16Y).                                 (15)

For 0<u<=1/2,

    G'/G=(3-4u)/(2u(1-u)) >= 2,
    c'/c=2u/(1+u^2) <= 4/5.                               (16)

Also 4x+G>2Y. Therefore F_u<0. Since

    F_x=cY+X>0,

implicit differentiation gives

    dL(t0,u)/du |_(u=u0) = -F_u/F_x > 0.                  (17)

Thus some u>u0 gives a strictly larger lower bound than L(t0,u0).
Combining (13) and (17),

    B_opt > L(t0,u0) > B_noise
          > B_both > B_int > B_tilt.                      (18)

This strict chain is analytic and does not rely on a numerical optimizer.

## 4. Numerical locator only

A high-precision scalar evaluation of (9), for orientation only, places the
two-parameter maximum near

    t ~= 0.9934103,
    u ~= 0.0974340,
    L ~= 0.3258669873.

These decimals are NOT the theorem and are not used in (18). The canonical
new constant is the exact supremum (11).

## 5. Scope

This advances the unconditional all-orders asymptotic lower bound. It does
not prove convergence, identify a limit, or supply a cross-order upper
comparison. The successor law remains the actual non-Gaussian Boolean law;
no Gaussian approximation is reapplied after the update.

The only new optimization is over two fixed scalar parameters already
admissible in the proved source and successor theorems. No family search,
finite signing census, or retired route is reopened.
