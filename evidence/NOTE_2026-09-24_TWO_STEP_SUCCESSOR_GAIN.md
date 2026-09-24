# A second actual update strictly improves the successor envelope

2026-09-24. Positive fixed-time repeated-update result.
The original convergence problem remains OPEN. No independent human or
proof-assistant review is claimed.

The optimized one-successor theorem gives a positive disagreement floor
for the ACTUAL non-Gaussian successor. Instead of stopping after converting
that floor into a deficit bound, apply the already-proved mean-update
inequality once more to the successor law. This extracts a further strict
gain.

## 1. First actual successor

Use the notation of
NOTE_2026-09-24_OPTIMIZED_SUCCESSOR_ENVELOPE.md. For fixed

    99/100 <= t <= 1,
    0 < u <= 1/2,

put

    H=H(t),
    e=e(t), f=f(t),
    A=(1-u)^2 e+u(1-u)f,
    d=2sqrt(u(1-u)/3),
    G=2u d H^(3/2),
    c=1+u^2.

After the first independent partial best-response update, the ACTUAL
successor laws satisfy

    e_1 >= A-u^2 alpha-o_L(1),                             (1)
    g_1=f_1-2e_1 >= G-o_L(1).                              (2)

No Gaussian successor law is used.

The sharp deficit closure gave the previous one-step envelope L(t,u),
characterized by

    [c L-A][L+G/2]=G^2/8.                                 (3)

## 2. Apply a second mean update to the actual successor

Choose a second fixed update probability

    0 <= v <= 1/2.

The cap-free mean-update inequality applies to arbitrary Boolean source
laws. Dropping only its nonnegative interior correction gives

    (1+v^2)alpha
      >=(1-v)^2 e_1+v(1-v)f_1-o_L(1)
       =(1-v^2)e_1+v(1-v)g_1-o_L(1).                     (4)

Insert (1)--(2). After n->infinity at fixed initial source cap,

    C(u,v) alpha >= D(t,u,v),                              (5)

where

    C(u,v)=1+v^2+(1-v^2)u^2
          =c+(1-u^2)v^2,                                  (6)

    D(t,u,v)=(1-v^2)A+v(1-v)G
            =A+Gv-(A+G)v^2.                               (7)

Therefore

    liminf alpha_n >= T(t,u,v):=D(t,u,v)/C(u,v).          (8)

The same-order regularization step removes the initial fixed cap exactly as
in the preceding source and successor theorems.

Define the exact two-update constant

    B_2 =
      sup_{99/100<=t<=1, 0<u<=1/2, 0<=v<=1/2}
          T(t,u,v).                                        (9)

Then liminf alpha_n>=B_2.

## 3. The second update is STRICTLY stronger than B_opt

The one-step function L(t,u) is continuous after its natural extension at
u=0, so on the compact parameter rectangle it attains B_opt. Since
B_opt>B_noise while u=0 gives only the pre-update source energy, a maximizing
pair (t_*,u_*) has u_*>0.

Fix such a pair and abbreviate

    x=L(t_*,u_*)=B_opt,
    A=A(t_*,u_*),
    G=G(t_*,u_*),
    c=1+u_*^2,
    X=cx-A>0.

Equation (3) gives

    X=G^2/[8(x+G/2)]
     =G^2/[4(2x+G)].                                      (10)

At the candidate alpha=x, subtract the right side of (5) from the left:

    C(u_*,v)x-D(t_*,u_*,v)
      =X-Gv+Kv^2,                                         (11)

with

    K=(1-u_*^2)x+A+G
     =2x+G-X.                                             (12)

Because A>0 and c<=5/4,

    X=cx-A<cx<=5x/4<2x,

so K>G. Hence

    v_*=G/(2K)

lies strictly between 0 and 1/2. Evaluating (11) there gives

    X-G^2/(4K).                                           (13)

But K=2x+G-X<2x+G, so by (10),

    G^2/(4K) > G^2/[4(2x+G)] = X.                         (14)

Thus (13) is strictly negative. Equivalently,

    T(t_*,u_*,v_*) > x=B_opt.                             (15)

Therefore

    B_2 > B_opt > B_noise > B_both > B_int > B_tilt.     (16)

The strict improvement is fully analytic. It does not rely on the
numerical optimizer or on the successor eta floor.

## 4. Numerical locator only

A scalar evaluation of (8), for orientation only, places the maximizing
triple near

    t ~= 0.9934102,
    u ~= 0.0974342,
    v ~= 0.0046322,
    T ~= 0.3258669876.

The theorem is the exact supremum (9), not this decimal.

## 5. What changed

This is the first use in the current September 24 chain where the positive
ACTUAL successor disagreement is fed into a further update and shown to
improve the unconditional constant strictly. It is still a fixed two-step
argument. It does not establish a nondecaying long-time gain, a summable
one-vertex excess, or convergence of alpha_n.
