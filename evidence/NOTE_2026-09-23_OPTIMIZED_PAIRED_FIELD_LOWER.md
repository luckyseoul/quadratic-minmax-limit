# Optimizing the paired-field mean update

2026-09-23. Positive quantitative improvement of the unconditional lower
bound. This uses the already-proved paired-field estimates and the same-order
regularization theorem; no new route premise is introduced. The original
convergence problem remains OPEN.

Put

    kappa = 2/pi,
    e0 = kappa/2,
    f0 = sqrt(kappa*(2/3+kappa/2)).

The September 17 paired-field theorem proves, after taking n->infinity at
each fixed operator cap and then removing the cap by same-order
regularization, that for every fixed p in [0,1],

    liminf alpha_n >= R(p),

where

    R(p) = [(1-p)^2 e0 + p(1-p) f0]/(1+p^2).              (1)

The previous published constant used p=1/10. There is no reason to freeze p.

Write a=f0-2e0>0. Expanding (1),

    R(p) = [e0 + a p + (e0-f0)p^2]/(1+p^2).

This is the Rayleigh quotient of

        [ e0       a/2   ]
        [ a/2    e0-f0   ]

at the vector (1,p). Therefore the maximum over real p is its largest
eigenvalue

    Bopt = e0 - f0/2 + (1/2)*sqrt(f0^2 + (f0-2e0)^2).     (2)

The maximizing p is

    pstar = (f0-2e0)/(f0 + sqrt(f0^2+(f0-2e0)^2))
          = 0.097102496812...                              (3)

and lies in [0,1], so it is admissible in the original mean-update
inequality. Numerically,

    Bopt = 0.3258474004377944...
    Bold = 0.3258407554555742...
    Bopt-Bold = 0.0000066449822202...

Thus the unconditional global bound improves to

    liminf_(n->infinity) alpha_n >= 0.3258474004377944...  (4)

The cap-removal argument is unchanged because pstar is a fixed scalar and
the limiting e0,f0 constants are independent of the fixed source cap.
No finite-order cutoff, convergence statement, or limit value is claimed.

Direct differentiation gives the same optimizer. The numerator of R'(p)
after multiplication by (1+p^2)^2 is

    (2e0-f0)p^2 - 2f0 p + (f0-2e0),

equivalently

    (f0-2e0)p^2 + 2f0 p - (f0-2e0)=0,

whose positive root is (3).
