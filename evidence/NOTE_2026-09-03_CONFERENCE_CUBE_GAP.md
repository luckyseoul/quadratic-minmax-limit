# Proposition 15.762: the integral conference cube gap

**Status:** proved for every symmetric conference matrix of order
`p^2+1`, with odd `p>=5`. This is a universal gap theorem and a sharp
counterexample-search criterion. It does not construct a nonregularizable
conference matrix and does not close residual (ii), E1, or the original
MathOverflow limit.

Let

\[
 C=C^t,\qquad C^2=p^2I,\qquad n=p^2+1,
 \qquad Q_C(x)={1\over2}x^tCx,
 \]

and let `x` be Boolean. Write

\[
 \Phi_{\rm sph}={pn\over2},\qquad
 \delta=\Phi_{\rm sph}-Q_C(x),\qquad
 z={C-pI\over2}x.                                      \tag{1}
\]

Every coordinate of `Cx` and `px` is odd, so `z` is integral. The
conference equation gives

\[
 Cz=-pz,\qquad \|z\|^2=p\delta.                         \tag{2}
\]

The second identity is the normalization guard that matters here:
`||z||^2=2p` means `Q=Phi_sph-2`, while `Q=Phi_sph-4` means
`||z||^2=4p`.

## The switched degree-parity lemma

Put `X=diag(x)`, `D=XCX`, and `w=Xz`. Then

\[
 D\mathbf1=p\mathbf1+2w,\qquad Dw=-pw,
 \qquad \sum_iw_i=-\delta,qquad \|w\|^2=p\delta.       \tag{3}
\]

Let `A=(J-I-D)/2` and let `d_i` be its degrees. For distinct `i,j`, the
off-diagonal conference equation is

\[
 0=(D^2)_{ij}=n-2-2(d_i+d_j)+4A_{ij}+4(A^2)_{ij}.       \tag{4}
\]

Since `n-2=p^2-1` is divisible by eight, (4) modulo four says that
`d_i+d_j` is even for every pair. Thus all degrees have one parity. From

\[
 d_i={p^2-p\over2}-w_i,                                 \tag{5}
\]

all coordinates of `w` also have one parity.

## Excluding the three gaps below eight

Both `Q_C(x)` and `Phi_sph` are odd, so `delta` is even.

If `w` is even, write `w=2v`; then `Dv=-pv`. At a coordinate of maximum
absolute value `m`,

\[
 pm\le\sum_{j\ne i}|v_j|=\|v\|_1-m,
 \qquad\text{so}\qquad \|v\|_1\ge(p+1)m.               \tag{6}
\]

For integral `v`, `||v||_1<=||v||^2`, and hence every nonzero integral
`-p` eigenvector has squared norm at least `p+1`. Also (3) makes `delta`
divisible by four in this branch. Therefore `delta=4` would give
`||v||^2=p`, contradicting (6), while `delta=2,6` have the wrong residue
modulo four.

If `w` is odd, then `||w||^2>=n`. This excludes `delta=2` for every odd
`p>=3`, and excludes `delta=6` for `p>=7`. At `p=5`, a sum of 26 odd
squares is `2 mod 8`, whereas `6p=30` is `6 mod 8`. Finally, an odd `w`
at `delta=4` is already impossible modulo four.

Thus, for every odd `p>=5`,

\[
 \boxed{\quad Cx=px\quad\hbox{or}\quad
 Q_C(x)\le {pn\over2}-8.\quad}                          \tag{7}
\]

Applying (7) to `-C` gives the two-sided criterion

\[
 \boxed{\quad \text{if }C\text{ has no Boolean }\pm p
 \text{ eigenvector, then }\Phi(C)\le {pn\over2}-8.\quad} \tag{8}
\]

No optimization at `Phi_sph-2`, `Phi_sph-4`, or `Phi_sph-6` is needed.
Absence of both Boolean eigenshells is already the complete certificate.

## Necessary form of the first possible shell

At `delta=8`, odd parity is impossible modulo eight, so `w=2v` and

\[
 Dv=-pv,\qquad \|v\|^2=2p,\qquad \sum_iv_i=-4.       \tag{9}
\]

If `||v||_infinity>=2`, (6) gives `||v||_1>=2p+2`, while integrality and
(9) give `||v||_1<=||v||^2=2p`, a contradiction. Hence
`v_i in {0,+-1}`, its support has size `2p`, and its sign counts are

\[
 \#\{v_i=1\}=p-2,\qquad \#\{v_i=-1\}=p+2.              \tag{10}
\]

Equations (9)--(10) are necessary, not sufficient. In particular this
proposition does not assert that the gap-eight shell is attained.

## Scope for residual (ii)

Proposition 15.762 is independent of the post-15.761 edge--Radon Boolean
fibre. It removes the proposed near-spectral MIQCP subproblem from any
conference-matrix counterexample route: a new candidate only needs exact
tests for Boolean `+p` and `-p` eigenvectors. None of the currently audited
order-122 Paley, Peisert, PN, or OA-derived families supplies a candidate
without such an eigenshell. No common residual graph or counterexample has
therefore been constructed. Residual (ii), E1, `L=1/2`, and the original
MathOverflow limit remain **OPEN**.
