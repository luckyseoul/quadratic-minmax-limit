# The scalar axioms allow the whole sandwich to oscillate

**Status.** Proved negative statement. The original limit stays OPEN.
Nothing here is a sequence of Seidel minima, and nothing was entered in
`src/original_mo_status.py`.

## 1. Statement

Let \(0<\lambda\le\Lambda\le 1/2\). There is an integer sequence \(m_n\) such that

- \(m_n\) is nondecreasing,
- \(m_{n+1}-m_n\le n-2\) for every large \(n\),
- \(\liminf m_n/n^{3/2}=\lambda\) and \(\limsup m_n/n^{3/2}=\Lambda\).

CORE §4 gives \(\liminf\alpha_n\ge 1/\pi\), and conference matrices give
\(\limsup\alpha_n\le 1/2\). CORE §2 gives the step ceiling \(n-2\) and
monotonicity of \(m_n\). Taking \(\lambda=1/\pi\) and \(\Lambda=1/2\) shows
that these axioms allow
\[
\limsup\alpha_n-\liminf\alpha_n=\tfrac12-\tfrac1\pi.
\]
A proof that the real sequence converges has to use something else.

`solution.md` §9 and §16 already say that some log-periodic integer
sequences fail to converge. Proposition 3.4 of
`evidence/NOTE_2026-10-11_FRAME_OBSTRUCTION_SINGLE_ROW_EXTENSION.md`
gives one example of amplitude \(0.03\) inside \([0.39,0.45]\). The
statement above is the sharp form: every pair of ends inside
\((0,1/2]\) occurs.

## 2. Proof

Write \(a=(\Lambda+\lambda)/2\) and \(b=(\Lambda-\lambda)/2\). When \(b>0\)
choose
\[
\omega=\min\bigl(1,\tfrac{3\lambda}{4b}\bigr),
\]
and when \(b=0\) take \(\omega=1\). Put \(g(t)=a+b\sin(\omega t)\) and
\(h(n)=n^{3/2}g(\ln n)\). Then
\[
h'(n)=n^{1/2}\bigl(\tfrac32 g(\ln n)+g'(\ln n)\bigr),
\]
and \(g'=b\omega\cos(\omega t)\), so
\[
c:=\tfrac32\lambda-b\omega\le\tfrac32 g+g'\le\tfrac32\Lambda+b\omega=:C.
\]
The choice of \(\omega\) gives \(b\omega\le\tfrac34\lambda\), hence
\(c\ge\tfrac34\lambda>0\). Also \(C\le\tfrac34+\tfrac34\lambda\le\tfrac98\).

Thus \(h(n+1)-h(n)\ge c\sqrt n\) and \(h(n+1)-h(n)\le C\sqrt{n+1}\).
Set \(m_n=\lfloor h(n)+1/2\rfloor\). If \(h\) increases by at least \(1\),
then \(m\) increases by at least \(1\). If \(h\) increases by at most \(D\),
then \(m\) increases by at most \(D+1\). Choose \(N\) so that
\(c\sqrt N\ge 1\) and \(C\sqrt{N+1}+1\le N-2\). Such an \(N\) exists
because the first bound grows like \(\sqrt N\) and the second is eventually
linear. For \(n\ge N\), \(m\) is nondecreasing and
\(m_{n+1}-m_n\le n-2\).

The steps of \(\omega\ln n\) tend to \(0\) and the sum diverges, so
\(\omega\ln n\bmod 2\pi\) comes arbitrarily close to every angle.
Therefore \(g(\ln n)\) has liminf \(a-b=\lambda\) and limsup \(a+b=\Lambda\).
Rounding changes \(\alpha_n\) by \(O(n^{-3/2})\), so \(m_n/n^{3/2}\) has
the same ends.

## 3. Check

`scripts/scalar_axiom_full_gap.py` builds the sequence for
\(\lambda=1/\pi\), \(\Lambda=1/2\). The derivative margin is positive,
the oscillation width equals \(1/2-1/\pi\), and on \(200000\) consecutive
orders from the first admissible index the integer sequence is strictly
increasing, respects \(n-2\), and comes within \(0.01\) of both ends.
`tests/test_scalar_axiom_full_gap.py` calls those functions.

## 4. What this does not say

The sequence is not realized by \(\min\Phi\). It does not prove that
\(\alpha_n\) diverges. It records that monotonicity, the step ceiling,
and the sandwich are compatible with every gap up to \(1/2-1/\pi\).
