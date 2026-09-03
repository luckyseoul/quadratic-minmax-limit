# Independent audit of the external exact `m_11=17` certificate

**Status:** the mandatory witness and exhaustive lower-bound replay passed.
This is a new finite datum for the calculator plot, not evidence that the
normalized sequence converges.

Source:
[CurtisAccelerate/antipodal-cut-code-k11](https://github.com/CurtisAccelerate/antipodal-cut-code-k11),
commit `977927bd16188eda7334b040601e9db6a321aaf4`.

The claimed theorem is

\[
 m_{11}=17,
 \qquad
 \rho(C_{\rm cut}(K_{11})+\langle\mathbf1\rangle)=19.
\]

The source was cloned fresh and its mandatory runner was executed unchanged:

```bash
git clone https://github.com/CurtisAccelerate/antipodal-cut-code-k11.git /tmp/antipodal-cut-code-k11
cd /tmp/antipodal-cut-code-k11
bash reproducibility/run_all.sh
```

The upper-bound verifier enumerated all `2^10=1024` projective spin states.
It found minimum `-15`, maximum `17`, maximum absolute value `17`, and 24
absolute extremizers.  The complete value distribution was

```text
-15: 29, -11: 95, -7: 142, -3: 206, 1: 209,
  5: 163,   9: 92, 13: 64, 17: 24.
```

The exact C++ lower-bound program then reproduced the advertised terminal
counter:

```text
TOTAL nodes=50778686 complete=936720 pass4=936720 passall=0
wall=48.73 sec maxrss=3648 KB
```

It switches one vertex positive, exhausts the 42 admissible sorted degree
sequences with `20<=e(H)<=22`, and checks all 1,024 subset inequalities on
every graph surviving its exact four-set pruning.  `passall=0` excludes
`Phi<=15`; parity of 55 signed terms then gives `m_11>=17`.  Together with
the witness this proves equality.

Audit hashes:

```text
15fa551ce2953a5a8bb350648879c9d478eb8a5f3def25bffdc32f9e86ec7299  m11_witness.json
da233ba73182f5738ef00d01e07b688c93ea25a5282e45bb579dafcf51025498  verify_witness.py
b98f8e55548b929aa45a89e84568e5e479bb9ef93dd76a7de1de519dac8d4c3a  lower_bound_enumerator.cpp
794763038eb40c834b5b8710b3425bb32958ca4d6bee89fb0fa6b97fd2e65541  witness_verification.log
478a88bcdde94c5d5f26a0561e128aae26b53540b76470eb21560e94df946a78  lower_bound_stdout.log
3f4c46fddda3119449efa492b4b24e1d12328cfc144f1f079a343f6e98838aa5  lower_bound_stderr.log
```

One packaging discrepancy was observed: the upstream README lists bundled
`logs/` and `SHA256SUMS`, but neither path exists at the audited commit.  The
runner creates the logs during replay, so this does not affect the theorem;
it does mean those advertised precomputed provenance files were unavailable.
