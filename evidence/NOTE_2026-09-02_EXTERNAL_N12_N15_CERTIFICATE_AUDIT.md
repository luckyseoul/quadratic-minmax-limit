# External exact values through order fifteen: provenance audit

**Status:** the public finite certificates support

\[
m_{12}=18,\qquad m_{13}=20,\qquad m_{14}=21,\qquad m_{15}=27.
\]

These are existing finite results, not new searches and not evidence of
convergence.  They were imported only to prevent the calculator plot and
future work from duplicating already completed exact-order calculations.

Source: [Robby955/mo-413935-research](https://github.com/Robby955/mo-413935-research),
commit `61ac2268b6f9234c3ea268e7a7072ac72141f36b`.

## Local checks

The repository was cloned fresh at the pinned commit.  The order-fifteen
default certifier completed with exit zero:

```bash
PYTHONDONTWRITEBYTECODE=1 python verification/research_order15_certify.py
```

It independently checked the `M=27` witness, the four threshold-23 extension
minima `27,29,29,29`, the two tower-receipt hashes, and the decompressed
1,313,164-record catalogue hash and count.  Its explicit trust boundary is
the completeness of the published nauty/labelg threshold tower; neither the
full top-level extension replay nor the roughly billion-record lower level
was regenerated here.

For orders thirteen and fourteen, the two published order-twelve threshold
survivors were evaluated directly in Python.  Both reproduced

```text
M=18, absolute maximizers=20, extension minimum=24,
optimal extension centers=772, deletion maxima=17x12.
```

The Paley conference witness reproduced `M(C14)=21`, and all fourteen
order-thirteen principal submatrices reproduced maximum `20`.  The lower
bound `m_13>=20` still trusts the published eight-shard nauty completeness
receipt over 1,018,997,864 residual graphs.  Given that result, monotonicity,
energy parity, and the checked witnesses give `m_14=21`.

The supplied C11 quick scanner did not compile unchanged under the local
GCC 14/15 toolchains because `-Wpedantic -Werror` promotes a pre-C23 array
qualifier diagnostic to an error.  This is an environmental portability
failure, not a passed check and not evidence against the certificate.  The
independent pure-Python survivor and conference evaluations above avoided
altering the upstream source.

Order twelve has a smaller independent deduction.  The separately replayed
certificate gives `m_11=17`; monotonicity makes `m_12>=17`, order-twelve
energies are even, and the checked order-twelve witness has maximum `18`.
Therefore `m_12=18` without trusting a second lower-bound enumeration.

## Local artifact hashes

```text
05f9ca599311f57e77f215d0a434eecc8484b77c13d6d6e75eee297e69f730b5  order13_full_receipt.txt
e99a33f27430eb55dd96fc2d1e9f3bbf17e786c81658ca2ab8faf465e20e7ae9  order15_full_receipt.txt
2eaa0bca7ce14a39188fec0d9979676c2978c7d3e5c638b93f3fbd84e2503624  order15_level14_catalogue.g6.gz
```

The exact-value table records the distinction between locally replayed
mathematics and externally trusted stream completeness.  No order-sixteen or
larger value is asserted.
