# Conditional-cross spectral bridge: order-six extremal formulation stopped

The exact Boolean-direction extremal formulation from
`scripts/conditional_spectral_bridge_extremal_probe.py` was dispatched for
the fixed order-six source (`Phi(A)=5`, conditional target `T=18`).  It was
configured to use 88 CP-SAT workers and a 600-second per-direction limit.

The first of 64 Boolean directions exceeded ten minutes of wall-clock work
without a certificate, a result file, or a solver error.  The process was
then stopped and the Soulkiller allocation was verified released.  Therefore
there is **no order-six extremal value**, witness, or numerical implication to
report.

This formulation is a dead end for extending the n=5 exhaustive screen: its
outer loop would require 64 such exact direction solves.  Any later order-six
work needs a materially stronger formulation or an analytic reduction, rather
than a longer run of this model.
