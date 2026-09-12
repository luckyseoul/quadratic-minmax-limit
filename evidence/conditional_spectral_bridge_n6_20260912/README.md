# Conditional-cross spectral bridge: exact order-six witness

Classification: exact finite minimax result for one source and one returned
conditional minimizer.  It is not an all-orders statement and does not
classify every order-six conditional minimizer.

For the recorded order-six `Phi(A)=5` source, Soulkiller ran the complete
36-bit CP-SAT cross-signing minimax model against all 4,096 ordered Boolean
pairs using 88 solver workers.  It proved conditional value `18`.  The
returned witness independently replays on the controller with conditional
objective `18`, `beta(B)=14`, and `||B||op=1+sqrt(5)`
(`3.23606797749979`).

The deterministic remote witness JSON has SHA-256
`eb74b9a668dc6a0647cbe72dc018d0c0098e3ebf21b8b2bb7cebe7fc6e3d132f`.
The shared small-order solver script SHA-256 is
`f8a70be97623086d28c666567398fb180f45ab741789b90d5f8dee8592ae5145`.

This point neither proves nor refutes the desired asymptotic
`||B_n||op=O(sqrt(n))` conditional spectral bridge.
