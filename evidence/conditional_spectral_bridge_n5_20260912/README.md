# Conditional-cross spectral bridge: exact order-five witness

Classification: exact finite minimax result for one source and one returned
conditional minimizer.  It is not an all-orders statement and does not
classify every order-five conditional minimizer.

For the negative-cycle order-five source (`Phi(A)=4`), CP-SAT optimizes the
complete 25-bit cross signing against all 1,024 ordered Boolean pairs.  On
Soulkiller with 88 solver workers it proved the conditional minimum is `13`.
The returned minimizer independently replays on the controller with
`Phi`-style conditional objective `13`, `beta(B)=13`, and
`||B||op=3.3722813232690143`.

The deterministic remote witness JSON has SHA-256
`e6a48f8bcd9d861e1fd62a0252afa2f1194d57942db9802f71c7b9d06eac1723`.
The original order-five run used the pre-generalization script SHA-256
`e651eaf58c5c12291387fa26735471d0d2f363bd89207e3d673396c8ea3e1174`.

This point neither proves nor refutes the desired asymptotic
`||B_n||op=O(sqrt(n))` conditional spectral bridge.
