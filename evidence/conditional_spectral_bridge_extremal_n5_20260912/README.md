# Conditional-cross spectral bridge: exact extremal order-five screen

This finite diagnostic asks a question not answered by a single conditional
minimizer: among all conditional minimizers for the fixed order-five source,
can a Boolean singular direction have unusually large energy?

For each of the 32 Boolean right vectors `y`, the CP-SAT model fixes the exact
conditional optimum `T=13` and maximizes `||B y||_2^2`.  Soulkiller completed
all 32 optimizations with 88 workers per optimization.  The exact maximum was
53, so the returned witness certifies only

`||B||_op^2 >= 53 / 5 = 10.6`.

The returned witness has `beta(B)=13`, numerical operator norm
`3.3722813232690143`, and independently replays to conditional value 13 and
Boolean-direction energy 53.  The result is a finite lower-bound screen.  It
does not optimize over arbitrary unit singular directions, bound all
conditional minimizers, establish an asymptotic rate, or prove the conditional
spectral bridge.

The solver is
`scripts/conditional_spectral_bridge_extremal_probe.py` at SHA-256
`dd813accc087854d7bd018b633fb2ae528ac19a911fe261330a6017b2fa20a9e`.
The result JSON SHA-256 is
`31d41f2cfd976844b7164a69fd58c3b1fb07e773ff305bc2c3948857b993de43`.
