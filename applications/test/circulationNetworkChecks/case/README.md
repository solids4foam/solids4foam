# circulationNetwork checks

`Test-circulationNetworkChecks` checks the lumped circulation network of
`src/circulationModels/circulation0D` on small circuits whose answers are
known in closed form. No mesh is needed: `regressionTest.sh` runs the
application in an empty directory, where it writes the waveform table it
reads.

The following are checked:

1. A coupling port feeding an inductor: at the end of a step, the port
   pressure is `P_C + L dQ/dt`, with `dQ/dt` taken analytically from the
   forcing polynomial, for steps from `1e-2` down to `1e-7` s.
2. The analytic derivatives of the `sine`, `halfSine` and `table` waveforms
   against central differences.
3. A flow source feeding an inductor through a valve: the pressure of the
   node the inductor holds is `P_C + L dQ/dt`, for a constant and a sine
   flow. This is a pressure row of the nodal solve, which its convergence
   test must include.
4. A flow source feeding a compliance through a hard valve (`epsilon 0`,
   `reverseFraction 0`), from an unpressurised and a pressurised
   compliance, with and without an inductor after the valve, in unit and in
   physiological SI values: the valve opens to `dP = (Q/CQ)^2`. A valve that
   should stay shut does, one that should open shares the flow correctly, and
   a circuit with no solution (a source pulling flow back through the valve)
   stops with an error.
5. A smooth valve (`epsilon 1`) with a regurgitant leak
   (`reverseFraction 0.05`): far from the switch its flow is that of the
   hard law, `CQ sqrt(dP)` forward and `-reverseFraction CQ sqrt(-dP)`
   backward, and with no leak a closed valve passes no flow.

The circulation models are built for OpenFOAM.com only, so the check is
skipped on the other forks.
