# mechanicalConstitutiveLaw law checks

`Test-mechanicalConstitutiveLawChecks` checks every mechanical constitutive
law on its own, at a single integration point of a one-cell mesh. Where
`Test-mechanicalConstitutiveLaw` exercises the manager on the laws a tutorial
happens to use, this covers each law whether or not a tutorial selects it.

The laws, with the constants of the tutorials that use them, are listed in
`constant/lawChecks`. For each, a fresh manager is built from a
`mechanicalProperties` holding that one law, and the following are checked:

1. The reference state, `F = I` or no strain, is stress free.
2. The tangent at the reference state, by central differences of the
   stress, has major symmetry.
3. For an isotropic law, the scalar tangent is `lambda + 2 mu` and the
   deviatoric scalar tangent `(4/3) mu`, both taken from that tangent.
4. For a finite-strain law, the Cauchy stress is objective under a
   superposed rotation, at a strain large enough to take the plastic laws
   past yield.
5. For a linear small-strain law, doubling the strain doubles the stress.
6. For a plastic law given `yieldBounds`, the equivalent stress past yield
   lies between the initial yield stress and the last value of the hardening
   table. The plasticity laws read `constant/plasticStrainVsYieldStress`,
   which has more than two rows, so that their nonlinear hardening branch is
   the one run.

Every law is also checked for a positive, finite stiffness and a nonzero
stress under load, so a law that is never reached cannot pass by returning
nothing. `electroMechanicalLaw` is listed twice, with its active tension off
at the reference state and on from the start, so that the objectivity check
sees the active term.

Run it with:

```bash
./regressionTest.sh
```

`tutorials/Alltest-regression` runs it with the tutorial regression tests.
Every law in the runtime selection table must have an entry in
`constant/lawChecks`: a law added without one fails the test rather than going
unchecked.
