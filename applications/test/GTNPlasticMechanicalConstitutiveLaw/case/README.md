# GTN material-point checks

`Test-GTNPlasticMechanicalConstitutiveLaw` drives `linearElasticGTNPlastic`
directly, at one point, over strain paths of many steps, so that its void
and plastic history can be read back. It checks that:

1. With no voids, the law is `linearElasticMisesPlastic`, step for step.
2. On a path of tensile mean stress, the returned state satisfies the
   discrete yield condition, the normality of the flow, the equivalence of
   plastic work and the void growth with nucleation, recomputed from the
   stress; with linear and nonlinear hardening, with the shear term on a
   path with shear, and with coalescence.
3. Under hydrostatic tension the voids grow, and under compression they
   close, the equations of 2 holding.
4. A point whose void fraction would pass the failed void fraction fails,
   stays failed, and keeps about 1% of its strength; and a step whose void
   growth equation has a root near the start of the step returns it, even
   where the equation changes sign again short of failure.
   `constant/perfectPlasticity` gives the constant yield stress of that
   case.
5. The void fraction converges at first order as the step is refined.
6. Under `logStrain`, on a path whose principal directions stay fixed, the
   Kirchhoff stress equals the small-strain law's at the logarithmic strain.
7. Inconsistent entries are refused.
8. In the non-local form, with the non-local void fraction supplied as an
   input: with it zero the stress is `linearElasticMisesPlastic`'s while the
   local voids nucleate; with it iterated to the local void fraction within
   each step, the local form is recovered; and the law declares the field.

The laws read `constant/linearHardening`, with two rows, and
`constant/nonlinearHardening`, a table of several.

Run it with:

```bash
./regressionTest.sh
```

`tutorials/Alltest-regression` runs it with the tutorial regression tests.
