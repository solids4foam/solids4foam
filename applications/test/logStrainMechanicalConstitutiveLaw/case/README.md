# logStrain material-point benchmarks

`Test-logStrainMechanicalConstitutiveLaw` drives the `logStrain` law at the
single integration point of a one-cell mesh, through the manager, over
deformation paths of many steps, so that the law's history is carried from
step to step as in a solid model. It checks that:

1. The tensor logarithm and exponential give their closed-form values, also
   for repeated and nearly repeated eigenvalues, and `exp(log(A)) = A`.
2. With a `linearElastic` sub-law the law is Hencky's hyperelastic law,
   `tau = 2 mu dev(ln V) + K tr(ln V) I`, after a path that stretches and
   rotates.
3. On a path whose principal directions stay fixed, the Kirchhoff stress of
   `logStrain` over `linearElasticMisesPlastic` equals the stress of
   `linearElasticMisesPlastic` itself driven by the logarithmic strain, step
   by step, through yield and a change of loading direction. This is exact
   for the multiplicative split.
4. A rigid rotation superposed after yield only rotates the stress.
5. In simple shear to a shear of 5, with perfect plasticity, the shear stress
   is within 0.1% of `sigmaY/sqrt(3)` and of
   `neoHookeanElasticMisesPlastic`, an independent multiplicative J2 law. The
   curves are written to `simpleShear.csv`.
6. A sub-law that cannot be driven by a trial elastic strain is refused.

The plastic sub-laws read `constant/linearHardening`, which has two rows, so
that the closed-form return mapping is run and check 3 is exact to round-off,
and `constant/perfectPlasticity`.

Run it with:

```bash
./regressionTest.sh
```

`tutorials/Alltest-regression` runs it with the tutorial regression tests.
